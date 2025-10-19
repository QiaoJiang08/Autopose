"""
HADDOCK 2.4 Adapter for Autopose
Parses HADDOCK protein-protein docking output files and prepares data for animation.
Fixed version addressing all bugs from Bugsv1Haddock.md
"""

import os
import re
import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa
from scipy.spatial.transform import Rotation as R


class HADDOCKAdapter:
    def __init__(self, mode='auto', enable_pedagogical_snap=False):
        """
        Initialize HADDOCK adapter with mode selection.
        
        Args:
            mode: 'auto', 'production', or 'demo'
                  'auto' - detect from input data
                  'production' - real HADDOCK data, exact geometry
                  'demo' - allow CA-level contacts and optional snapping
            enable_pedagogical_snap: Allow small nudge toward interface in demo mode
        """
        self.parser = PDBParser(QUIET=True)
        self.mode = mode
        self.enable_pedagogical_snap = enable_pedagogical_snap
        self.detected_mode = None  # Will be set after input validation
        
    def parse_cluster_energy(self, cluster_ener_path):
        """
        Parse cluster_ener.txt or cluster_haddock.txt for energy information.
        
        Expected format:
        Cluster Rank HADDOCKScore Evdw Eelec ...
        1       1    -123.4       -45.6 -78.9
        """
        cluster_data = []
        
        with open(cluster_ener_path, 'r') as f:
            lines = f.readlines()
        
        # Find header line
        header_idx = 0
        for i, line in enumerate(lines):
            if 'Cluster' in line and ('Rank' in line or 'HADDOCK' in line):
                header_idx = i
                break
        
        # Parse data lines
        for line in lines[header_idx + 1:]:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            parts = line.split()
            if len(parts) >= 3:
                try:
                    cluster_id = int(parts[0])
                    rank = int(parts[1])
                    haddock_score = float(parts[2])
                    
                    # Extract additional energy terms if available
                    evdw = float(parts[3]) if len(parts) > 3 else None
                    eelec = float(parts[4]) if len(parts) > 4 else None
                    edesol = float(parts[5]) if len(parts) > 5 else None
                    eair = float(parts[6]) if len(parts) > 6 else None
                    
                    cluster_data.append({
                        'cluster': cluster_id,
                        'rank': rank,
                        'haddock_score': haddock_score,
                        'evdw': evdw,
                        'eelec': eelec,
                        'edesol': edesol,
                        'eair': eair
                    })
                except (ValueError, IndexError):
                    continue
        
        return cluster_data
    
    def parse_cluster_rmsd(self, cluster_rmsd_path):
        """
        Parse cluster_rmsd.txt for RMSD information.
        
        Expected format:
        Cluster RMSD_from_overall_lowest
        1       0.0
        2       5.2
        """
        rmsd_data = {}
        
        with open(cluster_rmsd_path, 'r') as f:
            lines = f.readlines()
        
        # Find header line
        header_idx = 0
        for i, line in enumerate(lines):
            if 'Cluster' in line and 'RMSD' in line:
                header_idx = i
                break
        
        # Parse data lines
        for line in lines[header_idx + 1:]:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            parts = line.split()
            if len(parts) >= 2:
                try:
                    cluster_id = int(parts[0])
                    rmsd = float(parts[1])
                    rmsd_data[cluster_id] = rmsd
                except (ValueError, IndexError):
                    continue
        
        return rmsd_data
    
    def validate_model_quality(self, model_data):
        """
        Validate model quality and detect if it's sparse/CA-only or has proper interface.
        Returns dict with validation results.
        """
        validation = {
            'is_sparse': False,
            'is_ca_only': False,
            'is_interfacial': False,
            'atoms_per_chain_a': len(model_data['chain_a']),
            'atoms_per_chain_b': len(model_data['chain_b']),
            'min_interchain_distance': None,
            'heavy_atom_contacts': 0,
            'ca_contacts': 0,
            'suggested_mode': 'production'
        }
        
        # Check 1: Atom completeness
        avg_atoms_per_chain = (validation['atoms_per_chain_a'] + validation['atoms_per_chain_b']) / 2
        
        # CA-only detection: ~4-5 atoms per residue for CA-only
        if avg_atoms_per_chain < 300:
            validation['is_sparse'] = True
            if avg_atoms_per_chain < 100:
                validation['is_ca_only'] = True
        
        # Check 2: Interchain geometry
        min_dist = self.compute_interchain_distance(model_data['chain_a'], model_data['chain_b'])
        validation['min_interchain_distance'] = min_dist
        
        # Check if interfacial (min distance < 6-8 Å)
        if min_dist < 8.0:
            validation['is_interfacial'] = True
        
        # Check 3: Count contacts at different cutoffs
        # Heavy-atom contacts (4.5 Å)
        n_heavy_a, n_heavy_b = self.detect_interface_residues(
            model_data['chain_a'], model_data['chain_b'], cutoff=4.5
        )
        validation['heavy_atom_contacts'] = n_heavy_a + n_heavy_b
        
        # CA-CA contacts (8.0 Å) for demo mode
        if len(model_data['chain_a_ca']) > 0 and len(model_data['chain_b_ca']) > 0:
            n_ca_a, n_ca_b = self.detect_interface_residues(
                model_data['chain_a_ca'], model_data['chain_b_ca'], cutoff=8.0
            )
            validation['ca_contacts'] = n_ca_a + n_ca_b
        
        # Suggest mode
        if validation['is_sparse'] or not validation['is_interfacial']:
            validation['suggested_mode'] = 'demo'
        else:
            validation['suggested_mode'] = 'production'
        
        return validation
    
    def parse_haddock_model(self, model_pdb_path):
        """
        Parse a HADDOCK model PDB file and extract chain A and B coordinates.
        Returns both chains as separate structures.
        """
        structure = self.parser.get_structure('haddock_model', model_pdb_path)
        
        chain_a_atoms = []
        chain_b_atoms = []
        chain_a_residues = []
        chain_b_residues = []
        chain_a_ca_atoms = []
        chain_b_ca_atoms = []
        
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                
                for residue in chain:
                    if is_aa(residue):
                        residue_atoms = []
                        for atom in residue:
                            coords = atom.get_coord()
                            residue_atoms.append(coords)
                        
                        if chain_id == 'A':
                            chain_a_atoms.extend(residue_atoms)
                            chain_a_residues.append(residue)
                            if 'CA' in residue:
                                chain_a_ca_atoms.append(residue['CA'].get_coord())
                        elif chain_id == 'B':
                            chain_b_atoms.extend(residue_atoms)
                            chain_b_residues.append(residue)
                            if 'CA' in residue:
                                chain_b_ca_atoms.append(residue['CA'].get_coord())
        
        return {
            'chain_a': np.array(chain_a_atoms, dtype=float),
            'chain_b': np.array(chain_b_atoms, dtype=float),
            'chain_a_ca': np.array(chain_a_ca_atoms, dtype=float),
            'chain_b_ca': np.array(chain_b_ca_atoms, dtype=float),
            'chain_a_residues': chain_a_residues,
            'chain_b_residues': chain_b_residues,
            'structure': structure
        }
    
    def compute_kabsch_transform(self, ref_coords, mobile_coords):
        """
        Compute Kabsch transformation (rotation + translation) from mobile to ref.
        Returns rotation matrix and translation vector.
        
        Note: ref_coords and mobile_coords must have the same number of points.
        """
        # Ensure same number of points
        if len(ref_coords) != len(mobile_coords):
            # Take minimum number of points
            min_len = min(len(ref_coords), len(mobile_coords))
            print(f"Warning: Coordinate count mismatch. Ref: {len(ref_coords)}, Mobile: {len(mobile_coords)}. Using first {min_len} points.")
            ref_coords = ref_coords[:min_len]
            mobile_coords = mobile_coords[:min_len]
        
        # Center both coordinate sets
        ref_center = np.mean(ref_coords, axis=0)
        mobile_center = np.mean(mobile_coords, axis=0)
        
        ref_centered = ref_coords - ref_center
        mobile_centered = mobile_coords - mobile_center
        
        # Compute covariance matrix
        H = mobile_centered.T @ ref_centered
        
        # SVD to find rotation
        U, S, Vt = np.linalg.svd(H)
        
        # Ensure proper rotation (det = +1)
        d = np.linalg.det(Vt.T @ U.T)
        if d < 0:
            Vt[-1, :] *= -1
        
        rotation = Vt.T @ U.T
        translation = ref_center - rotation @ mobile_center
        
        return rotation, translation
    
    def apply_transform(self, coords, rotation, translation):
        """Apply rigid-body transformation to coordinates."""
        return coords @ rotation.T + translation
    
    def compute_rmsd(self, coords1, coords2):
        """Compute RMSD between two coordinate sets."""
        if len(coords1) != len(coords2):
            # Use shorter length
            min_len = min(len(coords1), len(coords2))
            coords1 = coords1[:min_len]
            coords2 = coords2[:min_len]
        
        if len(coords1) == 0:
            return None
            
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    def detect_interface_residues(self, chain_a_coords, chain_b_coords, cutoff=4.5):
        """
        Detect interface residues based on distance cutoff.
        Returns counts of atoms in interface.
        """
        interface_a = set()
        interface_b = set()
        
        for i, coord_a in enumerate(chain_a_coords):
            for j, coord_b in enumerate(chain_b_coords):
                dist = np.linalg.norm(coord_a - coord_b)
                if dist <= cutoff:
                    interface_a.add(i)
                    interface_b.add(j)
        
        return len(interface_a), len(interface_b)
    
    def compute_interchain_distance(self, chain_a_coords, chain_b_coords):
        """Compute minimum distance between chains."""
        min_dist = float('inf')
        for coord_a in chain_a_coords:
            for coord_b in chain_b_coords:
                dist = np.linalg.norm(coord_a - coord_b)
                if dist < min_dist:
                    min_dist = dist
        return min_dist
    
    def generate_docking_trajectory_with_rotation(self, chain_a_coords, chain_b_coords_target, 
                                                  chain_b_coords_start, num_frames=60,
                                                  enable_pedagogical_snap=False):
        """
        Generate a realistic docking trajectory with proper rigid-body rotation.
        Now includes deterministic initial rotation for visible motion.
        
        Args:
            chain_a_coords: Fixed chain A coordinates
            chain_b_coords_target: Final docked position of chain B (after alignment)
            chain_b_coords_start: Start position of chain B (displaced from target)
            num_frames: Number of frames to generate (default 60 for longer path)
            enable_pedagogical_snap: Apply small nudge toward interface in final frames (demo mode)
        
        Returns:
            List of coordinate arrays for chain B trajectory
        """
        # Compute centers
        center_a = np.mean(chain_a_coords, axis=0)
        center_b_target = np.mean(chain_b_coords_target, axis=0)
        center_b_start = np.mean(chain_b_coords_start, axis=0)
        
        # FIX B1: Apply deterministic 35° rotation to start pose for visible rotation
        approach_vec = center_b_target - center_a
        if np.linalg.norm(approach_vec) > 0:
            approach_vec = approach_vec / np.linalg.norm(approach_vec)
        else:
            approach_vec = np.array([1.0, 0.0, 0.0])
        
        # Create deterministic rotation (35 degrees around perpendicular axis)
        perpendicular_axis = np.cross(approach_vec, np.array([0, 0, 1]))
        if np.linalg.norm(perpendicular_axis) < 0.1:
            perpendicular_axis = np.cross(approach_vec, np.array([1, 0, 0]))
        perpendicular_axis = perpendicular_axis / np.linalg.norm(perpendicular_axis)
        
        initial_rotation_angle = np.radians(35)  # 35 degrees
        initial_rotation = R.from_rotvec(initial_rotation_angle * perpendicular_axis)
        
        # Apply initial rotation to start pose
        chain_b_start_rotated = initial_rotation.apply(chain_b_coords_start - center_b_start) + center_b_start
        
        # Compute rigid-body transform from rotated start to target
        rotation_matrix, translation = self.compute_kabsch_transform(
            chain_b_coords_target, chain_b_start_rotated
        )
        
        # Convert rotation matrix to scipy Rotation for interpolation
        rot_start = R.from_matrix(np.eye(3))
        rot_target = R.from_matrix(rotation_matrix)
        
        trajectory = []
        
        for i in range(num_frames):
            t = i / (num_frames - 1)
            
            # Smooth interpolation (cubic ease-in-out)
            smooth_t = 3 * t**2 - 2 * t**3
            
            # Interpolate rotation using spherical linear interpolation (slerp)
            rot_current = R.from_quat(
                self._slerp(rot_start.as_quat(), rot_target.as_quat(), smooth_t)
            )
            
            # Interpolate translation
            translation_current = translation * smooth_t
            
            # Apply interpolated transform
            center_b_start_array = center_b_start.reshape(1, 3)
            current_coords = self.apply_transform(
                chain_b_start_rotated - center_b_start,
                rot_current.as_matrix(),
                center_b_start + translation_current
            )
            
            # FIX B5: Optional pedagogical snap (demo mode only, last 2-3 frames)
            if enable_pedagogical_snap and i >= num_frames - 3:
                # Apply small nudge toward interface (max 2 Å)
                snap_t = (i - (num_frames - 3)) / 2.0  # 0 to 1 over last 3 frames
                snap_distance = 2.0 * snap_t  # Up to 2 Å
                snap_direction = center_a - np.mean(current_coords, axis=0)
                if np.linalg.norm(snap_direction) > 0:
                    snap_direction = snap_direction / np.linalg.norm(snap_direction)
                    current_coords = current_coords + snap_direction * snap_distance
            
            trajectory.append(current_coords)
            
            # Diagnostic for final frame
            if i == num_frames - 1:
                rmsd_final = self.compute_rmsd(current_coords, chain_b_coords_target)
                print(f"DEBUG: Final frame RMSD to target: {rmsd_final:.6f} Å")
                min_dist = self.compute_interchain_distance(chain_a_coords, current_coords)
                print(f"DEBUG: Final frame min interchain distance: {min_dist:.2f} Å")
        
        return trajectory
    
    def _slerp(self, q0, q1, t):
        """Spherical linear interpolation between two quaternions."""
        dot = np.dot(q0, q1)
        
        # If quaternions are close, use linear interpolation
        if abs(dot) > 0.9995:
            result = q0 + t * (q1 - q0)
            return result / np.linalg.norm(result)
        
        # Clamp dot product
        dot = np.clip(dot, -1.0, 1.0)
        
        # If dot < 0, negate one quaternion to take shorter path
        if dot < 0.0:
            q1 = -q1
            dot = -dot
        
        theta_0 = np.arccos(dot)
        theta = theta_0 * t
        
        q2 = q1 - q0 * dot
        q2 = q2 / np.linalg.norm(q2)
        
        return q0 * np.cos(theta) + q2 * np.sin(theta)
    
    def create_animation_data(self, monomer_pdb, model_pdbs, cluster_ener_path, 
                            cluster_rmsd_path, top_k=5):
        """
        Create animation data from HADDOCK output files.
        
        Args:
            monomer_pdb: Path to reference monomer PDB
            model_pdbs: List of paths to docked model PDB files
            cluster_ener_path: Path to cluster_ener.txt
            cluster_rmsd_path: Path to cluster_rmsd.txt
            top_k: Number of top models to include
            
        Returns:
            Dictionary with animation data for all models
        """
        # Parse cluster data
        cluster_energy = self.parse_cluster_energy(cluster_ener_path)
        cluster_rmsd = self.parse_cluster_rmsd(cluster_rmsd_path)
        
        # Sort by rank and take top K
        cluster_energy.sort(key=lambda x: x['rank'])
        top_clusters = cluster_energy[:top_k]
        
        # Parse reference monomer
        monomer_data = self.parse_haddock_model(monomer_pdb)
        ref_chain_a_ca = monomer_data['chain_a_ca']
        
        # Parse models and store top-ranked model for RMSD comparison
        animation_data = {
            'monomer': monomer_pdb,
            'models': [],
            'top_model_b_aligned': None  # Will be set to first model's aligned chain B
        }
        
        for i, model_path in enumerate(model_pdbs[:top_k]):
            if not os.path.exists(model_path):
                print(f"Warning: Model file {model_path} not found, skipping...")
                continue
            
            # Parse model
            model_data = self.parse_haddock_model(model_path)
            
            # Validate model quality and auto-detect mode
            validation = self.validate_model_quality(model_data)
            
            # Auto-detect mode on first model
            if i == 0:
                if self.mode == 'auto':
                    self.detected_mode = validation['suggested_mode']
                else:
                    self.detected_mode = self.mode
                
                print(f"\n{'='*60}")
                print(f"MODE DETECTION:")
                print(f"  Detected mode: {self.detected_mode}")
                print(f"  Atoms/chain: A={validation['atoms_per_chain_a']}, B={validation['atoms_per_chain_b']}")
                print(f"  Min interchain distance: {validation['min_interchain_distance']:.2f} Å")
                print(f"  Heavy-atom contacts (4.5Å): {validation['heavy_atom_contacts']}")
                print(f"  CA-CA contacts (8.0Å): {validation['ca_contacts']}")
                print(f"  Is sparse/CA-only: {validation['is_sparse']}/{validation['is_ca_only']}")
                print(f"  Is interfacial: {validation['is_interfacial']}")
                print(f"{'='*60}\n")
            
            # Store validation info for animation
            animation_data.setdefault('validation', validation)
            animation_data.setdefault('mode', self.detected_mode)
            
            # FIX B1: Compute transform for chain A and apply to BOTH chains
            print(f"\nProcessing model {i+1}...")
            rotation, translation = self.compute_kabsch_transform(
                ref_chain_a_ca,
                model_data['chain_a_ca']
            )
            
            # Apply transform to both chains to keep dimer internally consistent
            chain_a_aligned = self.apply_transform(model_data['chain_a'], rotation, translation)
            chain_b_aligned = self.apply_transform(model_data['chain_b'], rotation, translation)
            
            # Verify alignment
            rmsd_a = self.compute_rmsd(chain_a_aligned, monomer_data['chain_a'])
            if rmsd_a is not None:
                print(f"DEBUG: Chain A alignment RMSD: {rmsd_a:.6f} Å")
            else:
                print(f"DEBUG: Chain A alignment RMSD: Cannot compute (length mismatch)")
            
            # Store top model's chain B for RMSD calculations (FIX B2)
            if i == 0:
                animation_data['top_model_b_aligned'] = chain_b_aligned.copy()
            
            # FIX B2: Compute proper RMSD to top-ranked model
            if animation_data['top_model_b_aligned'] is not None and i > 0:
                rmsd_to_top = self.compute_rmsd(chain_b_aligned, animation_data['top_model_b_aligned'])
            else:
                rmsd_to_top = 0.0  # Top model has 0 RMSD to itself
            
            # Generate start position (displaced from target)
            center_a = np.mean(chain_a_aligned, axis=0)
            center_b = np.mean(chain_b_aligned, axis=0)
            approach_vec = center_b - center_a
            approach_distance = np.linalg.norm(approach_vec)
            
            if approach_distance > 0:
                approach_vec = approach_vec / approach_distance
            else:
                approach_vec = np.array([1.0, 0.0, 0.0])
            
            # Start position: 35-40 Å away from final position (longer path)
            start_distance = 35.0
            center_b_start = center_b + start_distance * approach_vec
            translation_to_start = center_b_start - center_b
            chain_b_start = chain_b_aligned + translation_to_start
            
            # Determine if pedagogical snap should be enabled
            enable_snap = (self.detected_mode == 'demo' and 
                          self.enable_pedagogical_snap and 
                          not validation['is_interfacial'])
            
            if enable_snap:
                print(f"  Demo mode: Pedagogical snap enabled (non-interfacial model)")
            
            # Generate trajectory with proper rotation (FIX B4, 60-80 frames for visible rotation)
            trajectory = self.generate_docking_trajectory_with_rotation(
                chain_a_aligned,
                chain_b_aligned,  # Target position (aligned)
                chain_b_start,    # Start position (displaced)
                num_frames=60,
                enable_pedagogical_snap=enable_snap
            )
            
            # Get cluster info
            cluster_info = top_clusters[i] if i < len(top_clusters) else {
                'cluster': i + 1,
                'rank': i + 1,
                'haddock_score': None,
                'evdw': None,
                'eelec': None,
                'edesol': None,
                'eair': None
            }
            
            # Get RMSD from cluster file
            cluster_id = cluster_info['cluster']
            rmsd_from_file = cluster_rmsd.get(cluster_id, None)
            
            # FIX B5: Detect interface on final aligned pose (dual metrics)
            # Heavy-atom contacts (4.5 Å)
            n_heavy_a, n_heavy_b = self.detect_interface_residues(
                chain_a_aligned,
                chain_b_aligned,
                cutoff=4.5
            )
            
            # CA-CA contacts (8.0 Å) for demo mode
            n_ca_a, n_ca_b = 0, 0
            if len(model_data['chain_a_ca']) > 0 and len(model_data['chain_b_ca']) > 0:
                chain_a_ca_aligned = self.apply_transform(model_data['chain_a_ca'], rotation, translation)
                chain_b_ca_aligned = self.apply_transform(model_data['chain_b_ca'], rotation, translation)
                n_ca_a, n_ca_b = self.detect_interface_residues(
                    chain_a_ca_aligned,
                    chain_b_ca_aligned,
                    cutoff=8.0
                )
            
            animation_data['models'].append({
                'chain_a': chain_a_aligned,
                'chain_b': chain_b_aligned,  # Store aligned final position
                'trajectory': trajectory,
                'info': {
                    'cluster': cluster_info['cluster'],
                    'rank': cluster_info['rank'],
                    'haddock_score': cluster_info['haddock_score'],
                    'evdw': cluster_info['evdw'],
                    'eelec': cluster_info['eelec'],
                    'edesol': cluster_info['edesol'],
                    'eair': cluster_info['eair'],
                    'rmsd': rmsd_from_file,  # RMSD from file
                    'rmsd_to_top': rmsd_to_top,  # Computed RMSD to top model
                    'heavy_contacts': n_heavy_a + n_heavy_b,  # Heavy-atom contacts
                    'ca_contacts': n_ca_a + n_ca_b,  # CA-CA contacts (8Å)
                    'min_interchain_distance': validation['min_interchain_distance'],
                    'interface_contacts_a': n_heavy_a,  # For backward compatibility
                    'interface_contacts_b': n_heavy_b
                }
            })
        
        return animation_data


if __name__ == "__main__":
    # Test the adapter
    adapter = HADDOCKAdapter()
    print("HADDOCKAdapter initialized successfully!")
