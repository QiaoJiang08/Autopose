import os
import re
import csv
import numpy as np
from Bio.PDB import PDBParser, Selection
from Bio.PDB.Polypeptide import is_aa

class ClusProAdapter:
    """Adapter for ClusPro protein-protein docking results"""
    
    def __init__(self):
        self.parser = PDBParser(QUIET=True)
        
    def parse_cluspro_scores(self, csv_path):
        """Parse ClusPro scores CSV file"""
        scores = []
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    # Handle different CSV formats
                    if 'Rank' in row and 'Score' in row:
                        scores.append({
                            'rank': int(row['Rank']),
                            'cluster': int(row.get('Cluster', 0)),
                            'cluster_size': int(row.get('Size', 1)),
                            'score': float(row['Score']),
                            'vdw': float(row.get('vdW', 0)),
                            'electrostatic': float(row.get('Electrostatic', 0)),
                            'desolvation': float(row.get('Desolvation', 0))
                        })
                    elif 'cluster' in row and 'score' in row:
                        scores.append({
                            'rank': int(row.get('rank', 0)),
                            'cluster': int(row['cluster']),
                            'cluster_size': int(row.get('size', 1)),
                            'score': float(row['score']),
                            'vdw': float(row.get('vdw', 0)),
                            'electrostatic': float(row.get('electrostatic', 0)),
                            'desolvation': float(row.get('desolvation', 0))
                        })
        except Exception as e:
            raise Exception(f"Error parsing ClusPro scores CSV: {str(e)}")
        return scores
    
    def parse_cluspro_model(self, pdb_path, monomer_pdb_path):
        """Parse ClusPro model PDB and extract chain transformations"""
        try:
            # Parse the docked complex
            complex_structure = self.parser.get_structure('complex', pdb_path)
            monomer_structure = self.parser.get_structure('monomer', monomer_pdb_path)
            
            # Extract chains from complex
            complex_chains = list(complex_structure.get_chains())
            if len(complex_chains) < 2:
                raise Exception("Complex PDB must contain at least 2 chains")
            
            # Assume first chain is the reference (same as monomer)
            # Second chain is the moving chain that docked
            ref_chain = complex_chains[0]
            moving_chain = complex_chains[1]
            
            # Extract coordinates for both chains
            ref_coords = []
            moving_coords = []
            ref_elements = []
            moving_elements = []
            
            for residue in ref_chain:
                if is_aa(residue):
                    for atom in residue:
                        coords = atom.get_coord()
                        ref_coords.append(coords)
                        ref_elements.append(atom.element)
            
            for residue in moving_chain:
                if is_aa(residue):
                    for atom in residue:
                        coords = atom.get_coord()
                        moving_coords.append(coords)
                        moving_elements.append(atom.element)
            
            ref_coords = np.array(ref_coords, dtype=float)
            moving_coords = np.array(moving_coords, dtype=float)
            
            # Calculate transformation from monomer to docked position
            transformation = self._calculate_rigid_body_transform(ref_coords, moving_coords)
            
            return {
                'ref_chain': {
                    'coords': ref_coords,
                    'elements': ref_elements,
                    'chain_id': ref_chain.id
                },
                'moving_chain': {
                    'coords': moving_coords,
                    'elements': moving_elements,
                    'chain_id': moving_chain.id
                },
                'transformation': transformation
            }
            
        except Exception as e:
            raise Exception(f"Error parsing ClusPro model: {str(e)}")
    
    def _calculate_rigid_body_transform(self, coords1, coords2):
        """Calculate rigid body transformation between two sets of coordinates"""
        # Center both coordinate sets
        center1 = np.mean(coords1, axis=0)
        center2 = np.mean(coords2, axis=0)
        
        centered1 = coords1 - center1
        centered2 = coords2 - center2
        
        # Calculate rotation matrix using SVD
        H = centered1.T @ centered2
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        
        # Ensure proper rotation matrix
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        # Calculate translation
        t = center2 - center1
        
        return {
            'rotation': R,
            'translation': t,
            'center1': center1,
            'center2': center2
        }
    
    def detect_interface_residues(self, complex_data, cutoff=4.5):
        """Detect interface residues between chains"""
        ref_coords = complex_data['ref_chain']['coords']
        moving_coords = complex_data['moving_chain']['coords']
        
        # Calculate distances between all atoms
        interface_atoms = []
        interface_residues = set()
        
        for i, ref_coord in enumerate(ref_coords):
            for j, mov_coord in enumerate(moving_coords):
                distance = np.linalg.norm(ref_coord - mov_coord)
                if distance <= cutoff:
                    interface_atoms.append((i, j, distance))
                    # Add residue information if available
                    interface_residues.add(f"REF_{i}")
                    interface_residues.add(f"MOV_{j}")
        
        return {
            'interface_atoms': interface_atoms,
            'interface_residues': list(interface_residues),
            'contact_count': len(interface_atoms)
        }
    
    def generate_docking_trajectory(self, complex_data, num_frames=120, start_distance=20.0):
        """Generate docking trajectory from separated to docked state"""
        ref_coords = complex_data['ref_chain']['coords']
        moving_coords = complex_data['moving_chain']['coords']
        transformation = complex_data['transformation']
        
        # Calculate approach vector
        ref_center = np.mean(ref_coords, axis=0)
        mov_center = np.mean(moving_coords, axis=0)
        approach_vector = mov_center - ref_center
        
        if np.linalg.norm(approach_vector) == 0:
            approach_vector = np.array([1.0, 0.0, 0.0])
        approach_vector = approach_vector / np.linalg.norm(approach_vector)
        
        # Generate start position (far from reference)
        start_center = ref_center + start_distance * approach_vector
        
        # Create trajectory frames
        trajectory = []
        for frame in range(num_frames):
            t = frame / (num_frames - 1)
            
            # Smooth easing function (cubic ease-in-out)
            smooth_t = 3 * t**2 - 2 * t**3
            
            # Interpolate center position
            current_center = start_center * (1 - smooth_t) + mov_center * smooth_t
            
            # Interpolate rotation (simplified)
            if t < 0.5:
                # First half: approach with minimal rotation
                current_coords = moving_coords + (current_center - mov_center)
            else:
                # Second half: apply final transformation
                rotation_t = (t - 0.5) * 2  # 0 to 1 in second half
                smooth_rot_t = 3 * rotation_t**2 - 2 * rotation_t**3
                
                # Interpolate rotation
                R_interp = (1 - smooth_rot_t) * np.eye(3) + smooth_rot_t * transformation['rotation']
                
                # Apply interpolated transformation
                centered_coords = moving_coords - transformation['center1']
                rotated_coords = centered_coords @ R_interp.T
                current_coords = rotated_coords + current_center
            
            trajectory.append({
                'coords': current_coords,
                'elements': complex_data['moving_chain']['elements'],
                'chain_id': complex_data['moving_chain']['chain_id'],
                'frame': frame,
                'progress': t
            })
        
        return trajectory
    
    def create_animation_data(self, monomer_pdb_path, model_paths, scores_csv_path, top_k=3):
        """Create complete animation data from ClusPro outputs"""
        try:
            # Parse scores
            scores = self.parse_cluspro_scores(scores_csv_path)
            scores = scores[:top_k]  # Take top K models
            
            # Parse models
            models_data = []
            for i, score in enumerate(scores):
                if i < len(model_paths):
                    model_data = self.parse_cluspro_model(model_paths[i], monomer_pdb_path)
                    interface_data = self.detect_interface_residues(model_data)
                    
                    models_data.append({
                        'model_index': i,
                        'score_data': score,
                        'complex_data': model_data,
                        'interface_data': interface_data,
                        'trajectory': self.generate_docking_trajectory(model_data)
                    })
            
            return {
                'monomer_path': monomer_pdb_path,
                'models': models_data,
                'total_models': len(models_data)
            }
            
        except Exception as e:
            raise Exception(f"Error creating animation data: {str(e)}")

if __name__ == "__main__":
    # Test the adapter
    adapter = ClusProAdapter()
    print("ClusProAdapter initialized successfully!") 