import os
import re
import numpy as np
import matplotlib
import random
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.vectors import Vector
from Bio.PDB.DSSP import DSSP
import tempfile
import subprocess

class DockingAnimatorEnhanced:
    def __init__(self, representation='sticks'):
        self.parser = PDBParser(QUIET=True)
        self.representation = representation  # 'sticks', 'ribbon', 'lines', 'spheres', 'surface', 'cartoon'
        
    def parse_pdb(self, pdb_path):
        """Parse PDB file and extract protein structure with complete bond information"""
        try:
            structure = self.parser.get_structure('protein', pdb_path)
            
            # Extract atoms and bonds
            protein_atoms = []
            protein_bonds = []
            atom_colors = []
            atom_names = []
            residue_info = []
            secondary_structure = {}
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if is_aa(residue):
                            residue_atoms = []
                            residue_atom_names = []
                            residue_atom_elements = []
                            
                            for atom in residue:
                                # Get atom coordinates
                                coords = atom.get_coord()
                                residue_atoms.append(coords)
                                residue_atom_names.append(atom.get_name())
                                residue_atom_elements.append(atom.element)
                                
                                # Assign colors based on element (PyMOL standard)
                                element = atom.element
                                if element == 'C':
                                    color = '#909090'  # Gray for carbon (PyMOL default)
                                elif element == 'N':
                                    color = '#3030FF'  # Blue for nitrogen
                                elif element == 'O':
                                    color = '#FF3030'  # Red for oxygen
                                elif element == 'S':
                                    color = '#FFFF30'  # Yellow for sulfur
                                elif element == 'P':
                                    color = '#FFA500'  # Orange for phosphorus
                                elif element == 'F':
                                    color = '#90FF90'  # Light green for fluorine
                                elif element == 'Cl':
                                    color = '#30FF30'  # Green for chlorine
                                else:
                                    color = '#FF1493'  # Pink for other elements
                                
                                atom_colors.append(color)
                            
                            # Add atoms to protein
                            if residue_atoms:
                                start_idx = len(protein_atoms)
                                protein_atoms.extend(residue_atoms)
                                atom_names.extend(residue_atom_names)
                                
                                # Create complete bonds within residue
                                self._create_residue_bonds(residue_atom_names, start_idx, protein_bonds)
                                
                                # Store residue information for ribbon/cartoon
                                residue_info.append({
                                    'start_idx': start_idx,
                                    'end_idx': len(protein_atoms) - 1,
                                    'atoms': residue_atoms,
                                    'names': residue_atom_names,
                                    'residue': residue
                                })
            
            # Create inter-residue bonds (peptide bonds)
            self._create_peptide_bonds(protein_atoms, atom_names, protein_bonds)
            
            # Try to get secondary structure information
            try:
                secondary_structure = self._get_secondary_structure(structure)
            except:
                # If DSSP fails, assign default secondary structure
                secondary_structure = self._assign_default_secondary_structure(residue_info)
            
            return np.array(protein_atoms, dtype=float), protein_bonds, atom_colors, residue_info, secondary_structure
            
        except Exception as e:
            raise Exception(f"Error parsing PDB file: {str(e)}")
    
    def _get_secondary_structure(self, structure):
        """Get secondary structure information using DSSP"""
        secondary_structure = {}
        try:
            for model in structure:
                for chain in model:
                    try:
                        dssp = DSSP(model, chain, 'dssp')
                        if dssp:
                            for key in dssp.keys():
                                residue_id = key[1]
                                ss = dssp[key][2]
                                secondary_structure[residue_id] = ss
                    except:
                        pass
        except:
            pass
        return secondary_structure
    
    def _assign_default_secondary_structure(self, residue_info):
        """Assign default secondary structure based on residue position"""
        secondary_structure = {}
        for i, residue in enumerate(residue_info):
            # Simple heuristic: alternate between helix and sheet
            if i % 3 == 0:
                secondary_structure[i] = 'H'  # Helix
            elif i % 3 == 1:
                secondary_structure[i] = 'E'  # Sheet
            else:
                secondary_structure[i] = 'C'  # Coil
        return secondary_structure
    
    def _create_residue_bonds(self, atom_names, start_idx, protein_bonds):
        """Create bonds within a residue"""
        for i, atom_name in enumerate(atom_names):
            if atom_name == 'N' and i + 1 < len(atom_names):
                # N-CA bond
                protein_bonds.append((start_idx + i, start_idx + i + 1))
            elif atom_name == 'CA' and i + 1 < len(atom_names):
                # CA-C bond
                protein_bonds.append((start_idx + i, start_idx + i + 1))
            elif atom_name == 'C' and i + 1 < len(atom_names):
                # C-O bond
                protein_bonds.append((start_idx + i, start_idx + i + 1))
            elif atom_name == 'CA' and 'CB' in atom_names:
                # CA-CB bond
                cb_idx = atom_names.index('CB')
                protein_bonds.append((start_idx + i, start_idx + cb_idx))
    
    def _create_peptide_bonds(self, protein_atoms, atom_names, protein_bonds):
        """Create peptide bonds between residues"""
        for i in range(len(protein_atoms) - 1):
            if atom_names[i] == 'C' and atom_names[i + 1] == 'N':
                protein_bonds.append((i, i + 1))
    
    def parse_dlg(self, dlg_path):
        """Parse DLG file and extract ligand poses with detailed docking information"""
        try:
            poses = []
            
            with open(dlg_path, 'r') as f:
                content = f.read()
            
            # Normalize line endings
            content = content.replace('\r\n', '\n')
            
            # Extract global docking information from the entire file
            global_info = self._extract_global_docking_info(content)
            
            # Find all MODEL blocks first
            model_pattern = r'DOCKED: MODEL\s+(\d+)\s*\n(.*?)(?=DOCKED: MODEL|\Z)'
            model_matches = re.findall(model_pattern, content, re.DOTALL)
            
            # Get all REMARK lines at the beginning of the file (for first model)
            before_first_model = content[:content.find('DOCKED: MODEL')]
            initial_remarks = re.findall(r'REMARK.*?\n', before_first_model, re.DOTALL)
            initial_remark_data = ''.join(initial_remarks)
            
            for i, (model_num, model_content) in enumerate(model_matches):
                # For first model, use initial remarks; for others, look for remarks before this model
                if i == 0:
                    remark_data = initial_remark_data
                else:
                    # Find remarks after this model declaration but before ATOM lines
                    model_start = content.find(f'DOCKED: MODEL        {model_num}')
                    if model_start > 0:
                        # Find the next DOCKED: ATOM line to get the end of remarks
                        atom_start = content.find('DOCKED: ATOM', model_start)
                        if atom_start > 0:
                            # Extract the section between MODEL and ATOM lines
                            remark_section = content[model_start:atom_start]
                            remark_lines = re.findall(r'REMARK.*?\n', remark_section, re.DOTALL)
                            remark_data = ''.join(remark_lines)
                        else:
                            remark_data = ''
                    else:
                        remark_data = ''
                
                # Extract pose-specific information
                pose_info = self._extract_remark_info(remark_data)
                
                # Fallback to global info if pose_info is missing
                if not pose_info.get('binding_energy') and global_info.get('binding_energy') is not None:
                    pose_info['binding_energy'] = global_info['binding_energy']
                if not pose_info.get('inhibition_constant') and global_info.get('inhibition_constant') is not None:
                    pose_info['inhibition_constant'] = global_info['inhibition_constant']
                if not pose_info.get('cluster') and global_info.get('cluster') is not None:
                    pose_info['cluster'] = global_info['cluster']
                if not pose_info.get('cluster_size') and global_info.get('cluster_size') is not None:
                    pose_info['cluster_size'] = global_info['cluster_size']
                if not pose_info.get('rank') and global_info.get('rank') is not None:
                    pose_info['rank'] = global_info['rank']
                if not pose_info.get('rmsd_from_top') and global_info.get('rmsd_from_top') is not None:
                    pose_info['rmsd_from_top'] = global_info['rmsd_from_top']
                if not pose_info.get('rotatable_bonds') and global_info.get('rotatable_bonds') is not None:
                    pose_info['rotatable_bonds'] = global_info['rotatable_bonds']
                if not pose_info.get('torsional_penalty') and global_info.get('torsional_penalty') is not None:
                    pose_info['torsional_penalty'] = global_info['torsional_penalty']
                if not pose_info.get('hydrogen_bonds') and global_info.get('hydrogen_bonds') is not None:
                    pose_info['hydrogen_bonds'] = global_info['hydrogen_bonds']
                if not pose_info.get('key_residues') and global_info.get('key_residues') is not None:
                    pose_info['key_residues'] = global_info['key_residues']
                
                # Ensure we always have a pose_info dictionary, even if parsing fails
                if not pose_info or all(v is None for v in pose_info.values()):
                    pose_info = {
                        'binding_energy': None,
                        'inhibition_constant': 'N/A',
                        'cluster': None,
                        'cluster_size': None,
                        'rank': None,
                        'rmsd_from_top': None,
                        'rotatable_bonds': None,
                        'torsional_penalty': None,
                        'hydrogen_bonds': None,
                        'key_residues': 'N/A'
                    }
                
                # Parse coordinates
                atoms = []
                for line in model_content.split('\n'):
                    if line.startswith('DOCKED: ATOM'):
                        parts = line.split()
                        if len(parts) >= 9:
                            try:
                                atom_num = int(parts[2])
                                element = parts[3]
                                x, y, z = float(parts[6]), float(parts[7]), float(parts[8])
                                atoms.append((element, x, y, z))
                            except (ValueError, IndexError):
                                continue
                
                if atoms:
                    pose_coords = []
                    pose_elements = []
                    pose_atom_names = []
                    
                    for element, x, y, z in atoms:
                        pose_coords.append([x, y, z])
                        pose_elements.append(element)
                        pose_atom_names.append(f"{element}{len(pose_coords)}")
                    
                    # Create ligand bonds
                    ligand_bonds = self._create_ligand_bonds(pose_elements, pose_atom_names)
                    
                    poses.append({
                        'coords': np.array(pose_coords, dtype=float),
                        'elements': pose_elements,
                        'bonds': ligand_bonds,
                        'info': pose_info
                    })
            

            
            if not poses:
                # Try alternative parsing for different DLG formats
                self._parse_alternative_dlg(dlg_path, poses)
            
            return poses
            
        except Exception as e:
            raise Exception(f"Error parsing DLG file: {str(e)}")
    
    def _extract_remark_info(self, remark_data):
        """Extract detailed docking information from REMARK lines with robust regex"""
        info = {
            'binding_energy': None,
            'inhibition_constant': None,
            'cluster': None,
            'cluster_size': None,
            'rank': None,
            'rmsd_from_top': None,
            'rotatable_bonds': None,
            'torsional_penalty': None,
            'hydrogen_bonds': None,
            'key_residues': None
        }
        
        if remark_data:
            # Use more robust regex patterns
            energy_match = re.search(r'binding energy\s*=\s*([-\d\.]+)', remark_data, re.IGNORECASE)
            if energy_match:
                info['binding_energy'] = float(energy_match.group(1))
            
            ki_match = re.search(r'inhibition constant\s*=\s*([^\\n]+?)(?=\\n|$)', remark_data, re.IGNORECASE)
            if ki_match:
                info['inhibition_constant'] = ki_match.group(1).strip()
            else:
                # Try alternative pattern
                ki_match = re.search(r'inhibition constant\s*=\s*([\\w\\s\\.μ]+)', remark_data, re.IGNORECASE)
                if ki_match:
                    info['inhibition_constant'] = ki_match.group(1).strip()
            
            cluster_match = re.search(r'cluster\s*=\s*(\d+)', remark_data, re.IGNORECASE)
            if cluster_match:
                info['cluster'] = int(cluster_match.group(1))
            
            cluster_size_match = re.search(r'cluster size\s*=\s*(\d+)', remark_data, re.IGNORECASE)
            if cluster_size_match:
                info['cluster_size'] = int(cluster_size_match.group(1))
            
            rank_match = re.search(r'rank\s*=\s*(\d+)', remark_data, re.IGNORECASE)
            if rank_match:
                info['rank'] = int(rank_match.group(1))
            
            rmsd_match = re.search(r'rmsd from top\s*=\s*([\d\.]+)', remark_data, re.IGNORECASE)
            if rmsd_match:
                info['rmsd_from_top'] = float(rmsd_match.group(1))
            
            rotatable_match = re.search(r'rotatable bonds\s*=\s*(\d+)', remark_data, re.IGNORECASE)
            if rotatable_match:
                info['rotatable_bonds'] = int(rotatable_match.group(1))
            
            torsional_match = re.search(r'torsional penalty\s*=\s*([-\d\.]+)', remark_data, re.IGNORECASE)
            if torsional_match:
                info['torsional_penalty'] = float(torsional_match.group(1))
            
            hbonds_match = re.search(r'hydrogen bonds\s*=\s*(\d+)', remark_data, re.IGNORECASE)
            if hbonds_match:
                info['hydrogen_bonds'] = int(hbonds_match.group(1))
            
            residues_match = re.search(r'key residues\s*=\s*([\w,]+)', remark_data, re.IGNORECASE)
            if residues_match:
                info['key_residues'] = residues_match.group(1).strip()
        
        return info
    
    def _extract_global_docking_info(self, content):
        """Extract global docking information from the entire DLG file content"""
        info = {
            'binding_energy': None,
            'inhibition_constant': None,
            'cluster': None,
            'cluster_size': None,
            'rank': None,
            'rmsd_from_top': None,
            'rotatable_bonds': None,
            'torsional_penalty': None,
            'hydrogen_bonds': None,
            'key_residues': None
        }
        
        # Extract energy - try multiple formats
        energy_match = re.search(r'Estimated Free Energy.*?=\s*([-\d\.]+)', content, re.IGNORECASE)
        if not energy_match:
            energy_match = re.search(r'Free Energy of Binding\s*=\s*([-\d\.]+)', content, re.IGNORECASE)
        if not energy_match:
            energy_match = re.search(r'binding energy\s*=\s*([-\d\.]+)', content, re.IGNORECASE)
        
        if energy_match:
            info['binding_energy'] = float(energy_match.group(1))
        
        # Extract Ki - try multiple formats
        ki_match = re.search(r'Inhibition Constant.*?=\s*([\deE\.\-]+(?:\s*[munp]?M)?)', content, re.IGNORECASE)
        if not ki_match:
            ki_match = re.search(r'Estimated Inhibition Constant.*?=\s*([\deE\.\-]+(?:\s*[munp]?M)?)', content, re.IGNORECASE)
        if not ki_match:
            ki_match = re.search(r'inhibition constant\s*=\s*([^\\n]+?)(?=\\n|$)', content, re.IGNORECASE)
        
        if ki_match:
            ki_value = ki_match.group(1).strip()
            # Clean up any extra text
            if 'REMARK' in ki_value:
                ki_value = ki_value.split('REMARK')[0].strip()
            info['inhibition_constant'] = ki_value
        
        # Search for other fields in AutoDock 4 format
        cluster_match = re.search(r'cluster\s*=\s*(\d+)', content, re.IGNORECASE)
        if cluster_match:
            info['cluster'] = int(cluster_match.group(1))
        
        cluster_size_match = re.search(r'cluster size\s*=\s*(\d+)', content, re.IGNORECASE)
        if cluster_size_match:
            info['cluster_size'] = int(cluster_size_match.group(1))
        
        rank_match = re.search(r'rank\s*=\s*(\d+)', content, re.IGNORECASE)
        if rank_match:
            info['rank'] = int(rank_match.group(1))
        
        rmsd_match = re.search(r'rmsd from top\s*=\s*([\d\.]+)', content, re.IGNORECASE)
        if rmsd_match:
            info['rmsd_from_top'] = float(rmsd_match.group(1))
        
        rotatable_match = re.search(r'rotatable bonds\s*=\s*(\d+)', content, re.IGNORECASE)
        if rotatable_match:
            info['rotatable_bonds'] = int(rotatable_match.group(1))
        
        torsional_match = re.search(r'torsional penalty\s*=\s*([-\d\.]+)', content, re.IGNORECASE)
        if torsional_match:
            info['torsional_penalty'] = float(torsional_match.group(1))
        
        hbonds_match = re.search(r'hydrogen bonds\s*=\s*(\d+)', content, re.IGNORECASE)
        if hbonds_match:
            info['hydrogen_bonds'] = int(hbonds_match.group(1))
        
        residues_match = re.search(r'key residues\s*=\s*([\w,]+)', content, re.IGNORECASE)
        if residues_match:
            info['key_residues'] = residues_match.group(1).strip()
        
        return info

    def _create_ligand_bonds(self, elements, atom_names):
        """Create bonds for ligand based on typical molecular connectivity"""
        bonds = []
        n_atoms = len(elements)
        
        # Create bonds based on typical molecular patterns
        for i in range(n_atoms - 1):
            # Connect adjacent atoms (simplified bonding)
            bonds.append((i, i + 1))
            
            # Add cross-bonds for ring structures or branched molecules
            if i < n_atoms - 2 and elements[i] == 'C' and elements[i + 2] == 'C':
                # Potential ring closure
                bonds.append((i, i + 2))
        
        return bonds
    
    def _parse_alternative_dlg(self, dlg_path, poses):
        """Alternative DLG parsing for different formats"""
        try:
            with open(dlg_path, 'r') as f:
                lines = f.readlines()
            
            current_coords = []
            current_elements = []
            current_names = []
            current_info = {}
            in_docked_section = False
            
            for line in lines:
                if line.startswith('DOCKED: MODEL'):
                    if current_coords:
                        ligand_bonds = self._create_ligand_bonds(current_elements, current_names)
                        poses.append({
                            'coords': np.array(current_coords, dtype=float),
                            'elements': current_elements,
                            'bonds': ligand_bonds,
                            'info': current_info
                        })
                        current_coords = []
                        current_elements = []
                        current_names = []
                        current_info = {}
                    in_docked_section = True
                elif line.startswith('REMARK') and in_docked_section:
                    # Parse REMARK line for info
                    if '=' in line:
                        parts = line.split('=', 1)
                        if len(parts) == 2:
                            key_part = parts[0].replace('REMARK', '').strip()
                            value = parts[1].strip()
                            
                            if key_part == 'binding energy':
                                current_info['binding_energy'] = float(value.split()[0])
                            elif key_part == 'inhibition constant':
                                current_info['inhibition_constant'] = value
                            elif key_part == 'cluster':
                                current_info['cluster'] = int(value)
                            elif key_part == 'cluster size':
                                current_info['cluster_size'] = int(value)
                            elif key_part == 'rank':
                                current_info['rank'] = int(value)
                            elif key_part == 'rmsd from top':
                                current_info['rmsd_from_top'] = float(value)
                            elif key_part == 'rotatable bonds':
                                current_info['rotatable_bonds'] = int(value)
                            elif key_part == 'torsional penalty':
                                current_info['torsional_penalty'] = float(value.split()[0])
                            elif key_part == 'hydrogen bonds':
                                current_info['hydrogen_bonds'] = int(value)
                            elif key_part == 'key residues':
                                current_info['key_residues'] = value
                elif line.startswith('DOCKED: ATOM') and in_docked_section:
                    parts = line.split()
                    if len(parts) >= 9:
                        try:
                            x, y, z = float(parts[6]), float(parts[7]), float(parts[8])
                            element = parts[3]
                            current_coords.append([x, y, z])
                            current_elements.append(element)
                            current_names.append(f"{element}{len(current_coords)}")
                        except ValueError:
                            continue
            
            if current_coords:
                ligand_bonds = self._create_ligand_bonds(current_elements, current_names)
                poses.append({
                    'coords': np.array(current_coords, dtype=float),
                    'elements': current_elements,
                    'bonds': ligand_bonds,
                    'info': current_info
                })
                
        except Exception as e:
            print(f"Alternative DLG parsing failed: {str(e)}")
    
    def get_ligand_colors(self, pose):
        """Get colors for ligand atoms - use consistent pink for all atoms"""
        colors = []
        for element in pose['elements']:
            # Use consistent pink color for all ligand atoms
            color = '#FF69B4'  # Hot pink for all ligand atoms
            colors.append(color)
        return colors
    
    def _calculate_center_of_mass(self, coords):
        """Calculate center of mass of coordinates"""
        return np.mean(coords, axis=0)
    
    def _calculate_rmsd(self, coords1, coords2):
        """Calculate RMSD between two sets of coordinates"""
        if len(coords1) != len(coords2):
            return float('inf')
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    def _align_poses(self, pose1, pose2):
        """Align two poses using Kabsch algorithm (simplified)"""
        # Center both poses
        center1 = self._calculate_center_of_mass(pose1)
        center2 = self._calculate_center_of_mass(pose2)
        
        centered1 = pose1 - center1
        centered2 = pose2 - center2
        
        # Calculate rotation matrix (simplified Kabsch)
        H = centered1.T @ centered2
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        
        # Ensure proper rotation matrix
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        return R, center1, center2
    
    def _generate_realistic_docking_path(self, protein_atoms, ligand_poses):
        """Generate realistic docking path with proper approach vector and rotation"""
        if len(ligand_poses) < 2:
            return ligand_poses
        
        # Get the final docked pose
        final_pose = ligand_poses[-1]['coords']
        
        # Calculate protein center (approximate binding site)
        protein_center = self._calculate_center_of_mass(protein_atoms)
        
        # Calculate approach vector from protein center to ligand center
        ligand_center = self._calculate_center_of_mass(final_pose)
        approach_vector = ligand_center - protein_center
        approach_distance = np.linalg.norm(approach_vector)
        
        # Normalize approach vector
        if approach_distance > 0:
            approach_vector = approach_vector / approach_distance
        
        # Generate start position (10 Å away from final position along approach vector)
        start_distance = 10.0  # Å
        start_center = ligand_center + start_distance * approach_vector
        
        # Create start pose by translating the final pose
        start_pose = final_pose + (start_center - ligand_center)
        
        # Create intermediate poses with rotation and translation
        num_intermediates = len(ligand_poses) - 1
        enhanced_poses = []
        
        # Pick a single source of truth for overlays (final docked pose is safest)
        overlay_info = ligand_poses[-1].get('info', ligand_poses[0].get('info', {}))
        
        # Add start pose
        enhanced_poses.append({
            'coords': start_pose,
            'elements': ligand_poses[0]['elements'],
            'bonds': ligand_poses[0]['bonds'],
            'info': overlay_info  # ✅ keep overlays
        })
        
        # Generate intermediate poses with realistic docking motion
        for i in range(num_intermediates):
            t = (i + 1) / (num_intermediates + 1)
            
            # Smooth interpolation function (ease-in-out)
            smooth_t = 3 * t**2 - 2 * t**3
            
            # Interpolate center position
            current_center = start_center * (1 - smooth_t) + ligand_center * smooth_t
            
            # Add some rotation during approach (realistic docking motion)
            rotation_angle = smooth_t * np.pi / 4  # 45 degrees total rotation
            cos_rot = np.cos(rotation_angle)
            sin_rot = np.sin(rotation_angle)
            
            # Simple rotation matrix around approach vector
            rotation_matrix = np.array([
                [cos_rot + approach_vector[0]**2 * (1 - cos_rot), 
                 approach_vector[0] * approach_vector[1] * (1 - cos_rot) - approach_vector[2] * sin_rot,
                 approach_vector[0] * approach_vector[2] * (1 - cos_rot) + approach_vector[1] * sin_rot],
                [approach_vector[1] * approach_vector[0] * (1 - cos_rot) + approach_vector[2] * sin_rot,
                 cos_rot + approach_vector[1]**2 * (1 - cos_rot),
                 approach_vector[1] * approach_vector[2] * (1 - cos_rot) - approach_vector[0] * sin_rot],
                [approach_vector[2] * approach_vector[0] * (1 - cos_rot) - approach_vector[1] * sin_rot,
                 approach_vector[2] * approach_vector[1] * (1 - cos_rot) + approach_vector[0] * sin_rot,
                 cos_rot + approach_vector[2]**2 * (1 - cos_rot)]
            ])
            
            # Apply rotation and translation
            centered_coords = final_pose - ligand_center
            rotated_coords = centered_coords @ rotation_matrix.T
            current_coords = rotated_coords + current_center
            
            enhanced_poses.append({
                'coords': current_coords,
                'elements': ligand_poses[0]['elements'],
                'bonds': ligand_poses[0]['bonds'],
                'info': overlay_info  # ✅ keep overlays
            })
        
        # Final pose already has info, but make sure it exists
        final_with_info = dict(ligand_poses[-1])
        final_with_info.setdefault('info', overlay_info)
        enhanced_poses.append(final_with_info)
        
        return enhanced_poses
    
    def _validate_poses(self, protein_atoms, ligand_poses):
        """Validate that poses are reasonable and properly aligned"""
        if len(ligand_poses) < 2:
            return ligand_poses
        
        # Check final pose alignment
        final_pose = ligand_poses[-1]['coords']
        protein_center = self._calculate_center_of_mass(protein_atoms)
        ligand_center = self._calculate_center_of_mass(final_pose)
        
        # Calculate distance between protein and ligand centers
        distance = np.linalg.norm(ligand_center - protein_center)
        
        # If distance is too large (>20 Å), adjust the final pose
        if distance > 20.0:
            print(f"Warning: Final pose too far from protein ({distance:.2f} Å). Adjusting...")
            
            # Move ligand closer to protein
            direction = (protein_center - ligand_center) / np.linalg.norm(protein_center - ligand_center)
            new_center = protein_center + 5.0 * direction  # 5 Å from protein center
            translation = new_center - ligand_center
            
            # Update final pose
            ligand_poses[-1]['coords'] = final_pose + translation
        
        return ligand_poses
    
    def create_animation(self, pdb_path, dlg_path, output_path):
        """Create the docking animation and save as MP4"""
        try:
            # Parse input files
            protein_atoms, protein_bonds, atom_colors, residue_info, secondary_structure = self.parse_pdb(pdb_path)
            ligand_poses = self.parse_dlg(dlg_path)
            
            if len(ligand_poses) == 0:
                raise Exception("No ligand poses found in DLG file")
            
            # Validate and enhance poses
            ligand_poses = self._validate_poses(protein_atoms, ligand_poses)
            ligand_poses = self._generate_realistic_docking_path(protein_atoms, ligand_poses)
            
            # Create animation with detailed overlays
            self._generate_animation_frames_with_overlays(protein_atoms, protein_bonds, atom_colors, 
                                                        residue_info, secondary_structure, ligand_poses, output_path)
            
        except Exception as e:
            raise Exception(f"Error creating animation: {str(e)}")
    
    def _generate_animation_frames_with_overlays(self, protein_atoms, protein_bonds, atom_colors, 
                                               residue_info, secondary_structure, ligand_poses, output_path):
        """Generate animation frames with detailed AutoDock4 overlays"""
        try:
            # Calculate bounds for consistent view
            all_coords = [protein_atoms]
            for pose in ligand_poses:
                all_coords.append(pose['coords'])
            
            all_coords_array = np.vstack(all_coords)
            x_min, y_min, z_min = all_coords_array.min(axis=0)
            x_max, y_max, z_max = all_coords_array.max(axis=0)
            
            # Add padding
            padding = max(x_max - x_min, y_max - y_min, z_max - z_min) * 0.1
            
            # Create frames for animation
            frames = []
            num_frames = 120  # More frames for smoother animation
            
            for frame in range(num_frames):
                # Create a new figure for each frame
                fig_frame = plt.figure(figsize=(19.2, 10.8), dpi=100)
                ax_frame = fig_frame.add_subplot(111, projection='3d')
                
                # Set the same bounds and labels
                ax_frame.set_xlim(x_min - padding, x_max + padding)
                ax_frame.set_ylim(y_min - padding, y_max + padding)
                ax_frame.set_zlim(z_min - padding, z_max + padding)
                ax_frame.set_xlabel('X (Å)')
                ax_frame.set_ylabel('Y (Å)')
                ax_frame.set_zlabel('Z (Å)')
                ax_frame.set_title(f'Autopose - AutoDock4 Protein-Ligand Docking Animation ({self.representation.title()})', fontsize=16)
                
                # Plot protein based on representation
                if self.representation == 'sticks':
                    self._plot_protein_sticks(ax_frame, protein_atoms, protein_bonds, atom_colors)
                elif self.representation == 'ribbon':
                    self._plot_protein_ribbon(ax_frame, protein_atoms, residue_info, atom_colors)
                elif self.representation == 'lines':
                    self._plot_protein_lines(ax_frame, protein_atoms, protein_bonds, atom_colors)
                elif self.representation == 'spheres':
                    self._plot_protein_spheres(ax_frame, protein_atoms, atom_colors)
                elif self.representation == 'surface':
                    self._plot_protein_surface(ax_frame, protein_atoms, atom_colors)
                elif self.representation == 'cartoon':
                    self._plot_protein_cartoon(ax_frame, protein_atoms, residue_info, secondary_structure)
                
                # Get current pose using smooth interpolation
                if len(ligand_poses) == 1:
                    pose = ligand_poses[0]
                    pose_info = pose.get('info', {})
                else:
                    # Use smooth interpolation between poses
                    pose_idx = (frame / num_frames) * (len(ligand_poses) - 1)
                    if pose_idx >= len(ligand_poses) - 1:
                        pose = ligand_poses[-1]
                        pose_info = pose.get('info', {})
                    else:
                        # Interpolate between two poses with smooth transition
                        idx1 = int(pose_idx)
                        idx2 = min(idx1 + 1, len(ligand_poses) - 1)
                        t = pose_idx - idx1
                        
                        # Smooth interpolation function
                        smooth_t = 3 * t**2 - 2 * t**3
                        
                        pose1 = ligand_poses[idx1]
                        pose2 = ligand_poses[idx2]
                        
                        # Interpolate coordinates
                        coords1 = pose1['coords']
                        coords2 = pose2['coords']
                        interpolated_coords = coords1 * (1 - smooth_t) + coords2 * smooth_t
                        
                        pose = {
                            'coords': interpolated_coords,
                            'elements': pose1['elements'],
                            'bonds': pose1['bonds'],
                            'info': pose1.get('info', {})  # Keep info from first pose
                        }
                        
                        # Use info from the first pose for overlay
                        pose_info = pose1.get('info', {})
                
                # Belt-and-suspenders fallback: if pose_info is empty, find nearest pose with info
                if not pose_info:
                    # fallback to nearest pose that has info (usually the final one)
                    for p in reversed(ligand_poses):
                        if p.get('info'):
                            pose_info = p['info']
                            break
                
                # Plot ligand
                if self.representation == 'sticks':
                    self._plot_ligand_sticks(ax_frame, pose)
                elif self.representation == 'ribbon':
                    self._plot_ligand_ribbon(ax_frame, pose)
                elif self.representation == 'lines':
                    self._plot_ligand_lines(ax_frame, pose)
                elif self.representation == 'spheres':
                    self._plot_ligand_spheres(ax_frame, pose)
                elif self.representation == 'surface':
                    self._plot_ligand_surface(ax_frame, pose)
                elif self.representation == 'cartoon':
                    self._plot_ligand_cartoon(ax_frame, pose)
                
                # Add detailed overlays (always draw, even with fallback values)
                # Top overlays (binding energy, inhibition constant, cluster info)
                y_top = 0.95
                y_spacing = 0.05
                
                # Always show binding energy (or N/A)
                if pose_info.get('binding_energy') is not None:
                    ax_frame.text2D(0.05, y_top, f'ΔG = {pose_info["binding_energy"]:.1f} kcal/mol', 
                                 transform=ax_frame.transAxes, fontsize=14, color='blue', weight='bold')
                else:
                    ax_frame.text2D(0.05, y_top, f'ΔG = N/A', 
                                 transform=ax_frame.transAxes, fontsize=14, color='blue', weight='bold')
                
                # Always show inhibition constant (or N/A)
                if pose_info.get('inhibition_constant') is not None:
                    ax_frame.text2D(0.05, y_top - y_spacing, f'Ki ≈ {pose_info["inhibition_constant"]}', 
                                 transform=ax_frame.transAxes, fontsize=12, color='blue')
                else:
                    ax_frame.text2D(0.05, y_top - y_spacing, f'Ki ≈ N/A', 
                                 transform=ax_frame.transAxes, fontsize=12, color='blue')
                
                # Show cluster info if available
                if pose_info.get('cluster') is not None and pose_info.get('cluster_size') is not None:
                    ax_frame.text2D(0.05, y_top - 2*y_spacing, f'Cluster {pose_info["cluster"]} ({pose_info["cluster_size"]} poses)', 
                                 transform=ax_frame.transAxes, fontsize=12, color='green')
                else:
                    ax_frame.text2D(0.05, y_top - 2*y_spacing, f'Cluster: N/A', 
                                 transform=ax_frame.transAxes, fontsize=12, color='green')
                
                # Show rank if available
                if pose_info.get('rank') is not None:
                    ax_frame.text2D(0.05, y_top - 3*y_spacing, f'Rank: {pose_info["rank"]}', 
                                 transform=ax_frame.transAxes, fontsize=12, color='green')
                else:
                    ax_frame.text2D(0.05, y_top - 3*y_spacing, f'Rank: N/A', 
                                 transform=ax_frame.transAxes, fontsize=12, color='green')
                
                # Bottom overlays (RMSD, torsional freedom)
                y_bottom = 0.05
                
                # Show RMSD if available (positioned to avoid axis overlap)
                if pose_info.get('rmsd_from_top') is not None:
                    ax_frame.text2D(0.05, -0.12, f'RMSD from Top Pose: {pose_info["rmsd_from_top"]:.2f} Å', 
                                 transform=ax_frame.transAxes, fontsize=12, color='red')
                else:
                    ax_frame.text2D(0.05, -0.12, f'RMSD from Top Pose: N/A', 
                                 transform=ax_frame.transAxes, fontsize=12, color='red')
                
                # Show rotatable bonds if available
                if pose_info.get('rotatable_bonds') is not None:
                    ax_frame.text2D(0.05, -0.17, f'Rotatable Bonds: {pose_info["rotatable_bonds"]}', 
                                 transform=ax_frame.transAxes, fontsize=12, color='orange')
                else:
                    ax_frame.text2D(0.05, -0.17, f'Rotatable Bonds: N/A', 
                                 transform=ax_frame.transAxes, fontsize=12, color='orange')
                
                # Show torsional penalty if available
                if pose_info.get('torsional_penalty') is not None:
                    ax_frame.text2D(0.05, -0.22, f'Torsional Penalty: {pose_info["torsional_penalty"]:.1f} kcal/mol', 
                                 transform=ax_frame.transAxes, fontsize=12, color='orange')
                else:
                    ax_frame.text2D(0.05, -0.22, f'Torsional Penalty: N/A', 
                                 transform=ax_frame.transAxes, fontsize=12, color='orange')
                
                # Final frame overlays (hydrogen bonds and key residues)
                if frame == num_frames - 1:  # Final frame
                    if pose_info.get('hydrogen_bonds') is not None:
                        ax_frame.text2D(0.6, 0.95, f'H-bonds: {pose_info["hydrogen_bonds"]}', 
                                     transform=ax_frame.transAxes, fontsize=14, color='purple', weight='bold')
                    else:
                        ax_frame.text2D(0.6, 0.95, f'H-bonds: N/A', 
                                     transform=ax_frame.transAxes, fontsize=14, color='purple', weight='bold')
                    
                    if pose_info.get('key_residues') is not None:
                        ax_frame.text2D(0.6, 0.90, f'Key Residues: {pose_info["key_residues"]}', 
                                     transform=ax_frame.transAxes, fontsize=12, color='purple')
                    else:
                        ax_frame.text2D(0.6, 0.90, f'Key Residues: N/A', 
                                     transform=ax_frame.transAxes, fontsize=12, color='purple')
                
                # Add legend
                ax_frame.legend()
                
                # Save frame
                frame_path = f"/tmp/frame_{frame:03d}.png"
                fig_frame.savefig(frame_path, dpi=100, bbox_inches='tight')
                frames.append(frame_path)
                
                # Close the frame figure to free memory
                plt.close(fig_frame)
            
            # Create video from frames
            self._create_video_from_frames(frames, output_path)
            
            # Clean up frame files
            for frame_path in frames:
                if os.path.exists(frame_path):
                    os.remove(frame_path)
            
            plt.close()
            
        except Exception as e:
            raise Exception(f"Error generating animation frames with overlays: {str(e)}")
    
    def _plot_protein_sticks(self, ax, protein_atoms, protein_bonds, atom_colors):
        """Plot protein as sticks with complete bonds (PyMOL style)"""
        # Plot bonds
        for bond in protein_bonds:
            atom1_idx, atom2_idx = bond
            if atom1_idx < len(protein_atoms) and atom2_idx < len(protein_atoms):
                x1, y1, z1 = protein_atoms[atom1_idx]
                x2, y2, z2 = protein_atoms[atom2_idx]
                ax.plot([x1, x2], [y1, y2], [z1, z2], 'k-', linewidth=2, alpha=0.8)
        
        # Plot atoms
        ax.scatter(protein_atoms[:, 0], protein_atoms[:, 1], protein_atoms[:, 2], 
                  c=atom_colors, s=40, alpha=0.9, label='Protein')
    
    def _plot_protein_ribbon(self, ax, protein_atoms, residue_info, atom_colors):
        """Plot protein as ribbon (PyMOL style)"""
        for residue in residue_info:
            start_idx = residue['start_idx']
            end_idx = residue['end_idx']
            
            if end_idx >= start_idx:
                # Get CA atoms for backbone
                ca_atoms = []
                for i in range(start_idx, end_idx + 1):
                    if i < len(protein_atoms):
                        ca_atoms.append(protein_atoms[i])
                
                if len(ca_atoms) > 1:
                    ca_array = np.array(ca_atoms)
                    ax.plot(ca_array[:, 0], ca_array[:, 1], ca_array[:, 2], 
                           'b-', linewidth=4, alpha=0.8, label='Protein Backbone')
    
    def _plot_protein_lines(self, ax, protein_atoms, protein_bonds, atom_colors):
        """Plot protein as lines (PyMOL style)"""
        # Plot bonds as thin lines
        for bond in protein_bonds:
            atom1_idx, atom2_idx = bond
            if atom1_idx < len(protein_atoms) and atom2_idx < len(protein_atoms):
                x1, y1, z1 = protein_atoms[atom1_idx]
                x2, y2, z2 = protein_atoms[atom2_idx]
                ax.plot([x1, x2], [y1, y2], [z1, z2], 'k-', linewidth=1, alpha=0.6)
        
        # Plot atoms as small points
        ax.scatter(protein_atoms[:, 0], protein_atoms[:, 1], protein_atoms[:, 2], 
                  c=atom_colors, s=10, alpha=0.7, label='Protein')
    
    def _plot_protein_spheres(self, ax, protein_atoms, atom_colors):
        """Plot protein as spheres (PyMOL style)"""
        # Plot atoms as spheres
        ax.scatter(protein_atoms[:, 0], protein_atoms[:, 1], protein_atoms[:, 2], 
                  c=atom_colors, s=100, alpha=0.8, label='Protein')
    
    def _plot_protein_surface(self, ax, protein_atoms, atom_colors):
        """Plot protein as surface (simplified PyMOL style)"""
        # Create a simplified surface using convex hull
        try:
            import scipy
            from scipy.spatial import ConvexHull
            hull = ConvexHull(protein_atoms)
            
            # Get the vertices of the convex hull
            vertices = protein_atoms[hull.vertices]
            
            # Create triangles from the convex hull
            triangles = []
            for simplex in hull.simplices:
                triangles.append(vertices[simplex])
            
            # Create a 3D collection of triangles
            tri_collection = Poly3DCollection(triangles, alpha=0.3, facecolor='lightblue', edgecolor='none')
            ax.add_collection3d(tri_collection)
            
        except (ImportError, Exception):
            # Fallback to spheres if scipy is not available
            self._plot_protein_spheres(ax, protein_atoms, atom_colors)
    
    def _plot_protein_cartoon(self, ax, protein_atoms, residue_info, secondary_structure):
        """Plot protein as cartoon (PyMOL style)"""
        # Group residues by secondary structure
        helix_residues = []
        sheet_residues = []
        coil_residues = []
        
        for i, residue in enumerate(residue_info):
            ss = secondary_structure.get(i, 'C')
            if ss == 'H':
                helix_residues.append(residue)
            elif ss == 'E':
                sheet_residues.append(residue)
            else:
                coil_residues.append(residue)
        
        # Plot helices as cylinders
        for residue in helix_residues:
            start_idx = residue['start_idx']
            end_idx = residue['end_idx']
            if end_idx >= start_idx:
                ca_atoms = []
                for i in range(start_idx, end_idx + 1):
                    if i < len(protein_atoms):
                        ca_atoms.append(protein_atoms[i])
                
                if len(ca_atoms) > 1:
                    ca_array = np.array(ca_atoms)
                    ax.plot(ca_array[:, 0], ca_array[:, 1], ca_array[:, 2], 
                           'r-', linewidth=6, alpha=0.8, label='Helix')
        
        # Plot sheets as arrows
        for residue in sheet_residues:
            start_idx = residue['start_idx']
            end_idx = residue['end_idx']
            if end_idx >= start_idx:
                ca_atoms = []
                for i in range(start_idx, end_idx + 1):
                    if i < len(protein_atoms):
                        ca_atoms.append(protein_atoms[i])
                
                if len(ca_atoms) > 1:
                    ca_array = np.array(ca_atoms)
                    ax.plot(ca_array[:, 0], ca_array[:, 1], ca_array[:, 2], 
                           'y-', linewidth=6, alpha=0.8, label='Sheet')
        
        # Plot coils as thin lines
        for residue in coil_residues:
            start_idx = residue['start_idx']
            end_idx = residue['end_idx']
            if end_idx >= start_idx:
                ca_atoms = []
                for i in range(start_idx, end_idx + 1):
                    if i < len(protein_atoms):
                        ca_atoms.append(protein_atoms[i])
                
                if len(ca_atoms) > 1:
                    ca_array = np.array(ca_atoms)
                    ax.plot(ca_array[:, 0], ca_array[:, 1], ca_array[:, 2], 
                           'g-', linewidth=2, alpha=0.6, label='Coil')
    
    def _plot_ligand_sticks(self, ax, pose, label='Ligand'):
        """Plot ligand as sticks with consistent pink color"""
        ligand_colors = self.get_ligand_colors(pose)
        
        # Plot bonds in pink to match atoms
        for bond in pose['bonds']:
            atom1_idx, atom2_idx = bond
            if atom1_idx < len(pose['coords']) and atom2_idx < len(pose['coords']):
                x1, y1, z1 = pose['coords'][atom1_idx]
                x2, y2, z2 = pose['coords'][atom2_idx]
                ax.plot([x1, x2], [y1, y2], [z1, z2], color='#FF69B4', linewidth=3, alpha=0.9)
        
        # Plot atoms in consistent pink
        ax.scatter(pose['coords'][:, 0], pose['coords'][:, 1], pose['coords'][:, 2],
                  c=ligand_colors, s=60, alpha=0.9, label=label)
    
    def _plot_ligand_ribbon(self, ax, pose):
        """Plot ligand as ribbon"""
        if len(pose['coords']) > 1:
            ax.plot(pose['coords'][:, 0], pose['coords'][:, 1], pose['coords'][:, 2],
                   'r-', linewidth=6, alpha=0.9, label='Ligand')
    
    def _plot_ligand_lines(self, ax, pose):
        """Plot ligand as lines"""
        ligand_colors = self.get_ligand_colors(pose)
        
        # Plot bonds as thin lines
        for bond in pose['bonds']:
            atom1_idx, atom2_idx = bond
            if atom1_idx < len(pose['coords']) and atom2_idx < len(pose['coords']):
                x1, y1, z1 = pose['coords'][atom1_idx]
                x2, y2, z2 = pose['coords'][atom2_idx]
                ax.plot([x1, x2], [y1, y2], [z1, z2], 'r-', linewidth=1, alpha=0.7)
        
        # Plot atoms as small points
        ax.scatter(pose['coords'][:, 0], pose['coords'][:, 1], pose['coords'][:, 2],
                  c=ligand_colors, s=15, alpha=0.8, label='Ligand')
    
    def _plot_ligand_spheres(self, ax, pose):
        """Plot ligand as spheres"""
        ligand_colors = self.get_ligand_colors(pose)
        ax.scatter(pose['coords'][:, 0], pose['coords'][:, 1], pose['coords'][:, 2],
                  c=ligand_colors, s=150, alpha=0.8, label='Ligand')
    
    def _plot_ligand_surface(self, ax, pose):
        """Plot ligand as surface"""
        try:
            import scipy
            from scipy.spatial import ConvexHull
            hull = ConvexHull(pose['coords'])
            
            vertices = pose['coords'][hull.vertices]
            triangles = []
            for simplex in hull.simplices:
                triangles.append(vertices[simplex])
            
            tri_collection = Poly3DCollection(triangles, alpha=0.4, facecolor='red', edgecolor='none')
            ax.add_collection3d(tri_collection)
            
        except (ImportError, Exception):
            self._plot_ligand_spheres(ax, pose)
    
    def _plot_ligand_cartoon(self, ax, pose):
        """Plot ligand as cartoon"""
        if len(pose['coords']) > 1:
            ax.plot(pose['coords'][:, 0], pose['coords'][:, 1], pose['coords'][:, 2],
                   'r-', linewidth=8, alpha=0.9, label='Ligand')
    
    def _create_video_from_frames(self, frame_paths, output_path):
        """Create MP4 video from frame images"""
        try:
            # Use ffmpeg to create video
            frame_pattern = "/tmp/frame_%03d.png"
            cmd = [
                'ffmpeg', '-y',  # Overwrite output file
                '-framerate', '12',  # 12 FPS for smoother animation
                '-i', frame_pattern,
                '-c:v', 'libx264',
                '-pix_fmt', 'yuv420p',
                '-vf', 'scale=1920:1080',  # 1080p resolution
                output_path
            ]
            
            subprocess.run(cmd, check=True, capture_output=True)
            
        except subprocess.CalledProcessError as e:
            raise Exception(f"FFmpeg error: {e.stderr.decode()}")
        except Exception as e:
            raise Exception(f"Error creating video: {str(e)}")

    def create_animation_with_pose_data(self, pdb_path, ligand_pose, output_path, pose_data=None):
        """Create animation for a single ligand pose with detailed pose data overlay and improved positioning."""
        protein_atoms, protein_bonds, atom_colors, residue_info, secondary_structure = self.parse_pdb(pdb_path)
        final_coords = ligand_pose['coords']
        ligand_center = self._calculate_center_of_mass(final_coords)
        protein_center = self._calculate_center_of_mass(protein_atoms)
        
        # Improved positioning: start closer to protein and end closer to binding site
        approach_vec = ligand_center - protein_center
        if np.linalg.norm(approach_vec) == 0:
            approach_vec = np.array([1.0, 0.0, 0.0])
        approach_vec = approach_vec / np.linalg.norm(approach_vec)
        
        # Start closer (5Å instead of 10Å) and end closer to protein
        start_distance = 5.0
        start_center = ligand_center + start_distance * approach_vec
        
        # Move final pose closer to protein if it's too far
        final_distance = np.linalg.norm(ligand_center - protein_center)
        if final_distance > 8.0:  # If more than 8Å away, move closer
            final_center = protein_center + 3.0 * approach_vec  # 3Å from protein center
            translation = final_center - ligand_center
            final_coords = final_coords + translation
            ligand_center = self._calculate_center_of_mass(final_coords)
        
        # Small random rotation
        import random
        theta = random.uniform(-0.3, 0.3)
        phi = random.uniform(-0.3, 0.3)
        psi = random.uniform(-0.3, 0.3)
        Rx = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
        Ry = np.array([[np.cos(phi), 0, np.sin(phi)], [0, 1, 0], [-np.sin(phi), 0, np.cos(phi)]])
        Rz = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])
        R = Rz @ Ry @ Rx
        centered = final_coords - ligand_center
        rotated = centered @ R.T
        start_coords = rotated + start_center
        
        # Animation
        try:
            all_coords = [protein_atoms, final_coords, start_coords]
            all_coords_array = np.vstack(all_coords)
            x_min, y_min, z_min = all_coords_array.min(axis=0)
            x_max, y_max, z_max = all_coords_array.max(axis=0)
            padding = max(x_max - x_min, y_max - y_min, z_max - z_min) * 0.1
            frames = []
            num_frames = 120
            for frame in range(num_frames):
                fig_frame = plt.figure(figsize=(19.2, 10.8), dpi=100)
                ax_frame = fig_frame.add_subplot(111, projection='3d')
                ax_frame.set_xlim(x_min - padding, x_max + padding)
                ax_frame.set_ylim(y_min - padding, y_max + padding)
                ax_frame.set_zlim(z_min - padding, z_max + padding)
                ax_frame.set_xlabel('X (Å)')
                ax_frame.set_ylabel('Y (Å)')
                ax_frame.set_zlabel('Z (Å)')
                ax_frame.set_title('Autopose - Protein-Ligand Docking Animation (Vina)', fontsize=16)
                self._plot_protein_sticks(ax_frame, protein_atoms, protein_bonds, atom_colors)
                # Interpolate ligand
                t = frame / (num_frames - 1)
                smooth_t = 3 * t**2 - 2 * t**3
                coords = start_coords * (1 - smooth_t) + final_coords * smooth_t
                pose_interp = {'coords': coords, 'elements': ligand_pose['elements'], 'bonds': ligand_pose['bonds']}
                self._plot_ligand_sticks(ax_frame, pose_interp, label='Ligand')
                
                # Overlay detailed pose data
                if pose_data is not None:
                    y_offset = 0.92
                    ax_frame.text2D(0.05, y_offset, f'Affinity: {pose_data["affinity"]:.2f} kcal/mol', 
                                   transform=ax_frame.transAxes, fontsize=14, color='red')
                    ax_frame.text2D(0.05, y_offset - 0.05, f'Pose Rank: {pose_data["rank"]}', 
                                   transform=ax_frame.transAxes, fontsize=14, color='red')
                    ax_frame.text2D(0.05, y_offset - 0.10, f'RMSD to Top: {pose_data["rmsd_lower"]:.2f} Å', 
                                   transform=ax_frame.transAxes, fontsize=14, color='red')
                
                ax_frame.legend()
                frame_path = f"/tmp/frame_{frame:03d}.png"
                fig_frame.savefig(frame_path, dpi=100, bbox_inches='tight')
                frames.append(frame_path)
                plt.close(fig_frame)
            self._create_video_from_frames(frames, output_path)
            for frame_path in frames:
                if os.path.exists(frame_path):
                    os.remove(frame_path)
            plt.close()
        except Exception as e:
            raise Exception(f"Error generating animation frames with pose data: {str(e)}")

    def _plot_ligand_sticks(self, ax, pose, label='Ligand'):
        ligand_colors = self.get_ligand_colors(pose)
        # Plot bonds in pink to match atoms
        for bond in pose['bonds']:
            atom1_idx, atom2_idx = bond
            if atom1_idx < len(pose['coords']) and atom2_idx < len(pose['coords']):
                x1, y1, z1 = pose['coords'][atom1_idx]
                x2, y2, z2 = pose['coords'][atom2_idx]
                ax.plot([x1, x2], [y1, y2], [z1, z2], color='#FF69B4', linewidth=3, alpha=0.9)
        # Plot atoms
        ax.scatter(pose['coords'][:, 0], pose['coords'][:, 1], pose['coords'][:, 2],
                  c=ligand_colors, s=60, alpha=0.9, label=label)

    def create_overlay_animation_with_pose_data(self, pdb_path, ligand_poses, output_path, pose_data_list=None):
        """Create overlay animation for multiple ligand poses with detailed pose data and improved positioning."""
        protein_atoms, protein_bonds, atom_colors, residue_info, secondary_structure = self.parse_pdb(pdb_path)
        try:
            # Improved positioning: ensure ligands start within frame and end closer to protein
            pose_finals = []
            pose_starts = []
            pose_centers = []
            protein_center = self._calculate_center_of_mass(protein_atoms)
            
            for i, pose in enumerate(ligand_poses):
                coords = pose['coords']
                center = self._calculate_center_of_mass(coords)
                
                # Move final pose closer to protein if it's too far
                final_distance = np.linalg.norm(center - protein_center)
                if final_distance > 8.0:  # If more than 8Å away, move closer
                    approach_vec = center - protein_center
                    if np.linalg.norm(approach_vec) == 0:
                        approach_vec = np.array([1.0, 0.0, 0.0])
                    approach_vec = approach_vec / np.linalg.norm(approach_vec)
                    final_center = protein_center + 3.0 * approach_vec  # 3Å from protein center
                    translation = final_center - center
                    coords = coords + translation
                    center = self._calculate_center_of_mass(coords)
                
                pose_finals.append(coords)
                pose_centers.append(center)
                
                # Generate start position closer to protein (5Å instead of 10Å)
                approach_vec = center - protein_center
                if np.linalg.norm(approach_vec) == 0:
                    approach_vec = np.array([1.0, 0.0, 0.0])
                approach_vec = approach_vec / np.linalg.norm(approach_vec)
                start_center = center + 5.0 * approach_vec
                
                # Apply small random rotation
                theta = random.uniform(-0.3, 0.3)
                phi = random.uniform(-0.3, 0.3)
                psi = random.uniform(-0.3, 0.3)
                Rx = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
                Ry = np.array([[np.cos(phi), 0, np.sin(phi)], [0, 1, 0], [-np.sin(phi), 0, np.cos(phi)]])
                Rz = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])
                R = Rz @ Ry @ Rx
                centered = coords - center
                rotated = centered @ R.T
                start_coords = rotated + start_center
                pose_starts.append(start_coords)
            
            # Calculate bounds including all poses
            all_coords = [protein_atoms] + pose_finals + pose_starts
            all_coords_array = np.vstack(all_coords)
            x_min, y_min, z_min = all_coords_array.min(axis=0)
            x_max, y_max, z_max = all_coords_array.max(axis=0)
            padding = max(x_max - x_min, y_max - y_min, z_max - z_min) * 0.1
            
            frames = []
            num_frames = 120
            base_colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta', 'yellow', 'brown', 'pink']
            
            # Animate all poses in parallel
            for frame in range(num_frames):
                fig_frame = plt.figure(figsize=(19.2, 10.8), dpi=100)
                ax_frame = fig_frame.add_subplot(111, projection='3d')
                ax_frame.set_xlim(x_min - padding, x_max + padding)
                ax_frame.set_ylim(y_min - padding, y_max + padding)
                ax_frame.set_zlim(z_min - padding, z_max + padding)
                ax_frame.set_xlabel('X (Å)')
                ax_frame.set_ylabel('Y (Å)')
                ax_frame.set_zlabel('Z (Å)')
                ax_frame.set_title('Autopose - Overlay of Multiple Vina Poses', fontsize=16)
                self._plot_protein_sticks(ax_frame, protein_atoms, protein_bonds, atom_colors)
                
                for i, pose in enumerate(ligand_poses):
                    color = base_colors[i % len(base_colors)]
                    alpha = 0.5 if len(ligand_poses) > 1 else 0.9
                    
                    # Interpolate this pose
                    t = frame / (num_frames - 1)
                    smooth_t = 3 * t**2 - 2 * t**3
                    coords = pose_starts[i] * (1 - smooth_t) + pose_finals[i] * smooth_t
                    pose_interp = {'coords': coords, 'elements': pose['elements'], 'bonds': pose['bonds']}
                    self._plot_ligand_sticks_overlay(ax_frame, pose_interp, color=color, alpha=alpha, label=f'Pose {i+1}')
                    
                    # Overlay detailed pose data
                    if pose_data_list and i < len(pose_data_list):
                        pose_data = pose_data_list[i]
                        y_offset = 0.92 - i * 0.08
                        ax_frame.text2D(0.05, y_offset, f'Pose {i+1} - Affinity: {pose_data["affinity"]:.2f} kcal/mol', 
                                       transform=ax_frame.transAxes, fontsize=12, color=color)
                        ax_frame.text2D(0.05, y_offset - 0.04, f'Rank: {pose_data["rank"]} | RMSD: {pose_data["rmsd_lower"]:.2f} Å', 
                                       transform=ax_frame.transAxes, fontsize=12, color=color)
                
                ax_frame.legend()
                frame_path = f"/tmp/frame_{frame:03d}.png"
                fig_frame.savefig(frame_path, dpi=100, bbox_inches='tight')
                frames.append(frame_path)
                plt.close(fig_frame)
            
            self._create_video_from_frames(frames, output_path)
            for frame_path in frames:
                if os.path.exists(frame_path):
                    os.remove(frame_path)
            plt.close()
        except Exception as e:
            raise Exception(f"Error generating overlay animation: {str(e)}")

    def _add_autodock4_overlays(self, ax, pose_info, frame, num_frames):
        """Add detailed AutoDock4 overlays to the animation (legacy method - now implemented directly in animation loop)"""
        # This method is kept for backward compatibility but overlays are now implemented directly in _generate_animation_frames_with_overlays
        pass

    def _plot_ligand_sticks_overlay(self, ax, pose, color='red', alpha=0.5, label='Ligand (overlay)'):
        # Plot bonds
        for bond in pose['bonds']:
            atom1_idx, atom2_idx = bond
            if atom1_idx < len(pose['coords']) and atom2_idx < len(pose['coords']):
                x1, y1, z1 = pose['coords'][atom1_idx]
                x2, y2, z2 = pose['coords'][atom2_idx]
                ax.plot([x1, x2], [y1, y2], [z1, z2], color=color, linewidth=3, alpha=alpha)
        # Plot atoms
        ax.scatter(pose['coords'][:, 0], pose['coords'][:, 1], pose['coords'][:, 2],
                  c=color, s=60, alpha=alpha, label=label)

if __name__ == "__main__":
    # Test the animator
    animator = DockingAnimatorEnhanced()
    print("DockingAnimatorEnhanced initialized successfully!") 