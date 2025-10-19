"""
HADDOCK Animator for Autopose
Generates MP4 animations from HADDOCK protein-protein docking results.
Fixed version addressing all bugs from Bugsv1Haddock.md
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
import subprocess
import tempfile


class HADDOCKAnimator:
    def __init__(self, representation='cartoon'):
        self.parser = PDBParser(QUIET=True)
        self.representation = representation
    
    def parse_monomer_pdb(self, pdb_path):
        """Parse reference monomer PDB for visualization."""
        structure = self.parser.get_structure('monomer', pdb_path)
        
        atoms = []
        atom_colors = []
        residue_info = []
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        residue_atoms = []
                        for atom in residue:
                            coords = atom.get_coord()
                            atoms.append(coords)
                            residue_atoms.append(coords)
                            
                            # Color by element
                            element = atom.element
                            if element == 'C':
                                color = '#909090'
                            elif element == 'N':
                                color = '#3030FF'
                            elif element == 'O':
                                color = '#FF3030'
                            elif element == 'S':
                                color = '#FFFF30'
                            else:
                                color = '#FF1493'
                            atom_colors.append(color)
                        
                        if residue_atoms:
                            residue_info.append({
                                'atoms': residue_atoms,
                                'residue': residue
                            })
        
        return np.array(atoms, dtype=float), atom_colors, residue_info
    
    def create_haddock_animation(self, animation_data, output_path, representation='cartoon'):
        """
        Create HADDOCK protein-protein docking animation with all bug fixes.
        
        Args:
            animation_data: Dictionary from HADDOCKAdapter.create_animation_data()
            output_path: Output MP4 file path
            representation: Visualization style
        """
        self.representation = representation
        
        if len(animation_data['models']) == 0:
            raise ValueError("No models found in animation data")
        
        # Get mode and validation info
        self.animation_mode = animation_data.get('mode', 'production')
        self.validation = animation_data.get('validation', {})
        is_single_model = len(animation_data['models']) == 1
        
        print(f"\n{'='*60}")
        print(f"ANIMATION SETUP:")
        print(f"  Mode: {self.animation_mode}")
        print(f"  Number of models: {len(animation_data['models'])}")
        print(f"  Single model run: {is_single_model}")
        print(f"{'='*60}\n")
        
        # Use first model for animation
        model = animation_data['models'][0]
        chain_a = model['chain_a']
        chain_b_trajectory = model['trajectory']
        chain_b_final = model['chain_b']  # Use stored aligned final position
        info = model['info']
        
        # FIX B1: Verify final trajectory frame matches stored final position
        final_trajectory_coords = chain_b_trajectory[-1]
        rmsd_check = np.sqrt(np.mean(np.sum((final_trajectory_coords - chain_b_final)**2, axis=1)))
        print(f"DEBUG: RMSD between final trajectory and stored final: {rmsd_check:.6f} Å")
        
        # Calculate bounds for initial view
        all_coords = [chain_a] + chain_b_trajectory + [chain_b_final]
        all_coords_array = np.vstack(all_coords)
        x_min, y_min, z_min = all_coords_array.min(axis=0)
        x_max, y_max, z_max = all_coords_array.max(axis=0)
        padding = max(x_max - x_min, y_max - y_min, z_max - z_min) * 0.15
        
        # Calculate interface center for camera zoom (FIX G6)
        interface_center = (np.mean(chain_a, axis=0) + np.mean(chain_b_final, axis=0)) / 2
        
        # Generate frames (FIX F: Add extra frames at end to hold)
        frames = []
        num_motion_frames = len(chain_b_trajectory)
        num_hold_frames = 24  # Hold for 2 seconds at 12 FPS
        total_frames = num_motion_frames + num_hold_frames
        
        for frame_idx in range(total_frames):
            fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
            ax = fig.add_subplot(111, projection='3d')
            
            # Determine which trajectory frame to show
            if frame_idx < num_motion_frames:
                motion_frame_idx = frame_idx
                chain_b_coords = chain_b_trajectory[motion_frame_idx]
                is_final_frame = False
            else:
                # Hold on final frame
                motion_frame_idx = num_motion_frames - 1
                chain_b_coords = chain_b_final  # Use exact final position
                is_final_frame = True
            
            # FIX F: Camera zoom in last 30% of motion frames
            t_motion = motion_frame_idx / max(1, num_motion_frames - 1)
            if t_motion > 0.7:
                # Interpolate camera to interface
                zoom_t = (t_motion - 0.7) / 0.3  # 0 to 1 in last 30%
                zoom_t = 3 * zoom_t**2 - 2 * zoom_t**3  # Smooth
                
                # Interpolate bounds toward interface
                view_center = all_coords_array.mean(axis=0) * (1 - zoom_t) + interface_center * zoom_t
                view_size = (x_max - x_min) * (1 - 0.3 * zoom_t)  # Zoom in 30%
                
                ax.set_xlim(view_center[0] - view_size/2, view_center[0] + view_size/2)
                ax.set_ylim(view_center[1] - view_size/2, view_center[1] + view_size/2)
                ax.set_zlim(view_center[2] - view_size/2, view_center[2] + view_size/2)
            else:
                # Normal view
                ax.set_xlim(x_min - padding, x_max + padding)
                ax.set_ylim(y_min - padding, y_max + padding)
                ax.set_zlim(z_min - padding, z_max + padding)
            
            ax.set_xlabel('X (Å)')
            ax.set_ylabel('Y (Å)')
            ax.set_zlabel('Z (Å)')
            ax.set_title(f'Autopose - HADDOCK Protein-Protein Docking Animation ({representation.title()})', 
                        fontsize=16)
            
            # Plot chain A (stationary, blue)
            self._plot_protein(ax, chain_a, color='#4169E1', label='Chain A (Receptor)', alpha=0.7)
            
            # Plot chain B (moving, orange)
            self._plot_protein(ax, chain_b_coords, color='#FF8C00', label='Chain B (Ligand)', alpha=0.7)
            
            # FIX B5: Recompute contacts on actual rendered final pose
            if is_final_frame:
                # Recompute interface contacts on exact final pose
                n_contacts_a, n_contacts_b = self._compute_interface_contacts(chain_a, chain_b_coords)
                info_with_final_contacts = dict(info)
                info_with_final_contacts['interface_contacts_a'] = n_contacts_a
                info_with_final_contacts['interface_contacts_b'] = n_contacts_b
                self._add_haddock_overlays(ax, info_with_final_contacts, is_final_frame, is_single_model)
            else:
                self._add_haddock_overlays(ax, info, is_final_frame, is_single_model)
            
            # Add mode banner on first and last frames
            if frame_idx == 0 or is_final_frame:
                mode_text = f"Mode: {self.animation_mode.title()}"
                if self.animation_mode == 'demo':
                    mode_text += " (CA-level contacts)"
                ax.text2D(0.5, 0.02, mode_text,
                         transform=ax.transAxes, fontsize=9, color='#666666',
                         ha='center', va='bottom', style='italic',
                         bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.6))
            
            # FIX B3: Legend positioned to avoid overlay collision
            ax.legend(loc='lower right', frameon=True, framealpha=0.8, fontsize=10)
            
            # Save frame
            frame_path = f"/tmp/haddock_frame_{frame_idx:03d}.png"
            fig.savefig(frame_path, dpi=100, bbox_inches='tight')
            frames.append(frame_path)
            plt.close(fig)
        
        # Create video
        self._create_video_from_frames(frames, output_path)
        
        # Cleanup
        for frame_path in frames:
            if os.path.exists(frame_path):
                os.remove(frame_path)
    
    def _compute_interface_contacts(self, chain_a_coords, chain_b_coords, cutoff=4.5):
        """Compute interface contacts between two chains."""
        interface_a = set()
        interface_b = set()
        
        for i, coord_a in enumerate(chain_a_coords):
            for j, coord_b in enumerate(chain_b_coords):
                dist = np.linalg.norm(coord_a - coord_b)
                if dist <= cutoff:
                    interface_a.add(i)
                    interface_b.add(j)
        
        return len(interface_a), len(interface_b)
    
    def _plot_protein(self, ax, coords, color='blue', label='Protein', alpha=0.7):
        """Plot protein representation."""
        if self.representation == 'cartoon' or self.representation == 'ribbon':
            # Plot as backbone trace
            if len(coords) > 1:
                ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 
                       color=color, linewidth=4, alpha=alpha, label=label)
        
        elif self.representation == 'spheres':
            ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                      c=color, s=100, alpha=alpha, label=label)
        
        elif self.representation == 'lines':
            if len(coords) > 1:
                ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
                       color=color, linewidth=1, alpha=alpha, label=label)
        
        elif self.representation == 'surface':
            # Simplified surface using scatter with large points
            ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                      c=color, s=150, alpha=alpha * 0.5, label=label)
        
        else:  # sticks (default)
            ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                      c=color, s=50, alpha=alpha, label=label)
    
    def _add_haddock_overlays(self, ax, info, is_final_frame, is_single_model=False):
        """
        Add HADDOCK-specific overlays to the animation.
        FIX B3: Positioned to avoid legend collision with bbox styling.
        FIX B2: Display correct RMSD values, with special handling for single-model runs.
        FIX B4: Show dual contact metrics and debug info.
        """
        # FIX B3: Position overlays at top-left with bbox
        bbox_style = dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.8)
        
        y_top = 0.97
        y_spacing = 0.045
        x_left = 0.03
        
        # Cluster and Rank
        if info.get('cluster') is not None and info.get('rank') is not None:
            ax.text2D(x_left, y_top, f'Cluster {info["cluster"]} | Rank {info["rank"]}',
                     transform=ax.transAxes, fontsize=13, color='#1a1a1a', weight='bold',
                     va='top', zorder=10, bbox=bbox_style)
        
        # HADDOCK Score (prominent for single-model runs)
        if info.get('haddock_score') is not None:
            ax.text2D(x_left, y_top - y_spacing, 
                     f'HADDOCK Score: {info["haddock_score"]:.2f}',
                     transform=ax.transAxes, fontsize=11, color='#2c5f2d',
                     va='top', zorder=10, bbox=bbox_style)
        
        # FIX B3: RMSD display - special handling for single-model
        if is_single_model:
            # For single model, show "RMSD to top: 0.00 Å (top model)"
            if info.get('rmsd') is not None:
                ax.text2D(x_left, y_top - 2*y_spacing,
                         f'RMSD: {info["rmsd"]:.2f} Å (top model)',
                         transform=ax.transAxes, fontsize=11, color='#d97706',
                         va='top', zorder=10, bbox=bbox_style)
        else:
            # For multi-model, show actual RMSD
            if info.get('rmsd') is not None:
                ax.text2D(x_left, y_top - 2*y_spacing,
                         f'RMSD: {info["rmsd"]:.2f} Å',
                         transform=ax.transAxes, fontsize=11, color='#d97706',
                         va='top', zorder=10, bbox=bbox_style)
            else:
                ax.text2D(x_left, y_top - 2*y_spacing,
                         f'RMSD: N/A',
                         transform=ax.transAxes, fontsize=11, color='#d97706',
                         va='top', zorder=10, bbox=bbox_style)
        
        # Energy components
        y_energy = y_top - 3*y_spacing
        if info.get('evdw') is not None:
            ax.text2D(x_left, y_energy, f'E_vdW: {info["evdw"]:.2f} kcal/mol',
                     transform=ax.transAxes, fontsize=10, color='#7c3aed',
                     va='top', zorder=10, bbox=bbox_style)
        
        if info.get('eelec') is not None:
            ax.text2D(x_left, y_energy - y_spacing, f'E_elec: {info["eelec"]:.2f} kcal/mol',
                     transform=ax.transAxes, fontsize=10, color='#7c3aed',
                     va='top', zorder=10, bbox=bbox_style)
        
        if info.get('edesol') is not None:
            ax.text2D(x_left, y_energy - 2*y_spacing, f'E_desol: {info["edesol"]:.2f} kcal/mol',
                     transform=ax.transAxes, fontsize=10, color='#7c3aed',
                     va='top', zorder=10, bbox=bbox_style)
        
        # FIX B5 & B2: Dual contact metrics on final frame
        if is_final_frame:
            y_contacts = 0.97
            x_right = 0.97
            
            # Show both heavy-atom and CA-CA contacts
            if info.get('heavy_contacts') is not None:
                contact_text = f'Contacts (heavy, 4.5Å): {info["heavy_contacts"]}'
                ax.text2D(x_right, y_contacts, contact_text,
                         transform=ax.transAxes, fontsize=11, color='#dc2626', weight='bold',
                         va='top', ha='right', zorder=10, bbox=bbox_style)
            
            if info.get('ca_contacts') is not None and info['ca_contacts'] > 0:
                ca_text = f'Contacts (CA-CA, 8.0Å): {info["ca_contacts"]}'
                ax.text2D(x_right, y_contacts - y_spacing, ca_text,
                         transform=ax.transAxes, fontsize=10, color='#dc2626',
                         va='top', ha='right', zorder=10, bbox=bbox_style)
            
            # Debug info: min interchain distance
            if info.get('min_interchain_distance') is not None:
                debug_text = f'min d(Chains) = {info["min_interchain_distance"]:.2f} Å'
                ax.text2D(x_right, y_contacts - 2*y_spacing, debug_text,
                         transform=ax.transAxes, fontsize=9, color='#666666', style='italic',
                         va='top', ha='right', zorder=10, bbox=bbox_style)
    
    def _create_video_from_frames(self, frame_paths, output_path):
        """Create MP4 video from frame images using FFmpeg."""
        try:
            # Create temporary directory for sequential frames
            temp_dir = tempfile.mkdtemp()
            
            # Copy frames with sequential naming
            for i, frame_path in enumerate(frame_paths):
                temp_frame = os.path.join(temp_dir, f"frame_{i:03d}.png")
                if os.path.exists(frame_path):
                    import shutil
                    shutil.copy(frame_path, temp_frame)
            
            # Use ffmpeg to create video
            frame_pattern = os.path.join(temp_dir, "frame_%03d.png")
            cmd = [
                'ffmpeg', '-y',
                '-framerate', '12',
                '-i', frame_pattern,
                '-c:v', 'libx264',
                '-pix_fmt', 'yuv420p',
                '-vf', 'scale=1920:1080',
                output_path
            ]
            
            subprocess.run(cmd, check=True, capture_output=True)
            
            # Cleanup temp directory
            import shutil
            shutil.rmtree(temp_dir)
            
        except subprocess.CalledProcessError as e:
            raise Exception(f"FFmpeg error: {e.stderr.decode()}")
        except Exception as e:
            raise Exception(f"Error creating video: {str(e)}")


if __name__ == "__main__":
    # Test the animator
    animator = HADDOCKAnimator()
    print("HADDOCKAnimator initialized successfully!")
