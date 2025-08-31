import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
import subprocess
import tempfile

matplotlib.use('Agg')  # Use non-interactive backend

class ClusProAnimator:
    """Specialized animator for ClusPro protein-protein docking results"""
    
    def __init__(self, representation='cartoon'):
        self.parser = PDBParser(QUIET=True)
        self.representation = representation
        
    def parse_monomer_pdb(self, pdb_path):
        """Parse monomer PDB file for reference structure"""
        try:
            structure = self.parser.get_structure('monomer', pdb_path)
            
            # Extract atoms and basic structure
            protein_atoms = []
            protein_elements = []
            residue_info = []
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if is_aa(residue):
                            residue_atoms = []
                            residue_elements = []
                            
                            for atom in residue:
                                coords = atom.get_coord()
                                residue_atoms.append(coords)
                                residue_elements.append(atom.element)
                            
                            if residue_atoms:
                                start_idx = len(protein_atoms)
                                protein_atoms.extend(residue_atoms)
                                protein_elements.extend(residue_elements)
                                
                                residue_info.append({
                                    'start_idx': start_idx,
                                    'end_idx': len(protein_atoms) - 1,
                                    'residue': residue
                                })
            
            return np.array(protein_atoms, dtype=float), protein_elements, residue_info
            
        except Exception as e:
            raise Exception(f"Error parsing monomer PDB: {str(e)}")
    
    def create_homodimer_animation(self, animation_data, output_path, representation='cartoon'):
        """Create homodimer animation from ClusPro data"""
        try:
            monomer_path = animation_data['monomer_path']
            models = animation_data['models']
            
            # Parse monomer structure
            monomer_atoms, monomer_elements, monomer_residues = self.parse_monomer_pdb(monomer_path)
            
            # Create animation frames
            frames = []
            num_frames = 120
            
            # Calculate bounds for consistent view
            all_coords = [monomer_atoms]
            for model in models:
                ref_coords = model['complex_data']['ref_chain']['coords']
                mov_coords = model['complex_data']['moving_chain']['coords']
                all_coords.extend([ref_coords, mov_coords])
            
            all_coords_array = np.vstack(all_coords)
            x_min, y_min, z_min = all_coords_array.min(axis=0)
            x_max, y_max, z_max = all_coords_array.max(axis=0)
            
            # Add padding
            padding = max(x_max - x_min, y_max - y_min, z_max - z_min) * 0.1
            
            # Create frames for each model
            for model_idx, model in enumerate(models):
                score_data = model['score_data']
                complex_data = model['complex_data']
                interface_data = model['interface_data']
                trajectory = model['trajectory']
                
                # Create frames for this model's trajectory
                for frame_idx, trajectory_frame in enumerate(trajectory):
                    # Create a new figure for each frame
                    fig_frame = plt.figure(figsize=(19.2, 10.8), dpi=100)
                    ax_frame = fig_frame.add_subplot(111, projection='3d')
                    
                    # Set consistent bounds
                    ax_frame.set_xlim(x_min - padding, x_max + padding)
                    ax_frame.set_ylim(y_min - padding, y_max + padding)
                    ax_frame.set_zlim(z_min - padding, z_max + padding)
                    ax_frame.set_xlabel('X (Å)')
                    ax_frame.set_ylabel('Y (Å)')
                    ax_frame.set_zlabel('Z (Å)')
                    ax_frame.set_title(f'Autopose - ClusPro Homodimer Docking - Model {model_idx + 1} ({representation.title()})', fontsize=16)
                    
                    # Plot reference chain (monomer)
                    self._plot_protein_representation(ax_frame, monomer_atoms, monomer_elements, 
                                                   monomer_residues, 'blue', 'Reference Chain', representation)
                    
                    # Plot moving chain at current trajectory position
                    current_coords = trajectory_frame['coords']
                    self._plot_protein_representation(ax_frame, current_coords, complex_data['moving_chain']['elements'],
                                                   [], 'red', 'Moving Chain', representation)
                    
                    # Add detailed overlays
                    self._add_cluspro_overlays(ax_frame, score_data, interface_data, model_idx + 1, len(models))
                    
                    # Add legend
                    ax_frame.legend()
                    
                    # Save frame
                    frame_path = f"cluspro_frame_{model_idx:02d}_{frame_idx:03d}.png"
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
            raise Exception(f"Error creating homodimer animation: {str(e)}")
    
    def _plot_protein_representation(self, ax, coords, elements, residues, color, label, representation):
        """Plot protein using specified representation"""
        if representation == 'cartoon':
            self._plot_protein_cartoon(ax, coords, residues, color, label)
        elif representation == 'sticks':
            self._plot_protein_sticks(ax, coords, color, label)
        elif representation == 'surface':
            self._plot_protein_surface(ax, coords, color, label)
        elif representation == 'ribbon':
            self._plot_protein_ribbon(ax, coords, residues, color, label)
        else:
            # Default to spheres
            self._plot_protein_spheres(ax, coords, color, label)
    
    def _plot_protein_cartoon(self, ax, coords, residues, color, label):
        """Plot protein as cartoon representation"""
        if len(coords) > 1:
            # Plot backbone as thick line
            ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 
                   color=color, linewidth=6, alpha=0.8, label=label)
            
            # Add some spheres at key points
            if len(coords) > 10:
                step = len(coords) // 10
                key_points = coords[::step]
                ax.scatter(key_points[:, 0], key_points[:, 1], key_points[:, 2],
                          c=color, s=100, alpha=0.9)
    
    def _plot_protein_sticks(self, ax, coords, color, label):
        """Plot protein as sticks representation"""
        # Plot atoms as spheres
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], 
                  c=color, s=60, alpha=0.8, label=label)
        
        # Add simple bonds between adjacent atoms
        for i in range(len(coords) - 1):
            x1, y1, z1 = coords[i]
            x2, y2, z2 = coords[i + 1]
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=color, linewidth=2, alpha=0.6)
    
    def _plot_protein_surface(self, ax, coords, color, label):
        """Plot protein as surface representation"""
        try:
            from scipy.spatial import ConvexHull
            hull = ConvexHull(coords)
            
            # Get the vertices of the convex hull
            vertices = coords[hull.vertices]
            
            # Create triangles from the convex hull
            triangles = []
            for simplex in hull.simplices:
                triangles.append(vertices[simplex])
            
            # Create a 3D collection of triangles
            tri_collection = Poly3DCollection(triangles, alpha=0.3, facecolor=color, edgecolor='none')
            ax.add_collection3d(tri_collection)
            
            # Add label
            ax.scatter(coords[0, 0], coords[0, 1], coords[0, 2], 
                      c=color, s=100, alpha=0.9, label=label)
            
        except ImportError:
            # Fallback to spheres if scipy is not available
            self._plot_protein_spheres(ax, coords, color, label)
    
    def _plot_protein_ribbon(self, ax, coords, residues, color, label):
        """Plot protein as ribbon representation"""
        if len(coords) > 1:
            # Plot backbone as ribbon
            ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 
                   color=color, linewidth=4, alpha=0.8, label=label)
    
    def _plot_protein_spheres(self, ax, coords, color, label):
        """Plot protein as spheres representation"""
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], 
                  c=color, s=80, alpha=0.8, label=label)
    
    def _add_cluspro_overlays(self, ax, score_data, interface_data, model_num, total_models):
        """Add ClusPro-specific overlays to the animation"""
        # Top overlays
        y_top = 0.95
        y_spacing = 0.05
        
        # Model information
        ax.text2D(0.05, y_top, f'Model {model_num}/{total_models}', 
                 transform=ax.transAxes, fontsize=14, color='blue', weight='bold')
        
        # Score information
        if score_data.get('score') is not None:
            ax.text2D(0.05, y_top - y_spacing, f'Total Score: {score_data["score"]:.1f} kcal/mol', 
                     transform=ax.transAxes, fontsize=12, color='blue')
        
        # Cluster information
        if score_data.get('cluster') is not None and score_data.get('cluster_size') is not None:
            ax.text2D(0.05, y_top - 2*y_spacing, f'Cluster {score_data["cluster"]} ({score_data["cluster_size"]} poses)', 
                     transform=ax.transAxes, fontsize=12, color='green')
        
        # Rank information
        if score_data.get('rank') is not None:
            ax.text2D(0.05, y_top - 3*y_spacing, f'Rank: {score_data["rank"]}', 
                     transform=ax.transAxes, fontsize=12, color='green')
        
        # Bottom overlays
        y_bottom = 0.05
        
        # Interface information
        if interface_data.get('contact_count') is not None:
            ax.text2D(0.05, -0.12, f'Interface Contacts: {interface_data["contact_count"]}', 
                     transform=ax.transAxes, fontsize=12, color='red')
        
        # Energy breakdown
        if score_data.get('vdw') is not None:
            ax.text2D(0.05, -0.17, f'vdW: {score_data["vdw"]:.1f} | Electrostatic: {score_data["electrostatic"]:.1f}', 
                     transform=ax.transAxes, fontsize=12, color='orange')
        
        if score_data.get('desolvation') is not None:
            ax.text2D(0.05, -0.22, f'Desolvation: {score_data["desolvation"]:.1f} kcal/mol', 
                     transform=ax.transAxes, fontsize=12, color='orange')
    
    def _create_video_from_frames(self, frame_paths, output_path):
        """Create MP4 video from frame images"""
        try:
            # Since we have non-sequential frame names, we need to use a different approach
            # Create a temporary directory and copy frames with sequential names
            import tempfile
            import shutil
            
            temp_dir = tempfile.mkdtemp()
            sequential_frames = []
            
            # Copy frames to temp directory with sequential names
            for i, frame_path in enumerate(frame_paths):
                if os.path.exists(frame_path):
                    new_name = f"frame_{i:04d}.png"
                    new_path = os.path.join(temp_dir, new_name)
                    shutil.copy2(frame_path, new_path)
                    sequential_frames.append(new_path)
            
            if not sequential_frames:
                raise Exception("No frame files found")
            
            # Use ffmpeg with sequential frame pattern
            frame_pattern = os.path.join(temp_dir, "frame_%04d.png")
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
            
            # Clean up temp directory
            shutil.rmtree(temp_dir)
            
        except subprocess.CalledProcessError as e:
            raise Exception(f"FFmpeg error: {e.stderr.decode()}")
        except Exception as e:
            raise Exception(f"Error creating video: {str(e)}")

if __name__ == "__main__":
    # Test the animator
    animator = ClusProAnimator()
    print("ClusProAnimator initialized successfully!")
