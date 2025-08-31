import os
import tempfile
from flask import Flask, request, render_template, send_file, jsonify
from werkzeug.utils import secure_filename
from docking_animator_enhanced import DockingAnimatorEnhanced
from cluspro_adapter import ClusProAdapter
from cluspro_animator import ClusProAnimator
import numpy as np

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['OUTPUT_FOLDER'] = 'outputs'

# Ensure directories exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['OUTPUT_FOLDER'], exist_ok=True)

ALLOWED_EXTENSIONS = {'pdb', 'dlg'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# --- Vina parser ---
def parse_vina_output_pdbqt(pdbqt_path):
    """Parse Vina output PDBQT and return list of poses (coords, elements, bonds)"""
    poses = []
    with open(pdbqt_path, 'r') as f:
        lines = f.readlines()
    pose_coords = []
    pose_elements = []
    pose_atom_names = []
    in_model = False
    for line in lines:
        if line.startswith('MODEL'):
            in_model = True
            pose_coords = []
            pose_elements = []
            pose_atom_names = []
        elif line.startswith('ENDMDL'):
            in_model = False
            if pose_coords:
                poses.append({
                    'coords': np.array(pose_coords, dtype=float),
                    'elements': pose_elements,
                    'bonds': [(i, i+1) for i in range(len(pose_coords)-1)]
                })
        elif in_model and (line.startswith('HETATM') or line.startswith('ATOM')):
            # Parse atom line
            parts = line.split()
            if len(parts) >= 9:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                element = line[76:78].strip() or line[12:14].strip()
                pose_coords.append([x, y, z])
                pose_elements.append(element)
                pose_atom_names.append(element + str(len(pose_coords)))
    # If only one model (no MODEL/ENDMDL), parse all HETATM as one pose
    if not poses and pose_coords:
        poses.append({
            'coords': np.array(pose_coords, dtype=float),
            'elements': pose_elements,
            'bonds': [(i, i+1) for i in range(len(pose_coords)-1)]
        })
    return poses

def parse_vina_log(log_path):
    """Parse Vina log file for detailed pose information"""
    pose_data = []
    with open(log_path, 'r') as f:
        for line in f:
            if line.strip().startswith('-----+'):  # Table header
                break
        for line in f:
            if line.strip() == '' or line.startswith('Writing output'): break
            parts = line.split()
            if len(parts) >= 4:
                try:
                    pose_rank = int(parts[0])
                    affinity = float(parts[1])
                    rmsd_lower = float(parts[2])
                    rmsd_upper = float(parts[3])
                    pose_data.append({
                        'rank': pose_rank,
                        'affinity': affinity,
                        'rmsd_lower': rmsd_lower,
                        'rmsd_upper': rmsd_upper
                    })
                except Exception:
                    continue
    return pose_data

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/generate_animation', methods=['POST'])
def generate_animation():
    try:
        docking_engine = request.form.get('docking_engine', 'autodock4')
        representation = request.form.get('representation', 'sticks')
        if representation not in ['sticks', 'ribbon', 'lines', 'spheres', 'surface', 'cartoon']:
            representation = 'sticks'

        if docking_engine == 'vina':
            vina_pdb_file = request.files.get('vina_pdb_file')
            vina_ligand_file = request.files.get('vina_ligand_file')
            vina_output_file = request.files.get('vina_output_file')
            vina_log_file = request.files.get('vina_log_file')
            vina_num_poses = int(request.form.get('vina_num_poses', 2))
            vina_pose_mode = request.form.get('vina_pose_mode', 'separate')
            if not (vina_pdb_file and vina_ligand_file and vina_output_file and vina_log_file):
                return jsonify({'error': 'All Vina files are required'}), 400
            vina_pdb_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(vina_pdb_file.filename or 'vina_receptor.pdb'))
            vina_ligand_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(vina_ligand_file.filename or 'vina_ligand.pdbqt'))
            vina_output_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(vina_output_file.filename or 'vina_output.pdbqt'))
            vina_log_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(vina_log_file.filename or 'vina_log.txt'))
            vina_pdb_file.save(vina_pdb_path)
            vina_ligand_file.save(vina_ligand_path)
            vina_output_file.save(vina_output_path)
            vina_log_file.save(vina_log_path)
            # Parse poses and detailed pose data
            ligand_poses = parse_vina_output_pdbqt(vina_output_path)
            vina_pose_data = parse_vina_log(vina_log_path)
            N = min(vina_num_poses, len(ligand_poses))
            ligand_poses = ligand_poses[:N]
            vina_pose_data = vina_pose_data[:N]
            protein_path = vina_pdb_path
            animator = DockingAnimatorEnhanced(representation=representation)
            output_files = []
            if vina_pose_mode == 'separate':
                for i, pose in enumerate(ligand_poses):
                    output_filename = f"docking_animation_vina_pose{i+1}.mp4"
                    output_path = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)
                    # Animate single pose, overlay detailed data
                    pose_data = vina_pose_data[i] if i < len(vina_pose_data) else None
                    animator.create_animation_with_pose_data(protein_path, pose, output_path, pose_data)
                    output_files.append({'filename': output_filename, 'download_url': f'/download/{output_filename}', 'pose_data': pose_data})
                return jsonify({
                    'success': True,
                    'message': f'{N} Vina pose animations generated (separate videos)',
                    'results': output_files
                })
            else:  # overlay mode
                output_filename = f"docking_animation_vina_overlay_{N}.mp4"
                output_path = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)
                animator.create_overlay_animation_with_pose_data(protein_path, ligand_poses, output_path, vina_pose_data)
                return jsonify({
                    'success': True,
                    'message': f'Overlay animation for {N} Vina poses generated',
                    'filename': output_filename,
                    'download_url': f'/download/{output_filename}'
                })
        elif docking_engine == 'cluspro':
            # ClusPro mode for protein-protein docking
            cluspro_monomer_file = request.files.get('cluspro_monomer_file')
            cluspro_models_files = request.files.getlist('cluspro_models_files')
            cluspro_scores_file = request.files.get('cluspro_scores_file')
            cluspro_num_models = int(request.form.get('cluspro_num_models', 2))
            cluspro_representation = request.form.get('cluspro_representation', 'cartoon')
            
            if not (cluspro_monomer_file and cluspro_models_files and cluspro_scores_file):
                return jsonify({'error': 'All ClusPro files are required'}), 400
            
            # Save uploaded files
            monomer_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(cluspro_monomer_file.filename or 'monomer.pdb'))
            cluspro_monomer_file.save(monomer_path)
            
            model_paths = []
            for model_file in cluspro_models_files:
                model_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(model_file.filename or 'model.pdb'))
                model_file.save(model_path)
                model_paths.append(model_path)
            
            scores_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(cluspro_scores_file.filename or 'scores.csv'))
            cluspro_scores_file.save(scores_path)
            
            # Create ClusPro animation
            try:
                adapter = ClusProAdapter()
                animation_data = adapter.create_animation_data(monomer_path, model_paths, scores_path, top_k=cluspro_num_models)
                
                if len(animation_data['models']) == 0:
                    return jsonify({'error': 'No models found in ClusPro data'}), 500
                
                # Generate animation
                output_filename = f"cluspro_homodimer_{cluspro_representation}.mp4"
                output_path = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)
                
                animator = ClusProAnimator(representation=cluspro_representation)
                animator.create_homodimer_animation(animation_data, output_path, cluspro_representation)
                
                if not os.path.exists(output_path):
                    return jsonify({'error': 'Failed to generate ClusPro animation'}), 500
                
                return jsonify({
                    'success': True,
                    'message': f'Autopose - ClusPro homodimer animation generated with {len(animation_data["models"])} models',
                    'filename': output_filename,
                    'download_url': f'/download/{output_filename}'
                })
                
            except Exception as e:
                return jsonify({'error': f'Error generating ClusPro animation: {str(e)}'}), 500
        else:
            # AutoDock4 mode (default)
            # Check if files were uploaded
            if 'pdb_file' not in request.files or 'dlg_file' not in request.files:
                return jsonify({'error': 'Both PDB and DLG files are required'}), 400
            
            pdb_file = request.files['pdb_file']
            dlg_file = request.files['dlg_file']
            
            # Check if files are selected
            if pdb_file.filename == '' or dlg_file.filename == '':
                return jsonify({'error': 'Please select both PDB and DLG files'}), 400
            
            # Check file extensions
            if not allowed_file(pdb_file.filename) or not allowed_file(dlg_file.filename):
                return jsonify({'error': 'Invalid file format. Please upload .pdb and .dlg files'}), 400
            
            # Get representation type
            # representation = request.form.get('representation', 'sticks') # This line is now redundant as representation is set above
            # if representation not in ['sticks', 'ribbon', 'lines', 'spheres', 'surface', 'cartoon']:
            #     representation = 'sticks'
            
            # Save uploaded files
            pdb_filename = secure_filename(pdb_file.filename or 'protein.pdb')
            dlg_filename = secure_filename(dlg_file.filename or 'docking.dlg')
            
            pdb_path = os.path.join(app.config['UPLOAD_FOLDER'], pdb_filename)
            dlg_path = os.path.join(app.config['UPLOAD_FOLDER'], dlg_filename)
            
            pdb_file.save(pdb_path)
            dlg_file.save(dlg_path)
            
            # Generate output filename
            output_filename = f"docking_animation_{representation}_{os.path.splitext(pdb_filename)[0]}.mp4"
            output_path = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)
            
            # Create animation
            animator = DockingAnimatorEnhanced(representation=representation)
            animator.create_animation(pdb_path, dlg_path, output_path)
            
            # Check if output file was created
            if not os.path.exists(output_path):
                return jsonify({'error': 'Failed to generate animation'}), 500
            
            return jsonify({
                'success': True,
                'message': f'Animation generated successfully with {representation} representation',
                'filename': output_filename,
                'download_url': f'/download/{output_filename}'
            })
    except Exception as e:
        return jsonify({'error': f'Error generating animation: {str(e)}'}), 500

@app.route('/download/<filename>')
def download_file(filename):
    try:
        file_path = os.path.join(app.config['OUTPUT_FOLDER'], filename)
        if os.path.exists(file_path):
            return send_file(file_path, as_attachment=True)
        else:
            return jsonify({'error': 'File not found'}), 404
    except Exception as e:
        return jsonify({'error': f'Error downloading file: {str(e)}'}), 500

@app.route('/test_sample')
def test_sample():
    """Test endpoint using sample files"""
    try:
        # Use sample files
        pdb_path = 'sample_protein_complete.pdb'
        dlg_path = 'sample_docking_detailed.dlg'  # Use detailed DLG for better testing
        
        if not os.path.exists(pdb_path) or not os.path.exists(dlg_path):
            return jsonify({'error': 'Sample files not found'}), 404
        
        # Test all representations
        representations = ['sticks', 'ribbon', 'lines', 'spheres', 'surface', 'cartoon']
        results = []
        
        for representation in representations:
            output_filename = f"sample_animation_{representation}.mp4"
            output_path = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)
            
            try:
                animator = DockingAnimatorEnhanced(representation=representation)
                animator.create_animation(pdb_path, dlg_path, output_path)
                
                if os.path.exists(output_path):
                    results.append({
                        'representation': representation,
                        'filename': output_filename,
                        'download_url': f'/download/{output_filename}',
                        'status': 'success'
                    })
                else:
                    results.append({
                        'representation': representation,
                        'status': 'failed',
                        'error': 'Output file not created'
                    })
            except Exception as e:
                results.append({
                    'representation': representation,
                    'status': 'failed',
                    'error': str(e)
                })
        
        return jsonify({
            'success': True,
            'message': 'Sample animations generated',
            'results': results
        })
        
    except Exception as e:
        return jsonify({'error': f'Error testing sample: {str(e)}'}), 500

@app.route('/test_sample_vina')
def test_sample_vina():
    try:
        vina_pdb_path = 'sample_vina_receptor.pdb'
        vina_ligand_path = 'sample_vina_ligand.pdbqt'
        vina_output_path = 'sample_vina_output.pdbqt'
        vina_log_path = 'sample_vina_log.txt'
        if not (os.path.exists(vina_pdb_path) and os.path.exists(vina_ligand_path) and os.path.exists(vina_output_path) and os.path.exists(vina_log_path)):
            return jsonify({'error': 'Sample Vina files are not found'}), 404
        
        # Parse poses and detailed pose data
        ligand_poses = parse_vina_output_pdbqt(vina_output_path)
        vina_pose_data = parse_vina_log(vina_log_path)
        
        if len(ligand_poses) == 0:
            return jsonify({'error': 'No poses found in Vina output file'}), 500
        
        representation = 'sticks'
        output_filename = 'sample_vina_animation.mp4'
        output_path = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)
        animator = DockingAnimatorEnhanced(representation=representation)
        
        # Use the first pose for testing
        pose_data = vina_pose_data[0] if vina_pose_data else None
        animator.create_animation_with_pose_data(vina_pdb_path, ligand_poses[0], output_path, pose_data)
        
        if not os.path.exists(output_path):
            return jsonify({'error': 'Failed to generate animation'}), 500
        return jsonify({
            'success': True,
            'message': f'Sample Vina animation generated with {len(ligand_poses)} poses available',
            'filename': output_filename,
            'download_url': f'/download/{output_filename}'
        })
    except Exception as e:
        return jsonify({'error': f'Error testing Vina sample: {str(e)}'}), 500

@app.route('/test_sample_cluspro')
def test_sample_cluspro():
    """Test endpoint using ClusPro sample files"""
    try:
        monomer_path = 'sample_protein_complete.pdb'
        model_paths = ['sample_cluspro_model_001.pdb', 'sample_cluspro_model_002.pdb']
        scores_path = 'sample_cluspro_scores.csv'
        
        if not (os.path.exists(monomer_path) and all(os.path.exists(p) for p in model_paths) and os.path.exists(scores_path)):
            return jsonify({'error': 'Sample ClusPro files not found'}), 404
        
        # Create ClusPro animation data
        adapter = ClusProAdapter()
        animation_data = adapter.create_animation_data(monomer_path, model_paths, scores_path, top_k=2)
        
        if len(animation_data['models']) == 0:
            return jsonify({'error': 'No models found in ClusPro data'}), 500
        
        # Generate animation
        representation = 'cartoon'
        output_filename = 'sample_cluspro_homodimer.mp4'
        output_path = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)
        
        animator = ClusProAnimator(representation=representation)
        animator.create_homodimer_animation(animation_data, output_path, representation)
        
        if not os.path.exists(output_path):
            return jsonify({'error': 'Failed to generate ClusPro animation'}), 500
        
        return jsonify({
            'success': True,
            'message': f'Autopose - ClusPro homodimer animation generated with {len(animation_data["models"])} models',
            'filename': output_filename,
            'download_url': f'/download/{output_filename}'
        })
        
    except Exception as e:
        return jsonify({'error': f'Error testing ClusPro sample: {str(e)}'}), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5002)  # Changed port to 5002 to avoid conflict 