# ClusPro Extension for Autopose

## Overview
Adds support for ClusPro protein-protein docking animations to Autopose.

## Features
- **Homodimer Animation**: Visualize protein-protein docking from separated to docked state
- **Multiple Models**: Animate top-K ClusPro models with different configurations
- **Interface Analysis**: Detect and display interface contacts between chains
- **Score Overlays**: Display ClusPro scoring (energy, vdW, electrostatic, desolvation)
- **Visualization Styles**: Cartoon, sticks, surface, ribbon, spheres

## Input Files Required
1. **Monomer PDB**: Reference protein structure (single chain)
2. **ClusPro Models**: Multiple PDB files with docked complexes
3. **Scores CSV**: ClusPro scoring information

## Usage

### Web Interface
1. Select "ClusPro (Protein-Protein)" as docking engine
2. Upload monomer PDB, model PDBs, and scores CSV
3. Choose number of models and visualization style
4. Generate and download MP4 animation

### Python API
```python
from cluspro_adapter import ClusProAdapter
from cluspro_animator import ClusProAnimator

# Create animation data
adapter = ClusProAdapter()
animation_data = adapter.create_animation_data(
    'monomer.pdb', ['model1.pdb', 'model2.pdb'], 'scores.csv', top_k=2
)

# Generate animation
animator = ClusProAnimator(representation='cartoon')
animator.create_homodimer_animation(animation_data, 'output.mp4', 'cartoon')
```

## Sample Files
- `sample_protein_complete.pdb`: Reference monomer
- `sample_cluspro_model_001.pdb`: Sample complex 1
- `sample_cluspro_model_002.pdb`: Sample complex 2
- `sample_cluspro_scores.csv`: Sample scores

## Test
```bash
python3 test_cluspro.py
```

## Output
- **Format**: MP4 (1080p, 12 FPS)
- **Content**: Animated homodimer formation with score overlays
- **Duration**: Configurable (default: 120 frames)

## Requirements
- Python 3.8+
- Biopython
- Matplotlib
- FFmpeg
- NumPy

The extension is now fully functional and integrated into Autopose!
