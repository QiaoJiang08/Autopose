# HADDOCK 2.4 Support in Autopose

Autopose now supports HADDOCK 2.4 protein-protein docking results, enabling professional visualization of homodimer and heterodimer formation.

## 📋 Overview

HADDOCK (High Ambiguity Driven protein-protein DOCKing) is a widely-used information-driven flexible docking approach for biomolecular complexes. Autopose integrates HADDOCK results to create high-quality MP4 animations showing:

- **Protein-Protein Docking**: Realistic approach trajectories for protein complex formation
- **HADDOCK Scores**: Display weighted energy scores and components
- **RMSD Analysis**: Show structural similarity to top-ranked models
- **Energy Breakdown**: Visualize van der Waals, electrostatic, desolvation, and AIR energies
- **Interface Contacts**: Highlight protein-protein interface residues

## 🗂️ Required Files

To generate a HADDOCK animation, you need:

1. **Monomer PDB** (`monomer.pdb`): Reference structure for chain A alignment
2. **Docked Model PDB(s)** (`best_1.pdb`, `best_2.pdb`, etc.): Docked complex structures with chains A and B
3. **Cluster Energy File** (`cluster_ener.txt`): HADDOCK scores and energy components
4. **Cluster RMSD File** (`cluster_rmsd.txt`): RMSD values from overall lowest energy structure

### File Format Details

#### cluster_ener.txt Format
```
# Cluster  Rank  HADDOCK-Score  Evdw      Eelec     Edesol    Eair
  1        1     -125.3         -45.2     -78.9     -12.5     -8.7
  2        2     -118.5         -42.8     -72.1     -10.3     -7.3
```

#### cluster_rmsd.txt Format
```
# Cluster  RMSD_from_overall_lowest
  1        0.0
  2        2.3
```

#### Model PDB Requirements
- Must contain chain A (receptor/reference)
- Must contain chain B (ligand/mobile chain)
- Standard PDB format with ATOM records

## 🎬 Animation Features

### Overlay Information Displayed

**During Animation:**
- Cluster ID and Rank
- HADDOCK Score (weighted energy)
- RMSD from top pose
- Van der Waals energy (E_vdW)
- Electrostatic energy (E_elec)
- Desolvation energy (E_desol)

**Final Frame:**
- Interface contact counts
- Number of residues in binding interface

### Visualization Styles

Choose from multiple representation styles:
- **Cartoon** (default): Smooth ribbon representation
- **Sticks**: Ball-and-stick model
- **Surface**: Molecular surface rendering
- **Ribbon**: Backbone trace
- **Spheres**: Space-filling model

## 🚀 Usage

### Web Interface

1. Navigate to `http://localhost:5002`
2. Select **HADDOCK 2.4** as docking engine
3. Upload required files:
   - Monomer PDB file
   - One or more docked model PDB files
   - cluster_ener.txt
   - cluster_rmsd.txt
4. Select number of models to animate (1-5)
5. Choose visualization style
6. Click "Generate Animation"

### Test with Sample Files

```bash
# Test HADDOCK implementation
python3 test_haddock.py

# Or use the web interface test button
# Click "🧪 Test HADDOCK Sample" button
```

### Programmatic Usage

```python
from haddock_adapter import HADDOCKAdapter
from haddock_animator import HADDOCKAnimator

# Parse HADDOCK results
adapter = HADDOCKAdapter()
animation_data = adapter.create_animation_data(
    monomer_pdb='monomer.pdb',
    model_pdbs=['best_1.pdb', 'best_2.pdb'],
    cluster_ener_path='cluster_ener.txt',
    cluster_rmsd_path='cluster_rmsd.txt',
    top_k=2
)

# Generate animation
animator = HADDOCKAnimator(representation='cartoon')
animator.create_haddock_animation(
    animation_data,
    output_path='haddock_docking.mp4',
    representation='cartoon'
)
```

## 📊 Algorithm Details

### Trajectory Generation

1. **Reference Alignment**: Chain A from docked model is superposed to reference monomer
2. **Coordinate Transformation**: Same transformation applied to chain B
3. **Start Position Generation**: Chain B positioned 25 Å away from final position
4. **Smooth Interpolation**: Cubic ease-in-out interpolation for realistic motion
5. **Rotation During Approach**: 30° rotation around approach vector for dynamic visualization

### Interface Detection

Interface residues are detected using a 4.5 Å distance cutoff between chain A and chain B atoms.

### Energy Parsing

The adapter extracts:
- **HADDOCK Score**: Weighted sum of all energy terms
- **E_vdW**: Van der Waals contribution
- **E_elec**: Electrostatic interaction energy
- **E_desol**: Desolvation energy
- **E_air**: Ambiguous interaction restraint energy

## 🔧 Technical Implementation

### Core Components

- **`haddock_adapter.py`**: Parses HADDOCK output files and prepares animation data
- **`haddock_animator.py`**: Generates MP4 animations from parsed data
- **Integration**: Seamlessly integrated into existing Autopose pipeline

### Key Features

- **Rigid-Body Transformations**: Proper rotation and translation calculations
- **BioPython Integration**: Uses Bio.PDB for structure parsing and superposition
- **FFmpeg Rendering**: High-quality 1080p MP4 output
- **Memory Efficient**: Frame-by-frame rendering with automatic cleanup

## 📝 Example Workflow

1. **Run HADDOCK Docking**: Complete your HADDOCK run and collect output files
2. **Identify Top Models**: Review `cluster_ener.txt` to select models
3. **Prepare Files**: Gather monomer PDB, model PDBs, and cluster files
4. **Generate Animation**: Upload to Autopose and create visualization
5. **Download and Share**: Get publication-ready MP4 animation

## 🎯 Best Practices

- **Model Selection**: Animate top 1-3 models for clarity
- **Visualization Style**: Use cartoon for protein-protein, sticks for detailed views
- **File Preparation**: Ensure PDB files are properly formatted with correct chain IDs
- **Quality Control**: Verify HADDOCK scores and RMSD values in overlay match your data

## 🐛 Troubleshooting

### Common Issues

**Issue**: Chain A/B not found in model PDB
- **Solution**: Ensure your PDB file has proper chain identifiers (A and B)

**Issue**: No animation generated
- **Solution**: Check that cluster files are properly formatted with correct headers

**Issue**: Overlay shows "N/A" values
- **Solution**: Verify cluster_ener.txt and cluster_rmsd.txt are correctly formatted

### File Format Requirements

- PDB files must follow standard PDB format
- Cluster files should have tab or space-separated values
- Headers should contain "Cluster" keyword for proper parsing

## 📚 References

- HADDOCK Web Server: https://wenmr.science.uu.nl/haddock2.4/
- HADDOCK Documentation: https://www.bonvinlab.org/software/haddock2.4/
- Autopose GitHub: https://github.com/QiaoJiang08/Autopose

## 🤝 Contributing

Contributions to improve HADDOCK support are welcome! Please submit issues or pull requests on GitHub.

---

**Ready to visualize your HADDOCK results?** Start with the sample files to test the implementation! ⭐

