# Autopose v2.1.0 - HADDOCK 2.4 Support

**Release Date**: October 19, 2025  
**Version**: 2.1.0  
**Type**: Major Feature Release

---

## 🎉 Major New Feature: HADDOCK 2.4 Support

Autopose v2.1.0 introduces comprehensive support for **HADDOCK 2.4** protein-protein docking results, making it the first tool to provide intelligent, publication-ready visualizations with automatic quality detection.

---

## ✨ What's New

### HADDOCK 2.4 Integration
- **Complete Pipeline**: Parse cluster_ener.txt, cluster_rmsd.txt, and model PDB files
- **Intelligent Mode Detection**: Automatically detects sparse/CA-only/full-atom models
- **Dual Contact Metrics**: Both heavy-atom (4.5Å) and CA-CA (8.0Å) contacts
- **Production Mode**: Exact geometry for real HADDOCK data
- **Demo Mode**: Enhanced CA-level contacts for educational purposes
- **Comprehensive Overlays**: HADDOCK scores, RMSD, energy components, interface analysis

### Enhanced Animation Quality
- **Visible Rotation**: Deterministic 35° rotation for clear molecular motion
- **Longer Path**: 35Å start distance (increased from 25Å)
- **More Frames**: 60 frames (increased from 40) for smoother animation
- **Smooth Interpolation**: SLERP for rigid-body rotation
- **Camera Zoom**: Automatic zoom to interface in final 30%
- **Frame Hold**: 2-second hold on final frame

### Intelligent Features
- **Auto-Mode Selection**: Detects model quality and selects appropriate mode
- **Self-Explanatory Output**: Clear logging of detection and quality metrics
- **Context-Aware Overlays**: Special handling for single-model runs
- **Debug Information**: Min interchain distance and contact counts
- **Mode Banner**: Clear labeling on first and last frames

### Optional Enhancements (Demo Mode)
- **Pedagogical Snap**: Optional ≤2Å nudge for non-interfacial demo models
- **Safety First**: Default OFF, requires explicit opt-in
- **Demo Only**: Never applied to production mode

---

## 🔧 Technical Improvements

### New Dependencies
```txt
scipy>=1.10.0  # For Rotation and SLERP interpolation
```

### New Components
- **haddock_adapter.py**: Parses HADDOCK output files with quality validation
- **haddock_animator.py**: Generates MP4 animations with mode-aware overlays
- **test_haddock.py**: Comprehensive test suite

### API Enhancements
```python
# Intelligent mode detection (recommended)
adapter = HADDOCKAdapter(mode='auto')

# Force specific modes
adapter = HADDOCKAdapter(mode='production')  # Real data
adapter = HADDOCKAdapter(mode='demo')        # Educational

# Optional pedagogical enhancement
adapter = HADDOCKAdapter(mode='demo', enable_pedagogical_snap=True)
```

---

## 📊 Quality Metrics

### Input Validation
Autopose now automatically detects:
- **Sparse models**: < 300 atoms per chain
- **CA-only models**: < 100 atoms per chain
- **Interfacial models**: Min distance < 8Å
- **Contact quality**: Both heavy-atom and CA-CA contacts

### Mode Detection Output
```
============================================================
MODE DETECTION:
  Detected mode: production
  Atoms/chain: A=1245, B=1183
  Min interchain distance: 2.8 Å
  Heavy-atom contacts (4.5Å): 156
  CA-CA contacts (8.0Å): 23
  Is sparse/CA-only: False/False
  Is interfacial: True
============================================================
```

---

## 🐛 Bug Fixes (from v2.0.0)

### Critical Fixes
- **Final Frame Gap**: Chains now exactly match docked coordinates (< 1e-6 Å RMSD)
- **RMSD Display**: Correct values with proper top-model comparison
- **Overlay Collision**: Fixed text overlap between legend and overlays
- **No Rotation**: Added proper rigid-body rotation with SLERP
- **Contact Calculation**: Fixed to use exact final rendered frame

### Improvements
- **Overlay Positioning**: Professional layout with bbox styling
- **Camera Behavior**: Smooth zoom and proper frame hold
- **Single-Model Handling**: Clear labeling for single-model runs
- **Error Messages**: More informative validation errors

---

## 📦 Supported Docking Engines

Autopose v2.1.0 now supports **four major docking platforms**:

| Engine | Type | Features |
|--------|------|----------|
| **AutoDock 4** | Protein-Ligand | DLG parsing, binding energies, cluster info |
| **AutoDock Vina** | Protein-Ligand | Multi-pose, overlay mode, score analysis |
| **ClusPro** | Protein-Protein | Homodimers, interface analysis, multiple models |
| **HADDOCK 2.4** | Protein-Protein | Intelligent modes, dual metrics, quality detection |

---

## 🚀 Getting Started

### Installation
```bash
git clone https://github.com/QiaoJiang08/Autopose.git
cd Autopose
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python3 app.py
```

### Quick Test
```bash
# Test HADDOCK implementation
python3 test_haddock.py

# Or use web interface
# Open http://localhost:5002
# Click "Test HADDOCK Sample" button
```

### Usage Example
```python
from haddock_adapter import HADDOCKAdapter
from haddock_animator import HADDOCKAnimator

# Create animation data (auto-detect mode)
adapter = HADDOCKAdapter(mode='auto')
animation_data = adapter.create_animation_data(
    monomer_pdb='monomer.pdb',
    model_pdbs=['model_1.pdb'],
    cluster_ener_path='cluster_ener.txt',
    cluster_rmsd_path='cluster_rmsd.txt',
    top_k=1
)

# Generate animation
animator = HADDOCKAnimator(representation='cartoon')
animator.create_haddock_animation(
    animation_data,
    output_path='haddock_docking.mp4',
    representation='cartoon'
)
```

---

## 📚 Documentation

- **[README.md](README.md)**: Main documentation and quick start
- **[HADDOCK_README.md](HADDOCK_README.md)**: Detailed HADDOCK guide
- **[CHANGELOG.md](CHANGELOG.md)**: Complete version history
- **[Bugsv1Haddock.md](Bugsv1Haddock.md)**: V1 bug fixes reference
- **[Bugsv2Haddock.md](Bugsv2Haddock.md)**: V2 improvements reference

---

## 🎯 Use Cases

### Research & Publications
- **Production Mode**: Exact geometry, heavy-atom contacts, professional output
- **Multiple Models**: Compare top-ranked HADDOCK poses
- **Publication Quality**: 1080p MP4 ready for presentations

### Education & Teaching
- **Demo Mode**: CA-level contacts for simplified models
- **Pedagogical Enhancements**: Optional cosmetic improvements
- **Clear Labeling**: Students understand when using placeholder data

### Quality Control
- **Automatic Validation**: Detects sparse or non-interfacial models
- **Comprehensive Metrics**: Heavy-atom and CA-CA contacts
- **Debug Information**: Min distances and quality indicators

---

## 🔄 Breaking Changes

**None.** This release is fully backward compatible with v2.0.0.

All new features are opt-in via new parameters:
- Default `mode='auto'` works with any input
- `enable_pedagogical_snap=False` by default
- Existing AutoDock4, Vina, and ClusPro functionality unchanged

---

## 🙏 Acknowledgments

Special thanks to the HADDOCK development team and the scientific community for feedback and suggestions that made this release possible.

---

## 📥 Download

**Source Code**: [v2.1.0.tar.gz](https://github.com/QiaoJiang08/Autopose/archive/v2.1.0.tar.gz)  
**Repository**: [https://github.com/QiaoJiang08/Autopose](https://github.com/QiaoJiang08/Autopose)  
**Issues**: [Report Bugs](https://github.com/QiaoJiang08/Autopose/issues)

---

## 📈 Statistics

- **New Files**: 3 (haddock_adapter.py, haddock_animator.py, test_haddock.py)
- **Modified Files**: 4 (app.py, templates/index.html, README.md, requirements.txt)
- **Lines Added**: ~2,500
- **Test Coverage**: Comprehensive HADDOCK test suite
- **Documentation**: 200+ lines of HADDOCK-specific docs

---

**Ready to visualize your HADDOCK results with intelligent, publication-quality animations?**  
Get started with Autopose v2.1.0! ⭐

---

*For previous releases, see [CHANGELOG.md](CHANGELOG.md)*

