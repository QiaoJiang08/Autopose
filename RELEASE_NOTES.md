# Release Notes - Autopose v2.0.0

## 🎉 Autopose v2.0.0 - ClusPro Extension

**Release Date**: January 20, 2025  
**Version**: 2.0.0  
**Type**: Major Release  

## 🚀 Major Release: Protein-Protein Docking Support

Autopose v2.0.0 introduces comprehensive support for ClusPro protein-protein docking animations, making it the first tool to provide professional-quality visualizations of homodimer formation.

## ✨ New Features

### **ClusPro Integration**
- **Protein-Protein Docking**: Animate homodimer formation from separated monomers
- **Interface Analysis**: Automatic detection and visualization of protein-protein contacts
- **Multiple Models**: Support for animating top-K ClusPro models with clustering information
- **Energy Breakdown**: Comprehensive overlays showing vdW, electrostatic, and desolvation energies
- **Professional Visualization**: High-quality animations perfect for publications and presentations

### **Enhanced Capabilities**
- **Improved AutoDock4**: Better overlay parsing and display
- **Enhanced Vina**: Optimized multi-pose handling and positioning
- **Better Error Handling**: Robust validation and user-friendly messages
- **Memory Optimization**: Efficient frame generation and cleanup

## 🐛 Bug Fixes

- Fixed AutoDock4 overlay display issues
- Resolved Vina pose positioning problems
- Corrected ClusPro frame generation
- Fixed file upload validation

## 🔄 Changes

- **Project Renamed**: From "Protein Docking Animation Generator" to "Autopose"
- **Repository URL**: Changed to `github.com/QiaoJiang08/Autopose`
- **Console Command**: Now available as `autopose` command
- **Animation Titles**: All animations now include "Autopose -" prefix

## 📦 Installation

```bash
git clone https://github.com/QiaoJiang08/Autopose.git
cd Autopose
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
brew install ffmpeg  # macOS
python3 app.py
```

## 🎯 What's New

- **ClusPro Support**: Complete protein-protein docking animation pipeline
- **Interface Visualization**: Automatic contact detection and highlighting
- **Multiple Model Animation**: Animate several docking models simultaneously
- **Professional Output**: Publication-ready MP4 animations
- **Comprehensive Documentation**: Detailed guides and examples

## 🔧 Technical Details

- **Python 3.8+** support
- **FFmpeg** integration for video generation
- **Matplotlib** 3D visualization
- **BioPython** molecular structure parsing
- **Flask** web interface

## 📚 Documentation

- [README.md](README.md) - Main documentation
- [ClusPro Guide](CLUSPRO_README.md) - ClusPro-specific instructions
- [Project Summary](PROJECT_SUMMARY.md) - Comprehensive overview
- [Changelog](CHANGELOG.md) - Detailed change history

## 🎬 Sample Data

Includes complete sample files for:
- **AutoDock4**: Protein + DLG files
- **AutoDock Vina**: Receptor + Ligand + Output + Log files
- **ClusPro**: Monomer + Models + Scores files

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.

## 🎉 Acknowledgments

Special thanks to the scientific community for feedback and suggestions that made this release possible.

---

**Ready to revolutionize molecular visualization?** Get started with Autopose v2.0.0! ⭐

**Download**: [Source Code](https://github.com/QiaoJiang08/Autopose/archive/v2.0.0.zip)  
**Documentation**: [Full Guide](README.md)  
**Issues**: [Report Bugs](https://github.com/QiaoJiang08/Autopose/issues)  
**Discussions**: [Community](https://github.com/QiaoJiang08/Autopose/discussions)
