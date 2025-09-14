# Changelog

All notable changes to Autopose will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-01-20

### Added
- **ClusPro Support**: Complete protein-protein docking animation support
- **Interface Analysis**: Automatic detection and visualization of protein-protein interfaces
- **Multiple Model Animation**: Support for animating multiple ClusPro models with clustering information
- **Enhanced Overlays**: Comprehensive energy breakdown (vdW, electrostatic, desolvation)
- **Professional Visualization**: High-quality homodimer formation animations
- **Sample Data**: Complete set of ClusPro test files and examples

### Enhanced
- **AutoDock4 Overlays**: Improved parsing and display of binding energies, inhibition constants, and cluster information
- **Vina Multi-Pose**: Better handling of multiple poses with separate and overlay modes
- **Error Handling**: Robust validation and user-friendly error messages
- **Memory Management**: Optimized frame generation and cleanup
- **Documentation**: Comprehensive guides and project summaries

### Fixed
- **AutoDock4 Overlay Display**: Resolved issues with overlay information not appearing
- **Vina Pose Positioning**: Fixed ligand positioning and movement issues
- **ClusPro Frame Generation**: Corrected frame naming and FFmpeg integration
- **File Upload Validation**: Improved file format checking and error reporting
- **Import Issues**: Removed unused moviepy dependency causing launch failures

### Changed
- **Project Name**: Renamed from "Protein Docking Animation Generator" to "Autopose"
- **Repository Structure**: Organized files for better maintainability
- **Installation Process**: Streamlined setup with improved dependency management

## [1.0.0] - 2024-12-01

### Added
- **AutoDock4 Support**: Complete DLG file parsing and animation generation
- **AutoDock Vina Support**: PDBQT file handling with multi-pose capabilities
- **Web Interface**: User-friendly Flask-based web application
- **Multiple Visualization Styles**: Sticks, ribbon, lines, spheres, surface, cartoon
- **Professional Output**: 1080p MP4 videos with H.264 encoding
- **Sample Data**: Test files for AutoDock4 and Vina workflows

### Technical
- **Core Animation Engine**: Matplotlib-based 3D molecular visualization
- **Video Generation**: FFmpeg integration for high-quality MP4 output
- **File Parsing**: Robust PDB, DLG, and PDBQT file handling
- **Trajectory Generation**: Smooth interpolation between docking poses
- **Overlay System**: Real-time display of docking metrics and scores

---

## Release Notes Format

### Version Numbering
- **MAJOR** (2.0.0): Breaking changes, major new features
- **MINOR** (2.1.0): New features, backward compatible
- **PATCH** (2.0.1): Bug fixes, minor improvements

### Release Types
- **Feature Release**: New functionality (minor version bump)
- **Bug Fix Release**: Bug fixes only (patch version bump)
- **Major Release**: Significant changes or breaking changes (major version bump)

### Release Checklist
- [ ] Update version numbers in all files
- [ ] Update CHANGELOG.md with new changes
- [ ] Create GitHub release with proper tags
- [ ] Update documentation if needed
- [ ] Test all functionality
- [ ] Verify sample data works correctly
