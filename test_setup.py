#!/usr/bin/env python3
"""
Test script for AutoDock Animation Generator
This script verifies the installation and creates sample files for testing.
"""

import os
import sys
import subprocess

def test_imports():
    """Test if all required packages can be imported"""
    print("Testing package imports...")
    
    try:
        import flask
        print("✓ Flask imported successfully")
    except ImportError as e:
        print(f"✗ Flask import failed: {e}")
        return False
    
    try:
        import Bio
        print("✓ Biopython imported successfully")
    except ImportError as e:
        print(f"✗ Biopython import failed: {e}")
        return False
    
    try:
        import matplotlib
        print("✓ Matplotlib imported successfully")
    except ImportError as e:
        print(f"✗ Matplotlib import failed: {e}")
        return False
    
    try:
        import moviepy
        print("✓ MoviePy imported successfully")
    except ImportError as e:
        print(f"✗ MoviePy import failed: {e}")
        return False
    
    try:
        import numpy
        print("✓ NumPy imported successfully")
    except ImportError as e:
        print(f"✗ NumPy import failed: {e}")
        return False
    
    return True

def test_ffmpeg():
    """Test if FFmpeg is available"""
    print("\nTesting FFmpeg installation...")
    
    try:
        result = subprocess.run(['ffmpeg', '-version'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print("✓ FFmpeg is available")
            return True
        else:
            print("✗ FFmpeg returned non-zero exit code")
            return False
    except FileNotFoundError:
        print("✗ FFmpeg not found in PATH")
        print("  Please install FFmpeg:")
        print("  macOS: brew install ffmpeg")
        print("  Ubuntu: sudo apt install ffmpeg")
        print("  Windows: Download from https://ffmpeg.org/download.html")
        return False
    except subprocess.TimeoutExpired:
        print("✗ FFmpeg test timed out")
        return False

def create_sample_files():
    """Create sample PDB and DLG files for testing"""
    print("\nCreating sample files...")
    
    # Sample PDB file (simple protein structure)
    sample_pdb = """ATOM      1  N   ALA A   1      27.451  14.105   5.420  1.00 20.00           N
ATOM      2  CA  ALA A   1      26.325  13.134   5.520  1.00 20.00           C
ATOM      3  C   ALA A   1      25.085  13.850   6.100  1.00 20.00           C
ATOM      4  O   ALA A   1      24.000  13.300   6.200  1.00 20.00           O
ATOM      5  CB  ALA A   1      26.700  12.200   4.350  1.00 20.00           C
ATOM      6  N   ALA A   2      25.200  15.100   6.600  1.00 20.00           N
ATOM      7  CA  ALA A   2      24.100  15.900   7.200  1.00 20.00           C
ATOM      8  C   ALA A   2      23.800  17.200   6.400  1.00 20.00           C
ATOM      9  O   ALA A   2      22.700  17.800   6.500  1.00 20.00           O
ATOM     10  CB  ALA A   2      24.500  16.200   8.650  1.00 20.00           C
TER
END"""
    
    # Sample DLG file (improved docking results with clear movement)
    sample_dlg = """DOCKED: MODEL        1
DOCKED: ATOM      1  C   LIG     1      30.000  20.000  10.000
DOCKED: ATOM      2  N   LIG     1      30.500  20.500  10.500
DOCKED: ATOM      3  O   LIG     1      29.500  19.500   9.500
DOCKED: MODEL        2
DOCKED: ATOM      1  C   LIG     1      28.500  18.500   8.500
DOCKED: ATOM      2  N   LIG     1      29.000  19.000   9.000
DOCKED: ATOM      3  O   LIG     1      28.000  18.000   8.000
DOCKED: MODEL        3
DOCKED: ATOM      1  C   LIG     1      27.000  17.000   7.000
DOCKED: ATOM      2  N   LIG     1      27.500  17.500   7.500
DOCKED: ATOM      3  O   LIG     1      26.500  16.500   6.500
DOCKED: MODEL        4
DOCKED: ATOM      1  C   LIG     1      26.000  16.000   6.000
DOCKED: ATOM      2  N   LIG     1      26.500  16.500   6.500
DOCKED: ATOM      3  O   LIG     1      25.500  15.500   5.500
DOCKED: MODEL        5
DOCKED: ATOM      1  C   LIG     1      25.500  15.500   5.500
DOCKED: ATOM      2  N   LIG     1      26.000  16.000   6.000
DOCKED: ATOM      3  O   LIG     1      25.000  15.000   5.000"""
    
    # Write sample files
    with open('sample_protein.pdb', 'w') as f:
        f.write(sample_pdb)
    
    with open('sample_docking.dlg', 'w') as f:
        f.write(sample_dlg)
    
    print("✓ Sample files created:")
    print("  - sample_protein.pdb")
    print("  - sample_docking.dlg")

def test_docking_animator():
    """Test the DockingAnimator class"""
    print("\nTesting DockingAnimator...")
    
    try:
        from docking_animator import DockingAnimator
        
        # Create animator instance
        animator = DockingAnimator()
        print("✓ DockingAnimator created successfully")
        
        # Test PDB parsing
        protein_atoms, protein_colors = animator.parse_pdb('sample_protein.pdb')
        print(f"✓ PDB parsed: {len(protein_atoms)} atoms found")
        
        # Test DLG parsing
        ligand_poses = animator.parse_dlg('sample_docking.dlg')
        print(f"✓ DLG parsed: {len(ligand_poses)} poses found")
        
        return True
        
    except Exception as e:
        print(f"✗ DockingAnimator test failed: {e}")
        return False

def main():
    """Run all tests"""
    print("AutoDock Animation Generator - Setup Test")
    print("=" * 50)
    
    # Test imports
    if not test_imports():
        print("\n❌ Import test failed. Please install missing packages:")
        print("pip install -r requirements.txt")
        return False
    
    # Test FFmpeg
    if not test_ffmpeg():
        print("\n❌ FFmpeg test failed. Please install FFmpeg.")
        return False
    
    # Create sample files
    create_sample_files()
    
    # Test DockingAnimator
    if not test_docking_animator():
        print("\n❌ DockingAnimator test failed.")
        return False
    
    print("\n✅ All tests passed! The application is ready to use.")
    print("\nTo start the web application:")
    print("python app.py")
    print("\nThen open http://localhost:5000 in your browser")
    print("\nYou can use the sample files (sample_protein.pdb, sample_docking.dlg) for testing.")
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 