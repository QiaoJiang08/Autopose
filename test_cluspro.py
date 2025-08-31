#!/usr/bin/env python3
"""Test script for ClusPro functionality"""

import os
import sys
from cluspro_adapter import ClusProAdapter
from cluspro_animator import ClusProAnimator

def test_cluspro_pipeline():
    """Test the complete ClusPro pipeline"""
    print("🧪 Testing ClusPro Pipeline...")
    
    try:
        # Test 1: Adapter functionality
        print("\n1. Testing ClusPro Adapter...")
        adapter = ClusProAdapter()
        
        # Parse scores
        scores = adapter.parse_cluspro_scores('sample_cluspro_scores.csv')
        print(f"   ✓ Parsed {len(scores)} scores")
        print(f"   ✓ Top score: {scores[0]['score']:.1f} kcal/mol")
        
        # Parse models
        model_data = adapter.parse_cluspro_model('sample_cluspro_model_001.pdb', 'sample_protein_complete.pdb')
        print(f"   ✓ Parsed model with {len(model_data['ref_chain']['coords'])} reference atoms")
        print(f"   ✓ Parsed model with {len(model_data['moving_chain']['coords'])} moving atoms")
        
        # Detect interface
        interface_data = adapter.detect_interface_residues(model_data)
        print(f"   ✓ Detected {interface_data['contact_count']} interface contacts")
        
        # Create animation data
        animation_data = adapter.create_animation_data(
            'sample_protein_complete.pdb',
            ['sample_cluspro_model_001.pdb', 'sample_cluspro_model_002.pdb'],
            'sample_cluspro_scores.csv',
            top_k=2
        )
        print(f"   ✓ Created animation data with {len(animation_data['models'])} models")
        
        # Test 2: Animator functionality
        print("\n2. Testing ClusPro Animator...")
        animator = ClusProAnimator(representation='cartoon')
        print(f"   ✓ Animator initialized with {animator.representation} representation")
        
        # Test monomer parsing
        monomer_atoms, monomer_elements, monomer_residues = animator.parse_monomer_pdb('sample_protein_complete.pdb')
        print(f"   ✓ Parsed monomer with {len(monomer_atoms)} atoms")
        
        # Test trajectory generation
        trajectory = adapter.generate_docking_trajectory(model_data, num_frames=10)
        print(f"   ✓ Generated trajectory with {len(trajectory)} frames")
        
        print("\n✅ All ClusPro tests passed!")
        return True
        
    except Exception as e:
        print(f"\n❌ ClusPro test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def test_cluspro_animation_generation():
    """Test actual animation generation"""
    print("\n🎬 Testing ClusPro Animation Generation...")
    
    try:
        # Create animation data
        adapter = ClusProAdapter()
        animation_data = adapter.create_animation_data(
            'sample_protein_complete.pdb',
            ['sample_cluspro_model_001.pdb'],
            'sample_cluspro_scores.csv',
            top_k=1
        )
        
        # Generate animation
        animator = ClusProAnimator(representation='cartoon')
        output_path = 'test_cluspro_animation.mp4'
        
        print(f"   Generating animation to {output_path}...")
        animator.create_homodimer_animation(animation_data, output_path, 'cartoon')
        
        if os.path.exists(output_path):
            print(f"   ✅ Animation generated successfully: {output_path}")
            print(f"   File size: {os.path.getsize(output_path)} bytes")
            return True
        else:
            print("   ❌ Animation file not created")
            return False
            
    except Exception as e:
        print(f"   ❌ Animation generation failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("🚀 ClusPro Extension Test Suite")
    print("=" * 40)
    
    # Test basic functionality
    if test_cluspro_pipeline():
        print("\n🎯 Basic functionality tests passed!")
        
        # Test animation generation
        if test_cluspro_animation_generation():
            print("\n🎉 All tests passed! ClusPro extension is working correctly.")
        else:
            print("\n⚠️  Animation generation test failed, but basic functionality works.")
    else:
        print("\n💥 Basic functionality tests failed. Please check the implementation.")
        sys.exit(1)
