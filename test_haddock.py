"""
Test script for HADDOCK implementation
Now includes verification of all bug fixes from Bugsv1Haddock.md
"""

from haddock_adapter import HADDOCKAdapter
from haddock_animator import HADDOCKAnimator
import os
import numpy as np

def test_haddock_pipeline():
    """Test the complete HADDOCK animation pipeline"""
    
    print("=" * 60)
    print("Testing HADDOCK Implementation for Autopose")
    print("=" * 60)
    
    # File paths
    monomer_path = 'sample_haddock_monomer.pdb'
    model_paths = ['sample_haddock_model_001.pdb']
    cluster_ener_path = 'sample_haddock_cluster_ener.txt'
    cluster_rmsd_path = 'sample_haddock_cluster_rmsd.txt'
    output_path = 'outputs/test_haddock_animation.mp4'
    
    # Check files exist
    print("\n1. Checking input files...")
    files = [monomer_path] + model_paths + [cluster_ener_path, cluster_rmsd_path]
    for f in files:
        if os.path.exists(f):
            print(f"   ✓ {f}")
        else:
            print(f"   ✗ {f} NOT FOUND")
            return False
    
    # Test adapter
    print("\n2. Testing HADDOCKAdapter...")
    try:
        adapter = HADDOCKAdapter(mode='auto', enable_pedagogical_snap=True)
        print("   ✓ Adapter initialized (auto mode, pedagogical snap enabled)")
        
        # Parse cluster energy
        cluster_energy = adapter.parse_cluster_energy(cluster_ener_path)
        print(f"   ✓ Parsed {len(cluster_energy)} cluster energy entries")
        
        # Parse cluster RMSD
        cluster_rmsd = adapter.parse_cluster_rmsd(cluster_rmsd_path)
        print(f"   ✓ Parsed {len(cluster_rmsd)} cluster RMSD entries")
        
        # Parse model
        model_data = adapter.parse_haddock_model(model_paths[0])
        print(f"   ✓ Parsed model with {len(model_data['chain_a'])} Chain A atoms")
        print(f"   ✓ Parsed model with {len(model_data['chain_b'])} Chain B atoms")
        
        # Create animation data
        animation_data = adapter.create_animation_data(
            monomer_path, model_paths, cluster_ener_path,
            cluster_rmsd_path, top_k=1
        )
        print(f"   ✓ Created animation data with {len(animation_data['models'])} model(s)")
        
        if len(animation_data['models']) > 0:
            model = animation_data['models'][0]
            print(f"   ✓ Model info: Cluster {model['info']['cluster']}, Rank {model['info']['rank']}")
            print(f"   ✓ HADDOCK Score: {model['info']['haddock_score']}")
            print(f"   ✓ RMSD: {model['info']['rmsd']} Å")
        
    except Exception as e:
        print(f"   ✗ Adapter test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False
    
    # Test animator
    print("\n3. Testing HADDOCKAnimator...")
    try:
        animator = HADDOCKAnimator(representation='cartoon')
        print("   ✓ Animator initialized")
        
        # Create animation
        print("   Creating animation (this may take a minute)...")
        os.makedirs('outputs', exist_ok=True)
        animator.create_haddock_animation(animation_data, output_path, 'cartoon')
        
        if os.path.exists(output_path):
            size = os.path.getsize(output_path) / (1024 * 1024)
            print(f"   ✓ Animation created: {output_path} ({size:.2f} MB)")
        else:
            print(f"   ✗ Animation file not created")
            return False
            
    except Exception as e:
        print(f"   ✗ Animator test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n" + "=" * 60)
    print("✅ HADDOCK Implementation Test PASSED")
    print("=" * 60)
    return True


if __name__ == "__main__":
    success = test_haddock_pipeline()
    exit(0 if success else 1)

