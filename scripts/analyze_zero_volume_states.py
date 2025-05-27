"""
Analyze zero-volume states in LQG for valence-4 nodes.

This script analyzes the mathematical structure of zero-volume states
in Loop Quantum Gravity for 4-valent nodes. It extends the previous 
work by investigating patterns and properties of these states.
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from collections import defaultdict, Counter
from lqg_utils import CF12j_numeric, get_allowed_j12_values

def build_volume_matrix(j1, j2, j3, j4):
    """
    Build the V^2 matrix for the given spin configuration using the correct
    intersection J = J₁₂ ∩ J₃₄ for the 4-valent intertwiner space.
    
    Args:
        j1, j2, j3, j4: The four spins of the valence-4 vertex
        
    Returns:
        V2 matrix and the list of intermediate coupling values (the intersection)
    """
    # Get allowed intermediate coupling values for each pair
    J12 = get_allowed_j12_values(j1, j2)
    J34 = get_allowed_j12_values(j3, j4)
    
    # The true basis for the 4-valent intertwiner space is the intersection
    J = sorted(set(J12) & set(J34))
    
    # If intersection is empty, this is a trivial zero-volume state
    if not J:
        # Return empty matrix to indicate no valid intertwiner space
        return np.array([]), []
    
    # Build V2 matrix using the intersection values
    V2 = np.zeros((len(J), len(J)))
    for a, ja in enumerate(J):
        for b, jb in enumerate(J):
            # Ordering: (j1,j2,j12), (j3,j4,j12), (j1,j2,j34), (j3,j4,j34)
            V2[a, b] = CF12j_numeric(j1, j2, ja, j3, j4, ja, j1, j2, jb, j3, j4, jb)
    
    return V2, J

def analyze_matrix_properties(matrix):
    """
    Analyze mathematical properties of a matrix.
    
    Args:
        matrix: The matrix to analyze
        
    Returns:
        Dictionary with various matrix properties
    """
    # Calculate basic properties
    det = np.linalg.det(matrix)
    eigenvals = np.linalg.eigvalsh(matrix)
    rank = np.linalg.matrix_rank(matrix)
    size = matrix.shape[0]
    
    # Count zero and non-zero eigenvalues (within numerical precision)
    zero_threshold = 1e-10
    zero_eigenvals = sum(1 for ev in eigenvals if abs(ev) < zero_threshold)
    nonzero_eigenvals = size - zero_eigenvals
    
    # Calculate the kernel dimension
    kernel_dim = size - rank
    
    # Check if matrix is singular
    is_singular = abs(det) < zero_threshold
    
    return {
        'determinant': det,
        'eigenvalues': eigenvals,
        'min_eigenvalue': min(abs(val) for val in eigenvals),
        'rank': rank,
        'size': size,
        'zero_eigenvals': zero_eigenvals,
        'nonzero_eigenvals': nonzero_eigenvals,
        'kernel_dimension': kernel_dim,
        'is_singular': is_singular
    }

def analyze_zero_volume_states():
    """
    Analyze patterns in zero-volume states.
    """
    # Define spin range to check
    spin_vals = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    
    # Dictionary to store results by kernel dimension
    results_by_kernel_dim = defaultdict(list)
    has_spin_half_count = 0
    total_states = 0
    kernel_dim_counts = Counter()
    
    # Track matrix sizes
    matrix_sizes = Counter()
    
    # Track progress
    total = len(spin_vals)**4
    count = 0
      # Track different types of zero-volume states
    trivial_zero_volume_count = 0  # When J₁₂ ∩ J₃₄ = ∅
    non_trivial_zero_volume_count = 0  # When intersection exists but kernel is non-trivial
    no_kernel_count = 0  # When intersection exists but no kernel
    
    print(f"Analyzing {total} spin configurations in range {min(spin_vals)}-{max(spin_vals)}...")
    
    for (j1, j2, j3, j4) in product(spin_vals, repeat=4):
        count += 1
        if count % (total // 20) == 0:
            print(f"Progress: {count}/{total} ({count/total*100:.1f}%)")
        
        # Get the coupling ranges for analysis
        J12 = get_allowed_j12_values(j1, j2)
        J34 = get_allowed_j12_values(j3, j4)
        
        # Build volume matrix using the intersection
        V2, J = build_volume_matrix(j1, j2, j3, j4)
        
        # Case 1: Trivial zero-volume state (no intersection)
        if len(J) == 0:
            trivial_zero_volume_count += 1
            # Check if this satisfies the Diophantine condition
            max_diff = max(abs(j1 - j2), abs(j3 - j4))
            min_sum = min(j1 + j2, j3 + j4)
            is_diophantine = max_diff > min_sum
            
            results_by_kernel_dim['trivial_no_intersection'].append({
                'spins': (j1, j2, j3, j4),
                'J12': J12,
                'J34': J34,
                'intersection': J,
                'diophantine_condition': is_diophantine,
                'max_diff': max_diff,
                'min_sum': min_sum
            })
            continue
        
        matrix_sizes[V2.shape[0]] += 1
        
        # Skip cases where intertwiner space is 1 dimensional (no kernel possible)
        if V2.shape[0] <= 1:
            no_kernel_count += 1
            continue
              # Analyze matrix properties
        props = analyze_matrix_properties(V2)
        
        # Only process singular matrices (zero volume states)
        if props['is_singular']:
            total_states += 1
            non_trivial_zero_volume_count += 1
            
            # Add to appropriate category based on kernel dimension
            kernel_dim = props['kernel_dimension']
            kernel_dim_counts[kernel_dim] += 1
            
            # Check if this state has spin 1/2 on any edge
            has_spin_half = any(j == 0.5 for j in [j1, j2, j3, j4])
            if has_spin_half:
                has_spin_half_count += 1
                
            # Store detailed information for analysis
            result = {
                'spins': (j1, j2, j3, j4),
                'J12': J12,
                'J34': J34,
                'intersection': J,
                'matrix_props': props,
                'has_spin_half': has_spin_half
            }
            
            results_by_kernel_dim[kernel_dim].append(result)
        else:
            no_kernel_count += 1
      # Calculate statistics and print results
    total_analyzed = total_states + trivial_zero_volume_count + no_kernel_count
    print(f"\n=== ZERO-VOLUME STATE ANALYSIS RESULTS ===")
    print(f"Total configurations analyzed: {total_analyzed}")
    print(f"Trivial zero-volume states (J₁₂ ∩ J₃₄ = ∅): {trivial_zero_volume_count}")
    print(f"Non-trivial zero-volume states (kernel exists): {total_states}")  
    print(f"No kernel states (full rank matrices): {no_kernel_count}")
    
    if total_states > 0:
        print(f"Percentage with at least one spin-1/2 edge: {has_spin_half_count/total_states*100:.1f}%")
        
        print("\nKernel dimension statistics (non-trivial cases):")
        for dim, count in sorted(kernel_dim_counts.items()):
            print(f"  Kernel dimension {dim}: {count} states ({count/total_states*100:.1f}%)")
    
    # Analyze the trivial zero-volume states to validate Diophantine condition
    if trivial_zero_volume_count > 0:
        print(f"\nAnalyzing {trivial_zero_volume_count} trivial zero-volume states:")
        diophantine_satisfied = 0
        for case in results_by_kernel_dim['trivial_no_intersection']:
            if case['diophantine_condition']:
                diophantine_satisfied += 1
        print(f"States satisfying Diophantine condition: {diophantine_satisfied}/{trivial_zero_volume_count}")
        print(f"Percentage: {diophantine_satisfied/trivial_zero_volume_count*100:.1f}%")
    
    print("\nMatrix size statistics:")
    for size, count in sorted(matrix_sizes.items()):
        print(f"  {size}x{size} matrices: {count} cases")
      # Analyze some specific interesting cases
    if results_by_kernel_dim:
        # Create detailed catalog
        create_zero_volume_catalog(results_by_kernel_dim)
        
        # Find the largest kernel dimension (excluding trivial cases)
        numeric_dims = [k for k in results_by_kernel_dim.keys() if isinstance(k, int)]
        if numeric_dims:
            max_kernel_dim = max(numeric_dims)
            if max_kernel_dim > 0:
                print(f"\nAnalyzing cases with largest kernel dimension ({max_kernel_dim}):")
                for i, case in enumerate(results_by_kernel_dim[max_kernel_dim][:3]):
                    spins = case['spins']
                    props = case['matrix_props']
                    print(f"\nCase {i+1}: Spins {spins}")
                    print(f"  Matrix size: {props['size']}x{props['size']}")
                    print(f"  Matrix rank: {props['rank']} (kernel dimension: {props['kernel_dimension']})")
                    print(f"  Determinant: {props['determinant']:.10e}")
                print(f"  Eigenvalues: {props['eigenvalues']}")
          # Create visualization
        if kernel_dim_counts:  # Only create if we have non-trivial states
            visualize_kernel_dimensions(kernel_dim_counts)
            visualize_spin_half_correlation(results_by_kernel_dim)
    
    return results_by_kernel_dim

def visualize_kernel_dimensions(kernel_dim_counts):
    """
    Create a visualization of kernel dimension distribution.
    """
    plt.figure(figsize=(10, 6))
    dims = sorted(kernel_dim_counts.keys())
    counts = [kernel_dim_counts[dim] for dim in dims]
    
    plt.bar(dims, counts)
    plt.xlabel('Kernel Dimension')
    plt.ylabel('Number of States')
    plt.title('Distribution of Kernel Dimensions for Zero-Volume States')
    plt.xticks(dims)
    plt.grid(axis='y', alpha=0.3)
    
    # Add percentage labels
    total = sum(counts)
    for i, count in enumerate(counts):
        plt.text(dims[i], count, f'{count/total*100:.1f}%', 
                 ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig('kernel_dimension_distribution.png')
    print("\nSaved visualization to 'kernel_dimension_distribution.png'")

def visualize_spin_half_correlation(results_by_kernel_dim):
    """
    Visualize correlation between spin-1/2 presence and kernel dimension.
    """
    # Only consider numeric kernel dimensions (exclude 'trivial_no_intersection')
    numeric_dims = [k for k in results_by_kernel_dim.keys() if isinstance(k, int)]
    if not numeric_dims:
        print("No non-trivial kernel dimensions to visualize.")
        return
        
    dims = sorted(numeric_dims)
    with_spin_half = []
    without_spin_half = []
    
    for dim in dims:
        cases = results_by_kernel_dim[dim]
        has_half = sum(1 for case in cases if case.get('has_spin_half', False))
        no_half = len(cases) - has_half
        with_spin_half.append(has_half)
        without_spin_half.append(no_half)
    
    plt.figure(figsize=(12, 7))
    
    bar_width = 0.35
    index = np.arange(len(dims))
    
    plt.bar(index, with_spin_half, bar_width, label='With spin-1/2')
    plt.bar(index + bar_width, without_spin_half, bar_width, label='Without spin-1/2')
    
    plt.xlabel('Kernel Dimension')
    plt.ylabel('Number of States')
    plt.title('Presence of Spin-1/2 Edge vs. Kernel Dimension')
    plt.xticks(index + bar_width/2, dims)
    plt.legend()
    
    # Add percentage labels for each group
    for i, (has_half, no_half) in enumerate(zip(with_spin_half, without_spin_half)):
        total = has_half + no_half
        if total > 0:
            plt.text(i, has_half/2, f'{has_half/total*100:.1f}%', 
                    ha='center', va='center')
            plt.text(i + bar_width, no_half/2, f'{no_half/total*100:.1f}%', 
                    ha='center', va='center')
    
    plt.tight_layout()
    plt.savefig('spin_half_correlation.png')
    print("Saved visualization to 'spin_half_correlation.png'")

def detailed_case_study(j1, j2, j3, j4):
    """
    Perform a detailed case study of a specific spin configuration.
    """
    print(f"\n=== DETAILED CASE STUDY: Spins ({j1}, {j2}, {j3}, {j4}) ===")
    
    # Get the coupling ranges for analysis
    J12 = get_allowed_j12_values(j1, j2)
    J34 = get_allowed_j12_values(j3, j4)
    
    print(f"J12 values (j1={j1}, j2={j2}): {J12}")
    print(f"J34 values (j3={j3}, j4={j4}): {J34}")
    
    # Build volume matrix using the corrected intersection
    V2, J = build_volume_matrix(j1, j2, j3, j4)
    
    print(f"Intersection J = J12 ∩ J34: {J}")
    
    if len(J) == 0:
        print("TRIVIAL ZERO-VOLUME STATE: No intersection between coupling ranges")
        
        # Check Diophantine condition
        max_diff = max(abs(j1 - j2), abs(j3 - j4))
        min_sum = min(j1 + j2, j3 + j4)
        is_diophantine = max_diff > min_sum
        
        print(f"Diophantine condition: max(|j₁-j₂|,|j₃-j₄|) > min(j₁+j₂,j₃+j₄)")
        print(f"  max({abs(j1-j2)}, {abs(j3-j4)}) > min({j1+j2}, {j3+j4})")
        print(f"  {max_diff} > {min_sum} → {is_diophantine}")
        
        return None, J, None
    
    # Analyze matrix properties
    props = analyze_matrix_properties(V2)
    
    print(f"Matrix size: {props['size']}x{props['size']}")
    print(f"Matrix rank: {props['rank']} (kernel dimension: {props['kernel_dimension']})")
    print(f"Determinant: {props['determinant']:.10e}")
    print(f"Eigenvalues: {props['eigenvalues']}")
    
    # If the kernel is non-trivial, compute basis vectors
    if props['kernel_dimension'] > 0:
        # Use SVD to find the kernel (null space)
        U, S, Vh = np.linalg.svd(V2)
        null_mask = (abs(S) < 1e-10)
        null_space = Vh[null_mask]
        
        print(f"\nKernel basis vectors ({null_space.shape[0]} vectors):")
        for i, vec in enumerate(null_space):
            print(f"Vector {i+1}: {vec}")
          # Check if the kernel vectors correspond to specific intertwiner combinations
        print("\nIntertwiner interpretation of kernel vectors:")
        for i, vec in enumerate(null_space):
            print(f"Vector {i+1} coefficients correspond to intermediate couplings:")
            for j, coef in enumerate(vec):
                if abs(coef) > 1e-6:  # Only print significant contributions
                    print(f"  J = {J[j]}: {coef:.6f}")
    
    # Print the full matrix for inspection
    print("\nFull matrix:")
    np.set_printoptions(precision=6, suppress=True, linewidth=120)
    print(V2)
    
    return V2, J, props

def create_zero_volume_catalog(results_by_kernel_dim):
    """
    Create a detailed catalog of all zero-volume states found.
    """
    print("\n" + "="*60)
    print("DETAILED ZERO-VOLUME STATE CATALOG")
    print("="*60)
    
    # First, catalog trivial zero-volume states
    if 'trivial_no_intersection' in results_by_kernel_dim:
        trivial_cases = results_by_kernel_dim['trivial_no_intersection']
        print(f"\n1. TRIVIAL ZERO-VOLUME STATES (J₁₂ ∩ J₃₄ = ∅): {len(trivial_cases)} cases")
        print("-" * 60)
        
        for i, case in enumerate(trivial_cases[:10]):  # Show first 10
            spins = case['spins']
            J12 = case['J12']
            J34 = case['J34']
            is_diophantine = case['diophantine_condition']
            max_diff = case['max_diff']
            min_sum = case['min_sum']
            
            print(f"\nCase {i+1}: j₁={spins[0]}, j₂={spins[1]}, j₃={spins[2]}, j₄={spins[3]}")
            print(f"  J₁₂ = {J12}")
            print(f"  J₃₄ = {J34}")
            print(f"  Intersection: ∅")
            print(f"  Diophantine condition: max(|j₁-j₂|,|j₃-j₄|) > min(j₁+j₂,j₃+j₄)")
            print(f"    max({abs(spins[0]-spins[1])}, {abs(spins[2]-spins[3])}) > min({spins[0]+spins[1]}, {spins[2]+spins[3]})")
            print(f"    {max_diff} > {min_sum} → {is_diophantine}")
            
        if len(trivial_cases) > 10:
            print(f"\n... and {len(trivial_cases) - 10} more trivial cases")
      # Then, catalog any non-trivial zero-volume states
    non_trivial_found = False
    for kernel_dim in sorted(k for k in results_by_kernel_dim.keys() if isinstance(k, int)):
        if kernel_dim > 0:
            non_trivial_found = True
            cases = results_by_kernel_dim[kernel_dim]
            print(f"\n2. NON-TRIVIAL ZERO-VOLUME STATES (kernel dimension {kernel_dim}): {len(cases)} cases")
            print("-" * 60)
            
            for i, case in enumerate(cases):
                spins = case['spins']
                J12 = case['J12']
                J34 = case['J34']
                J = case['intersection']
                props = case['matrix_props']
                
                print(f"\nCase {i+1}: j₁={spins[0]}, j₂={spins[1]}, j₃={spins[2]}, j₄={spins[3]}")
                print(f"  J₁₂ = {J12}")
                print(f"  J₃₄ = {J34}")
                print(f"  Intersection J = {J}")
                print(f"  Matrix size: {props['size']}×{props['size']}")
                print(f"  Rank: {props['rank']}, Kernel dimension: {props['kernel_dimension']}")
                print(f"  Determinant: {props['determinant']:.2e}")
                print(f"  Min eigenvalue: {props['min_eigenvalue']:.2e}")
    
    if not non_trivial_found:
        print(f"\n2. NON-TRIVIAL ZERO-VOLUME STATES: NONE FOUND! ✓")
        print("-" * 60)
        print("This confirms the Diophantine root catalog assertion:")
        print("There are NO non-trivial 4-valent kernel states in the spin range 0.5-3.0")

if __name__ == "__main__":
    # Analyze all zero-volume states
    results = analyze_zero_volume_states()
    
    # Perform detailed case studies of interesting configurations
    # Case 1: A high-dimensional matrix with large kernel
    detailed_case_study(2.5, 3.0, 0.5, 0.5)
    
    # Case 2: A state without any spin-1/2
    detailed_case_study(3.0, 3.0, 3.0, 3.0)
    
    # Case 3: A state with typical behavior
    detailed_case_study(1.5, 1.5, 1.5, 1.5)
    
    # Create a catalog of all zero-volume states found
    create_zero_volume_catalog(results)
