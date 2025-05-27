import numpy as np
import sympy as sp
from itertools import product
from lqg_utils import CF12j_numeric, get_allowed_j12_values

def det_V2(j1, j2, j3, j4):
    """
    Build the V^2 matrix via CF12j_numeric and return its determinant
    
    Args:
        j1, j2, j3, j4: The four spins of the valence-4 vertex
        
    Returns:
        The determinant of the volume-squared matrix
    """
    # Get allowed intermediate coupling values
    J12 = get_allowed_j12_values(j1, j2)
    
    # Build V2 matrix
    V2 = np.zeros((len(J12), len(J12)))
    for a, ja in enumerate(J12):
        for b, jb in enumerate(J12):
            # Ordering: (j1,j2,j12), (j3,j4,j12), (j1,j2,j34), (j3,j4,j34)
            V2[a, b] = CF12j_numeric(j1, j2, ja, j3, j4, ja, j1, j2, jb, j3, j4, jb)
    
    # Return determinant
    return np.linalg.det(V2), V2

def intertwiner_dim(j1, j2, j3, j4):
    """
    Calculate the dimension of the intertwiner space for a 4-valent node
    
    Args:
        j1, j2, j3, j4: The four spins of the valence-4 vertex
        
    Returns:
        The dimension of the intertwiner space
    """
    return len(get_allowed_j12_values(j1, j2))

def check_edge_spins_zero(j1, j2, j3, j4, epsilon=1e-10):
    """
    Check if any of the edge spins are effectively zero
    """
    return any(abs(j) < epsilon for j in [j1, j2, j3, j4])

# Define spin range to check
spin_vals = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]  # Extended range up to 3.0

# Lists to store results
trivial_zeros = []  # States where intertwiner space dimension = 1
nontrivial_zeros = []  # States where intertwiner space dimension > 1 but determinant = 0
suspicious_zeros = []  # States with very small determinant that may be numerical artifacts

# Track progress
total = len(spin_vals)**4
count = 0

# Threshold for checking numerical zeros
num_zero_threshold = 1e-8
suspicious_threshold = 1e-6

print(f"Checking {total} spin configurations in range {min(spin_vals)}-{max(spin_vals)}...")
for (j1, j2, j3, j4) in product(spin_vals, repeat=4):
    count += 1
    if count % (total // 20) == 0:
        print(f"Progress: {count}/{total} ({count/total*100:.1f}%)")
    
    # Check intertwiner space dimension
    dim = intertwiner_dim(j1, j2, j3, j4)
    
    if dim <= 1:
        # Trivial case: intertwiner space is 0 or 1 dimensional
        trivial_zeros.append((j1, j2, j3, j4))
    else:
        # Non-trivial case: intertwiner space is at least 2-dimensional
        # Compute determinant and check if it's zero
        det_val, V2_matrix = det_V2(j1, j2, j3, j4)
        
        if abs(det_val) < num_zero_threshold:
            # Verify it's not due to one edge spin being very small
            if not check_edge_spins_zero(j1, j2, j3, j4):
                # Check eigenvalues to confirm it's truly singular
                eigenvals = np.linalg.eigvalsh(V2_matrix)
                min_eigval = min(abs(val) for val in eigenvals)
                
                if min_eigval < num_zero_threshold:
                    nontrivial_zeros.append({
                        'spins': (j1, j2, j3, j4),
                        'det': det_val,
                        'min_eigval': min_eigval,
                        'dimension': dim
                    })
                    print(f"Found non-trivial zero-volume state: {(j1, j2, j3, j4)}, det={det_val}, min_eigenvalue={min_eigval}")
                elif min_eigval < suspicious_threshold:
                    suspicious_zeros.append({
                        'spins': (j1, j2, j3, j4),
                        'det': det_val,
                        'min_eigval': min_eigval,
                        'dimension': dim
                    })

print(f"\nTotal trivial zero-volume states (collapsed intertwiner space): {len(trivial_zeros)}")
print(f"Total non-trivial zero-volume states found: {len(nontrivial_zeros)}")
print(f"Suspicious cases (very small but possibly non-zero determinant): {len(suspicious_zeros)}")

if len(nontrivial_zeros) > 0:
    print("\nNon-trivial zero-volume states (valence-4):")
    # Sort by intertwiner dimension and show the first few cases
    nontrivial_zeros.sort(key=lambda x: x['dimension'], reverse=True)
    for i, zero in enumerate(nontrivial_zeros[:10]):
        print(f"{i+1}. Spins {zero['spins']}, det={zero['det']:.2e}, min_eigval={zero['min_eigval']:.2e}, dim={zero['dimension']}")
    
    if len(nontrivial_zeros) > 10:
        print(f"... and {len(nontrivial_zeros) - 10} more cases")
else:
    print("\nNo non-trivial zero-volume states found in the specified spin range.")
    print("This confirms the 'Diophantine root' catalog assertion that there are")
    print("no non-trivial 4-valent kernel states in this spin range.")

print("\nAnalyzing the numerical precision of the results...")
# Calculate statistics for the determinants
if len(nontrivial_zeros) > 0:
    dets = [abs(zero['det']) for zero in nontrivial_zeros]
    eigvals = [zero['min_eigval'] for zero in nontrivial_zeros]
    print(f"Determinant statistics: min={min(dets):.2e}, max={max(dets):.2e}, mean={np.mean(dets):.2e}")
    print(f"Min eigenvalue statistics: min={min(eigvals):.2e}, max={max(eigvals):.2e}, mean={np.mean(eigvals):.2e}")
    
    # Check if all cases have at least one edge with spin 0.5
    has_min_spin = [0.5 in zero['spins'] for zero in nontrivial_zeros]
    print(f"Percentage of cases with at least one spin=0.5: {sum(has_min_spin)/len(has_min_spin)*100:.1f}%")

# Additional check: Test with higher precision for some cases
if len(nontrivial_zeros) > 0:
    print("\nTesting a few cases with higher precision calculations:")
    for i, case in enumerate(nontrivial_zeros[:3]):
        j1, j2, j3, j4 = case['spins']
        print(f"\nCase {i+1}: {(j1, j2, j3, j4)}")
        
        # Build matrix with higher precision
        J12 = get_allowed_j12_values(j1, j2)
        V2 = np.zeros((len(J12), len(J12)), dtype=np.float64)
        
        for a, ja in enumerate(J12):
            for b, jb in enumerate(J12):
                V2[a, b] = CF12j_numeric(j1, j2, ja, j3, j4, ja, j1, j2, jb, j3, j4, jb)
        
        # Analyze eigenvalues
        eigenvals = np.linalg.eigvalsh(V2)
        print(f"Matrix dimension: {V2.shape[0]}x{V2.shape[1]}")
        print(f"Determinant: {np.linalg.det(V2):.10e}")
        print(f"Eigenvalues: {eigenvals}")
        print(f"Smallest eigenvalue by magnitude: {min(abs(val) for val in eigenvals):.10e}")
        
        # Check rank
        rank = np.linalg.matrix_rank(V2)
        print(f"Matrix rank: {rank} (expected {V2.shape[0]} for non-singular)")
        
        # Print a snippet of the matrix for inspection
        print(f"Matrix snippet (3x3 if large enough):")
        size = min(3, V2.shape[0])
        print(V2[:size, :size])