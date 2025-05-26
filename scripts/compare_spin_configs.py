import sympy as sp
import pandas as pd
import numpy as np
import os
from sympy import Rational, factorial, sqrt, prod
import matplotlib.pyplot as plt
from itertools import product

def CF12j_num(j1, j2, j12, j3, j4, j23, j5, j6, j34, j7, j8, j45):
    # convert to rational internally for exact prefactor
    a, b, c = map(Rational, (j1, j2, j12))
    d, e, f = map(Rational, (j3, j4, j23))
    g, h, i = map(Rational, (j5, j6, j34))
    k, l, m = map(Rational, (j7, j8, j45))
    triples = [(a, b, c), (d, e, f), (g, h, i), (k, l, m)]
    
    # Check if any triple doesn't satisfy the triangle inequality
    for (x, y, z) in triples:
        if x + y < z or x + z < y or y + z < x:
            return 0.0
            
    Delta = sqrt(prod(
        factorial(-x + y + z) * factorial(x - y + z) * factorial(x + y - z)
        / factorial(x + y + z + 1)
        for (x, y, z) in triples
    ))
    num = [sp.Rational(1, 2), -c, c + 1, -f, f + 1]
    den = [a + b - c + 1, d + e - f + 1, g + h - i + 1, k + l - m + 1]
    
    try:
        val = Delta * sp.N(sp.hyper(num, den, 1))
        # Handle potential complex values
        if val.is_real:
            return float(val)
        else:
            print(f"Warning: Complex value encountered: {val}")
            return float(abs(val))
    except Exception as e:
        print(f"Warning: Error computing CF12j: {e}")
        return 0.0

def get_allowed_j12_values(j1, j2):
    """Get allowed intermediate coupling values for given j1, j2"""
    j_min = abs(j1 - j2)
    j_max = j1 + j2
    # Generate all half-integer values between j_min and j_max
    step = 0.5 if (j_min % 1 != 0 or j_max % 1 != 0) else 1.0
    return np.arange(j_min, j_max + step, step).tolist()

def compute_volume_spectrum(j1, j2, j3, j4):
    """Compute the volume spectrum for given spins"""
    # Get allowed intermediate coupling values
    J12 = get_allowed_j12_values(j1, j2)
    J34 = get_allowed_j12_values(j3, j4)
    
    print(f"  Intermediate couplings for (j1,j2)=({j1},{j2}): {J12}")
    print(f"  Intermediate couplings for (j3,j4)=({j3},{j4}): {J34}")
    
    # Get the unique allowed intermediate couplings for matrix size
    all_j_intermediate = sorted(set(J12 + J34))
    n = len(all_j_intermediate)
    
    print(f"  Matrix size: {n}x{n}")
    
    # Build numeric V2
    V2_num = np.zeros((n, n))
      # Map from intermediate coupling values to matrix indices
    index_map = {j: i for i, j in enumerate(all_j_intermediate)}
    
    # Fill the matrix
    for ja in J12:
        for jb in J34:
            a = index_map[ja]
            b = index_map[jb]
            # Corrected spin mirroring with proper intermediate spin ordering for 12j symbol
            # Ordering: (j1,j2,j12), (j3,j4,j12), (j1,j2,j34), (j3,j4,j34)
            val = CF12j_num(j1, j2, ja, j3, j4, ja, j1, j2, jb, j3, j4, jb)
            V2_num[a, b] = val
            # Mirror for symmetry
            V2_num[b, a] = val
    
    # Print the matrix for debugging
    print(f"  V2 matrix:\n{V2_num}")
    
    # Compute eigenvalues
    eigenvalues = np.linalg.eigvals(V2_num)
    
    # Filter out eigenvalues that are effectively zero
    significant_eigenvalues = np.array([ev for ev in eigenvalues if abs(ev) > 1e-10])
    
    # Take absolute value before sqrt to handle any potential numerical issues
    result = sorted(np.sqrt(abs(significant_eigenvalues)), reverse=True)
    print(f"  Found {len(result)} significant eigenvalues")
    
    return result

def compute_multiple_spectra():
    """Compute spectra for different spin configurations"""
    # Define spin configurations to test
    spin_configs = [
        (0.5, 0.5, 0.5, 0.5, "j₁=j₂=j₃=j₄=1/2"),
        (1.0, 1.0, 0.5, 1.5, "j₁=j₂=1, j₃=1/2, j₄=3/2"),
        (1.0, 1.0, 1.0, 1.0, "j₁=j₂=j₃=j₄=1"),
        (1.5, 1.5, 1.5, 1.5, "j₁=j₂=j₃=j₄=3/2"),
        (0.5, 1.0, 1.5, 2.0, "j₁=1/2, j₂=1, j₃=3/2, j₄=2"),
    ]
    
    # Create results dataframe
    results = []
    
    for j1, j2, j3, j4, label in spin_configs:
        print(f"Computing spectrum for {label}...")
        spectrum = compute_volume_spectrum(j1, j2, j3, j4)
        
        # Create a flat data structure for the results
        for i, val in enumerate(spectrum):
            results.append({
                'configuration': label,
                'j1': j1,
                'j2': j2,
                'j3': j3,
                'j4': j4,
                'eigenvalue_index': i,
                'eigenvalue': val
            })
        
        print(f"Found {len(spectrum)} non-zero eigenvalues")
    
    # Create the dataframe and save to CSV
    os.makedirs('data', exist_ok=True)
    df = pd.DataFrame(results)
    csv_path = 'data/volume_spectra_comparison.csv'
    df.to_csv(csv_path, index=False)
    
    print(f"Spectra for different spin configurations saved to {csv_path}")
    
    # Create a table view that's easy to read
    table_df = df.pivot_table(
        index='configuration', 
        columns='eigenvalue_index',
        values='eigenvalue',
        aggfunc='first'
    )
    
    # Save the table to CSV as well
    table_csv_path = 'data/volume_spectra_table.csv'
    table_df.to_csv(table_csv_path)
    print(f"Formatted table saved to {table_csv_path}")
    
    # Visualization
    plt.figure(figsize=(12, 8))
    
    for config_name in df['configuration'].unique():
        config_data = df[df['configuration'] == config_name]
        plt.plot(
            config_data['eigenvalue_index'], 
            config_data['eigenvalue'], 
            'o-', 
            label=config_name
        )
    
    plt.xlabel('Eigenvalue Index')
    plt.ylabel('Volume Eigenvalue')
    plt.title('Volume Spectra for Different Spin Configurations')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('data/volume_spectra_comparison.png')
    print("Plot saved to data/volume_spectra_comparison.png")

if __name__ == "__main__":
    compute_multiple_spectra()
