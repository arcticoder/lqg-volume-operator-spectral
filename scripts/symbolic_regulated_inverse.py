"""
Symbolic calculation of the regulated inverse volume operator in LQG.
This script performs symbolic analysis of the volume operator and its inverse.
"""

import sympy as sp
from sympy import Matrix, simplify, latex, symbols, Rational, Integer, zeros
from pathlib import Path
from lqg_utils import CF12j_symbolic, get_allowed_j12_values, create_symbolic_spins

# Define epsilon regulator
epsilon = symbols('epsilon', positive=True)

def compute_regulated_inverse_operator(j_values=None, symbolic=False, numeric_check=False):
    """
    Compute the regulated inverse operator V_inv^(epsilon) and
    verify its boundedness and self-adjointness properties.
    
    Args:
        j_values: Tuple of (j1, j2, j3, j4) spin values, or None to use symbolic variables
        symbolic: If True, keep j_values symbolic, otherwise convert to numeric
        numeric_check: If True, attempt to evaluate symbolic expressions numerically
    
    Returns:
        Dictionary containing the analysis results
    """
    # If no j_values are provided or symbolic=True, use symbolic variables
    if j_values is None or symbolic:
        j_values = create_symbolic_spins()
        symbolic_mode = True
    else:
        j_values = tuple(map(Rational, j_values))  # Convert to Rational for exact calculations
        symbolic_mode = False
        
    # Disable numeric_check for fully symbolic expressions
    if symbolic_mode:
        numeric_check = False
        
    # Unpack spin values
    j1_val, j2_val, j3_val, j4_val = j_values
    
    # Create readable label for this configuration
    if symbolic_mode:
        config_label = "symbolic"
    else:
        config_label = f"j1={float(j1_val)},j2={float(j2_val)},j3={float(j3_val)},j4={float(j4_val)}"
    print(f"Computing regulated inverse for {config_label} configuration")
    
    # Define the allowed intermediate couplings
    J12 = get_allowed_j12_values(j1_val, j2_val)
    J34 = get_allowed_j12_values(j3_val, j4_val)
    
    # Get the unique allowed intermediate couplings for matrix size
    if symbolic_mode:
        all_j_intermediate = J12 + J34
        # For symbolic mode, we can't sort, just eliminate duplicates
        all_j_intermediate = list(dict.fromkeys(all_j_intermediate))  # preserve order but remove duplicates
    else:
        # Convert to Rational for exact arithmetic
        all_j_intermediate = sorted(set([j for j in J12 + J34]))
    
    n = len(all_j_intermediate)
    
    # Create index mapping
    index_map = {j: i for i, j in enumerate(all_j_intermediate)}
    
    V2_symbolic = zeros(n, n)
    for ja in J12:
        for jb in J34:
            a = index_map[ja]
            b = index_map[jb]
            # Calculate the matrix element symbolically
            val = CF12j_symbolic(j1_val, j2_val, ja, j3_val, j4_val, ja, 
                                j1_val, j2_val, jb, j3_val, j4_val, jb)
            V2_symbolic[a, b] = val
            # Mirror for symmetry
            V2_symbolic[b, a] = val
    
    print("Symbolic V² matrix:")
    print(V2_symbolic)
    
    # Convert to Matrix for eigenvalue calculation
    V2_matrix = Matrix(V2_symbolic)
    
    try:
        # Compute eigenvalues and eigenvectors
        if symbolic_mode:
            # For symbolic mode, use symbolic computation
            eigenvalues = V2_matrix.eigenvals()
            eigenvectors = V2_matrix.eigenvects()
            
            print("\nSymbolic eigenvalues:")
            for val, mult in eigenvalues.items():
                # For symbolic expressions, avoid expansion which can cause problems
                print(f"lambda = {val} (multiplicity: {mult})")
        else:
            # For numeric mode with exact Rational values
            # First convert to a numeric matrix for numerical stability
            V2_numeric = V2_matrix.evalf()
            eigenvalues_numeric = V2_numeric.eigenvals()
            eigenvectors_numeric = V2_numeric.eigenvects()
            
            # Convert back to dictionary format
            eigenvalues = eigenvalues_numeric
            eigenvectors = eigenvectors_numeric
            
            print("\nNumeric eigenvalues:")
            for val, mult in eigenvalues.items():
                print(f"lambda = {val.evalf()} (multiplicity: {mult})")
    except Exception as e:
        print(f"Error computing eigenvalues: {e}")
        print("Falling back to numerical approximation")
        
        # Convert to float matrix and use numpy for eigenvalues
        import numpy as np
        try:
            # First try with SymPy's N function for better conversion
            V2_float = np.array(V2_matrix.evalf()).astype(float)
        except Exception:
            # If that fails, try direct conversion (might lose some precision)
            print("Warning: Direct numerical conversion may lose precision")
            V2_float = np.array([[float(V2_matrix[i, j].evalf()) if not V2_matrix[i, j].has(sp.Symbol) else 0.0 
                               for j in range(V2_matrix.cols)] 
                              for i in range(V2_matrix.rows)])
        
        evals, evecs = np.linalg.eigh(V2_float)
        
        # Convert back to SymPy format
        eigenvalues = {sp.Float(e): 1 for e in evals}
        eigenvectors = [(sp.Float(evals[i]), 1, [Matrix(evecs[:, i])]) for i in range(len(evals))]
        
        print("\nNumerical eigenvalues (numpy):")
        for val, mult in eigenvalues.items():
            print(f"lambda = {val} (multiplicity: {mult})")
    
    # Construct the projectors and check for any zero eigenvalues
    projectors = []
    zero_eigenvalues = []
    nonzero_eigenvalues = []
    
    print("\nConstructing projectors and checking for zero eigenvalues...")
    for val, mult, basis in eigenvectors:
        # Check if eigenvalue is close to zero (handle symbolic case separately)
        is_zero = False
        if symbolic_mode:
            # For symbolic mode, we can't numerically check if it's zero
            # We rely on SymPy's symbolic zero detection
            if val.is_zero or val == 0:
                is_zero = True
            # Don't try numerical evaluation for symbolic expressions
        else:
            # For numeric values, check if they're close to zero
            try:
                if abs(complex(val.evalf())) < 1e-10:
                    is_zero = True
            except Exception as e:
                # If conversion to complex fails, assume non-zero
                print(f"  Note: Could not check if eigenvalue is zero numerically: {e}")
                is_zero = False
        
        if is_zero:
            print(f"Zero eigenvalue found: {val}")
            zero_eigenvalues.append((val, basis))
        else:
            nonzero_eigenvalues.append((val, basis))
        
        # For each eigenvector, construct the projector
        for vec in basis:
            # Normalize the eigenvector
            norm = sp.sqrt(sum(v**2 for v in vec))
            normalized_vec = vec / norm
            # Calculate the projector P = |v⟩⟨v|
            projector = Matrix(normalized_vec) * Matrix(normalized_vec).transpose()
            projectors.append((val, projector))
    
    # Define the kernel projector (if any zero eigenvalues exist)
    P_kernel = zeros(n, n)
    if zero_eigenvalues:
        print("Kernel projector found.")
        for val, basis in zero_eigenvalues:
            for vec in basis:
                # Normalize the eigenvector
                norm = sp.sqrt(sum(v**2 for v in vec))
                normalized_vec = vec / norm
                # Add to the kernel projector
                P_kernel += Matrix(normalized_vec) * Matrix(normalized_vec).transpose()
    else:
        print("No kernel projector (all eigenvalues are non-zero).")
    
    # Define the regulated inverse operator
    V_inv_reg = zeros(n, n)
    
    # Add the non-zero eigenvalue contributions
    for val, basis in nonzero_eigenvalues:
        for vec in basis:
            # Normalize the eigenvector
            norm = sp.sqrt(sum(v**2 for v in vec))
            normalized_vec = vec / norm
            # Add (1/lambda) * P to the inverse
            projector = Matrix(normalized_vec) * Matrix(normalized_vec).transpose()
            V_inv_reg += (1/val) * projector
    
    # Add the kernel contribution with epsilon regulator
    if not P_kernel.is_zero_matrix:
        V_inv_reg += (1/epsilon) * P_kernel
    
    # Verify boundedness: check that all 1/lambda_k are finite for non-zero eigenvalues
    print("\nVerifying boundedness:")
    is_bounded = True
    
    for val, _ in nonzero_eigenvalues:
        inv_val = 1/val
        print(f"1/lambda = {inv_val}")
        
        # Only attempt numerical evaluation if requested and not in fully symbolic mode
        if numeric_check and not symbolic_mode:
            try:
                num_val = complex(val.evalf())
                num_inv_val = 1.0/num_val
                print(f"  λ = {num_val:.10f}")
                print(f"  1/λ = {num_inv_val:.10f}")
                
                if abs(num_val) < 1e-10:  # Close to zero
                    is_bounded = False
                    print(f"  Warning: Eigenvalue too close to zero: {num_val}")
            except Exception as e:
                print(f"  Could not evaluate numerically: {e}")
                # For partially symbolic expressions, we'll attempt to substitute values
                if not symbolic_mode:
                    print("  Attempting to evaluate with simplified expressions...")
                    try:
                        # Try with some simplifications first
                        simpler_val = val.simplify()
                        # Avoid expansion of hypergeometric functions
                        if 'hyper' in str(simpler_val):
                            print("  Contains hypergeometric functions, skipping numerical check")
                        else:
                            num_val = complex(simpler_val.evalf())
                            print(f"  Simplified λ = {num_val:.10f}")
                    except Exception as e2:
                        print(f"  Still could not evaluate: {e2}")
        elif symbolic_mode:
            print("  Symbolic mode: Skipping numerical evaluation")
        
        # Symbolic boundedness check
        if inv_val.has(sp.zoo) or inv_val.has(sp.oo) or inv_val.has(-sp.oo):
            is_bounded = False
            print(f"  Warning: Unbounded term 1/{val} = {inv_val}")
    
    if is_bounded:
        print("  The operator is bounded (all 1/lambda_k are finite).")
    else:
        print("  The operator is unbounded (some 1/lambda_k are infinite).")
        
    # Verify self-adjointness
    print("\nVerifying self-adjointness:")
    # For a matrix to be self-adjoint (Hermitian), it should equal its conjugate transpose
    V_inv_adj = V_inv_reg.conjugate().transpose()
    diff = simplify(V_inv_reg - V_inv_adj)
    is_self_adjoint = diff.is_zero_matrix
    
    if is_self_adjoint:
        print("  The operator is self-adjoint (V_inv = V_inv†).")
    else:
        print("  The operator is not self-adjoint.")
        print("  V_inv - V_inv† = ")
        print(diff)
    
    # Display the explicit form of the regulated inverse operator
    print("\nRegulated Inverse Operator V_inv^(epsilon):")
    print(simplify(V_inv_reg))
    
    # LaTeX representation for the regulated inverse operator
    latex_V_inv = latex(simplify(V_inv_reg))
    print("\nLaTeX form of V_inv^(epsilon):")
    print(latex_V_inv)
    
    # Create a spectral representation expression
    spectral_form = ""
    for val, _ in nonzero_eigenvalues:
        spectral_form += f"\\frac{{1}}{{{latex(val)}}}P_{{\\lambda_{{{latex(val)}}}}} + "
    
    # Add the kernel term if present
    if not P_kernel.is_zero_matrix:
        spectral_form += f"\\frac{{1}}{{\\epsilon}}P_{{\\ker V}}"
    else:
        # Remove the trailing " + "
        spectral_form = spectral_form[:-3]
    
    print("\nSpectral form of V_inv^(epsilon):")
    print(f"V_inv^(epsilon) = {spectral_form}")
    
    # Create a simpler form with numerical indices for the projectors
    simple_spectral_form = ""
    for i, (val, _) in enumerate(nonzero_eigenvalues):
        simple_spectral_form += f"\\frac{{1}}{{{latex(val)}}}P_{{{i+1}}} + "
    
    # Add the kernel term if present
    if not P_kernel.is_zero_matrix:
        simple_spectral_form += f"\\frac{{1}}{{\\epsilon}}P_{{0}}"
    else:
        # Remove the trailing " + "
        simple_spectral_form = simple_spectral_form[:-3]
    
    print("\nSimplified spectral form of V_inv^(epsilon):")
    print(f"V_inv^(epsilon) = {simple_spectral_form}")
    
    # Collect results in a dictionary
    results = {
        'config_label': config_label,
        'j_values': j_values,
        'V2_matrix': V2_symbolic,
        'eigenvalues': eigenvalues,
        'is_bounded': is_bounded,
        'is_self_adjoint': is_self_adjoint,
        'V_inv_reg': V_inv_reg,
        'latex_V_inv': latex_V_inv,
        'spectral_form': spectral_form,
        'simple_spectral_form': simple_spectral_form,
        'nonzero_eigenvalues': nonzero_eigenvalues,
    }
    
    # Save results to file
    data_dir = Path('data')
    data_dir.mkdir(exist_ok=True)
    
    # Generate filename based on configuration
    if symbolic_mode:
        output_file = data_dir / 'regulated_inverse_operator_symbolic.txt'
    else:
        spins_str = f"j{float(j1_val)}-{float(j2_val)}-{float(j3_val)}-{float(j4_val)}".replace(".", "_")
        output_file = data_dir / f'regulated_inverse_operator_{spins_str}.txt'
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(f"Regulated Inverse Volume Operator Analysis for {config_label}\n")
        f.write("=======================================\n\n")
        f.write("V² Matrix:\n")
        f.write(str(V2_symbolic) + "\n\n")
        
        f.write("Eigenvalues of V²:\n")
        for val, mult in eigenvalues.items():
            f.write(f"lambda = {val} (multiplicity: {mult})\n")
        f.write("\n")
        
        f.write("Boundedness Analysis:\n")
        for val, _ in nonzero_eigenvalues:
            inv_val = 1/val
            f.write(f"1/lambda = {inv_val}\n")
        if is_bounded:            
            f.write("Conclusion: The operator is bounded (all 1/lambda_k are finite).\n\n")
        else:
            f.write("Conclusion: The operator is unbounded (some 1/lambda_k are infinite).\n\n")
        
        f.write("Self-adjointness Analysis:\n")
        if is_self_adjoint:
            f.write("Conclusion: The operator is self-adjoint (V_inv = V_inv†).\n\n")
        else:
            f.write("Conclusion: The operator is not self-adjoint.\n")
            f.write(f"V_inv - V_inv† = {diff}\n\n")
        
        f.write("Regulated Inverse Operator V_inv^(epsilon):\n")
        f.write(str(simplify(V_inv_reg)) + "\n\n")
        
        f.write("LaTeX form of V_inv^(epsilon):\n")
        f.write(latex_V_inv + "\n\n")
        
        f.write("Spectral form of V_inv^(epsilon):\n")
        f.write(f"V_inv^(epsilon) = {spectral_form}\n\n")
        
        f.write("Simplified spectral form of V_inv^(epsilon):\n")
        f.write(f"V_inv^(epsilon) = {simple_spectral_form}\n")
    
    print(f"Results saved to {output_file}")
    
    return results

if __name__ == "__main__":
    # Analyze with symbolic parameters - explicitly set numeric_check=False for symbolic mode
    compute_regulated_inverse_operator(symbolic=True, numeric_check=False)
    print("Regulated inverse operator analysis with symbolic parameters complete.")
    
    # Test with specific numeric configurations
    spin_configs = [
        (0.5, 0.5, 0.5, 0.5, "j₁=j₂=j₃=j₄=1/2"),
        # Uncomment to test other configurations
        # (1.0, 1.0, 1.0, 1.0, "j₁=j₂=j₃=j₄=1"),
        # (1.5, 1.5, 1.5, 1.5, "j₁=j₂=j₃=j₄=3/2"),
        # (0.5, 1.0, 1.5, 2.0, "j₁=1/2, j₂=1, j₃=3/2, j₄=2"),
    ]
    
    for j1_val, j2_val, j3_val, j4_val, label in spin_configs:
        print(f"\nAnalyzing configuration: {label}")
        results = compute_regulated_inverse_operator((j1_val, j2_val, j3_val, j4_val), numeric_check=True)
        
        # For the j=1/2 configuration, try explicit substitutions for J_{1/21/2}
        if j1_val == 0.5 and j2_val == 0.5 and j3_val == 0.5 and j4_val == 0.5:
            print("\nAttempting explicit evaluation with J_{1/21/2} = 0:")
            try:
                # Try substituting J_{1/21/2} = 0
                J_sym = sp.symbols('J_{1/21/2}')
                eigen_expr = results['eigenvalues']
                for val in eigen_expr:
                    numeric_val = val.subs(J_sym, 0).evalf()
                    print(f"Eigenvalue = {numeric_val}")
                    inv_val = 1 / numeric_val
                    print(f"1/Eigenvalue = {inv_val}")
            except Exception as e:
                print(f"Error with substitution: {e}")
                
            print("\nAttempting explicit evaluation with J_{1/21/2} = 1:")
            try:
                # Try substituting J_{1/21/2} = 1
                numeric_val = val.subs(J_sym, 1).evalf()
                print(f"Eigenvalue = {numeric_val}")
                inv_val = 1 / numeric_val
                print(f"1/Eigenvalue = {inv_val}")
            except Exception as e:
                print(f"Error with substitution: {e}")
    
    print("Regulated inverse operator analysis complete.")
