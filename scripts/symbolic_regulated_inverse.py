import sympy as sp
import numpy as np
import os
from sympy import Rational, factorial, sqrt, prod, Matrix, simplify, latex, symbols

# Define symbolic spin variables and epsilon regulator
j1, j2, j3, j4, epsilon = symbols('j1 j2 j3 j4 epsilon', positive=True)

def CF12j_symbolic(j1, j2, j12, j3, j4, j23, j5, j6, j34, j7, j8, j45):
    """
    Symbolic implementation of the 12j coefficients for the volume operator.
    Returns the expression in symbolic form.
    """
    # convert to rational for exact computations
    a, b, c = j1, j2, j12
    d, e, f = j3, j4, j23
    g, h, i = j5, j6, j34
    k, l, m = j7, j8, j45
    triples = [(a, b, c), (d, e, f), (g, h, i), (k, l, m)]
    
    # Define the symbolic Delta factor
    Delta_expr = sp.sqrt(prod(
        sp.factorial(-x + y + z) * sp.factorial(x - y + z) * sp.factorial(x + y - z)
        / sp.factorial(x + y + z + 1)
        for (x, y, z) in triples
    ))
    
    # Define hypergeometric function parameters
    num = [sp.Rational(1, 2), -c, c + 1, -f, f + 1]
    den = [a + b - c + 1, d + e - f + 1, g + h - i + 1, k + l - m + 1]
    
    # Return the symbolic expression
    return Delta_expr * sp.hyper(num, den, 1)

def compute_regulated_inverse_operator():
    """
    Compute the regulated inverse operator V_inv^(epsilon) for j=1/2 case
    and verify its boundedness and self-adjointness properties.
    """
    # For symbolic analysis, we use j=1/2 case
    j_val = Rational(1, 2)
    
    # Define the allowed intermediate couplings
    # For j1=j2=j3=j4=1/2, J12 can be 0 or 1
    J12 = [sp.Integer(0), sp.Integer(1)]
    
    # Build symbolic V2 matrix
    V2_symbolic = sp.zeros(2, 2)
    for a, ja in enumerate(J12):
        for b, jb in enumerate(J12):
            # Calculate the matrix element symbolically
            val = CF12j_symbolic(j_val, j_val, ja, j_val, j_val, ja, j_val, j_val, jb, j_val, j_val, jb)
            V2_symbolic[a, b] = val
    
    print("Symbolic V² matrix:")
    print(V2_symbolic)
    
    # Convert to Matrix for eigenvalue calculation
    V2_matrix = Matrix(V2_symbolic)
    
    # Compute eigenvalues and eigenvectors symbolically
    eigenvalues = V2_matrix.eigenvals()
    eigenvectors = V2_matrix.eigenvects()
    
    print("\nSymbolic eigenvalues:")
    for val, mult in eigenvalues.items():
        # Simplify the expression
        simplified_val = simplify(val)
        print(f"lambda = {simplified_val} (multiplicity: {mult})")
        
    # Construct the projectors and check for any zero eigenvalues
    projectors = []
    zero_eigenvalues = []
    nonzero_eigenvalues = []
    
    print("\nConstructing projectors and checking for zero eigenvalues...")
    for val, mult, basis in eigenvectors:
        simplified_val = simplify(val)
        # Store eigenvalue and check if it's zero
        if simplified_val.is_zero:
            zero_eigenvalues.append((simplified_val, basis))
        else:
            nonzero_eigenvalues.append((simplified_val, basis))
        
        # For each eigenvector, construct the projector
        for vec in basis:
            # Normalize the eigenvector
            norm = sp.sqrt(sum(v**2 for v in vec))
            normalized_vec = vec / norm
            # Calculate the projector P = |v⟩⟨v|
            projector = Matrix(normalized_vec) * Matrix(normalized_vec).transpose()
            projectors.append((simplified_val, projector))
    
    # Define the kernel projector (if any zero eigenvalues exist)
    P_kernel = sp.zeros(2, 2)
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
    V_inv_reg = sp.zeros(2, 2)
    
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
    
    # Verify boundedness: check that all 1/lambda_k are finite for non-zero eigenvalues    print("\nVerifying boundedness:")
    is_bounded = True
    for val, _ in nonzero_eigenvalues:
        inv_val = 1/val
        print(f"1/lambda = {simplify(inv_val)}")
        
        # Try to evaluate numerically to check boundedness more explicitly
        try:
            num_val = float(val.evalf())
            num_inv_val = float(inv_val.evalf())
            print(f"  λ = {num_val:.10f}")
            print(f"  1/λ = {num_inv_val:.10f}")
            
            if abs(num_val) < 1e-10:  # Close to zero
                is_bounded = False
                print(f"  Warning: Eigenvalue too close to zero: {num_val}")
        except Exception as e:
            print(f"  Could not evaluate numerically: {e}")
            
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
    
    # Save results to file
    os.makedirs('data', exist_ok=True)
    with open('data/regulated_inverse_operator.txt', 'w', encoding='utf-8') as f:
        f.write("Regulated Inverse Volume Operator Analysis\n")
        f.write("=======================================\n\n")
        f.write("V² Matrix:\n")
        f.write(str(V2_symbolic) + "\n\n")
        
        f.write("Eigenvalues of V²:\n")
        for val, mult in eigenvalues.items():
            f.write(f"lambda = {simplify(val)} (multiplicity: {mult})\n")
        f.write("\n")
        
        f.write("Boundedness Analysis:\n")
        for val, _ in nonzero_eigenvalues:
            inv_val = 1/val
            f.write(f"1/lambda = {simplify(inv_val)}\n")
        if is_bounded:            f.write("Conclusion: The operator is bounded (all 1/lambda_k are finite).\n\n")
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

if __name__ == "__main__":
    compute_regulated_inverse_operator()
    print("Regulated inverse operator analysis complete. Results saved to data/regulated_inverse_operator.txt")
