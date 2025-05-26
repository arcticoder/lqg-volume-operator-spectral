import sympy as sp
import pandas as pd
import numpy as np
import os
from sympy import Rational, factorial, sqrt, prod, Matrix, simplify, latex, symbols

# Define symbolic spin variables
j1, j2, j3, j4 = symbols('j1 j2 j3 j4', positive=True)

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

def compute_symbolic_spectrum():
    """
    Compute symbolic expressions for eigenvalues and projectors of the volume operator
    for the j=1/2 case.
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
        print(f"LaTeX: {latex(simplified_val)}")
        
        # Calculate the volume eigenvalue (square root)
        volume_val = sp.sqrt(abs(simplified_val))
        print(f"Volume eigenvalue: {simplify(volume_val)}")
        print(f"Volume eigenvalue LaTeX: {latex(simplify(volume_val))}")
        print()
    
    print("\nSymbolic eigenvectors and projectors:")
    for val, mult, basis in eigenvectors:
        print(f"For eigenvalue lambda = {simplify(val)}:")
        for i, vec in enumerate(basis):
            # Normalize the eigenvector
            norm = sp.sqrt(sum(v**2 for v in vec))
            normalized_vec = vec / norm
            print(f"Eigenvector {i+1}: {normalized_vec}")
            
            # Calculate the projector P = |v⟩⟨v|
            projector = Matrix(normalized_vec) * Matrix(normalized_vec).transpose()
            print(f"Projector {i+1}:")
            print(projector)
            print()
    
    # Calculate general form for arbitrary j values
    print("\nGeneral formula for eigenvalues lambda_k(j_i):")
    print("For arbitrary spins j1, j2, j3, j4, the eigenvalues are given by the roots of the characteristic polynomial of V²")
    print("where V² is built from the 12j symbols with the appropriate spin configurations.")
    
    # Save results to file using UTF-8 encoding
    os.makedirs('data', exist_ok=True)
    with open('data/symbolic_volume_expressions.txt', 'w', encoding='utf-8') as f:
        f.write("Symbolic Volume Operator Expressions\n")
        f.write("===================================\n\n")
        f.write("V² Matrix:\n")
        f.write(str(V2_symbolic) + "\n\n")
        
        f.write("Eigenvalues:\n")
        for val, mult in eigenvalues.items():
            f.write(f"lambda = {simplify(val)} (multiplicity: {mult})\n")
            f.write(f"LaTeX: {latex(simplify(val))}\n")
            f.write(f"Volume eigenvalue: {simplify(sp.sqrt(abs(val)))}\n")
            f.write(f"Volume eigenvalue LaTeX: {latex(simplify(sp.sqrt(abs(val))))}\n\n")
        
        f.write("Eigenvectors and Projectors:\n")
        for val, mult, basis in eigenvectors:
            f.write(f"For eigenvalue lambda = {simplify(val)}:\n")
            for i, vec in enumerate(basis):
                norm = sp.sqrt(sum(v**2 for v in vec))
                normalized_vec = vec / norm
                f.write(f"Eigenvector {i+1}: {normalized_vec}\n")
                
                projector = Matrix(normalized_vec) * Matrix(normalized_vec).transpose()
                f.write(f"Projector {i+1}:\n{projector}\n\n")
        
        # Add a section for the spectral form of V²
        f.write("\nSpectral form of V²:\n")
        f.write("V² = sum(lambda_k * P_k)\n\n")
        
        # Write eigenvalues as lambda_k(j_i)
        f.write("For j1=j2=j3=j4=1/2:\n")
        for val, _ in eigenvalues.items():
            simplified_val = simplify(val)
            f.write(f"lambda_k = {simplified_val}\n")
            
        # Also save the latex format for the volume eigenvalues (square root of abs of eigenvalues)
        f.write("\nVolume Eigenvalues in LaTeX:\n")
        for val, _ in eigenvalues.items():
            volume_val = sp.sqrt(abs(val))
            f.write(f"lambda_k = {latex(simplify(volume_val))}\n")

if __name__ == "__main__":
    os.makedirs('data', exist_ok=True)
    compute_symbolic_spectrum()
    print("Symbolic analysis complete. Results saved to data/symbolic_volume_expressions.txt")
