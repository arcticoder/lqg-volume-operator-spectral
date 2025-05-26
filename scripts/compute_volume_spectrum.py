import sympy as sp
import pandas as pd
import numpy as np
import os
from sympy import Rational, factorial, sqrt, prod

# 1. Numeric spins
j1 = j2 = j3 = j4 = 0.5
J12 = [0.0, 1.0]

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

# 2. Build numeric V2
V2_num = np.zeros((2, 2))
for a, ja in enumerate(J12):
    for b, jb in enumerate(J12):
        # Corrected spin mirroring with proper intermediate spin ordering for 12j symbol
        # Ordering: (j1,j2,j12), (j3,j4,j12), (j1,j2,j34), (j3,j4,j34)
        V2_num[a, b] = CF12j_num(j1, j2, ja, j3, j4, ja, j1, j2, jb, j3, j4, jb)

# 3. Compute numeric eigenvalues
eigenvalues = np.linalg.eigvals(V2_num)

# 4. Save to CSV
os.makedirs('data', exist_ok=True)
df = pd.DataFrame({'eigenvalue': eigenvalues})
csv_path = 'data/volume_spectrum.csv'
df.to_csv(csv_path, index=False)

print(f"Numeric eigenvalues written to {csv_path}")
