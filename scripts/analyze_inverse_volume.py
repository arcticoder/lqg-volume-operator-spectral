import numpy as np
import pandas as pd
import os
import sympy as sp
from sympy import Rational, factorial, sqrt, prod

# Define the spins and J12 values (same as in compute_volume_spectrum.py)
j1 = j2 = j3 = j4 = 0.5
J12 = [0.0, 1.0]

def CF12j_num(j1, j2, j12, j3, j4, j23, j5, j6, j34, j7, j8, j45):
    # convert to rational internally for exact prefactor
    a, b, c = map(Rational, (j1, j2, j12))
    d, e, f = map(Rational, (j3, j4, j23))
    g, h, i = map(Rational, (j5, j6, j34))
    k, l, m = map(Rational, (j7, j8, j45))
    triples = [(a, b, c), (d, e, f), (g, h, i), (k, l, m)]
    Delta = sqrt(prod(
        factorial(-x + y + z) * factorial(x - y + z) * factorial(x + y - z)
        / factorial(x + y + z + 1)
        for (x, y, z) in triples
    ))
    num = [sp.Rational(1, 2), -c, c + 1, -f, f + 1]
    den = [a + b - c + 1, d + e - f + 1, g + h - i + 1, k + l - m + 1]
    val = Delta * sp.N(sp.hyper(num, den, 1))
    return float(val)

# Build V2_num
V2_num = np.zeros((2, 2))
for a, ja in enumerate(J12):
    for b, jb in enumerate(J12):
        V2_num[a, b] = CF12j_num(j1, j2, ja, j3, j4, ja, j1, j2, jb, j3, j4, jb)

# Compute eigenvalues and (orthonormal) eigenvectors
eigvals, eigvecs = np.linalg.eigh(V2_num)

ε = 1e-6  # small regulator
P_nonzero = []
P_kernel = None  # Initialize P_kernel
for val, vec in zip(eigvals, eigvecs.T):
    P = np.outer(vec, vec)
    if abs(val) > 1e-12:
        P_nonzero.append((val, P))
    else:
        P_kernel = P

# If we didn't find any kernel vectors, set P_kernel to zero matrix
if P_kernel is None:
    P_kernel = np.zeros_like(V2_num)

# regulated inverse (as a matrix)
V_inv = sum((1/val)*P for val, P in P_nonzero) + (1/ε)*P_kernel

# operator norm of the non-kernel part
V_inv_nonker = sum((1/val)*P for val, P in P_nonzero)
norm = np.linalg.norm(V_inv_nonker, 2)
print("Non-kernel operator norm:", norm)
# confirm norm < ∞ as ε→0

# 4. Save to CSV
os.makedirs('data', exist_ok=True)
df = pd.DataFrame({'eigenvalue': eigvals, 'norm': norm})
csv_path = 'data/inverse_volume_spectrum.csv'
df.to_csv(csv_path, index=False)