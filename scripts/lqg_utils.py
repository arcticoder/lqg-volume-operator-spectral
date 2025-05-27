"""
Common utilities for Loop Quantum Gravity volume operator calculations.
This module contains shared functions for both symbolic and numeric calculations.
"""

import sympy as sp
import numpy as np
from sympy import Rational, factorial, sqrt, prod, symbols, gamma

def CF12j_symbolic(j1, j2, j12, j3, j4, j23, j5, j6, j34, j7, j8, j45):
    """
    Symbolic implementation of the 12j coefficients for the volume operator.
    
    Args:
        j1, j2, j12: First triple of angular momenta
        j3, j4, j23: Second triple of angular momenta
        j5, j6, j34: Third triple of angular momenta
        j7, j8, j45: Fourth triple of angular momenta
    
    Returns:
        The symbolic expression for the 12j coefficient
    """
    # Define the spin variables
    a, b, c = j1, j2, j12
    d, e, f = j3, j4, j23
    g, h, i = j5, j6, j34
    k, l, m = j7, j8, j45
    triples = [(a, b, c), (d, e, f), (g, h, i), (k, l, m)]
    
    # Check triangle inequalities symbolically if possible
    for (x, y, z) in triples:
        # Skip symbolic check as it will be handled at evaluation time
        if hasattr(x, 'is_Symbol') or hasattr(y, 'is_Symbol') or hasattr(z, 'is_Symbol'):
            continue
        if x + y < z or x + z < y or y + z < x:
            return sp.Integer(0)
    
    # Define the symbolic Delta factor using gamma functions for better generalization
    # with symbolic parameters (gamma is a generalization of factorial)
    Delta_expr = sqrt(prod(
        gamma(-x + y + z + 1) * gamma(x - y + z + 1) * gamma(x + y - z + 1)
        / gamma(x + y + z + 2)
        for (x, y, z) in triples
    ))
    
    # Define hypergeometric function parameters
    num = [Rational(1, 2), -c, c + 1, -f, f + 1]
    den = [a + b - c + 1, d + e - f + 1, g + h - i + 1, k + l - m + 1]
    
    # Return the symbolic expression
    return Delta_expr * sp.hyper(num, den, 1)

def CF12j_numeric(j1, j2, j12, j3, j4, j23, j5, j6, j34, j7, j8, j45):
    """
    Numeric implementation of the 12j coefficients for the volume operator.
    
    Args:
        j1, j2, j12: First triple of angular momenta
        j3, j4, j23: Second triple of angular momenta
        j5, j6, j34: Third triple of angular momenta
        j7, j8, j45: Fourth triple of angular momenta
    
    Returns:
        A float value for the 12j coefficient
    """
    # Convert to rational internally for exact prefactor
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
    num = [Rational(1, 2), -c, c + 1, -f, f + 1]
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
    """
    Get allowed intermediate coupling values for given j1, j2
    
    Args:
        j1, j2: Angular momenta values
        
    Returns:
        List of allowed intermediate coupling values
    """
    # Handle symbolic case
    if hasattr(j1, 'is_Symbol') or hasattr(j2, 'is_Symbol'):
        return [symbols(f'J_{{{j1}{j2}}}')]
    
    j_min = abs(float(j1) - float(j2))
    j_max = float(j1) + float(j2)
    # Generate all half-integer values between j_min and j_max
    step = 0.5 if (j_min % 1 != 0 or j_max % 1 != 0) else 1.0
    return [j for j in np.arange(j_min, j_max + step, step)]

def create_symbolic_spins(j_values=None):
    """
    Create symbolic spin variables with optional numeric values.
    
    Args:
        j_values: Optional tuple of (j1, j2, j3, j4) numeric values
        
    Returns:
        Tuple of (j1, j2, j3, j4) symbolic variables or numeric values
    """
    if j_values is None:
        # Create symbolic variables
        return symbols('j1 j2 j3 j4', positive=True)
    else:
        # Use the provided numeric values
        return j_values
