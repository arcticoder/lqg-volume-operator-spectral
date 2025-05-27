**find_zero_volume_valence4.py**
A companion script that
-   builds the V² matrix **only** from J₁₂,    
-   checks its determinant for zeros,    
-   flags any non-trivial singularities (and “suspicious” near-zeros),    
-   and prints out any found cases.
Even though it uses the old “J₁₂-only” logic, it can still serve as a regression test: we can run it before/after our fix to make sure nothing unexpected pops up.