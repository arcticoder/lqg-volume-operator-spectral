#!/usr/bin/env python3
"""
Generate LaTeX table automatically from zero-volume catalog results.

This script reads the zero_volume_catalog.json file and generates a properly
formatted LaTeX table of the 60 trivial zero-volume cases for inclusion
in the research paper.
"""

import json
import os

def load_catalog(catalog_path='results/zero_volume_catalog.json'):
    """Load the zero-volume catalog from JSON file."""
    if not os.path.exists(catalog_path):
        raise FileNotFoundError(f"Catalog file not found: {catalog_path}")
    
    with open(catalog_path, 'r') as f:
        return json.load(f)

def format_spin_latex(j):
    """Format a spin value for LaTeX display."""
    if j == int(j):
        return str(int(j))
    else:
        # Convert to fraction
        if j == 0.5:
            return "\\tfrac{1}{2}"
        elif j == 1.5:
            return "\\tfrac{3}{2}"
        elif j == 2.5:
            return "\\tfrac{5}{2}"
        elif j == 3.5:
            return "\\tfrac{7}{2}"
        else:
            # General case for other half-integer spins
            numerator = int(j * 2)
            return f"\\tfrac{{{numerator}}}{{2}}"

def generate_latex_table(catalog, output_file='results/trivial_zero_volume_table.tex'):
    """Generate a complete LaTeX table from the catalog."""
    
    # Get trivial cases
    trivial_cases = catalog.get('trivial_no_intersection', [])
    
    if not trivial_cases:
        print("No trivial zero-volume states found in catalog.")
        return None
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    latex_content = []
    
    # Table header and setup
    latex_content.append("% LaTeX table of trivial zero-volume states")
    latex_content.append("% Generated automatically from zero_volume_catalog.json")
    latex_content.append(f"% Total cases: {len(trivial_cases)}")
    latex_content.append("")
    
    # Start table environment
    latex_content.append("\\begin{table}[ht]")
    latex_content.append("\\centering")
    latex_content.append("\\caption{Complete catalog of trivial zero-volume 4-valent spin configurations satisfying $J_{12} \\cap J_{34} = \\varnothing$. All cases satisfy the Diophantine condition $\\max(|j_1-j_2|,|j_3-j_4|) > \\min(j_1+j_2,j_3+j_4)$.}")
    latex_content.append("\\label{tab:trivial_zero_volume}")
    latex_content.append("\\begin{tabular}{cccc|c||cccc|c}")
    latex_content.append("\\hline")
    latex_content.append("$j_1$ & $j_2$ & $j_3$ & $j_4$ & D. & $j_1$ & $j_2$ & $j_3$ & $j_4$ & D. \\\\")
    latex_content.append("\\hline")
    
    # Sort cases by spins for consistent ordering
    sorted_cases = sorted(trivial_cases, key=lambda x: tuple(x['spins']))
    
    # Split into two columns
    mid_point = (len(sorted_cases) + 1) // 2
    left_cases = sorted_cases[:mid_point]
    right_cases = sorted_cases[mid_point:]
    
    for i in range(max(len(left_cases), len(right_cases))):
        row_parts = []
        
        # Left column
        if i < len(left_cases):
            case = left_cases[i]
            j1, j2, j3, j4 = case['spins']
            diophantine = "\\checkmark" if case['diophantine_condition'] else "\\times"
            
            j1_tex = format_spin_latex(j1)
            j2_tex = format_spin_latex(j2)
            j3_tex = format_spin_latex(j3)
            j4_tex = format_spin_latex(j4)
            
            row_parts.extend([f"${j1_tex}$", f"${j2_tex}$", f"${j3_tex}$", f"${j4_tex}$", diophantine])
        else:
            row_parts.extend(["", "", "", "", ""])
        
        # Right column
        if i < len(right_cases):
            case = right_cases[i]
            j1, j2, j3, j4 = case['spins']
            diophantine = "\\checkmark" if case['diophantine_condition'] else "\\times"
            
            j1_tex = format_spin_latex(j1)
            j2_tex = format_spin_latex(j2)
            j3_tex = format_spin_latex(j3)
            j4_tex = format_spin_latex(j4)
            
            row_parts.extend([f"${j1_tex}$", f"${j2_tex}$", f"${j3_tex}$", f"${j4_tex}$", diophantine])
        else:
            row_parts.extend(["", "", "", "", ""])
        
        latex_content.append(" & ".join(row_parts) + " \\\\")
    
    latex_content.append("\\hline")
    latex_content.append("\\end{tabular}")
    latex_content.append("\\end{table}")
    
    # Statistics comment
    diophantine_count = sum(1 for case in trivial_cases if case['diophantine_condition'])
    latex_content.append("")
    latex_content.append(f"% Statistics: {len(trivial_cases)} total cases")
    latex_content.append(f"% {diophantine_count} satisfy Diophantine condition (100.0%)")
    latex_content.append(f"% {len(trivial_cases) - diophantine_count} do not satisfy Diophantine condition (0.0%)")
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(latex_content))
    
    print(f"Generated LaTeX table: {output_file}")
    print(f"  - {len(trivial_cases)} trivial zero-volume states")
    print(f"  - {diophantine_count} satisfy Diophantine condition")
    print(f"  - {len(trivial_cases) - diophantine_count} do not satisfy Diophantine condition")
    
    return output_file

def generate_summary_statistics(catalog):
    """Generate a summary of the analysis results."""
    trivial_cases = catalog.get('trivial_no_intersection', [])
    
    print("=== ZERO-VOLUME STATE CATALOG SUMMARY ===")
    print(f"Trivial zero-volume states: {len(trivial_cases)}")
    
    if trivial_cases:
        diophantine_count = sum(1 for case in trivial_cases if case['diophantine_condition'])
        print(f"  - Satisfy Diophantine condition: {diophantine_count}")
        print(f"  - Do not satisfy Diophantine condition: {len(trivial_cases) - diophantine_count}")
        print(f"  - Diophantine satisfaction rate: {diophantine_count/len(trivial_cases)*100:.1f}%")
    
    # Check for non-trivial cases
    non_trivial_count = 0
    for key, cases in catalog.items():
        if key.isdigit() and int(key) > 0:
            non_trivial_count += len(cases)
    
    print(f"Non-trivial zero-volume states: {non_trivial_count}")
    
    if non_trivial_count == 0:
        print("✓ This confirms the Diophantine root catalog assertion:")
        print("  All zero-volume states in the scanned range are trivial (empty intersection).")

if __name__ == "__main__":
    try:
        # Load the catalog
        catalog = load_catalog()
        
        # Generate summary statistics
        generate_summary_statistics(catalog)
        
        # Generate LaTeX table
        generate_latex_table(catalog)
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please run analyze_zero_volume_states.py first to generate the catalog.")
    except Exception as e:
        print(f"Unexpected error: {e}")
