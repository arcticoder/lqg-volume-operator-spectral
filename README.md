# LQG Volume Operator Spectral Analysis

This repository contains scripts for computing and analyzing the spectrum of the Loop Quantum Gravity (LQG) volume operator and its inverse.

## Overview

The volume operator in Loop Quantum Gravity has a discrete spectrum, which we compute using the 12j-symbol formulation. This project provides tools to:

1. Compute the volume spectrum for different spin network configurations
2. Analyze the inverse volume operator, which is important for regularization schemes

## Requirements

See `requirements.txt` for the necessary Python packages. Install them with:

```bash
pip install -r requirements.txt
```

## Usage

### Computing Volume Spectrum

The script `scripts/compute_volume_spectrum.py` calculates the eigenvalues of the volume-squared operator for a given set of spins:

```bash
python scripts/compute_volume_spectrum.py
```

The resulting spectrum is saved to `data/volume_spectrum.csv`.

By default, it uses j₁ = j₂ = j₃ = j₄ = 1/2, but you can modify the spin values in the script to explore other configurations.

### Analyzing Inverse Volume

The script `scripts/analyze_inverse_volume.py` computes the inverse volume operator spectrum and operator norms:

```bash
python scripts/analyze_inverse_volume.py
```

The output is saved to `data/inverse_volume_spectrum.csv`.

## Mathematical Background

The scripts implement the volume operator using 12j-symbols with specific spin configurations. The volume-squared operator is constructed by computing matrix elements:

V²_{ab} = CF12j(j₁, j₂, j₁₂ᵃ, j₃, j₄, j₁₂ᵃ, j₁, j₂, j₁₂ᵇ, j₃, j₄, j₁₂ᵇ)

where j₁₂ᵃ and j₁₂ᵇ run through the allowed intermediate coupling values.

## Extending the Analysis

### Comparing Different Spin Configurations

The script `scripts/compare_spin_configs.py` automatically computes spectra for multiple spin configurations:

```bash
python scripts/compare_spin_configs.py
```

This generates:
- `data/volume_spectra_comparison.csv`: Detailed eigenvalue data
- `data/volume_spectra_table.csv`: A formatted table comparing spectra
- `data/volume_spectra_comparison.png`: A visualization of different spectra

Pre-configured spin configurations include:
- j₁ = j₂ = j₃ = j₄ = 1/2
- j₁ = j₂ = 1, j₃ = 1/2, j₄ = 3/2
- j₁ = j₂ = j₃ = j₄ = 1
- j₁ = j₂ = j₃ = j₄ = 3/2
- j₁ = 1/2, j₂ = 1, j₃ = 3/2, j₄ = 2

To add more configurations, modify the `spin_configs` list in the script.

### Custom Configurations

For custom analysis, you can also modify the `j1`, `j2`, `j3`, `j4` values and the `J12` array in the original scripts to include all possible intermediate coupling values allowed by recoupling theory.
