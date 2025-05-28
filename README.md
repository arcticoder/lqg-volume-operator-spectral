# LQG Volume Operator Spectral Analysis

This repository contains scripts and tools for computing and analyzing the spectrum of the Loop Quantum Gravity (LQG) volume operator and its kernel (zero-volume) states.

**📖 [View the complete research documentation on GitHub Pages](https://arcticoder.github.io/lqg-volume-kernel-catalog/)**

## Repository Name

**lqg-volume-kernel-catalog**

## Description

Scripts to:
- Compute the LQG volume operator spectrum via SU(2) 12j-symbols.
- Identify trivial and non-trivial zero-volume states for 4-valent nodes.
- Analyze kernel dimensions of the volume-squared matrix.
- Generate LaTeX tables and figures summarizing results.
- Validate the Diophantine characterization of zero-volume states.

## Topics

`loop-quantum-gravity` `volume-operator` `quantum-gravity` `su2` `recoupling-theory` `12j-symbols` `latex`

## Requirements

See `requirements.txt` for dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### 1. Analyze Zero-Volume States (4-valent)

```bash
python scripts/analyze_zero_volume_states.py
```
- Scans spin configurations (j=0.5 to 3.0).
- Identifies trivial zero-volume states (J₁₂ ∩ J₃₄ = ∅).
- Confirms no non-trivial zero-volume states in the tested range.
- Outputs summary statistics and generates:
  - `results/kernel_dimension_distribution.png`
  - `results/zero_volume_catalog.json`

### 2. Generate LaTeX Table of Trivial Zero-Volume States

```bash
python scripts/generate_latex_table.py
```
- Produces `results/trivial_zero_volume_table.tex` with a single-column summary of trivial zero-volume cases.

### 3. Create Spin-1/2 Correlation Placeholder Figure

```bash
python scripts/generate_spin_half_correlation.py
```
- Generates `results/figures/spin_half_correlation.png` showing the absence of non-trivial kernel states.

### 4. Generate Volume Spectrum (Optional)

```bash
python scripts/compute_volume_spectrum.py
```
- Computes eigenvalues of the volume-squared operator for specified spin configurations.
- Outputs to `data/volume_spectrum.csv`.

## Results

- **Total configurations scanned**: 1,296 (j_i ∈ {0.5, 1.0, …, 3.0})
- **Trivial zero-volume states**: 60 (4.6%)
- **Non-trivial zero-volume states**: 0 (0.0%)
- **Full-rank matrices (non-zero volume)**: 674 (52.0%)
- **Other kernel matrices (not zero-volume)**: 562 (43.4%)

All trivial zero-volume cases satisfy the Diophantine condition:

```
max(|j₁−j₂|, |j₃−j₄|) > min(j₁+j₂, j₃+j₄)
```

## Directory Structure

```
├── scripts/
│   ├── analyze_zero_volume_states.py
│   ├── generate_latex_table.py
│   ├── generate_spin_half_correlation.py
│   ├── compute_volume_spectrum.py
│   └── ...
├── results/
│   ├── trivial_zero_volume_table.tex
│   ├── zero_volume_catalog.json
│   ├── kernel_dimension_distribution.png
│   └── figures/
│       └── spin_half_correlation.png
├── data/
│   └── volume_spectrum.csv
├── README.md
└── requirements.txt
```
