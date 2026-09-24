# Linear System Resolution using the Gauss-Seidel Method

This numerical analysis project was conducted as part of the engineering curriculum in the GM department at **INSA Rouen Normandie** (Academic Year 2024-2025). Supervised by **Mr. J.G. Caputo**, this work was carried out by **Mouad Sheradj Drissi** and **Amr Zaki Salih**.

---

## Overview

The goal of this project is to study, implement, and benchmark the iterative **Gauss-Seidel method** for solving linear systems of the form $Ax = b$.

Two distinct algorithmic formulations are implemented and evaluated across scientific programming languages (**Python**, **Matlab**, **Fortran**)[cite: 8]:
1. **Program A (Matrix Decomposition Formulation):** Leverages the splitting $A = D + L + U$ expressed iteratively as $X_{k+1} = FX_k + G$, where $F = -(D+L)^{-1}U$ and $G = (D+L)^{-1}b$.
2. **Program B (Component-wise / Direct Iterative Formulation):** Directly computes each updated entry $x_i^{(k+1)}$ sequentially without explicit matrix inversion.

---

## Key Numerical Findings

- **Matrix vs. Component-wise Performance:**
  - The vectorized matrix approach (Program A) leverages low-level BLAS/NumPy routines and significantly outperforms explicit scalar loop formulations (Program B) as matrix dimensions grow.
- **Diagonal Dominance & Conditioning:**
  - Regularized Hilbert matrices ($A = H_n + \omega I_n$) guarantee strict diagonal dominance and ensure rapid convergence.
  - Matrices lacking diagonal dominance or displaying poor conditioning diverge or require substantially more iterations to reach tolerance.
- **Single vs. Double Precision (Fortran):**
  - Single precision (32-bit float, ~7 significant digits) plateaus due to early round-off errors.
  - Double precision (64-bit float, ~15 significant digits) achieves residual convergence down to $10^{-9}$ while avoiding numerical stagnation.

---

## Repository Structure

    
    ├── Rapport_projet_1_MMSN.pdf    # Full technical report (INSA Rouen)
    ├── codesA/                       # Matrix implementations (X_{k+1} = F*X_k + G)
    │   ├── code gauss seidel pythonA.py
    │   ├── GaussSeidelA.f            # Fortran implementation (Single Precision)
    │   ├── GaussSeidelADouble.f      # Fortran implementation (Double Precision)
    │   └── Makefile
    ├── codesB/                       # Direct scalar component-wise implementations
    │   ├── code gauss seidel pythonB.py
    │   ├── GaussSeidelB.f            # Fortran implementation
    │   └── Makefile
    └── README.md
---
## Compilation & Execution

### Python Setup

*Requires Python 3 and scientific computing packages*

**To run the matrix implementation**:
    
    pip install numpy scipy matplotlib

### Fortran Setup

*Requires a standard Fortran compiler such as gfortran*

---

## Compile the programs:

### Compile Program A (Single Precision)
    gfortran -g GaussSeidelA.f -o GS_A

### Compile Program A (Double Precision)
    gfortran -g GaussSeidelADouble.f -o GS_A_Double

### Compile Program B (Component-wise)
    gfortran -g GaussSeidelB.f -o GS_B

**Run an executable**:
    
    ./GS_A
---
## Authors & Acknowledgments

**Mouad Sheradj Drissi**

**Amr Zaki Salih**

Academic Advisor: **Mr. J.G. Caputo** (INSA Rouen Normandie)

--- 
## Institution

INSA Rouen Normandie
Department of Mathematical Engineering

