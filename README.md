# Surface-Tension-by-FD

# 2D-ADER: Finite Difference Young-Laplace Method & Experimental Droplet Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

This repository contains the official code, numerical simulations, and experimental datasets for the **2D-ADER Finite Difference (FD) Young-Laplace method**. The suite provides an automated pipeline to determine interfacial/surface tension ($\Gamma$) and drop geometrical parameters from axisymmetric pendant and sessile drop profiles without requiring shooting algorithms or non-linear ODE integrations.

---

## 1. Overview & Theoretical Background

Unlike conventional Axisymmetric Drop Shape Analysis (ADSA) algorithms that iteratively solve the non-linear Young-Laplace differential equations from an estimated apex curvature, the 2D-ADER framework computes interface curvature directly from experimental edge coordinates using high-order central finite difference stencils.

For an axisymmetric drop in mechanical and gravitational equilibrium, the hydrostatic pressure jump across the liquid-fluid interface obeys the Young-Laplace equation:

$$\Delta P(y) = \Delta P_0 + \Delta \rho g y = 2 H(y) \Gamma$$

Where:
* $\Delta \rho = \rho_{\text{liquid}} - \rho_{\text{fluid}}$ is the density difference ($\text{kg/m}^3$).
* $g$ is the local acceleration due to gravity ($\text{m/s}^2$).
* $y$ is the vertical elevation oriented along the gravity vector ($\text{m}$).
* $2H(y) = \kappa_1 + \kappa_2$ is the total mean curvature ($\text{m}^{-1}$).
* $\Gamma$ is the interfacial/surface tension ($\text{N/m}$).

### Dimensionless Formulation
Scaling spatial coordinates by the equatorial radius $$R_0 = r_{\text{pixels}} \cdot \text{pixel\_size}$$:

$$\bar{x} = \frac{x}{R_0}, \quad \bar{y} = \frac{y}{R_0}$$

The dimensionless mean curvature $2\bar{H}$ is evaluated as:

$$2\bar{H}(\bar{y}) = -\frac{\frac{d^2\bar{x}}{d\bar{y}^2}}{\left[1 + \left(\frac{d\bar{x}}{d\bar{y}}\right)^2\right]^{3/2}} + \frac{1}{\bar{x}\left[1 + \left(\frac{d\bar{x}}{d\bar{y}}\right)^2\right]^{1/2}}$$

Substituting into the balance equation yields a linear relationship over the hydrostatic region:

$$2\bar{H}(\bar{y}) = B \bar{y} + \text{Constant}$$

Where $B$ is the dimensionless Bond number / shape factor:

$$B = \frac{\Delta \rho g R_0^2}{\Gamma} \implies \Gamma = \frac{\Delta \rho g R_0^2}{B}$$

---

## 2. Repository Structure

```text
.
├── LICENSE                          # MIT License for software and analysis code
├── CITATION.cff                     # Citation metadata for automated Zenodo indexing
├── README.md                        # Root documentation (this file)
│
├── src/                             # Core Python computational pipeline
│   ├── FD_method.py                 # Multi-core processing, batch evaluation, and reporting
│   ├── fd_utils.py                  # High-order FD stencils, PBL sweeps, geometry routines
│   ├── requirements.txt             # Python environment dependencies
│   └── README.md                    # Detailed software documentation and execution flags
│
├── experimental_data/               # Experimental edge profiles and summary reports
│   ├── README.md                    # Data dictionary and experimental metadata
│   ├── small_drops/                 # Coordinate files categorized by droplet volume/regime
│   ├── intermediate_drops/
│   ├── large_drops/
│   └── Summary_Results.xlsx         # Comparative metrics (2D-ADER vs. OpenDrop)
│
└── numerical_simulations/           # Theoretical and synthetic reference profiles
    ├── profiles/                    # Ground-truth profiles generated via ODE integration
    └── images/                      # Rendered synthetic drop profiles for noise validation
