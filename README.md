# 2D-ADER: Finite Difference Young-Laplace Method & Experimental Droplet Analysis

[![DOI](https://doi.org/10.5281/zenodo.23405107)](https://doi.org/10.5281/zenodo.23405107)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

This repository contains the complete computational implementation, experimental datasets, and numerical benchmarks for the **2D-ADER Finite Difference (FD) Young-Laplace method**. The framework enables high-precision determination of interfacial/surface tension ($\Gamma$) and geometrical properties of axisymmetric pendant and sessile drops without requiring non-linear shooting algorithms or ODE integration from the apex.

---

## 1. Overview & Theoretical Background

Standard Axisymmetric Drop Shape Analysis (ADSA) algorithms iteratively integrate non-linear Young-Laplace differential equations from an estimated apex curvature. In contrast, the 2D-ADER method calculates interface curvature directly from experimental edge coordinates using high-order finite difference stencils.

For an axisymmetric liquid droplet in gravitational and mechanical equilibrium, the hydrostatic pressure jump across the interface is described by the Young-Laplace equation:

$$\Delta P(y) = \Delta P_0 + \Delta \rho g y = 2 H(y) \Gamma$$

Where:
* $\Delta \rho = \rho_{\text{liquid}} - \rho_{\text{fluid}}$ is the density difference between phases ($\text{kg/m}^3$).
* $g$ is the local acceleration due to gravity ($\text{m/s}^2$).
* $y$ is the vertical elevation coordinate oriented along the gravity vector ($\text{m}$).
* $2H(y) = \kappa_1 + \kappa_2$ is the total interface mean curvature ($\text{m}^{-1}$).
* $\Gamma$ is the interfacial/surface tension ($\text{N/m}$).

### Dimensionless Formulation
Spatial coordinates are non-dimensionalized by the equatorial radius $R_0 = r_{\text{pixels}} \cdot \text{pixel\_size}$:

$$\bar{x} = \frac{x}{R_0}, \quad \bar{y} = \frac{y}{R_0}$$

The dimensionless mean curvature $2\bar{H}$ is evaluated analytically along the profile as:

$$2\bar{H}(\bar{y}) = - \frac{\frac{d^2 \bar{x}}{d \bar{y}^2}}{\left[1 + \left(\frac{d \bar{x}}{d \bar{y}}\right)^2\right]^{3/2}} + \frac{1}{\bar{x} \left[1 + \left(\frac{d \bar{x}}{d \bar{y}}\right)^2\right]^{1/2}}$$

Substituting this expression into the hydrostatic equilibrium yields a linear relationship across the hydrostatic region:

$$2\bar{H}(\bar{y}) = B \bar{y} + \text{Constant}$$

Where $B$ is the dimensionless Bond number / shape factor:

$$B = \frac{\Delta \rho g R_0^2}{\Gamma} \implies \Gamma = \frac{\Delta \rho g R_0^2}{B}$$

Linear regression over the hydrostatic profile yields $B$ and directly determines $\Gamma$.

---

## 2. Repository Structure

```text
.
├── LICENSE                          # MIT License for code and analysis scripts
├── CITATION.cff                     # Citation metadata for automatic Zenodo indexing
├── README.md                        # Project documentation (this file)
│
├── src/                             # Core computational code
│   ├── FD_method.py                 # Multiprocessing analysis pipeline and batch processing
│   ├── fd_utils.py                  # High-order FD stencils, PBL sweeps, and geometry routines
│   └── requirements.txt             # Environment dependencies
│
├── experimental_data/               # Experimental drop profiles and summary reports
│   ├── README.md                    # Detailed data dictionary and measurement protocols
│   ├── small_drops/                 # Edge coordinate text files (Profile_{id}.txt)
│   ├── intermediate_drops/
│   ├── large_drops/
│   └── Summary_Results.xlsx         # Comparative metrics (2D-ADER vs. OpenDrop)
│
└── numerical_simulations/           # Theoretical benchmark profiles and verification
    ├── profiles/                    # Ground-truth profiles obtained via ODE integration
    └── images/                      # Rendered synthetic drop images for noise robustness tests
```

---

## 3. Installation & Requirements

The software runs on Linux, macOS, and Windows with **Python 3.8+**.

1. **Clone the repository:**
   ```bash
   git clone https://github.com/robert31626/Surface-Tension-by-FD
   cd your-repo-name
   ```

2. **Create a virtual environment (optional but recommended):**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies:**
   ```bash
   pip install numpy scipy pandas matplotlib openpyxl
   ```

---

## 4. Software Usage & Execution (`src/`)

### Input Data Format
Edge coordinates must be provided as tab- or space-delimited text files named either `Profile_{index}.txt` or `Perfil_{index}.txt`, structured in three columns:
* **Column 0:** Right branch $X$ coordinate (pixels)
* **Column 1:** Left branch $X$ coordinate (pixels)
* **Column 2:** Vertical $Y$ coordinate (pixels)

### Configuration
All parameters are managed in the `CONFIG` dictionary located at the top of `src/FD_method.py`:

| Parameter | Default Value | Description |
| :--- | :--- | :--- |
| `INPUT_DIR` | `r".../Video {video}"` | Source directory containing raw profile coordinate files |
| `OUTPUT_DIR` | `r".../Video {video}"` | Destination directory for export tables and diagnostic plots |
| `VIDEOS` | `[2]` | List of video/experiment sequence identifiers |
| `PROFILE_INDEXES` | `np.arange(0, 2)` | Range or array of profile indices to process |
| `PIXEL_SIZE` | `1.3661e-6` | Optical magnification calibration factor ($\text{m/pixel}$) |
| `FLUID_DENSITY` | `997.5` | Working fluid density ($\text{kg/m}^3$) |
| `GRAVITY` | `9.78` | Local gravitational acceleration ($\text{m/s}^2$) |
| `NEEDLE_DIAMETER` | `0.0015` | Outer capillary needle diameter ($\text{m}$) |
| `SURFACE_TENSION_REF`| `0.07278` | Reference surface tension for Worthington number ($\text{N/m}$) |
| `PBLS` | `np.arange(100, 101)` | Sampling stride values (Point Boundary Length) to evaluate |
| `PBL_EVAL_MIN` | `100` | Lower bound of the evaluation plateau window |
| `PBL_EVAL_MAX` | `100` | Upper bound of the evaluation plateau window |
| `TAIL_PERCENTAGE` | `90.0` | Upper profile percentage included in the linear fit |
| `NUM_WORKERS` | `4` | Number of CPU cores for parallel execution (`multiprocessing`) |
| `SAVE_DIAGNOSTIC_PLOTS`| `True` | Flag to generate and save diagnostic PNG plots |

### Running the Analysis
To execute the pipeline:
```bash
python src/FD_method.py
```

### Generated Outputs
The pipeline outputs structured datasets and graphical diagnostics into `OUTPUT_DIR`:
1. `FD_Analysis_Results[_Video_{id}].xlsx` & `.txt`: Primary summary table reporting $R_0$, $B$, $\Gamma$, contact angle, volume, area, and Worthington number.
2. `FD_Standard_Deviations[_Video_{id}].xlsx`: Linear regression standard errors across evaluated strides.
3. `FD_B_Branches[_Video_{id}].xlsx`: Branch-resolved $B$ values for left and right drop profiles.
4. `FD_Gamma_vs_PBLs[_Video_{id}].xlsx` & `FD_B_vs_PBLs[_Video_{id}].xlsx`: Global stability matrices as a function of stride length.
5. Diagnostic Plots (`.png`): Stride sensitivity plots ($B$, $\Gamma$, standard error vs. PBL) and curvature vs. elevation fit curves.

---

## 5. Experimental Data Dictionary

The file `experimental_data/Summary_Results.xlsx` compiles experimental results across drop size regimes, contrasting the 2D-ADER FD method with standard **OpenDrop** (OP) fits:

| Column | Unit | Description |
| :--- | :---: | :--- |
| `OP_theta` | $\text{deg}$ | Symmetry axis rotation angle identified by OpenDrop |
| `R0` | $\text{mm}$ | Mean apex radius of curvature ($R_0$) |
| `OP_Bo` | — | Mean Bond number relative to apex curvature (OpenDrop) |
| `OP_S_Bo` | — | Standard deviation of Bond number (OpenDrop) |
| `OP_gamma` | $\text{N/m}$ | Mean surface tension determined by OpenDrop |
| `OP_S_gamma` | $\text{N/m}$ | Standard deviation (repeatability) of surface tension (OpenDrop) |
| `d` | $\text{mm}$ | Mean equatorial droplet diameter |
| `B_FD` | — | Mean Bond number / shape parameter from 2D-ADER FD method |
| `S_B_FD` | — | Standard deviation of $B$ across the sequence (FD method) |
| `<Sigma_B>_FD` | — | Mean standard error ($\sigma_B$) of the linear regression slope |
| `Gamma_FD` | $\text{N/m}$ | Mean surface tension computed via 2D-ADER FD method |
| `S_Gama_FD` | $\text{N/m}$ | Standard deviation of surface tension across sequence (FD method) |
| `Shape_parameter` | — | Mean drop deviation parameter from spherical cap geometry |
| `Contact_Angle` | $\text{deg}$ | Mean contact angle ($\theta$) |
| `Volume` | $\mu\text{L}$ | Axisymmetric volume calculated via numerical disk slicing integration |
| `Area` | $\text{mm}^2$ | Projected 2D profile area computed via trapezoidal integration |
| `Neumann_number` | — | Dimensionless Neumann number |
| `Worthington_number`| — | Dimensionless Worthington number ($Wo = \frac{\Delta \rho g V}{\pi \Gamma_{\text{ref}} d_{\text{needle}}}$) |
| `Fitted_points` | integer | Mean number of spatial points utilized in the linear regression |

---

## 6. Citation & Attribution

When using this code, methodology, or dataset in academic publications, please cite the software archive:

```bibtex
@software{bahena_2026_zenodo,
  author       = {Bahena, Roberto U.},{Escovar, Juan V.}
  title        = {Finite Difference Young-Laplace Method and Experimental Droplet Datasets},
  month        = oct,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.23405107},
  url          = {https://doi.org/10.5281/zenodo.23405107}
}
```

---

## 7. License

* The software and analytical scripts in `src/` are distributed under the [MIT License](LICENSE).
* The experimental and numerical datasets in `experimental_data/` and `numerical_simulations/` are licensed under [Creative Commons Attribution 4.0 International (CC-BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
