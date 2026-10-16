================================================================================
FINITE DIFFERENCE YOUNG-LAPLACE METHOD (2D-ADER)
SUPPLEMENTAL MATERIAL SOFTWARE DOCUMENTATION
================================================================================

1. OVERVIEW
-----------
This package provides a robust, automated Python implementation of the 2D-ADER 
Finite Difference (FD) Young-Laplace method for surface tension determination 
and geometrical analysis of axisymmetric pendant or sessile drops.

Unlike traditional numerical methods that rely on shooting algorithms to solve 
the non-linear Young-Laplace Ordinary Differential Equations (ODEs) from an 
estimated apex curvature, this method directly calculates the interface curvature 
from experimental edge coordinates using high-order finite difference stencils. 
Linear regression across the hydrostatic pressure region then yields the 
dimensionless shape factor (B) and the interfacial tension (Gamma), providing 
high stability, reproducibility, and computational efficiency.


2. FILE STRUCTURE
-----------------
The package is organized into two primary Python scripts and this documentation:

- FD_method.py:
    Main analysis pipeline. Manages multi-core parallel processing (via 
    multiprocessing.Pool), batch execution across image/profile series, 
    diagnostic visualization generation, and automated export of structured 
    Excel (.xlsx) and text (.txt) reports.

- fd_utils.py:
    Dedicated mathematical and geometrical utility module. Contains 
    high-order central finite difference stencils, multi-stride interleaved 
    curvature evaluation (PBL sweeps), apex localization, contact angle 
    estimation, drop volume integration, and projected cross-sectional area 
    calculation.

- README.txt:
    This documentation file, explaining the methodology, configuration, 
    and output formats.


3. THEORETICAL BACKGROUND AND MATHEMATICAL FORMULATION
------------------------------------------------------
For an axisymmetric liquid droplet suspended from a capillary needle under 
gravitational equilibrium, the pressure difference across the interface satisfies 
the Young-Laplace equation:

    Delta_P(y) = Delta_P_0 + Delta_rho * g * y = 2 * H(y) * Gamma

where:
    Delta_rho : Density difference between liquid and fluid phases [kg/m^3]
    g         : Local gravitational acceleration [m/s^2]
    y         : Vertical coordinate oriented along the gravity vector [m]
    2*H(y)    : Total mean curvature (sum of principal curvatures kappa_1 + kappa_2) [1/m]
    Gamma     : Interfacial / surface tension [N/m]

By scaling spatial dimensions by the equatorial radius R_0 = r_pixels * pixel_size,
the dimensionless coordinates are:
    x_bar = x / R_0
    y_bar = y / R_0

The dimensionless mean curvature 2*H_bar is given analytically by:
    2*H_bar = - (d^2 x_bar / dy_bar^2) / [1 + (dx_bar / dy_bar)^2]^(3/2)
              + 1 / (x_bar * [1 + (dx_bar / dy_bar)^2]^(1/2))

Substituting into the Young-Laplace equation yields a linear relation:
    2*H_bar(y_bar) = B * y_bar + Constant

where B is the dimensionless Bond number / shape factor:
    B = (Delta_rho * g * R_0^2) / Gamma

Consequently, the surface tension is obtained directly from the fitted slope B:
    Gamma = (Delta_rho * g * R_0^2) / B


4. GEOMETRICAL PARAMETERS
-------------------------
In addition to surface tension, the pipeline computes key geometrical indicators:

- Equatorial Radius (R_0):
    Half-width of the widest horizontal section of the droplet:
    R_0 = r_pixels * pixel_size [m]

- Projected Profile Area:
    Cross-sectional projected 2D area calculated via numerical trapezoidal rule:
    Area = 0.5 * sum((x_i + x_{i-1}) * (y_i - y_{i-1})) * R_0^2 * 10^6 [mm^2]

- Drop Volume:
    Axisymmetric volume computed by numerical disk slicing integration:
    Volume = sum(pi * ((x_i - x_{-i})/2)^2 * Delta_y_i) * R_0^3 * 10^9 [uL]

- Contact Angle (theta) & Shape Parameter:
    Estimated using the 2016 spherical cap deviation protocol. A linear fit to 
    the upper 1% of points near the capillary contact determines the contact angle, 
    while the shape parameter quantifies departure from spherical geometry.

- Worthington Number (Wo):
    Dimensionless ratio of drop weight to capillary retention force:
    Wo = (Delta_rho * g * Volume_m3) / (pi * Gamma_ref * Needle_diameter)


5. INPUT DATA FORMAT
--------------------
Profile coordinates should be provided as tab- or space-delimited text files 
named either:
    Profile_{index}.txt  or  Perfil_{index}.txt

Each file contains the detected edge coordinates with columns:
    Column 0: Right branch X coordinate (pixels)
    Column 1: Left branch X coordinate (pixels)
    Column 2: Vertical Y coordinate (pixels)


6. CONFIGURATION GUIDE (FD_method.py)
-------------------------------------
All user-configurable parameters are defined in the 'CONFIG' dictionary at the 
top of FD_method.py:

Parameter               Default Value               Description
--------------------------------------------------------------------------------
INPUT_DIR               r".../Video {video}"        Directory containing coordinate files.
OUTPUT_DIR              r".../Video {video}"        Destination directory for outputs.
VIDEOS                  [2]                         List of video/experiment IDs to iterate.
PROFILE_INDEXES         np.arange(0, 2)             Range or list of profile indices.
PIXEL_SIZE              1.3661e-6                   Optical needle calibration [m/pixel].
FLUID_DENSITY           997.5                       Fluid density [kg/m^3].
GRAVITY                 9.78                        Local gravitational acceleration [m/s^2].
NEEDLE_DIAMETER         0.0015                      Outer needle diameter [m] (1.5 mm).
SURFACE_TENSION_REF     0.07278                     Reference tension for Worthington num [N/m].
PBLS                    np.arange(100, 101)         Stride values (PBL) to evaluate.
PBL_EVAL_MIN            100                         Start of evaluation plateau window [PBL].
PBL_EVAL_MAX            100                         End of evaluation plateau window [PBL].
TAIL_PERCENTAGE         90.0                        % of upper profile used for regression.
NUM_WORKERS             4                           CPU cores for parallel execution.
SAVE_DIAGNOSTIC_PLOTS   True                        Save individual profile PNG plots.


7. OUTPUT FILES
---------------
Upon execution, the script generates the following files in the specified OUTPUT_DIR:

1. FD_Analysis_Results[_Video_{id}].xlsx & .txt:
   Main summary report containing:
   - Profile               : Profile index number
   - Radius_m              : Equatorial radius R_0 in meters [m]
   - B_Laplace             : Dimensionless Young-Laplace slope (B)
   - Gamma_Laplace_N_per_m : Computed surface tension [N/m]
   - Shape_Parameter       : Dimensionless drop shape deviation parameter
   - Contact_Angle_deg     : Estimated contact angle [deg]
   - Drop_Volume_uL        : Calculated drop volume in microliters [uL]
   - Area_mm2              : Projected 2D profile area in square millimeters [mm^2]
   - Worthington_Number    : Dimensionless Worthington number
   - Fitted_Points         : Number of data points used in linear regression

2. FD_Standard_Deviations[_Video_{id}].xlsx:
   Standard error of the linear regression slope across evaluated PBL strides, 
   separated into 'Left_Branch' and 'Right_Branch' sheets.

3. FD_B_Branches[_Video_{id}].xlsx:
   Branch-resolved B values across all evaluated PBL strides ('Left_Branch' 
   and 'Right_Branch' sheets).

4. FD_Gamma_vs_PBLs[_Video_{id}].xlsx & FD_B_vs_PBLs[_Video_{id}].xlsx:
   Global overview tables of surface tension (N/m) and B as a function of PBL.

5. Diagnostic Figures (PNG format):
   - B_vs_PBL_Profile_{i}.png          : Slope consistency vs sampling stride.
   - Gamma_vs_PBL_Profile_{i}.png      : Surface tension vs sampling stride.
   - StdDev_vs_PBL_Profile_{i}.png     : Linear regression standard error vs stride.
   - Points_vs_PBL_Profile_{i}.png     : Number of regression points vs stride.
   - curvature_vs_Y_pbl_{stride}_...   : Direct curvature vs elevation profile fits.


8. SYSTEM REQUIREMENTS AND INSTALLATION
---------------------------------------
The software runs on Windows, Linux, and macOS with Python 3.8 or higher.

Required dependencies:
    pip install numpy scipy pandas matplotlib openpyxl

To run the analysis:
    python FD_method.py


9. CITATION AND ATTRIBUTION
---------------------------
When using this code or method in scientific publications, please cite the 
corresponding journal article and this supplemental material repository.
================================================================================
