"""
fd_utils.py - Finite Difference Young-Laplace Drop Analysis Utilities

This module provides mathematical and geometrical utilities for analyzing
axisymmetric pendant and sessile drop profiles using the 2D-ADER Finite Difference (FD)
formulation of the Young-Laplace equation.

Main capabilities:
- High-order (5-point) numerical differentiation for profile coordinates.
- Curvature calculation via Young-Laplace formulation.
- Multi-stride (PBL: Pixels Between Lines) sweeping and linear regression.
- Drop volume integration via numerical disk slicing.
- Shape parameter and contact angle estimation.
"""

import math
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
from scipy.optimize import curve_fit, OptimizeWarning


def first_derivative_5p(x, y):
    """
    Computes the first derivative dx/dy using a 5-point central difference stencil:
    dx/dy ≈ (-x[i+2] + 8*x[i+1] - 8*x[i-1] + x[i-2]) / (12 * h)

    Parameters:
    -----------
    x : array-like
        Dependent coordinate values (horizontal position).
    y : array-like
        Evenly spaced independent coordinate values (vertical position).

    Returns:
    --------
    np.ndarray
        Approximation of dx/dy for internal points (length len(x) - 4).
    """
    x = np.asarray(x, dtype=float)
    if len(x) < 5:
        return np.array([], dtype=float)
    h = y[1] - y[0]
    return (-x[4:] + 8.0 * x[3:-1] - 8.0 * x[1:-3] + x[:-4]) / (12.0 * h)


def second_derivative_5p(x, y):
    """
    Computes the second derivative d^2x/dy^2 using a 5-point central difference stencil:
    d^2x/dy^2 ≈ (-x[i+2] + 16*x[i+1] - 30*x[i] + 16*x[i-1] - x[i-2]) / (12 * h^2)

    Parameters:
    -----------
    x : array-like
        Dependent coordinate values (horizontal position).
    y : array-like
        Evenly spaced independent coordinate values (vertical position).

    Returns:
    --------
    np.ndarray
        Approximation of d^2x/dy^2 for internal points (length len(x) - 4).
    """
    x = np.asarray(x, dtype=float)
    if len(x) < 5:
        return np.array([], dtype=float)
    h = y[1] - y[0]
    return (-x[4:] + 16.0 * x[3:-1] - 30.0 * x[2:-2] + 16.0 * x[1:-3] - x[:-4]) / (12.0 * (h**2))


def compute_laplace_curvature(x, y, pbl, offset):
    """
    Computes the mean curvature (2H) of an axisymmetric drop profile for a given
    subsampling stride (PBL) and starting index offset.

    Curvature formula:
        2H = - (x'') / (1 + x'^2)^(3/2) + 1 / (x * (1 + x'^2)^(1/2))

    Parameters:
    -----------
    x : array-like
        Normalized horizontal coordinates.
    y : array-like
        Normalized vertical coordinates.
    pbl : int
        Pixels Between Lines (stride step size for finite differences).
    offset : int
        Starting index for the interleaved subgrid (0 <= offset < pbl).

    Returns:
    --------
    laplace_curvature : np.ndarray
        Computed curvature values along the profile slice.
    h_coords : np.ndarray
        Corresponding vertical coordinates (y values).
    """
    x_sub = np.asarray(x[offset::pbl], dtype=float)
    y_sub = np.asarray(y[offset::pbl], dtype=float)

    if len(x_sub) < 5:
        return np.array([], dtype=float), np.array([], dtype=float)

    dxdy = first_derivative_5p(x_sub, y_sub)
    ddxdy = second_derivative_5p(x_sub, y_sub)

    h_coords = y_sub[2:-2]
    x_inner = x_sub[2:-2]

    denom1 = (1.0 + dxdy**2)**1.5
    denom2 = x_inner * ((1.0 + dxdy**2)**0.5)

    lap = -(ddxdy) / denom1 + 1.0 / denom2

    # Remove the boundary endpoint to mitigate edge boundary artifacts
    h_coords = h_coords[:-1]
    lap = lap[:-1]

    return lap, h_coords


def sweep_laplace(x, y, pbl):
    """
    Performs an interleaved subgrid sweep across all starting offsets (0 to pbl-1).
    Pools and sorts all resulting curvature measurements by vertical position Y.

    Parameters:
    -----------
    x : array-like
        Normalized horizontal coordinates.
    y : array-like
        Normalized vertical coordinates.
    pbl : int
        Stride size (PBL).

    Returns:
    --------
    laplace_sorted : np.ndarray
        Pooled curvature values sorted by vertical position.
    y_sorted : np.ndarray
        Sorted vertical positions corresponding to laplace_sorted.
    """
    y_list = []
    laplace_list = []

    for offset in range(pbl):
        lap, h = compute_laplace_curvature(x, y, pbl, offset)
        if len(lap) > 0:
            y_list.append(h)
            laplace_list.append(lap)

    if len(laplace_list) > 0:
        laplace_all = np.concatenate(laplace_list)
        y_all = np.concatenate(y_list)

        sort_indices = np.argsort(y_all)
        return laplace_all[sort_indices], y_all[sort_indices]

    return np.array([], dtype=float), np.array([], dtype=float)


def multi_pbls_laplace(
    X,
    Y,
    pbl_range,
    tail_percentage=90,
    output_dir=None,
    branch_name="Branch",
    profile_index=None,
    save_plots=False
):
    """
    Evaluates profile curvature across multiple PBL strides, isolates the
    linear hydrostatic region (upper tail), and extracts the slope (B).

    Parameters:
    -----------
    X : array-like
        Normalized horizontal coordinates (r / R0).
    Y : array-like
        Normalized vertical coordinates (y / R0).
    pbl_range : array-like
        Array or list of PBL (stride) values.
    tail_percentage : float, optional
        Percentage of points from the tail to retain for linear fit (default: 90).
    output_dir : str, optional
        Directory where diagnostic plots will be stored if save_plots is True.
    branch_name : str, optional
        Identifier of the profile branch ('Right' or 'Left').
    profile_index : int or str, optional
        Profile index for tracking in plot titles and filenames.
    save_plots : bool, optional
        Whether to generate diagnostic Curvature vs Y figures (default: False).

    Returns:
    --------
    slopes : list of float
        Estimated Young-Laplace slope (B) for each PBL value.
    std_devs : list of float
        Standard error of the linear regression for each PBL value.
    points : list of int
        Number of points included in the linear fit for each PBL value.
    """
    slopes = []
    std_devs = []
    points = []

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)

    for pbl in pbl_range:
        laplace_vals, yp_vals = sweep_laplace(X, Y, pbl)

        if len(yp_vals) == 0:
            slopes.append(np.nan)
            std_devs.append(np.nan)
            points.append(0)
            continue

        # Select the upper percentage of points (drop tail toward needle)
        tail_count = int(len(yp_vals) * tail_percentage / 100.0)
        start = max(0, len(yp_vals) - tail_count)
        end = len(yp_vals)

        yp_sel = yp_vals[start:end]
        laplace_sel = laplace_vals[start:end]
        num_pts = len(yp_sel)
        points.append(num_pts)

        # Minimum points required for linear regression
        if num_pts >= 3:
            slope, intercept, r_val, p_val, slope_std_err = stats.linregress(yp_sel, laplace_sel)
            fitted_vals = slope * yp_sel + intercept
        else:
            slope = np.nan
            slope_std_err = np.nan
            fitted_vals = np.full_like(yp_sel, np.nan, dtype=float)

        slopes.append(slope)
        std_devs.append(slope_std_err)

        if save_plots and output_dir is not None:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.scatter(yp_vals, laplace_vals, color='blue', s=2, alpha=0.5, label='All Curvature Points')
            ax.scatter(yp_sel, laplace_sel, color='red', s=4, label=f'Fitted Points ({tail_percentage}% tail)')
            if not np.isnan(slope):
                ax.plot(yp_sel, fitted_vals, color='green', linewidth=1.5, label=f'Linear Fit (Slope = {slope:.4f})')
            ax.plot(Y, X + 1.0, color='black', linestyle='--', alpha=0.7, label='Profile Contour (Shifted +1)')
            ax.set_xlabel('Normalized Vertical Position (Y)')
            ax.set_ylabel('Dimensionless Curvature (2H)')
            prof_info = f" - Profile {profile_index}" if profile_index is not None else ""
            ax.set_title(f'Curvature vs Y (PBL = {pbl}, {branch_name} Branch{prof_info})')
            ax.legend(loc='best')
            ax.set_ylim(1.0, 2.75)

            prof_str = f"_profile_{profile_index}" if profile_index is not None else ""
            plot_name = f'curvature_vs_Y_pbl_{pbl}_{branch_name.lower()}{prof_str}.png'
            fig.savefig(os.path.join(output_dir, plot_name), dpi=100, bbox_inches='tight')
            plt.close(fig)

    return slopes, std_devs, points


def fit_quadratic_apex(x_data, y_data, a_fixed=1.0):
    """
    Fits a parabolic model with a fixed peak/vertex value A:
        F(x) = A_fixed + B * (x - x0)^2

    Used to accurately determine the translation offset (apex coordinate x0)
    of a profile or reference circle.

    Parameters:
    -----------
    x_data : array-like
        Independent coordinate variable (Y values near apex).
    y_data : array-like
        Dependent coordinate variable (X values near apex).
    a_fixed : float, optional
        Fixed peak height (default: 1.0).

    Returns:
    --------
    B : float or None
        Curvature coefficient of the parabola.
    x0 : float or None
        Apex position (translation coordinate).
    popt : tuple or None
        Optimized parameter vector [B, x0].
    """
    def model_fixed_a(x, B, x0):
        return a_fixed + B * (x - x0)**2

    x_data = np.asarray(x_data, dtype=float)
    y_data = np.asarray(y_data, dtype=float)

    if len(x_data) < 3:
        return None, None, None

    # Estimate initial guess using quadratic polynomial fit
    try:
        poly_coeffs = np.polyfit(x_data, y_data, 2)
        c2, c1, _ = poly_coeffs
        b_init = c2
        x0_init = -c1 / (2.0 * c2) if abs(c2) > 1e-10 else float(x_data[np.argmin(np.abs(y_data - np.mean(y_data)))])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            popt, _ = curve_fit(model_fixed_a, x_data, y_data, p0=[b_init, x0_init], maxfev=10000)
            return popt[0], popt[1], popt
    except Exception:
        return None, None, None


def calculate_shape_parameter_2016(x, y, output_dir=None, profile_index=None, save_plot=False):
    """
    Computes the drop shape parameter and contact angle based on the 2016
    spherical cap comparison methodology.

    Parameters:
    -----------
    x : array-like
        Full concatenated horizontal profile coordinates (left[::-1] and right).
    y : array-like
        Full concatenated vertical profile coordinates.
    output_dir : str, optional
        Directory to save diagnostic figure if save_plot is True.
    profile_index : int or str, optional
        Profile index for tracking.
    save_plot : bool, optional
        Whether to generate contact angle comparison plot (default: False).

    Returns:
    --------
    shape_param : float
        Calculated shape parameter.
    contact_angle_deg : float
        Estimated contact angle in degrees.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Translate vertical coordinate to the apex height
    y = y - y[np.argmax(x)]

    half_len = len(x) // 2
    x_left = x[:half_len][::-1]
    y_left = y[:half_len][::-1]
    x_right = x[half_len:]
    y_right = y[half_len:]

    # Threshold: First 1% of points near the needle contact base
    threshold = 0.01
    fit_indices_left = [i for i in range(len(x_left)) if (i / len(x_left)) < threshold]
    fit_Y_left = [y_left[i] for i in fit_indices_left]
    fit_X_left = [x_left[i] for i in fit_indices_left]

    fit_indices_right = [i for i in range(len(x_right)) if (i / len(x_right)) < threshold]
    fit_Y_right = [y_right[i] for i in fit_indices_right]
    fit_X_right = [x_right[i] for i in fit_indices_right]

    coeffs_left = np.polyfit(fit_X_left, fit_Y_left, 1)
    coeffs_right = np.polyfit(fit_X_right, fit_Y_right, 1)

    slope_left = coeffs_left[0]
    slope_right = coeffs_right[0]
    contact_angle_left = math.atan(slope_left)
    contact_angle_right = math.atan(slope_right)

    contact_angle = math.pi - (abs(contact_angle_left) + abs(contact_angle_right)) / 2.0
    kx = 2.0 * math.pi * (1.0 - math.cos(contact_angle))

    # Reference unit circle
    theta = np.linspace(0, 2 * math.pi, 3000)
    x_circle_ref = np.cos(theta)
    y_circle_ref = np.sin(theta)

    if save_plot and output_dir is not None:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(x, y, label='Experimental Profile', color='red')
        ax.plot(fit_X_left, np.polyval(coeffs_left, fit_X_left), label='Left Linear Fit', color='green')
        ax.plot(fit_X_right, np.polyval(coeffs_right, fit_X_right), label='Right Linear Fit', color='blue')
        ax.plot(x_circle_ref, y_circle_ref, label='Reference Circle', color='black', linestyle='--')
        ax.set_aspect('equal')
        prof_info = f" (Profile {profile_index})" if profile_index is not None else ""
        ax.set_title(f'Profile Contact Angle Estimation{prof_info}')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.legend(loc='best')

        prof_str = f"_profile_{profile_index}" if profile_index is not None else ""
        plot_name = f'contact_angle_comparison{prof_str}.png'
        fig.savefig(os.path.join(output_dir, plot_name), dpi=100, bbox_inches='tight')
        plt.close(fig)

    min_x = min(np.min(y_left), np.min(y_right))
    min_c = min(y_circle_ref)

    if min_x < min_c:
        kn = 4.0 * math.pi
        shape_param = abs(kx - kn)
    else:
        # Approximate apex height of reference circle
        fit_indices_circ = [i for i in range(len(x_circle_ref)) if abs(x_circle_ref[i] - 1.0) < 0.0003]
        fit_Y_circ = [y_circle_ref[i] for i in fit_indices_circ]
        fit_X_circ = [x_circle_ref[i] for i in fit_indices_circ]

        _, Y0_theo, _ = fit_quadratic_apex(fit_Y_circ, fit_X_circ, a_fixed=1.0)
        Yt_apex = Y0_theo if Y0_theo is not None else 0.0

        # Approximate apex height of experimental profile
        fit_indices_exp = [i for i in range(len(x)) if abs(x[i] - 1.0) < 0.0003]
        fit_Y_exp = [y[i] for i in fit_indices_exp]
        fit_X_exp = [x[i] for i in fit_indices_exp]

        _, Y0_exp, _ = fit_quadratic_apex(fit_Y_exp, fit_X_exp, a_fixed=1.0)
        Yx_apex = Y0_exp if Y0_exp is not None else 0.0

        y_circle_shifted = y_circle_ref - Yt_apex + Yx_apex
        y_min = min(y)
        y_max = max(y)

        y_circle_cut = [yc for yc in y_circle_shifted if y_min <= yc <= y_max]
        if len(y_circle_cut) > 0:
            h = abs(max(y_circle_cut) - min(y_circle_cut))
        else:
            h = abs(max(y_circle_shifted) - min(y_circle_shifted))

        kn = 2.0 * math.pi * h
        shape_param = abs(kx - kn)

    return float(shape_param), float(np.degrees(contact_angle))


def calculate_drop_volume(x, y):
    """
    Computes the dimensionless volume of the axisymmetric drop using numerical
    disk integration:
        V = \\sum \\pi * r_i^2 * \\Delta y_i

    Parameters:
    -----------
    x : array-like
        Normalized horizontal coordinates of the drop contour (left[::-1] + right).
    y : array-like
        Normalized vertical coordinates of the drop contour.

    Returns:
    --------
    volume : float
        Dimensionless volume of the drop profile.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    volume = 0.0
    half_len = len(x) // 2
    for i in range(1, half_len):
        radius = (x[i] - x[-i]) / 2.0
        dy = y[i + 1] - y[i]
        volume += math.pi * (radius**2) * dy

    volume = -volume  # Ensures positive sign given contour orientation
    return float(volume)


def calculate_profile_area(x, y):
    """
    Computes the dimensionless cross-sectional projection area of the drop profile
    using the trapezoidal numerical integration rule:
        Area = 0.5 * \\sum |(x[i] + x[i-1]) * (y[i] - y[i-1])|

    Parameters:
    -----------
    x : array-like
        Normalized horizontal coordinates of the drop contour (left[::-1] + right).
    y : array-like
        Normalized vertical coordinates of the drop contour.

    Returns:
    --------
    area : float
        Dimensionless cross-sectional projection area of the drop.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    area = 0.5 * np.sum((x[1:] + x[:-1]) * (y[1:] - y[:-1]))
    return float(abs(area))


# Backward-compatibility aliases matching original naming
Profile_volume = calculate_drop_volume
Profile_area = calculate_profile_area
shape_parameter_2016 = calculate_shape_parameter_2016
barrido_laplace = sweep_laplace
laplace = compute_laplace_curvature
first_deriv = first_derivative_5p
sec_deriv = second_derivative_5p
fit_quadratic_fixed_A = fit_quadratic_apex
