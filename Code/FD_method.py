"""
FD_method.py - Finite Difference Young-Laplace Drop Analysis (2D-ADER)
======================================================================

Description:
------------
This script performs automated surface tension, shape parameter, contact angle,
and volume analysis of axisymmetric pendant or sessile drops from experimental
boundary coordinate profiles.

The mathematical approach uses the 2D-ADER Finite Difference (FD) Young-Laplace
method to evaluate drop surface curvature without iterative non-linear ODE shooting,
providing robust, high-accuracy surface tension measurements across user-defined
stride samplings (PBL: Pixels Between Lines).

Designed for supplemental material repositories and scientific publication:
- Fully self-documenting with comprehensive English comments and type annotations.
- Modularized architecture using the companion 'fd_utils.py' module.
- Centralized configuration block for directories, physical constants, and parameters.
- Multi-core parallel processing (multiprocessing.Pool) for rapid batch evaluation.
- Outputs structured Excel (.xlsx) and text (.txt) reports along with publication-quality plots.

Dependencies:
-------------
numpy, pandas, scipy, matplotlib, openpyxl, fd_utils
"""

import os
import time
from multiprocessing import Pool, Lock
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fd_utils import (
    multi_pbls_laplace,
    calculate_drop_volume,
    calculate_profile_area,
    calculate_shape_parameter_2016,
)

# ==============================================================================
# CONFIGURATION AND EXPERIMENTAL PARAMETERS
# ==============================================================================
CONFIG = {
    # -------------------------------------------------------------------------
    # Directory and File Settings
    # -------------------------------------------------------------------------
    # If the input directory path contains '{video}', it will be formatted with
    # each entry in the VIDEOS list. Otherwise, a single directory is processed.
    "INPUT_DIR": r"G:\Mi unidad\Artículo\Septiembre\Video {video}",
    "OUTPUT_DIR": r"G:\Mi unidad\Artículo\Septiembre\Video {video}",

    # Video or experimental series indices to process.
    # Set to [None] or a list of specific video IDs (e.g., [2], [0, 1, 2]).
    "VIDEOS": [2],

    # Range of profile indices to analyze (e.g., np.arange(0, 2) or np.arange(0, 100))
    "PROFILE_INDEXES": np.arange(0, 2),

    # -------------------------------------------------------------------------
    # Physical Constants and Optical Calibration
    # -------------------------------------------------------------------------
    "PIXEL_SIZE": 1.3661e-6,      # [m/pixel] Optical calibration factor (needle scale)
    "FLUID_DENSITY": 997.5,       # [kg/m^3] Density of liquid phase (e.g., pure water at 20°C)
    "GRAVITY": 9.78,              # [m/s^2] Local gravitational acceleration
    "NEEDLE_DIAMETER": 0.0015,    # [m] Outer needle diameter (1.5 mm, used for Worthington number)
    "SURFACE_TENSION_REF": 0.07278, # [N/m] Reference surface tension (used for Worthington number)

    # -------------------------------------------------------------------------
    # Finite Difference Analysis Settings
    # -------------------------------------------------------------------------
    # Array of PBL (stride) values to evaluate.
    # Typical multi-PBL range: np.arange(10, 151, 5) or single stride: np.arange(100, 101)
    "PBLS": np.arange(100, 101),

    # Plateau evaluation window: the PBL range averaged to report final B & Gamma
    "PBL_EVAL_MIN": 100,
    "PBL_EVAL_MAX": 100,

    # Percentage of profile points (from apex toward needle base) used for linear fit
    "TAIL_PERCENTAGE": 90.0,

    # -------------------------------------------------------------------------
    # Execution and Visualization Settings
    # -------------------------------------------------------------------------
    "NUM_WORKERS": 4,             # Number of parallel CPU workers
    "SAVE_DIAGNOSTIC_PLOTS": True # Generate and save diagnostic PNG figures
}

# Thread/process lock for synchronized file writing
write_lock = Lock()


def process_single_profile(profile_index, directory, out_directory, config):
    """
    Processes a single coordinate profile using the 2D-ADER Finite Difference Young-Laplace method.

    Workflow:
    1. Loads coordinate data from 'Profile_{profile_index}.txt' or 'Perfil_{profile_index}.txt'.
    2. Identifies equatorial radius (R0) and normalizes coordinate axes.
    3. Evaluates curvature across specified PBL values for left and right drop branches.
    4. Computes dimensionless shape parameter B and surface tension Gamma (N/m).
    5. Calculates drop volume, contact angle, shape parameter, and Worthington number.
    6. Generates diagnostic plots (B, Gamma, StdDev, Points vs PBL).

    Parameters:
    -----------
    profile_index : int
        Index number of the drop profile.
    directory : str
        Path to input directory containing profile coordinate files.
    out_directory : str
        Path to output directory for saving plots and results.
    config : dict
        Configuration dictionary containing physical parameters and settings.

    Returns:
    --------
    dict or None
        Structured analysis results dictionary, or None if processing failed.
    """
    # Locate profile file (supports English or Spanish naming conventions)
    coords_file = None
    for candidate_name in [f"Profile_{profile_index}.txt", f"Perfil_{profile_index}.txt"]:
        candidate_path = os.path.join(directory, candidate_name)
        if os.path.exists(candidate_path):
            coords_file = candidate_path
            break

    if coords_file is None:
        print(f"[Warning] Profile file for index {profile_index} not found in {directory}")
        return None

    try:
        # Load profile coordinates
        # Expected format: [xr, xl, yr, yl] or [xr, xl, y]
        coordinates = np.loadtxt(coords_file)
        xr = coordinates[:, 0]
        xl = coordinates[:, 1]
        yr = coordinates[:, 2]
        yl = coordinates[:, 2]

        # Calculate equatorial radius in pixels
        r_pixels = (xr[np.argmax(xr)] - xl[np.argmin(xl)]) / 2.0

        # Center horizontal coordinates and normalize by equatorial radius
        center = (np.max(xr) + np.min(xl)) / 2.0
        xr_norm = (xr - center) / r_pixels
        xl_norm = (xl - center) / r_pixels
        yr_norm = (yr - np.min(yr)) / r_pixels
        yl_norm = (yl - np.min(yl)) / r_pixels

        # Extract configuration constants
        pbls = config["PBLS"]
        pixel_size = config["PIXEL_SIZE"]
        density = config["FLUID_DENSITY"]
        gravity = config["GRAVITY"]
        tail_pct = config["TAIL_PERCENTAGE"]
        save_plots = config["SAVE_DIAGNOSTIC_PLOTS"]
        pbl_min = config["PBL_EVAL_MIN"]
        pbl_max = config["PBL_EVAL_MAX"]
        needle_d = config["NEEDLE_DIAMETER"]
        gamma_ref = config["SURFACE_TENSION_REF"]

        # Capillary scaling factor for surface tension:
        # Gamma = (Delta_rho * g * R0^2) / B
        capillary_factor = density * gravity * ((r_pixels * pixel_size) ** 2)

        # ---------------------------------------------------------------------
        # Curvature evaluation for Right and Left profile branches
        # ---------------------------------------------------------------------
        slopes_r, std_devs_r, points_r = multi_pbls_laplace(
            xr_norm, yr_norm, pbls,
            tail_percentage=tail_pct,
            output_dir=out_directory,
            branch_name="Right",
            profile_index=profile_index,
            save_plots=save_plots
        )

        slopes_l, std_devs_l, points_l = multi_pbls_laplace(
            -xl_norm, yl_norm, pbls,
            tail_percentage=tail_pct,
            output_dir=out_directory,
            branch_name="Left",
            profile_index=profile_index,
            save_plots=save_plots
        )

        # Calculate branch-specific and average surface tensions
        gamma_r = [capillary_factor / slope if (slope and not np.isnan(slope)) else np.nan for slope in slopes_r]
        gamma_l = [capillary_factor / slope if (slope and not np.isnan(slope)) else np.nan for slope in slopes_l]
        av_slope = [(sr + sl) / 2.0 for sr, sl in zip(slopes_r, slopes_l)]
        av_gamma = [capillary_factor / s if (s and not np.isnan(s)) else np.nan for s in av_slope]

        # Isolate the evaluation plateau window
        eval_indices = np.where((pbls >= pbl_min) & (pbls <= pbl_max))[0]
        if len(eval_indices) > 0:
            b_eval = [av_slope[i] for i in eval_indices]
            gamma_eval = [av_gamma[i] for i in eval_indices]
            points_eval = [(points_r[i] + points_l[i]) / 2.0 for i in eval_indices]
        else:
            b_eval = av_slope
            gamma_eval = av_gamma
            points_eval = [(pr + pl) / 2.0 for pr, pl in zip(points_r, points_l)]

        b_final = float(np.nanmean(b_eval))
        gamma_final = float(np.nanmean(gamma_eval))
        points_final = float(np.nanmean(points_eval))

        # ---------------------------------------------------------------------
        # Geometrical Analysis: Shape Parameter, Contact Angle, and Volume
        # ---------------------------------------------------------------------
        xx_full = np.concatenate((xl_norm[::-1], xr_norm))
        yx_full = np.concatenate((yl_norm[::-1], yr_norm))

        shape_param, contact_angle = calculate_shape_parameter_2016(
            xx_full, yx_full,
            output_dir=out_directory,
            profile_index=profile_index,
            save_plot=save_plots
        )

        # Volume calculation in cubic meters (scaled by R0^3)
        volume_m3 = calculate_drop_volume(xx_full, yx_full) * ((r_pixels * pixel_size) ** 3)
        volume_uL = volume_m3 * 1e9  # 1 m^3 = 10^9 uL

        # Worthington number: Wo = (Delta_rho * g * V) / (pi * gamma_ref * d_needle)
        worthington_num = (density * gravity * volume_m3) / (np.pi * gamma_ref * needle_d)

        # Profile projected cross-sectional area (scaled by R0^2)
        # Area in square meters: area_norm * (r_pixels * pixel_size)**2
        # Area in square millimeters (1 m^2 = 1e6 mm^2)
        area_m2 = calculate_profile_area(xx_full, yx_full) * ((r_pixels * pixel_size) ** 2)
        area_mm2 = area_m2 * 1e6

        # ---------------------------------------------------------------------
        # Diagnostic Figures Generation
        # ---------------------------------------------------------------------
        if save_plots:
            # 1. Slope (B) vs PBL
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(pbls, slopes_r, marker='o', color='black', label='Right Branch')
            ax.plot(pbls, slopes_l, marker='x', color='red', label='Left Branch')
            ax.plot(pbls, av_slope, marker='d', color='green', label='Average B')
            if len(eval_indices) > 0:
                ax.plot(pbls[eval_indices], b_eval, marker='s', color='blue',
                        label=f'Evaluation Window ({pbl_min}-{pbl_max} PBLs)')
            ax.set_xlabel('Pixels Between Lines (PBL)')
            ax.set_ylabel('Slope (B)')
            ax.set_xscale('log')
            ax.set_title(f'B vs PBL - Young-Laplace (Profile {profile_index})')
            ax.legend(loc='best')
            fig.savefig(os.path.join(out_directory, f'B_vs_PBL_Profile_{profile_index}.png'),
                        dpi=100, bbox_inches='tight')
            plt.close(fig)

            # 2. Surface Tension (Gamma) vs PBL
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.scatter(pbls, gamma_r, color='black', label='Right Branch', s=15)
            ax.scatter(pbls, gamma_l, color='red', label='Left Branch', s=15)
            ax.scatter(pbls, av_gamma, color='green', label='Average Gamma', s=15)
            if len(eval_indices) > 0:
                ax.scatter(pbls[eval_indices], gamma_eval, color='blue',
                           label=f'Evaluation Window ({pbl_min}-{pbl_max} PBLs)', s=20)
            ax.set_xlabel('Pixels Between Lines (PBL)')
            ax.set_ylabel('Surface Tension (N/m)')
            ax.set_xscale('log')
            ax.set_title(f'Surface Tension vs PBL (Profile {profile_index})')
            ax.legend(loc='best')
            fig.savefig(os.path.join(out_directory, f'Gamma_vs_PBL_Profile_{profile_index}.png'),
                        dpi=100, bbox_inches='tight')
            plt.close(fig)

            # 3. Standard Deviation vs PBL
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.scatter(pbls, std_devs_r, color='black', label='Right Branch', s=15)
            ax.scatter(pbls, std_devs_l, color='red', label='Left Branch', s=15)
            ax.set_xlabel('Pixels Between Lines (PBL)')
            ax.set_ylabel('Standard Deviation')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_title(f'Standard Deviation vs PBL (Profile {profile_index})')
            ax.legend(loc='best')
            fig.savefig(os.path.join(out_directory, f'StdDev_vs_PBL_Profile_{profile_index}.png'),
                        dpi=100, bbox_inches='tight')
            plt.close(fig)

            # 4. Number of Points vs PBL
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.scatter(pbls, points_r, color='black', label='Right Branch', s=15)
            ax.scatter(pbls, points_l, color='red', label='Left Branch', s=15)
            ax.scatter(pbls, [(pr + pl) / 2.0 for pr, pl in zip(points_r, points_l)],
                       color='green', label='Average Points', s=15)
            ax.set_xlabel('Pixels Between Lines (PBL)')
            ax.set_ylabel('Number of Points')
            ax.set_xscale('log')
            ax.set_title(f'Number of Points vs PBL (Profile {profile_index})')
            ax.legend(loc='best')
            fig.savefig(os.path.join(out_directory, f'Points_vs_PBL_Profile_{profile_index}.png'),
                        dpi=100, bbox_inches='tight')
            plt.close(fig)

        # ---------------------------------------------------------------------
        # Package and return structured profile results
        # ---------------------------------------------------------------------
        return {
            'result': {
                'Profile': profile_index,
                'Radius_m': r_pixels * pixel_size,
                'B_Laplace': b_final,
                'Gamma_Laplace_N_per_m': gamma_final,
                'Shape_Parameter': shape_param,
                'Contact_Angle_deg': contact_angle,
                'Drop_Volume_uL': volume_uL,
                'Area_mm2': area_mm2,
                'Worthington_Number': worthington_num,
                'Fitted_Points': points_final,
            },
            'std_devs_left': dict(zip(pbls, std_devs_l)),
            'std_devs_right': dict(zip(pbls, std_devs_r)),
            'b_left': dict(zip(pbls, slopes_l)),
            'b_right': dict(zip(pbls, slopes_r)),
            'pbls': pbls,
            'Gamma_average': av_gamma,
            'B_average': av_slope
        }

    except Exception as err:
        print(f"[Error] Failed to process profile {profile_index}: {err}")
        import traceback
        traceback.print_exc()
        return None


def _worker_wrapper(args):
    """
    Multiprocessing worker wrapper function.
    Unpacks arguments and invokes process_single_profile.
    """
    profile_index, directory, out_directory, config = args
    result = process_single_profile(profile_index, directory, out_directory, config)
    return result, profile_index


def run_pipeline(directory, out_directory, config, video_label=None):
    """
    Executes the batch processing pipeline across all designated profiles in parallel.
    Compiles summary tables and exports them into Excel and TXT formats.

    Parameters:
    -----------
    directory : str
        Input directory containing profile coordinate files.
    out_directory : str
        Output destination directory.
    config : dict
        Pipeline configuration dictionary.
    video_label : str or int, optional
        Label used in output filenames (e.g., video ID).
    """
    start_time = time.time()
    profile_indexes = config["PROFILE_INDEXES"]
    num_workers = config["NUM_WORKERS"]
    video_suffix = f"_Video_{video_label}" if video_label is not None else ""

    print("\n" + "=" * 70)
    print(f"FD Method Analysis (2D-ADER Young-Laplace){video_suffix}")
    print(f"Profiles to evaluate : {len(profile_indexes)}")
    print(f"Input directory      : {directory}")
    print(f"Output directory     : {out_directory}")
    print(f"Parallel workers     : {num_workers}")
    print("=" * 70 + "\n")

    os.makedirs(out_directory, exist_ok=True)
    worker_args = [(idx, directory, out_directory, config) for idx in profile_indexes]

    results_dict = {}
    std_devs_left_list = []
    std_devs_right_list = []
    b_left_list = []
    b_right_list = []
    gamma_avg_list = []
    b_avg_list = []
    ordered_indices = []
    pbls_global = None

    with Pool(num_workers) as pool:
        processed_count = 0
        total_profiles = len(profile_indexes)

        for result_data, profile_idx in pool.imap_unordered(_worker_wrapper, worker_args):
            processed_count += 1
            if result_data is not None:
                results_dict[profile_idx] = result_data['result']
                ordered_indices.append(profile_idx)
                std_devs_left_list.append(result_data['std_devs_left'])
                std_devs_right_list.append(result_data['std_devs_right'])
                b_left_list.append(result_data['b_left'])
                b_right_list.append(result_data['b_right'])
                gamma_avg_list.append(result_data['Gamma_average'])
                b_avg_list.append(result_data['B_average'])

                if pbls_global is None:
                    pbls_global = result_data['pbls']

                res = result_data['result']
                print(f"  [{processed_count}/{total_profiles}] [OK] Profile {profile_idx} completed")
                print(f"      Radius             : {res['Radius_m']:.6e} m")
                print(f"      B Laplace          : {res['B_Laplace']:.6f}")
                print(f"      Gamma Laplace      : {res['Gamma_Laplace_N_per_m']:.6f} N/m")
                print(f"      Shape Parameter    : {res['Shape_Parameter']:.3f}")
                print(f"      Contact Angle      : {res['Contact_Angle_deg']:.2f} deg")
                print(f"      Drop Volume        : {res['Drop_Volume_uL']:.6e} uL")
                print(f"      Profile Area       : {res['Area_mm2']:.4f} mm^2")
                print(f"      Worthington Number : {res['Worthington_Number']:.6f}")
                print(f"      Fitted Points      : {res['Fitted_Points']:.1f}")
            else:
                print(f"  [{processed_count}/{total_profiles}] [FAIL] Profile {profile_idx} - Processing error\n")

    # Sort results by profile index
    sorted_profiles = sorted(results_dict.keys())
    results_list = [results_dict[p] for p in sorted_profiles]

    if not results_list:
        print("[Warning] No profiles were successfully processed.")
        return

    # -------------------------------------------------------------------------
    # 1. Export Main Summary Results (Excel & TXT)
    # -------------------------------------------------------------------------
    with write_lock:
        df_results = pd.DataFrame(results_list)
        excel_out = os.path.join(out_directory, f"FD_Analysis_Results{video_suffix}.xlsx")
        txt_out = os.path.join(out_directory, f"FD_Analysis_Results{video_suffix}.txt")

        df_results.to_excel(excel_out, index=False)
        with open(txt_out, 'w') as f:
            f.write(df_results.to_string(index=False))

        print("\n" + "=" * 70)
        print(f"Main results saved to : {excel_out}")
        print(f"Text summary saved to : {txt_out}")
        print(f"Total profiles saved  : {len(results_list)}")
        print("=" * 70)
        print("\nSummary Table:")
        print(df_results.to_string(index=False))

    # -------------------------------------------------------------------------
    # 2. Export Standard Deviations per Branch across PBLs
    # -------------------------------------------------------------------------
    if std_devs_left_list and std_devs_right_list and pbls_global is not None:
        with write_lock:
            rows_sd_left = []
            rows_sd_right = []
            for p_idx in sorted_profiles:
                idx_pos = ordered_indices.index(p_idx)
                row_l = {'Profile': p_idx}
                row_r = {'Profile': p_idx}
                for pbl in pbls_global:
                    row_l[f"PBL_{pbl}"] = std_devs_left_list[idx_pos].get(pbl, np.nan)
                    row_r[f"PBL_{pbl}"] = std_devs_right_list[idx_pos].get(pbl, np.nan)
                rows_sd_left.append(row_l)
                rows_sd_right.append(row_r)

            df_sd_left = pd.DataFrame(rows_sd_left)
            df_sd_right = pd.DataFrame(rows_sd_right)
            sd_file = os.path.join(out_directory, f"FD_Standard_Deviations{video_suffix}.xlsx")

            with pd.ExcelWriter(sd_file, engine='openpyxl') as writer:
                df_sd_left.to_excel(writer, sheet_name='Left_Branch', index=False)
                df_sd_right.to_excel(writer, sheet_name='Right_Branch', index=False)

            print(f"\nStandard deviations saved to: {sd_file}")

    # -------------------------------------------------------------------------
    # 3. Export B per Branch across PBLs
    # -------------------------------------------------------------------------
    if b_left_list and b_right_list and pbls_global is not None:
        with write_lock:
            rows_b_left = []
            rows_b_right = []
            for p_idx in sorted_profiles:
                idx_pos = ordered_indices.index(p_idx)
                row_l = {'Profile': p_idx}
                row_r = {'Profile': p_idx}
                for pbl in pbls_global:
                    row_l[f"PBL_{pbl}"] = b_left_list[idx_pos].get(pbl, np.nan)
                    row_r[f"PBL_{pbl}"] = b_right_list[idx_pos].get(pbl, np.nan)
                rows_b_left.append(row_l)
                rows_b_right.append(row_r)

            df_b_left = pd.DataFrame(rows_b_left)
            df_b_right = pd.DataFrame(rows_b_right)
            b_file = os.path.join(out_directory, f"FD_B_Branches{video_suffix}.xlsx")

            with pd.ExcelWriter(b_file, engine='openpyxl') as writer:
                df_b_left.to_excel(writer, sheet_name='Left_Branch', index=False)
                df_b_right.to_excel(writer, sheet_name='Right_Branch', index=False)

            print(f"Branch-specific B values saved to: {b_file}")

    # -------------------------------------------------------------------------
    # 4. Export Gamma & B vs PBLs
    # -------------------------------------------------------------------------
    if gamma_avg_list and b_avg_list and pbls_global is not None:
        with write_lock:
            rows_gamma = []
            rows_b = []
            for p_idx in sorted_profiles:
                idx_pos = ordered_indices.index(p_idx)
                g_vals = gamma_avg_list[idx_pos]
                b_vals = b_avg_list[idx_pos]
                row_g = {'Profile': p_idx}
                row_b = {'Profile': p_idx}
                for pbl, gamma, b_val in zip(pbls_global, g_vals, b_vals):
                    row_g[f"Gamma_{pbl}"] = gamma
                    row_b[f"B_{pbl}"] = b_val
                rows_gamma.append(row_g)
                rows_b.append(row_b)

            df_gamma_pbls = pd.DataFrame(rows_gamma)
            gamma_file = os.path.join(out_directory, f"FD_Gamma_vs_PBLs{video_suffix}.xlsx")
            df_gamma_pbls.to_excel(gamma_file, index=False)

            df_b_pbls = pd.DataFrame(rows_b)
            b_pbls_file = os.path.join(out_directory, f"FD_B_vs_PBLs{video_suffix}.xlsx")
            df_b_pbls.to_excel(b_pbls_file, index=False)

            print(f"Gamma vs PBLs saved to: {gamma_file}")
            print(f"B vs PBLs saved to    : {b_pbls_file}")

    # Execution timing summary
    elapsed = time.time() - start_time
    hrs = int(elapsed // 3600)
    mins = int((elapsed % 3600) // 60)
    secs = elapsed % 60
    print("\n" + "=" * 70)
    print(f"Total Execution Time : {hrs}h {mins}m {secs:.2f}s")
    print("=" * 70 + "\n")


def main():
    """
    Main entry point for batch processing.
    """
    videos = CONFIG.get("VIDEOS")

    if videos and len(videos) > 0 and videos[0] is not None:
        for vid in videos:
            in_dir = CONFIG["INPUT_DIR"].format(video=vid)
            out_dir = CONFIG["OUTPUT_DIR"].format(video=vid)
            run_pipeline(in_dir, out_dir, CONFIG, video_label=vid)
    else:
        in_dir = CONFIG["INPUT_DIR"]
        out_dir = CONFIG["OUTPUT_DIR"]
        run_pipeline(in_dir, out_dir, CONFIG, video_label=None)


if __name__ == "__main__":
    main()
