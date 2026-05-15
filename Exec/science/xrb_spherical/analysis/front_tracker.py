#!/usr/bin/env python

import sys
import glob
import yt
import numpy as np
import pandas as pd
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from yt.frontends.boxlib.api import CastroDataset
from yt.units import cm
from dataclasses import dataclass


def get_flame_front(theta_1d, enuc, threshold=1e-5, percent=1e-2):
    """
    Procedure to determine flame front:
    1) Use enuc as metric to determine flame front position
    2) Determine the global max of that quantity (enuc)
    3) Determine the minimum value required to consider averaging zones based on global max
    4) Do a radial average of the data set to convert to 1D as a function of theta
    5) Determine flame front theta where the radially averaged quantity drops to
       percent * averaged_max of that quantity. Default percent is 1%
    """

    # Find the minimum threshold to include in radial average
    # Consider zones that are larger than minimum value
    # and determine the sum of the valid zone count in radial direction
    valid_enuc_zones = enuc > enuc.max() * threshold
    count_valid_enuc_zones  = valid_enuc_zones.sum(axis=0)

    # Find the radially averaged enuc where we only consider valid zones
    sum_enuc = np.where(valid_enuc_zones, enuc, 0.0).sum(axis=0)
    avg_enuc = np.where(count_valid_enuc_zones > 0, sum_enuc / count_valid_enuc_zones, 0.0)

    # Find theta where averaged field drops to percent * averaged_max
    avg_max_enuc = avg_enuc.max()
    avg_max_enuc_index = np.argmax(avg_enuc)

    # Now assuming flame moves forward in theta, find theta such that the field drops below some threshold of the averaged max
    flame_index = avg_enuc[avg_max_enuc_index:] <= percent * avg_max_enuc

    # Find the first theta that the averaged enuc drops below the threshold.
    flame_theta = theta_1d[avg_max_enuc_index:][flame_index][0]

    # Also find theta where max averaged enuc is
    avg_max_enuc_theta = theta_1d[avg_max_enuc_index]

    return avg_max_enuc, flame_theta, avg_max_enuc_theta

def get_flame_tail(r, theta_1d, abar, surface_height):
    '''
    Procedure to determine flame tail:
    1) Use abar as metric. Tail is when there are still fresh (He4) fuel at the bottom
    2) Select zones that have abar < abar_theshold and below the surface height + 2000cm
    3) Do radial average of the selected abar to convert 1D as a function of theta.
       Choose the first valid zone.
    '''
    # Select valid zones for abar
    abar_threshold = 4.2
    r_threshold = surface_height + 2000*cm
    valid_abar_zones = (abar < abar_threshold) & (r < r_threshold)
    count_valid_abar_zones  = valid_abar_zones.sum(axis=0)

    # Find the radially averaged quantity where we only consider valid zones
    sum_abar = np.where(valid_abar_zones, abar, 0.0).sum(axis=0)
    averaged_abar = np.where(count_valid_abar_zones > 0, sum_abar / count_valid_abar_zones, 56.0)

    # Find the last theta of the valid zone
    flame_tail_index = averaged_abar < abar_threshold
    flame_tail = theta_1d[flame_tail_index][0]

    return flame_tail

def get_ash_front(r, theta_1d, abar, surface_height, abar_threshold=5):
    '''
    Procedure to determine ash front:
    1) Use abar as metric to determine ash front positin
    2) Select zones above Ni56 layer, and choose zones that have 4 < abar < 56
    3) Do radial average of abar to convert to 1D as a function of theta.
       If there are no valid zones for a given theta, it will be assigned with abar=4
    4) Choose ash front with the first theta that have abar = 4
    '''

    # Select valid zones for abar
    valid_abar_zones = (abar > abar_threshold) & (abar < 55) & (r > surface_height)
    count_valid_abar_zones  = valid_abar_zones.sum(axis=0)

    # Find the radially averaged quantity where we only consider valid zones
    sum_abar = np.where(valid_abar_zones, abar, 0.0).sum(axis=0)
    averaged_abar = np.where(count_valid_abar_zones > 0, sum_abar / count_valid_abar_zones, 4.0)

    # Find the last theta of the valid zone
    ash_index = averaged_abar > abar_threshold
    ash_theta = theta_1d[ash_index][-1]

    return ash_theta


def get_ocean_height(r, theta_1d, abar, flame_tail, surface_height):
    """
    Procedure to determine burned ocean height:
    1) Use abar as metric
    2) Find valid zones that have 5 < abar < 55
    3) Along each theta, find the maximum r of the valid zone to
       be our ocean height so that we have H(theta).
    4) Average over H(theta) for theta < flame_tail.
       This gives us the averaged burned ocean height
    """

    # Select valid zones for abar
    valid_abar_zones = (abar > 5) & (abar < 55)
    r_valid = np.where(valid_abar_zones, r, 0.0)

    # Get the r of the ocean top as a function of theta
    ocean_r_1d= r_valid.max(axis=0)

    # Location of max r of the valid zone
    ocean_r = ocean_r_1d.max()
    ocean_theta = theta_1d[np.argmax(ocean_r_1d)]

    # Average over ocean_r_1d for theta < flame_tail
    valid_r_1d = ocean_r_1d[theta_1d < flame_tail]

    # Relative height compared to the surface layer.
    ocean_height = np.average(valid_r_1d) - surface_height.value

    return ocean_r, ocean_theta, ocean_height

def track_front(ds, threshold=1e-5, percent=1e-2):
    '''
    This function tracks multiple quantities of the xrb flame.
    1. Flame front and flame tail in [rad]
    2. Ash front in [rad]
    3. Depth and Height of the burned ocean in [cm]
    '''

    time = ds.current_time.in_units("ms")
    ds.force_periodicity()

    # get level 0 covering grid
    dims = ds.domain_dimensions.copy()
    dims[2] = 1
    cg = ds.smoothed_covering_grid(level=0, left_edge=ds.domain_left_edge, dims=dims)

    # Get enuc, abar ,and theta data
    enuc = cg["boxlib", "enuc"][:, :, 0].to_ndarray()
    abar = cg["boxlib", "abar"][:, :, 0].to_ndarray()
    r = cg["index", "r"][:, :, 0].to_ndarray() # in cm
    theta = cg["index", "theta"][:, :, 0].to_ndarray()
    theta_1d = theta[0, :]

    # This is the relative height where we reach isentropic atmosphere.
    star_height = ds.domain_left_edge[0].in_units("cm") + ds.parameters.get("problem.H_star")*cm
    surface_height = star_height + 1.5 * ds.parameters.get("problem.atm_delta")*cm

    avg_max_enuc, flame_theta, avg_max_enuc_theta = get_flame_front(theta_1d, enuc, threshold=threshold, percent=percent)
    flame_tail = get_flame_tail(r, theta_1d, abar, surface_height)
    ocean_r, ocean_theta, ocean_height = get_ocean_height(r, theta_1d, abar, flame_tail, surface_height)

    abar_threshold = 5
    ash_theta = get_ash_front(r, theta_1d, abar, surface_height, abar_threshold=abar_threshold)

    # Find ash velocity
    # Given this theta, find r such that abar is closest to the abar threshold
    # then use this coordinate to find ash velocity of the front
    ash_theta_idx = np.argmin(np.abs(theta_1d - ash_theta))
    abar_1d = abar[:, ash_theta_idx]
    ash_r_idx = np.argmin(np.abs(abar_1d - abar_threshold))
    ash_velocity = cg["boxlib", "y_velocity"][ash_r_idx, ash_theta_idx, 0].to_ndarray()

    # Find coriolis parameter
    omega = 2.0 * np.pi / ds.parameters.get("castro.rotational_period")
    coriolis_param = 2.0*omega*np.cos(flame_theta)

    # # Find burning timescale -- using for estimating flame speed
    # t_nuc = None
    # if ("boxlib", "eint_e") in ds.field_list:
    #     # Get a characteristic position in the flame front
    #     # This can be used to get characteristic values of different quantities
    #     e_theta_idx = np.argmin(np.abs(theta_1d - avg_max_enuc_theta))
    #     enuc_1d = enuc[:, e_theta_idx]
    #     e_r_idx = np.argmin(np.abs(enuc_1d - avg_max_enuc))

    #     eint = cg["boxlib", "eint_e"][e_r_idx, e_theta_idx, 0].to_ndarray()
    #     t_nuc = eint / avg_max_enuc

    # # Now find diffusion coefficient, D
    # D = None
    # if ("boxlib", "diff_coeff") in ds.field_list:
    #     diff_coeff = cg["boxlib", "diff_coeff"][:, :, 0].to_ndarray()
    #     density = cg["boxlib", "density"][:, :, 0].to_ndarray()
    #     max_diff_coeff = diff_coeff.max()
    #     diff_coeff_cond = diff_coeff > max_diff_coeff * 1e-3
    #     valid_diff_coeff = diff_coeff[diff_coeff_cond]
    #     percentile = 99
    #     D = np.percentile(valid_diff_coeff, percentile, method="inverted_cdf", weights=density[diff_coeff_cond])

    # Find burning timescale and diffusion coefficient.
    # Consider density average over region with valid enuc zones
    enuc_mask = enuc > enuc.max() * threshold
    density = cg["boxlib", "density"][:, :, 0].to_ndarray()[enuc_mask]

    t_nuc = None
    if ("boxlib", "eint_e") in ds.field_list:
        eint = cg["boxlib", "eint_e"][:, :, 0].to_ndarray()[enuc_mask]
        avg_eint = np.average(eint, weights=density)
        avg_enuc = np.average(enuc[enuc_mask], weights=density)
        t_nuc = avg_eint / avg_enuc

    D = None
    if ("boxlib", "diff_coeff") in ds.field_list:
        diff_coeff = cg["boxlib", "diff_coeff"][:, :, 0].to_ndarray()[enuc_mask]
        D = np.average(diff_coeff)

    # Returns 10 quantities
    # 1) file name of the dataset
    # 2) time in ms
    # 3) theta that corresponds to the flame front
    # 4) theta that corresponds to the flame tail
    # 5) theta that corresponds to the ash front
    # 6) theta velocity in ash [cm/s]
    # 7) averaged ocean height in [cm]
    # 8) characteristic burning timescale
    # 9) characteristic diffusion coefficient
    # 10) coriolis parameter at flame front theta

    tracking_data = {"fname":str(ds),
                     "time[ms]": float(time),
                     "flame_theta": flame_theta,
                     "flame_tail": flame_tail,
                     "ash_theta": ash_theta,
                     "ash_velocity": ash_velocity,
                     "ocean_height": ocean_height,
                     "t_nuc": t_nuc,
                     "D": D,
                     "coriolis_param": coriolis_param}

    return tracking_data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    This file tracks the flame front and ash front and writes them into a csv file.
    """)

    parser.add_argument('fnames', nargs='+', type=str,
                        help="Dataset file names for tracking flame front.")
    parser.add_argument('--percent', '-p', default=1e-2, type=float,
                        help="""Float number between (0, 1]. Representing the percent of
                        the averaged maximum of enuc used to track the flame.""")
    parser.add_argument('--threshold', '-t', default=1.e-5, type=float,
                        help="""Float number between (0, 1]. Representing the percent of
                        the global maximum of enuc used to select valid zones
                        for averaging""")
    parser.add_argument('--jobs', '-j', default=1, type=int,
                        help="""Number of workers to process plot files in parallel""")
    parser.add_argument('--out', '-o', default="front_tracking.csv", type=str,
                        help="""Output filename for the tracking information""")

    args = parser.parse_args()

    if args.percent <= 0.0 or args.percent > 1.0:
        parser.error("percent must be a float between (0, 1]")

    if args.threshold <= 0.0 or args.percent > 1.0:
        parser.error("threshold must be a float between (0, 1]")

    ts = [CastroDataset(fname) for fname in args.fnames]
    tracking_data_list = []

    ###
    ### Parallelize the loop. Copied from flame_wave/analysis/front_tracker.py
    ###
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        future_to_index = {
            executor.submit(track_front, ds,
                            threshold=args.threshold, percent=args.percent): i
            for i, ds in enumerate(ts)
        }
        try:
            for future in as_completed(future_to_index):
                i = future_to_index.pop(future)
                try:
                    tracking_data_list.append(future.result())
                except Exception as exc:
                    print(f"{args.fnames[i]} generated an exception: {exc}", file=sys.stderr, flush=True)
        except KeyboardInterrupt:
            print(
                "\n*** got ctrl-c, cancelling remaining tasks and waiting for existing ones to finish...\n",
                flush=True,
            )
            executor.shutdown(wait=True, cancel_futures=True)
            sys.exit(1)

    # Write to file
    df = pd.DataFrame(tracking_data_list)

    # Sort by time
    df = df.sort_values("time[ms]").reset_index(drop=True)

    # Write to csv
    df.to_csv(args.out, index=False)
