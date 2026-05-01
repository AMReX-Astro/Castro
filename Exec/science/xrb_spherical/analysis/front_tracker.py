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


def track_front(fname, threshold=1e-6, percent=1e-3):
    '''
    This function tracks the flame front and ash front position for a given dataset.

    Procedure to determine flame front:
    1) Use enuc as metric to determine flame front position
    2) Determine the global max of that quantity (enuc)
    3) Determine the minimum value required to consider averaging zones based on global max
    4) Do a radial average of the data set to convert to 1D as a function of theta
    5) Determine flame front theta where the radially averaged quantity drops to
       percent * averaged_max of that quantity.

    Procedure to determine ash front:
    1) Use abar as metric to determine ash front positin
    2) Select zones above Ni56 layer, and choose zones that have 4 < abar < 56
    3) Do radial average of abar to convert to 1D as a function of theta.
       If there are no valid zones for a given theta, it will be assigned with abar=4
    4) Choose ash front with the first theta that have abar = 4
    '''

    ds = CastroDataset(fname)
    time = ds.current_time.in_units("ms")
    ds.force_periodicity()

    # get level 0 covering grid
    dims = ds.domain_dimensions.copy()
    dims[2] = 1
    cg = ds.smoothed_covering_grid(level=0, left_edge=ds.domain_left_edge, dims=dims)

    # Get enuc, abar ,and theta data
    enuc  = cg["boxlib", "enuc"][:, :, 0].to_ndarray()
    abar  = cg["boxlib", "abar"][:, :, 0].to_ndarray()
    r     = cg["index", "r"][:, :, 0].to_ndarray()
    theta = cg["index", "theta"][:, :, 0].to_ndarray()
    theta_1d = theta[0, :]

    # Find the minimum threshold to include in radial average
    # Consider zones that are larger than minimum value
    # and determine the sum of the valid zone count in radial direction
    valid_enuc_zones = enuc > enuc.max() * threshold

    # This is the relative height where we reach isentropic atmosphere.
    rel_height = ds.parameters.get("problem.H_star") + 1.5 * ds.parameters.get("problem.atm_delta")
    surface_height = ds.domain_left_edge[0].in_units("cm") + rel_height*cm

    # Select valid zones for abar
    valid_abar_zones = (abar > 4.1) & (abar < 55) & (r > surface_height)
    count_valid_enuc_zones  = valid_enuc_zones.sum(axis=0)
    count_valid_abar_zones  = valid_abar_zones.sum(axis=0)

    # Find the radially averaged quantity where we only consider valid zones
    sum_enuc = np.where(valid_enuc_zones, enuc, 0.0).sum(axis=0)
    averaged_enuc = np.where(count_valid_enuc_zones > 0, sum_enuc / count_valid_enuc_zones, 0.0)

    sum_abar = np.where(valid_abar_zones, abar, 0.0).sum(axis=0)
    averaged_abar = np.where(count_valid_abar_zones > 0, sum_abar / count_valid_abar_zones, 4.0)

    # Find theta where averaged field drops to percent * averaged_max
    averaged_max_enuc = averaged_enuc.max()
    averaged_max_enuc_index = np.argmax(averaged_enuc)

    # Now assuming flame moves forward in theta, find theta such that the field drops below some threshold of the averaged max
    flame_index = averaged_enuc[averaged_max_enuc_index:] <= percent * averaged_max_enuc

    # Find the first theta that the averaged enuc drops below the threshold.
    flame_theta = theta_1d[averaged_max_enuc_index:][flame_index][0]

    # Find the last theta of the valid zone
    ash_index = averaged_abar == 4.0
    ash_theta = theta_1d[ash_index][0]

    # Returns 5 quantities
    # 1) file name of the dataset
    # 2) time in ms
    # 3) theta that corresponds to the flame front
    # 4) theta that corresponds to the maximum averaged value
    # 5) maximum averaged value
    # 6) theta that corresponds to the ash front

    tracking_data = {"fname":str(ds), "time[ms]": float(time),
                     "flame_theta": float(flame_theta),
                     "theta_max_avg":float(theta_1d[averaged_max_enuc_index]),
                     "max_avg_enuc": float(averaged_max_enuc),
                     "ash_theta": float(ash_theta)}

    return tracking_data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    This file tracks the flame front and ash front and writes them into a csv file.
    """)

    parser.add_argument('fnames', nargs='+', type=str,
                        help="Dataset file names for tracking flame front.")
    parser.add_argument('--percent', '-p', default=1e-3, type=float,
                        help="""Float number between (0, 1]. Representing the percent of
                        the averaged maximum of enuc used to track the flame.""")
    parser.add_argument('--threshold', '-t', default=1.e-6, type=float,
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

    tracking_data_list = []

    ###
    ### Parallelize the loop. Copied from flame_wave/analysis/front_tracker.py
    ###
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        future_to_index = {
            executor.submit(track_front, fname,
                            threshold=args.threshold, percent=args.percent): i
            for i, fname in enumerate(args.fnames)
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
