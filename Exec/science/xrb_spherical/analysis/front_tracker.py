#!/usr/bin/env python3

import sys
import glob
import yt
import numpy as np
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from yt.frontends.boxlib.api import CastroDataset
from yt.units import cm

def track_flame_front(ds, args):
    '''
    This function tracks the flame front position for a given dataset.
    It returns a list of the form: [Time (in ms), Theta, averaged_max_field, Theta_max, max_field]

    Procedure to determine flame front:
    1) User selects a quantity to use as a metric: enuc or Temp
    2) Determine the global max of that quantity
    3) Determine the minimum value required to consider averaging zones based on global max
    4) Do a radial average of the data set to convert to 1D as a function of theta
    5) Determine flame front at theta where the radially averaged quantity drops to
       percent * averaged_max of that quantity.
    '''

    time = ds.current_time.in_units("ms")
    rr = ds.domain_right_edge[0].in_units("cm")
    rl = ds.domain_left_edge[0].in_units("cm")

    thetar = ds.domain_right_edge[1]
    thetal = ds.domain_left_edge[1]
    dtheta = thetar - thetal

    max_level = ds.index.max_level
    ref_ratio = int(np.prod(ds.ref_factors[0:max_level]))

    # Assume default resolution using finest grid
    default_res = ds.domain_dimensions * ref_ratio

    ###
    ### Select data by choosing rays at different theta
    ###

    # Get different possible cell-centered thetas
    thetas = np.linspace(thetal, thetar, default_res[1], endpoint=False) + 0.5 * dtheta / default_res[1]

    # Container for the radially averaged field quantity, enuc or temp
    averaged_field = []

    # First determine the global max of field quantity
    max_val = ds.all_data()[args.field].max()

    # Determine a threshold of selecting zones for the average, i.e. minimum value allowed
    min_val = max_val * args.threshold

    # track the theta that has the maximum global value
    max_theta_loc = 0.0

    # Loop over different thetas
    for theta in thetas:
        # simply go from (rmin, theta) -> (rmax, theta). Doesn't need to convert to physical R-Z
        ray = ds.ray((rl, theta, 0*cm), (rr, theta, 0*cm))

        # sort by "t", which goes from 0 to 1 representing the spatial order.
        # isrt = np.argsort(ray["t"])

        # Do the tracking
        if any(ray[args.field) == max_val):
            max_theta_loc = theta

        # Consider zones that are larger than minimum value
        valid_zones = ray[args.field] > min_val
        valid_values = ray[args.field][valid_zones]

        if len(valid_values) > 0:
            averaged_field.append(valid_values.mean())
        else:
            averaged_field.append(0.0)

    # Now Determine the index of the maximum radially averaged field
    max_index = np.argmax(averaged_field)

    # Now assuming flame moves forward in theta, find theta such that the field drops below some threshold of the averaged max
    loc_index = averaged_field[max_index:] < args.percent * max(averaged_field)

    # Find the first theta that the field drops below the threshold.
    theta_loc = thetas[max_index:][loc_index][0]

    # Returns 5 quantities
    # 1) time in ms
    # 2) theta that corresponds to the maximum averaged value
    # 3) maximum averaged value
    # 4) theta that corresponds to the maximum global value
    # 5) maximum global value
    timeTheta = [time, theta_loc, max(averaged_field), max_theta_loc, max_val]

    return timeTheta


def process_dataset(fname, args):
    ds = CastroDataset(fname)

    # Returns a list [time, theta, max averaged value, theta_max, max value]
    return track_flame_front(ds, args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    This file tracks the flame front and writes them into a txt file.
    """)

    parser.add_argument('--fnames', nargs='+', type=str,
                        help="Dataset file names for tracking flame front.")
    parser.add_argument('--field', '-f', default="enuc", type=str,
                        metavar="{enuc, Temp}",
                        help="""field parameter used as metric to determine flame front.
                        Choose between {enuc, Temp}""")
    parser.add_argument('--percent', '-p', default=1e-3, type=float,
                        help="""Float number between (0, 1]. Representing the percent of
                        the averaged maximum of the field quantity used to track the flame.""")
    parser.add_argument('--threshold', '-t', default=1.e-6, type=float,
                        help="""Float number between (0, 1]. Representine the percent of
                        the global maximum of the field quantity used to select valid zones
                        for averaging""")
    parser.add_argument('--jobs', '-j', default=1, type=int,
                        help="""Number of workers to process plot files in parallel"""
    parser.add_argument('--out', '-o', default="front_tracking.dat", type=str,
                        help="""Output filename for the tracking information""")

    args = parser.parse_args()

    metric_quantities = ["enuc", "Temp"]
    if args.field not in metric_quantities:
        parser.error("field must be either enuc or Temp")

    if args.percent <= 0.0 or args.percent > 1.0:
        parser.error("percent must be a float between (0, 1]")

    if args.threshold <= 0.0 or args.percent > 1.0:
        parser.error("threshold must be a float between (0, 1]")

    timeThetaArray = []

    ###
    ### Parallelize the loop. Copied from flame_wave/analysis/front_tracker.py
    ###
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
        future_to_index = {
            executor.submit(process_dataset, fname, args): i
            for i, fname in enumerate(args.fnames)
        }
        try:
            for future in as_completed(future_to_index):
                i = future_to_index.pop(future)
                try:
                    timeThetaArray.append(future.result())
                except Exception as exc:
                    print(f"{args.fnames[i]} generated an exception: {exc}", file=sys.stderr, flush=True)
        except KeyboardInterrupt:
            print(
                "\n*** got ctrl-c, cancelling remaining tasks and waiting for existing ones to finish...\n",
                flush=True,
            )
            executor.shutdown(wait=True, cancel_futures=True)
            sys.exit(1)

    # Sort array by time and write to file
    timeThetaArray.sort(key=lambda x: x[0])
    timeThetaArray = np.array(timeThetaArray)

    np.savetxt(args.out, timeThetaArray, delimiter=',')
