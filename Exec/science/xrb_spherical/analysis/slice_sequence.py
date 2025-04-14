#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pandas as pd
from slice import slice
from concurrent.futures import ProcessPoolExecutor, as_completed

parser = argparse.ArgumentParser(description="""
This script uses the front_tracking.dat from front_tracker.py
and slice.py to create a sequence of slice plots along with
flame front position.""")

parser.add_argument('tracking_fname', type=str,
                   help="txt file generated from front_tracker.py to track flame front position.")
parser.add_argument('-f', '--fields', nargs='+', type=str,
                    help="field parameters for plotting, e.g. enuc abar.")
parser.add_argument('-w', '--width', default=4.0, type=float,
                    help="scaling for the domain width of the slice plot")
parser.add_argument('-r', '--dr', default=0.15, type=float,
                    help="""Distance between upper r and lower r shown in the SlicePlot.
                    Assumed in unit km. This is used to control center and width of the SlicePlot""")
parser.add_argument('--jobs', '-j', default=1, type=int,
                    help="""Number of workers to plot in parallel""")

args = parser.parse_args()

# data has columns: fname, time, front_theta, theta_max_avg, max_avg, theta_max, max_val.
# See front_tracker.py for more info
tracking_data = pd.read_csv(args.tracking_fname)

# Get file name and theta of flame front
fnames = tracking_data["fname"]
front_thetas = tracking_data["front_theta"]

# Parallelize the plotting
with ProcessPoolExecutor(max_workers=args.jobs) as executor:
    future_to_index = {
        executor.submit(slice, [fname], args.fields, widthScale=args.width
                        dr=args.dr, theta=front_thetas[i]): i
        for i, fname in enumerate(fnames)
    }
    try:
        for future in as_completed(future_to_index):
            i = future_to_index.pop(future)
            try:
                future.result()
            except Exception as exc:
                print(f"{fnames[i]} generated an exception: {exc}", file=sys.stderr, flush=True)
    except KeyboardInterrupt:
        print(
            "\n*** got ctrl-c, cancelling remaining tasks and waiting for existing ones to finish...\n",
            flush=True,
        )
        executor.shutdown(wait=True, cancel_futures=True)
        sys.exit(1)
