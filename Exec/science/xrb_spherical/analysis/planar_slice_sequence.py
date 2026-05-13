#!/usr/bin/env python

import sys
import argparse
import numpy as np
import pandas as pd
from yt.frontends.boxlib.api import CastroDataset
from planar_slice import planar_slice
from concurrent.futures import ProcessPoolExecutor, as_completed

parser = argparse.ArgumentParser(description="""
This script uses the planar_slice.py to create a
time sequence of planar slice plots""")

parser.add_argument('fnames', nargs='+', type=str,
                    help="""dataset file names for time sequence plotting.""")
parser.add_argument('-f', '--fields', nargs='+', type=str,
                    help="""field parameters for plotting. Accepts one or more datasets.
                    If multiple parameters are given, a grid of slice plots of different
                    field parameters will be plotted for a given fname.
                    Note that either fnames or fields must be single valued.
                    """)
parser.add_argument("--figsize", nargs=2, type=float, default=[16, 9],
                    metavar=("WIDTH", "HEIGHT"), help="Figure size in inches.")
parser.add_argument("--xmin", type=float, default=None, metavar="THETA",
                    help="Minimum theta for plot xlim")
parser.add_argument("--xmax", type=float, default=None, metavar="THETA",
                    help="Maximum theta for plot xlim")
parser.add_argument("--ymin", type=float, default=None, metavar="R",
                    help="Minimum r [km] for plot ylim")
parser.add_argument("--ymax", type=float, default=None, metavar="R",
                    help="Maximum r [km] for plot ylim")
parser.add_argument("--contour-field", default=None,
                    help="Field variable to use for overplotting contour lines (e.g. 'pressure').")
parser.add_argument("--overplot-fine-levels", action="store_true",
                    help="""Overplot finer AMR level data on top of level 0.
                    This can make plotting a lot slower""")
parser.add_argument("--annotate-front", action="store_true",
                    help="Overlay vertical lines to indicate flame front and ash front.")
parser.add_argument("--annotate-velocity-streamlines", action="store_true",
                    help="Overlay velocity streamlines.")
parser.add_argument("--annotate-acceleration-streamlines", action="store_true",
                    help="Overlay acceleration streamlines.")
parser.add_argument("--annotate-grids", action="store_true",
                    help="Overlay AMR grid boundaries colored by level.")
parser.add_argument("--time-interval", type=float, default=None,
                     help="Time interval [ms] at which the plot file are divisible by.")
parser.add_argument('--jobs', '-j', default=1, type=int,
                    help="""Number of workers to plot in parallel""")

args = parser.parse_args()

fnames = args.fnames
if args.time_interval is not None:
    fnames = [fname for fname in fnames
              if np.isclose(CastroDataset(fname).current_time.in_units("ms").value % args.time_interval, 0.0,
                            rtol=0, atol=1e-3)] # only use atol

# Parallelize the plotting
with ProcessPoolExecutor(max_workers=args.jobs) as executor:
    future_to_index = {
        executor.submit(planar_slice, [fname], args.fields,
                        figsize=args.figsize,
                        xmin=args.xmin, xmax=args.xmax,
                        ymin=args.ymin, ymax=args.ymax,
                        contour_field=args.contour_field,
                        overplot_fine_levels=args.overplot_fine_levels,
                        annotate_front=args.annotate_front,
                        annotate_velocity_streamlines=args.annotate_velocity_streamlines,
                        annotate_acceleration_streamlines=args.annotate_acceleration_streamlines,
                        annotate_grids=args.annotate_grids,
                        outName=fname+"_planar_slice.png"): i
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
