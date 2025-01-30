#!/usr/bin/env python3

import sys
import glob
import yt
import numpy as np
from yt.frontends.boxlib.api import CastroDataset


def track_flame_front(ds, metric):
    '''
    This function tracks the flame front position for a given dataset.
    It returns a tuple of the form: (Time (in ms), Theta)

    Procedure to determine flame front:
    1) User selects e a quantity to use as a metric: enuc or Temp
    2) Determine the global max of that quantity
    3) Determine the theta position of the cell that contains the global max
    4) Do a radial average of the data set to convert to 1D as a function of theta
    5) Determine flame front at theta where the metric quantity drops to
       metric_number * global_max of that quantity.
    '''

    time = ds.current_time.in_units("ms")


    ad = ds.all_data()


    # 1) Global max temperature: this can be used to track initial
    # detonation resulted from the initial temperature perturbation
    # 2) Use radially averaged enuc then find flame front is determined by
    # when the averaged enuc drops below some percentage of max global enuc


    return timeTheta


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
    parser.add_argument('--percent', '-p', default=0.001, type=float,
                        help="""Float number between (0, 1]. Representing the percent of
                        the global maximum of the field quantity used to track the flame.""")

    args = parser.parse_args()

    metric_quantities = ["enuc", "Temp"]
    if args.field not in metric_quantities:
        parser.error("field must be either enuc or Temp")

    if args.percent <= 0.0 or args.percent > 1.0:
        parser.error("metric must be a float between (0, 1]")

    metric = (args.field, args.percent)
    timeThetaArray = []

    for fname in args.fnames:
        ds = CastroDataset(fname)

        # Get tuple in form (theta, time)
        timeTheta= track_flame_front(ds, metric)
        timeThetaArray.append(timeTheta)

    # Sort array by time and write to file
    timeThetaArray.sort()
    timeThetaArray = np.array(timeThetaArray)

    np.savetxt('front_tracking.dat', timeThetaArray, delimiter=',')
