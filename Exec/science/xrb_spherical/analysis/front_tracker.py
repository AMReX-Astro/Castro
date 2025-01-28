#!/usr/bin/env python3

import sys
import glob
import yt
import numpy as np
from yt.frontends.boxlib.api import CastroDataset

'''
This file tracks the flame front and writes them into a txt file.

Usage: ./front_tracker.py plotFiles_*
'''

def track_flame_front(ds):
    '''
    This script tracks the flame front position for a given dataset.
    It returns a tuple of the form: (Time, Theta)
    Note: time is in milisecond.
    '''

    time = ds.current_time.in_units("ms")

    # How to determine the flame front?
    # 1) Global max temperature: this can be used to track initial
    # detonation resulted from the initial temperature perturbation
    # 2) Use vertically averaged max enuc to determine the actual flame front


    return timeTheta


if __name__ == "__main__":

    ts = sys.argv[1:]

    timeThetaArray = []
    for fname in ts:
        ds = CastroDataset(fname)

        # Get tuple in form (theta, time)
        timeTheta= track_flame_front(ds)
        timeThetaArray.append(timeTheta)

    # Sort array by time and write to file
    timeThetaArray.sort()
    timeThetaArray = np.array(timeThetaArray)

    np.savetxt('front_tracking.dat', timeThetaArray, delimiter=',')
