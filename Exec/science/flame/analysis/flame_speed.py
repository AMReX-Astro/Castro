#!/usr/bin/env python3

from __future__ import print_function

import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt
import sys
import numpy as np


def get_T_profile(plotfile):

    ds = yt.load(plotfile)

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot
    srt = np.argsort(ad['x'])
    x_coord = np.array(ad['x'][srt])
    temp = np.array(ad['Temp'][srt])

    return time, x_coord, temp


def find_x_for_T(x, T, T_0=1.e9):
    """ given a profile x(T), find the x_0 that corresponds to T_0 """
    
    # our strategy here assumes that the hot ash is in the early part
    # of the profile.  We then find the index of the first point where
    # T drops below T_0
    idx = np.where(T < T_0)[0][0]
    
    T1 = T[idx-1]
    x1 = x[idx-1]

    T2 = T[idx]
    x2 = x[idx]

    slope = (x2 - x1)/(T2 - T1)
    
    return x1 + slope*(T_0 - T1)


def find_flame_width(x, T):
    """ given a profile T(x), find the flame width """

    gradT = np.max(np.abs( (T[1:] - T[:-1])/(x[1:] - x[:-1]) ))
    dT = T.max() - T.min()
    
    return dT/gradT

if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--skip", type=int, default=10,
                   help="interval between plotfiles")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")

    args = p.parse_args()

    prefix = args.plotfiles[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.plotfiles], key=int)

    t = []
    v = []
    w = []

    xpos_old = None

    for n in range(0, len(plot_nums), args.skip):

        pfile = "{}{}".format(prefix, plot_nums[n])

        time, x, T = get_T_profile(pfile)
        xpos = find_x_for_T(x, T)
        width = find_flame_width(x, T)

        if xpos_old is not None:

            # difference with the previous file to find the det speed
            # note: the corresponding time is centered in the interval
            v.append((xpos - xpos_old)/(time - time_old))
            w.append(0.5*(width + width_old))
            t.append(0.5*(time_old + time))

        time_old = time
        xpos_old = xpos
        width_old = width

    f = plt.figure()
    f.set_size_inches(7.0, 9.0)

    ax_s = f.add_subplot(211)
    ax_w = f.add_subplot(212)

    ax_s.plot(t, v)
    ax_w.plot(t, w)

    ax_s.set_yscale("log")
    ax_w.set_yscale("log")

    ax_w.set_xlabel("time (s)")
    ax_s.set_ylabel("flame speed (cm/s)")
    ax_w.set_ylabel("flame width (cm)")

    plt.savefig("speed.png")

