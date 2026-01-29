#!/usr/bin/env python


import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt
import sys
import numpy as np
from scipy.optimize import brentq

yt.funcs.mylog.setLevel(50)

class Profile:
    """read a plotfile using yt and store the 1d profile for T and enuc"""

    def __init__(self, plotfile):
        ds = yt.load(plotfile, hint="castro")

        time = float(ds.current_time)
        ad = ds.all_data()

        # Sort the ray values by 'x' so there are no discontinuities
        # in the line plot
        srt = np.argsort(ad['x'])
        x_coord = np.array(ad['x'][srt])
        temp = np.array(ad['Temp'][srt])
        enuc = np.array(ad['enuc'][srt])

        self.time = time
        self.x = x_coord
        self.T = temp
        self.enuc = enuc

    def get_dx(self):
        """ return the grid dx -- assumes uniform spacing"""
        return self.x[1] - self.x[0]

    def find_x_for_T(self, T_0=1.e9):
        """ given a profile x(T), find the x_0 that corresponds to T_0 """

        # our strategy here assumes that the hot ash is in the early
        # part of the profile.  We then fit a quartic conservative interpolant
        # through the 4 zones surrounding the desired temperature and root find to
        # find the x that corresponds to the T

        # find the index of the first point where T drops below T_0
        idx = np.where(self.T < T_0)[0][0]

        # the 4 points for our quartic will then be: idx-2, idx-1, idx, and idx+1
        T1 = self.T[idx-2]
        x1 = self.x[idx-2]

        T2 = self.T[idx-1]
        x2 = self.x[idx-1]

        T3 = self.T[idx]
        x3 = self.x[idx]

        T4 = self.T[idx+1]
        x4 = self.x[idx+1]

        dx = x2 - x1

        s1 = (T4 - 3*T3 + 3*T2 - T1)/(6*dx**3)
        s2 = (T3 - 2*T2 + T1)/(2*dx**2)
        s3 = (-5*T4 + 27*T3 - 15*T2 - 7*T1)/(24*dx)
        s4 = (-T3 + 26*T2 - T1)/24.0

        x_flame = brentq(lambda x: s1*(x-x2)**3 + s2*(x-x2)**2 + s3*(x-x2) + s4 - T_0,
                         x1, x4)

        return x_flame

    def find_flame_width(self):
        """ given a profile T(x), find the flame width """

        gradT = np.max(np.abs( (self.T[1:] - self.T[:-1])/(self.x[1:] - self.x[:-1]) ))
        dT = self.T.max() - self.T.min()

        return dT/gradT


if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--skip", type=int, default=1,
                   help="interval between plotfiles")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")

    args = p.parse_args()

    prefix = args.plotfiles[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.plotfiles], key=int)

    t = []
    v = []
    vs = []
    w = []

    x1_old = None
    x2_old = None
    x3_old = None

    for n in range(0, len(plot_nums), args.skip):

        pfile = f"{prefix}{plot_nums[n]}"
        p = Profile(pfile)

        x1 = p.find_x_for_T(T_0=1.5e9)
        x2 = p.find_x_for_T(T_0=1.e9)
        x3 = p.find_x_for_T(T_0=2.e9)
        width = p.find_flame_width()

        if x1_old is not None:

            # difference with the previous file to find the det speed
            # note: the corresponding time is centered in the interval
            v1 = (x1 - x1_old)/(p.time - time_old)
            v2 = (x2 - x2_old)/(p.time - time_old)
            v3 = (x3 - x3_old)/(p.time - time_old)
            v_avg = (v1 + v2 + v3)/3.0
            v_sigma = np.std([v1, v2, v3])

            v.append(v_avg)
            vs.append(v_sigma)

            w.append(0.5*(width + width_old))
            t.append(0.5*(time_old + p.time))

        time_old = p.time
        x1_old = x1
        x2_old = x2
        x3_old = x3
        width_old = width

    f = plt.figure()
    f.set_size_inches(7.0, 9.0)

    ax_s = f.add_subplot(211)
    ax_w = f.add_subplot(212)

    ax_s.plot(t, v)
    ax_w.plot(t, w)

    ax_s.set_yscale("symlog")
    ax_w.set_yscale("log")

    ax_w.set_xlabel("time (s)")
    ax_s.set_ylabel("flame speed (cm/s)")
    ax_w.set_ylabel("flame width (cm)")

    plt.tight_layout()
    plt.savefig("speed.png")

    for to, vo, vso, wo in zip(t, v, vs, w):
        print(f"{to:10.3g} : {vo:15.8g} +/- {vso:15.8g}  |  {wo:15.8g}")
