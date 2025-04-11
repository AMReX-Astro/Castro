#!/usr/bin/env python3

from __future__ import print_function

import argparse
import sys

import numpy as np
from cycler import cycler
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import yt
from yt.frontends.boxlib.api import CastroDataset


def get_T_profile(plotfile):

    ds = CastroDataset(plotfile)

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot

    dimension = ds.dimensionality
    coords = {"cartesian":"x",
              "cylindrical":"z",
              "spherical":"r"}

    coord = coords[ds.geometry]

    problo = ds.domain_left_edge
    probhi = ds.domain_right_edge

    srt = np.argsort(ad[coord])
    x_coord = np.array(ad[coord][srt])
    temp = np.array(ad['Temp'][srt])
    analytic_temp = np.array(ad['analytic'][srt])

    return time, x_coord, temp, analytic_temp


def doit(prefix, nums, skip, limitlabels, xmin, xmax):

    f = plt.figure()
    f.set_size_inches(7.0, 7.0)

    ax_T = f.add_subplot(111)

    # Get colors
    numplots = len(range(0, len(nums), skip))
    colors = plt.cm.nipy_spectral(np.linspace(0, 1, numplots))

    if limitlabels > 1:
        skiplabels = int(numplots / limitlabels)
    elif limitlabels < 0:
        print("Illegal value for limitlabels: %.0f" % limitlabels)
        sys.exit()
    else:
        skiplabels = 1
    index = 0

    for n in range(0, len(nums), skip):
        pfile = "{}{}".format(prefix, nums[n])
        i = int(n / skip)
        time, x, T, analytic_T = get_T_profile(pfile)

        if index % skiplabels == 0:
            ax_T.plot(x, T, "-", color=colors[i], label="simulation: t = {:7.5f} s".format(time))
            ax_T.plot(x, analytic_T, "--", color=colors[i], label="analytic: t = {:7.5f} s".format(time))
        else:
            ax_T.plot(x, T, "-", color=colors[i])
            ax_T.plot(x, analytic_T, "--", color=colors[i])

        index += 1

    ax_T.legend(frameon=False)
    ax_T.set_ylabel("T (K)")

    if xmax > 0:
        ax_T.set_xlim(xmin, xmax)

    f.savefig("diffuse.png")


if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--skip", type=int, default=10,
                   help="interval between plotfiles")
    p.add_argument("--xmin", type=float, default=0,
                   help="minimum x-coordinate to show")
    p.add_argument("--xmax", type=float, default=-1,
                   help="maximum x-coordinate to show")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")
    p.add_argument("--limitlabels", type=float, default=1.,
                   help="Show all labels (default) or reduce to ~ given value")

    args = p.parse_args()

    plot_prefix = args.plotfiles[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.plotfiles], key=int)

    doit(plot_prefix, plot_nums, args.skip, args.limitlabels, args.xmin, args.xmax)
