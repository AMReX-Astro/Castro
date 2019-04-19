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

## Define RGBA to HEX
def rgba_to_hex(rgba):
    r = int(rgba[0]*255.0)
    g = int(rgba[1]*255.0)
    b = int(rgba[2]*255.0)
    return '#{:02X}{:02X}{:02X}'.format(r, g, b)

def get_Te_profile(plotfile):

    ds = yt.load(plotfile)

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot
    srt = np.argsort(ad['x'])
    x_coord = np.array(ad['x'][srt])
    temp = np.array(ad['Temp'][srt])
    enuc = np.array(ad['enuc'][srt])

    return time, x_coord, temp, enuc


def doit(prefix, nums, skip, limitlabels):

    f = plt.figure()
    f.set_size_inches(7.0, 9.0)

    ax_T = f.add_subplot(211)
    ax_e = f.add_subplot(212)

    # Get set of colors to use and apply to plot
    numplots = int(len(nums) / skip)
    cm = plt.get_cmap('nipy_spectral')
    clist = [cm(0.95*i/numplots) for i in range(numplots + 1)]
    hexclist = [rgba_to_hex(ci) for ci in clist]
    ax_T.set_prop_cycle(cycler('color', hexclist))
    ax_e.set_prop_cycle(cycler('color', hexclist))

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

        time, x, T, enuc = get_Te_profile(pfile)

        if index % skiplabels == 0:
            ax_T.plot(x, T, label="t = {:6.4g} s".format(time))
        else:
            ax_T.plot(x, T)

        ax_e.plot(x, enuc)

        index = index + 1

    ax_T.legend(frameon=False)
    ax_T.set_ylabel("T (K)")

    ax_e.set_yscale("log")
    ax_e.set_ylabel(r"$S_\mathrm{nuc}$ (erg/g/s)")
    ax_e.set_xlabel("x (cm)")

    f.savefig("flame.png")


if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--skip", type=int, default=10,
                   help="interval between plotfiles")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")
    p.add_argument("--limitlabels", type=float, default=1.,
                   help="Show all labels (default) or reduce to ~ given value")

    args = p.parse_args()

    plot_prefix = args.plotfiles[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.plotfiles], key=int)

    doit(plot_prefix, plot_nums, args.skip, args.limitlabels)
