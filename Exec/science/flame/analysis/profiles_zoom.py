#!/usr/bin/env python


import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt
import sys
import numpy as np


def get_Te_profile(plotfile):

    ds = yt.load(plotfile, hint="castro")

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot
    srt = np.argsort(ad['x'])
    x_coord = np.array(ad['x'][srt])
    temp = np.array(ad['Temp'][srt])
    enuc = np.array(ad['enuc'][srt])

    return time, x_coord, temp, enuc


def doit(pprefix, nums, skip):

    f = plt.figure()
    f.set_size_inches(7.0, 9.0)

    ax_T = f.add_subplot(211)
    ax_e = f.add_subplot(212)

    for n in range(0, len(nums), skip):

        pfile = f"{prefix}{nums[n]}"

        time, x, T, enuc = get_Te_profile(pfile)

        ax_T.plot(x, T, 'rx', label=f"t = {time:6.4g} s")
        ax_e.plot(x, enuc)

    ax_T.legend(frameon=False)
    ax_T.set_ylabel("T (K)")
    ax_T.set_xlim(22000, 28000)

    ax_e.set_yscale("log")
    ax_e.set_ylabel(r"$S_\mathrm{nuc}$ (erg/g/s)")
    ax_e.set_xlabel("x (cm)")

    ax_e.set_xlim(22000, 28000)

    f.savefig("flame.png")


if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--skip", type=int, default=10,
                   help="interval between plotfiles")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")

    args = p.parse_args()

    prefix = args.plotfiles[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.plotfiles], key=int)
    print(plot_nums)
    doit(prefix, plot_nums, args.skip)
