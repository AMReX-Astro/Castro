#!/usr/bin/env python


import argparse

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import yt

def get_Te_profile(plotfile):

    ds = yt.load(plotfile, hint="castro")

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot
    srt = np.argsort(ad['x'])
    x_coord = np.array(ad['x'][srt])
    temp = np.array(ad['Temp'][srt])
    pres = np.array(ad['pressure'][srt])
    dens = np.array(ad['density'][srt])
    vel = np.array(ad['x_velocity'][srt])
    enuc = np.array(ad['enuc'][srt])

    return time, x_coord, temp, pres, dens, vel, enuc


def doit(prefix, nums, xmax):

    f = plt.figure()
    f.set_size_inches(9.0, 9.0)

    ax_T = f.add_subplot(221)
    ax_r = f.add_subplot(222)
    ax_p = f.add_subplot(223)
    ax_v = f.add_subplot(224)

    for n in range(0, len(nums)):

        pfile = f"{prefix}{nums[n]}"

        time, x, T, p, rho, v, enuc = get_Te_profile(pfile)

        ax_T.plot(x, T, color=f"C{n}", label=f"t = {time:6.4g} s")
        ax_r.plot(x, rho, color=f"C{n}")
        ax_p.plot(x, p, color=f"C{n}")
        ax_v.plot(x, v, color=f"C{n}")

    ax_T.legend(frameon=False, fontsize="small")
    ax_T.set_ylabel("T (K)")
    if xmax > 0:
        ax_T.set_xlim(0, xmax)

    ax_r.set_ylabel(r"$\rho (g\ cm^{-3})$")
    if xmax > 0:
        ax_r.set_xlim(0, xmax)

    ax_p.set_ylabel(r"$P (erg\ cm^{-3})$")
    if xmax > 0:
        ax_p.set_xlim(0, xmax)

    ax_v.set_ylabel("v (cm/s)")
    if xmax > 0:
        ax_v.set_xlim(0, xmax)

    f.tight_layout()
    f.savefig("flame_start_finish.png")


if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--xmax", type=float, default=-1,
                   help="maximum x-coordinate to show")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")

    args = p.parse_args()

    plot_prefix = args.plotfiles[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.plotfiles], key=int)

    doit(plot_prefix, plot_nums, args.xmax)
