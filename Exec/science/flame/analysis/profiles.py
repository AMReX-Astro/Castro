#!/usr/bin/env python3

from __future__ import print_function

import argparse
import glob
import sys

import matplotlib
import numpy as np
from cycler import cycler

matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt

yt.set_log_level(50)

dt = 0.005
tmax = 0.05

times = np.arange(0.0, tmax + dt, dt)

def find_files(plist):

    mask = np.zeros(len(times))
    files_to_plot = []
    for pfile in plist:
        ds = yt.load(pfile, hint="castro")
        for k, t in enumerate(times):
            if mask[k]:
                continue
            if ds.current_time >= t:
                files_to_plot.append(pfile)
                print(f"plotting {pfile} at t = {t}")
                mask[k] = 1.0
                break

    return files_to_plot

## Define RGBA to HEX
def rgba_to_hex(rgba):
    r = int(rgba[0]*255.0)
    g = int(rgba[1]*255.0)
    b = int(rgba[2]*255.0)
    return '#{:02X}{:02X}{:02X}'.format(r, g, b)

def get_Tev_profile(plotfile):

    ds = yt.load(plotfile, hint="castro")

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot
    srt = np.argsort(ad['x'])
    x_coord = np.array(ad['x'][srt])
    temp = np.array(ad['Temp'][srt])
    enuc = np.array(ad['enuc'][srt])
    vx = np.array(ad['x_velocity'][srt])

    return time, x_coord, temp, enuc, vx


def doit(pfiles, limitlabels, xmin, xmax):

    f = plt.figure()
    f.set_size_inches(7.0, 10.0)

    ax_T = f.add_subplot(311)
    ax_e = f.add_subplot(312)
    ax_v = f.add_subplot(313)

    # Get set of colors to use and apply to plot
    numplots = int(len(pfiles))
    cm = plt.get_cmap('nipy_spectral')
    clist = [cm(i) for i in np.linspace(0.0, 0.95, numplots, endpoint=True)]
    hexclist = [rgba_to_hex(ci) for ci in clist]
    ax_T.set_prop_cycle(cycler('color', hexclist))
    ax_e.set_prop_cycle(cycler('color', hexclist))
    ax_v.set_prop_cycle(cycler('color', hexclist))

    if limitlabels > 1:
        skiplabels = int(numplots / limitlabels)
    elif limitlabels < 0:
        print("Illegal value for limitlabels: %.0f" % limitlabels)
        sys.exit()
    else:
        skiplabels = 1
    index = 0

    for pf in pfiles:

        time, x, T, enuc, velx = get_Tev_profile(pf)

        if index % skiplabels == 0:
            ax_T.plot(x, T, label="t = {:6.4g} s".format(time))
        else:
            ax_T.plot(x, T)

        ax_e.plot(x, enuc)

        ax_v.plot(x, velx)

        index = index + 1

    ax_T.legend(frameon=False)
    ax_T.set_ylabel("T (K)")

    if xmax > 0:
        ax_T.set_xlim(xmin, xmax)

    ax_e.set_yscale("log")
    ax_e.set_ylabel(r"$S_\mathrm{nuc}$ (erg/g/s)")

    cur_lims = ax_e.get_ylim()
    ax_e.set_ylim(1.e-10*cur_lims[-1], cur_lims[-1])

    if xmax > 0:
        ax_e.set_xlim(xmin, xmax)

    ax_v.set_ylabel("v (cm/s)")
    ax_v.set_xlabel("x (cm)")

    if xmax > 0:
        ax_v.set_xlim(xmin, xmax)

    f.tight_layout()
    f.savefig("flame.png")


if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--prefix", type=str, default="plt",
                   help="plotfile prefix")
    p.add_argument("--xmin", type=float, default=0,
                   help="minimum x-coordinate to show")
    p.add_argument("--xmax", type=float, default=-1,
                   help="maximum x-coordinate to show")
    p.add_argument("--limitlabels", type=float, default=1.,
                   help="Show all labels (default) or reduce to ~ given value")

    args = p.parse_args()

    pfiles = glob.glob(f"{args.prefix}?????")
    pfiles += glob.glob(f"{args.prefix}??????")
    pfiles += glob.glob(f"{args.prefix}???????")

    plot_nums = sorted([p.split(args.prefix)[1] for p in pfiles], key=int)

    psorted = [f"{args.prefix}{q}" for q in plot_nums]

    plot_want = find_files(psorted)
    print(plot_want)

    doit(plot_want, args.limitlabels, args.xmin, args.xmax)
