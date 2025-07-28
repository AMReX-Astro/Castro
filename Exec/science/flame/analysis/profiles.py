#!/usr/bin/env python3


import argparse
import glob
import sys

import matplotlib
import numpy as np
from cycler import cycler
from matplotlib import ticker

matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt

yt.set_log_level(50)

def find_files(plist, times):

    mask = np.zeros(len(times))
    files_to_plot = []
    for pfile in plist:
        ds = yt.load(pfile, hint="castro")
        for k, t in enumerate(times):
            if mask[k]:
                continue
            if ds.current_time >= t:
                files_to_plot.append(pfile)
                print(f"plotting {pfile} at t = {float(ds.current_time):7.3g}")
                mask[k] = 1.0
                break

    return files_to_plot

## Define RGBA to HEX
def rgba_to_hex(rgba):
    r = int(rgba[0]*255.0)
    g = int(rgba[1]*255.0)
    b = int(rgba[2]*255.0)
    return f'#{r:02X}{g:02X}{b:02X}'

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
            ax_T.plot(x, T, label=f"t = {time:6.4g} s")
        else:
            ax_T.plot(x, T)

        ax_e.plot(x, enuc)

        ax_v.plot(x, velx)

        index = index + 1

    ax_T.legend(frameon=False)
    ax_T.set_ylabel("T (K)")
    ax_T.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))

    if xmax > 0:
        ax_T.set_xlim(xmin, xmax)

    max_enuc = np.abs(enuc).max()
    ax_e.set_yscale("symlog", linthresh=1.e-6 * max_enuc)
    ax_e.set_ylabel(r"$S_\mathrm{nuc}$ (erg/g/s)")

    if xmax > 0:
        ax_e.set_xlim(xmin, xmax)

    ax_v.set_ylabel("v (cm/s)")
    ax_v.set_xlabel("x (cm)")
    ax_v.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))

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
    p.add_argument("--dt", type=float, default=0.005,
                   help="time spacing between plotfiles")
    p.add_argument("--tmax", type=float, default=-1.0,
                   help="maximum time to show")
    p.add_argument("--limitlabels", type=float, default=1.,
                   help="Show all labels (default) or reduce to ~ given value")

    args = p.parse_args()

    pfiles = glob.glob(f"{args.prefix}?????")
    pfiles += glob.glob(f"{args.prefix}??????")
    pfiles += glob.glob(f"{args.prefix}???????")

    plot_nums = sorted([p.split(args.prefix)[1] for p in pfiles], key=int)

    psorted = [f"{args.prefix}{q}" for q in plot_nums]

    if args.tmax == -1:
        ds = yt.load(psorted[-1])
        args.tmax = float(ds.current_time) * 1.01

    times = np.arange(0.0, args.tmax + args.dt, args.dt)

    plot_want = find_files(psorted, times)

    doit(plot_want, args.limitlabels, args.xmin, args.xmax)
