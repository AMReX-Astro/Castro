#!/usr/bin/env python3

# Take a sequence of plotfiles and plot T and enuc vs. position

import argparse
import sys

import numpy as np
from cycler import cycler
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import math
import yt

## Define RGBA to HEX
def rgba_to_hex(rgba):
    r = int(rgba[0]*255.0)
    g = int(rgba[1]*255.0)
    b = int(rgba[2]*255.0)
    return f'#{r:02X}{g:02X}{b:02X}'

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

def get_nuc_profile(plotfile):

    ds = yt.load(plotfile, hint="castro")

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot
    srt = np.argsort(ad['x'])
    x_coord = np.array(ad['x'][srt])

    nuc_fracs = []

    for f in ds.field_list:
        if f[1][0] == "X":
            nuc_fracs.append(np.array(ad[f[1]][srt]))

    return time, x_coord, nuc_fracs

def plot_Te(prefix, nums, skip, limitlabels, xmin, xmax):

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

        pfile = f"{prefix}{nums[n]}"

        time, x, T, enuc = get_Te_profile(pfile)

        if index % skiplabels == 0:
            ax_T.plot(x, T, label=f"t = {time:6.4g} s")
        else:
            ax_T.plot(x, T)

        ax_e.plot(x, enuc)

        index = index + 1

    ax_T.legend(frameon=False)
    ax_T.set_ylabel("T (K)")

    if xmax > 0:
        ax_T.set_xlim(xmin, xmax)

    ax_e.set_yscale("log")
    ax_e.set_ylabel(r"$S_\mathrm{nuc}$ (erg/g/s)")
    ax_e.set_xlabel("x (cm)")

    cur_lims = ax_e.get_ylim()
    ax_e.set_ylim(1.e-10*cur_lims[-1], cur_lims[-1])

    if xmax > 0:
        ax_e.set_xlim(xmin, xmax)

    f.savefig("det_Te.png")

def plot_nuc_frac(prefix, nums, skip, limitlabels, xmin, xmax):

    f = plt.figure()
    f.set_size_inches(16.0, 10.0)

    # Get set of colors to use and apply to plot
    numplots = int(len(nums) / skip)
    cm = plt.get_cmap('nipy_spectral')
    clist = [cm(0.95*i/numplots) for i in range(numplots + 1)]
    hexclist = [rgba_to_hex(ci) for ci in clist]
    
    if limitlabels > 1:
        skiplabels = int(numplots / limitlabels)
    elif limitlabels < 0:
        print("Illegal value for limitlabels: %.0f" % limitlabels)
        sys.exit()
    else:
        skiplabels = 1
        
    pfile = f"{prefix}{nums[1]}"
    ds = yt.load(pfile, hint="castro")
    
    nuc_fracs = [f[1] for f in ds.field_list if f[1][0] == "X"]
    print(nuc_fracs)
    N = len(nuc_fracs)
    nrows = math.ceil(math.sqrt(N))
    ncols = math.ceil(math.sqrt(N))
    for i in range(N):
        ax = f.add_subplot(nrows, ncols, i+1)

        index = 0
        for n in range(0, len(nums), skip):

            pfile = f"{prefix}{nums[n]}"

            time, x, nuc_prof = get_nuc_profile(pfile)

            if i == 0 and index % skiplabels == 0:
                ax.plot(x, nuc_prof[i], label=f"t = {time:6.4g} s")
            else:
                ax.plot(x, nuc_prof[i])

            index = index + 1

        ax.legend(frameon=False)
        ax.set_ylabel(nuc_fracs[i])        
        ax.set_yscale("log")

        if xmax > 0:
            ax.set_xlim(xmin, xmax)
    
    f.savefig("det_nuc.png")
    
def doit(prefix, nums, skip, limitlabels, xmin, xmax, do_nuc_fracs=False):

    if do_nuc_fracs:
        plot_nuc_frac(prefix, nums, skip, limitlabels, xmin, xmax)
    else:
        plot_Te(prefix, nums, skip, limitlabels, xmin, xmax)

if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--skip", type=int, default=1,
                   help="interval between plotfiles")
    p.add_argument("--xmin", type=float, default=0,
                   help="minimum x-coordinate to show")
    p.add_argument("--xmax", type=float, default=-1,
                   help="maximum x-coordinate to show")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")
    p.add_argument("--limitlabels", type=float, default=1.,
                   help="Show all labels (default) or reduce to ~ given value")
    p.add_argument("--do_nuc_fracs", type=bool, default=False,
                   help="Set to True if want to plot nuc fracs, otherwise Temp and enuc plot")
    
    args = p.parse_args()

    plot_prefix = args.plotfiles[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.plotfiles], key=int)

    doit(plot_prefix, plot_nums, args.skip, args.limitlabels, args.xmin, args.xmax, do_nuc_fracs=args.do_nuc_fracs)
