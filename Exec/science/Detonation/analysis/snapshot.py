#!/usr/bin/env python3

# Plot the structure of a single plotfile, including the most abundance nuclei

import argparse
import re
import sys

import matplotlib
import matplotlib.ticker as mticker

import numpy as np
from cycler import cycler

matplotlib.use('agg')
import math

import matplotlib.pyplot as plt
import yt

# important nuclei and their colors
main_nuclei = {"X(He4)": "C0",
               "X(Fe52)": "C1",
               "X(Fe54)": "C2",
               "X(Fe56)": "C3",
               "X(Ni56)": "C4",
               "X(Ni58)": "C6"}


# Extract number from nuc list
def nuc_list_filter(nuc):

    match = re.search(r'\d+', nuc)
    if match:
        return int(match.group())

    return 0


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

    nuc_names = [f[1] for f in ds.field_list if f[1][0] == "X"]
    nuc_names.sort(key=nuc_list_filter)

    nuc_fracs = [np.array(ad[nuc][srt]) for nuc in nuc_names]

    return time, x_coord, nuc_fracs, nuc_names


def doit(pf, xmin, xmax, nuc_thresh):

    f = plt.figure()

    # Get set of colors to use and apply to plot

    f.set_size_inches(7.0, 10.0)

    ax_T = f.add_subplot(311)
    ax_e = f.add_subplot(312)
    ax_X = f.add_subplot(313)

    time, x, T, enuc = get_Te_profile(pf)

    _, __, mass_fractions, nuc_names = get_nuc_profile(pf)

    ax_T.plot(x, T)
    ax_e.plot(x, enuc)

    for nuc, X in zip(nuc_names, mass_fractions):
        if nuc in main_nuclei and X.max() > nuc_thresh:
            ax_X.plot(x, X, label=nuc, color=main_nuclei[nuc], lw=2)
        elif X.max() > 0.1:
            ax_X.plot(x, X, label=nuc, ls="--")
        elif X.max() > nuc_thresh:
            ax_X.plot(x, X, label=nuc, ls=":")

    ax_T.set_ylabel("T (K)")
    ax_T.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    ax_T.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))

    if xmax > 0:
        ax_T.set_xlim(xmin, xmax)
        ax_e.set_xlim(xmin, xmax)
        ax_X.set_xlim(xmin, xmax)

    ax_e.set_yscale("log")
    ax_e.set_ylabel(r"$S_\mathrm{nuc}$ (erg/g/s)")
    cur_lims = ax_e.get_ylim()
    ax_e.set_ylim(1.e-10*cur_lims[-1], cur_lims[-1])
    ax_e.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))

    ax_X.set_yscale("log")
    ax_X.set_ylabel(r"X")
    ax_X.set_xlabel("x (cm)")
    ax_X.set_ylim(0.1 * nuc_thresh, 1.1)
    ax_X.legend(ncol=2, fontsize="small", loc="upper left")
    ax_X.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))

    f.text(0.02, 0.02, f"t = {float(time):8.3f} s", transform=f.transFigure)
    f.tight_layout()
    f.savefig(f"snapshot_{pf}.png")



if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--xmin", type=float, default=0,
                   help="minimum x-coordinate to show")
    p.add_argument("--xmax", type=float, default=-1,
                   help="maximum x-coordinate to show")
    p.add_argument("--nuc_thresh", type=float, default=1.e-2,
                   help="minimum value of X for a nuclei to be plotted")
    p.add_argument("plotfile", type=str, nargs=1,
                   help="plotfiles to plot")

    args = p.parse_args()

    doit(args.plotfile[0], args.xmin, args.xmax, args.nuc_thresh)
