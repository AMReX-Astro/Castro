#!/usr/bin/env python3

# Plot the structure of a single plotfile, including the most abundance nuclei

import argparse
import os
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

import pynucastro as pyna

# important nuclei and their colors
main_nuclei = {"X(He4)": "C0",
               "X(H1)": "C1",
               "X(Na22)": "C2",
               "X(Ca40)": "C3",
               "X(Ti44)": "C4",
               "X(Cr48)": "C5",
               "X(Fe52)": "C6",
               "X(Ni56)": "C7"}


# Extract number from nuc list
def nuc_list_filter(nuc):

    match = re.search(r'\d+', nuc)
    if match:
        return int(match.group())

    return 0


def get_rTe_profile(plotfile):

    ds = yt.load(plotfile, hint="castro")

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot
    srt = np.argsort(ad['x'])
    x_coord = np.array(ad['x'][srt])
    density = np.array(ad['density'][srt])
    temp = np.array(ad['Temp'][srt])
    enuc = np.array(ad['enuc'][srt])
    return time, x_coord, density, temp, enuc

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

    f = plt.figure(figsize=(7.0, 10.0), constrained_layout=True)

    # Get set of colors to use and apply to plot

    ax_T = f.add_subplot(311)
    ax_e = f.add_subplot(312)
    ax_X = f.add_subplot(313)

    time, x, rho, T, enuc = get_rTe_profile(pf)

    _, __, mass_fractions, nuc_names = get_nuc_profile(pf)

    ax_T.plot(x, T, color="C0")
    ax_rho = ax_T.twinx()
    ax_rho.plot(x, rho, color="C1")

    ax_e.plot(x, enuc)

    nplot = 0
    for nuc, X in zip(nuc_names, mass_fractions):
        _name = re.search(r'\(([^)]+)\)', nuc).group(1)
        nucleus = pyna.Nucleus(_name)
        if nuc in main_nuclei and X.max() > nuc_thresh:
            ax_X.plot(x, X, label=rf"${nucleus.pretty}$", color=main_nuclei[nuc], lw=2)
            nplot += 1
        elif X.max() > 0.1:
            ax_X.plot(x, X, label=rf"${nucleus.pretty}$", ls="--")
            nplot += 1
        elif X.max() > nuc_thresh:
            ax_X.plot(x, X, label=rf"${nucleus.pretty}$", ls=":")
            nplot += 1

    ax_T.set_ylabel("T (K)", color="C0")
    ax_T.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    ax_T.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))

    ax_rho.set_ylabel(r"$\rho\ (\mathrm{g/cm^3})$", color="C1")
    ax_rho.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))

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

    ax_T.set_xmargin(0)
    ax_e.set_xmargin(0)
    ax_X.set_xmargin(0)

    ax_T.grid(ls=":", axis="x")
    ax_e.grid(ls=":")
    ax_X.grid(ls=":")

    ncol = 1
    if nplot >= 5:
        ncol = 2

    ax_X.legend(ncol=ncol, fontsize="small", columnspacing=1.0) #, loc="upper left")
    ax_X.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))

    f.text(0.02, 0.005, f"t = {float(time):8.3f} s", transform=f.transFigure)

    prefix = os.getcwd().split("/")[-1]
    time_str = f"{time:05.3f}s"
    f.savefig(f"snapshot_{prefix}_{pf}_{time_str}.png")



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
