#!/usr/bin/env python


import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt
import sys
import numpy as np


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
    rho_He = np.array(ad['rho_he4'][srt])

    return time, x_coord, temp, enuc, rho_He


def doit(plotfiles, xmin, xmax):

    f = plt.figure()
    f.set_size_inches(7.0, 9.0)

    ax_T = f.add_subplot(311)
    ax_e = f.add_subplot(312)
    ax_he = f.add_subplot(313)

    for p in plotfiles:
        if plotfiles[p] is None:
            continue

        time, x, T, enuc, rho_he4 = get_Te_profile(plotfiles[p])

        ax_T.plot(x, T, label=p)
        ax_e.plot(x, enuc)
        ax_he.plot(x, rho_he4)

    ax_T.legend(frameon=False)

    ax_T.set_ylabel(r"$T$ (K)")
    ax_e.set_ylabel(r"$H_\mathrm{nuc}$ (erg/g/s)")
    ax_he.set_ylabel(r"$\rho X({}^4\mathrm{He})$ (g/cm${}^3$)")
    ax_he.set_xlabel("x (cm)")

    ax_e.set_yscale("log")
    ax_he.set_yscale("log")

    ax_T.set_xlim(left=xmin, right=xmax)
    ax_e.set_xlim(left=xmin, right=xmax)
    ax_he.set_xlim(left=xmin, right=xmax)

    f.tight_layout()

    f.savefig("det.png")


if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--ctu", help="plotfile with Strang CTU splitting", type=str, default=None)
    p.add_argument("--lobatto2", help="plotfile with Lobatto SDC-2", type=str, default=None)
    p.add_argument("--lobatto4", help="plotfile with Lobatto SDC-4", type=str, default=None)
    p.add_argument("--radau4", help="plotfile with Radau SDC-4", type=str, default=None)
    p.add_argument("-x", help="minimum x value", type=float, default=0.0)
    p.add_argument("-X", help="maximum x value", type=float, default=None)

    args = p.parse_args()

    pfiles = {}
    pfiles["Strang CTU"] = args.ctu
    pfiles["Gauss-Lobatto SDC-2"] = args.lobatto2
    pfiles["Gauss-Lobatto SDC-4"] = args.lobatto4
    pfiles["Radau SDC-4"] = args.radau4

    doit(pfiles, args.x, args.X)
