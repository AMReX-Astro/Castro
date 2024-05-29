#!/usr/bin/env python3

# Take a sequence of plotfiles and plot T and enuc vs. position

import matplotlib
import numpy as np

matplotlib.use('agg')

import matplotlib.pyplot as plt

import matplotlib.ticker as mticker

import detonation

def plot_Te(ddir):

    f = plt.figure()

    f.set_size_inches(7.5, 9.0)

    ax_T = f.add_subplot(211)
    ax_e = f.add_subplot(212)

    d = detonation.Detonation(ddir)
    profile = d.get_data(-1)

    ax_T.plot(profile.x, profile.T)

    ishk = profile.shock == 1.0
    ax_T.scatter(profile.x[ishk], profile.T[ishk])

    ax_e.plot(profile.x, profile.enuc)
    ax_e.scatter(profile.x[ishk], profile.enuc[ishk])

    ax_T.set_ylabel("T (K)")
    ax_T.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    ax_T.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))

    ax_e.set_yscale("log")
    ax_e.set_ylabel(r"$S_\mathrm{nuc}$ (erg/g/s)")
    ax_e.set_xlabel("x (cm)")
    ax_e.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    cur_lims = ax_e.get_ylim()
    ax_e.set_ylim(1.e-10*cur_lims[-1], cur_lims[-1])

    f.tight_layout()
    f.savefig("shock_flag.png")


if __name__ == "__main__":

    #ddir = "res12.288km"
    #ddir = "res1.536km"
    #ddir = "res0.192km"
    ddir = "res0.003km"

    plot_Te(ddir)
