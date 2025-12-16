#!/usr/bin/env python

# Take a sequence of plotfiles and plot T and enuc vs. position

import matplotlib
import numpy as np

matplotlib.use('agg')

import matplotlib.pyplot as plt

import matplotlib.ticker as mticker

import detonation

def plot_Te(data):

    f = plt.figure()

    f.set_size_inches(7.5, 9.0)

    ax_T = f.add_subplot(211)
    ax_e = f.add_subplot(212)

    for n, _d in enumerate(data):

        ddir, label = _d

        d = detonation.Detonation(ddir)
        profile = d.get_data(-1)

        ax_T.plot(profile.x, profile.T, label=label)
        ax_e.plot(profile.x, profile.enuc, label=label)

    ax_T.set_ylabel("T (K)")
    ax_T.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    ax_T.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    ax_T.legend()

    ax_e.set_yscale("log")
    ax_e.set_ylabel(r"$S_\mathrm{nuc}$ (erg/g/s)")
    ax_e.set_xlabel("x (cm)")
    ax_e.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    cur_lims = ax_e.get_ylim()
    ax_e.set_ylim(1.e-10*cur_lims[-1], cur_lims[-1])

    f.tight_layout()
    f.savefig("profile_compare.png")


if __name__ == "__main__":


    data = [("res0.003km", "300 cm"),
            ("res0.024km", "2400 cm"),
            ("res0.192km", "0.192 km"),
            ("res1.536km", "1.536 km"),
            ("res12.288km", "12.288 km")]

    plot_Te(data)
