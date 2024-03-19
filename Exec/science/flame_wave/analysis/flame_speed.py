#!/usr/bin/env python3

import sys
import pandas as pd
from scipy.stats import linregress

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

plt.rcParams.update \
(
    {
        "font.family": "serif"
    }
)

def measure_and_plot(dat, stab_ind):
    """
    Plots front_location vs. time on the current pyplot figure, as well as any
    reasonable linear fits. *stab_ind* should give an index for which the flame
    front has stabilized, and the slope obtained from linear regression will yield
    an accurate measurement of the front speed.
    """

    slopes = dict()

    radarr = dat["enuc[0.1%]"]
    tarr = dat["time"]

    plt.plot(tarr * 1000, radarr, linewidth=2)

    indx = tarr.argsort()
    tarr = tarr[indx][stab_ind:]
    radarr = radarr[indx][stab_ind:]

    m, b, r, _, sd = linregress(tarr, radarr)
    print("{:>20}\t{:.2e}\t{:.3f}\t{:.2e}".format("enuc[0.1%]", m, r, sd))

    if r > 0.8:

        plt.plot(dat["time"] * 1000, m * dat["time"] + b, "k--", linewidth=2)

    ax = plt.gca()
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.major.formatter._useMathText = True
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.xlabel("t [ms]")
    plt.ylabel("r [cm]")
    plt.title("Front Position vs. Time", fontsize=24)

    for item in (ax.xaxis.label, ax.yaxis.label):

        item.set_fontsize(20)
        item.set_color("dimgrey")

    for item in (ax.get_xticklabels() + ax.get_yticklabels()) + [ax.yaxis.offsetText]:

        item.set_fontsize(16)
        item.set_color("dimgrey")

    return m

if __name__ == "__main__":

    try:

        dat = pd.read_csv(sys.argv[1], sep="\s+")
        stab_ind = int(sys.argv[2])

        plt.style.use("ggplot")
        measure_and_plot(dat, stab_ind)
        plt.show()

    except (ValueError, IndexError):

        print("Usage: ./fit_speed.py data_file starting_index")
