#!/usr/bin/env python
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt
import argparse
from yt.frontends.boxlib.api import CastroDataset

SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 20
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)
plt.rc('xtick.major', size=7, width=2)
plt.rc('xtick.minor', size=5, width=1)
plt.rc('ytick.major', size=7, width=2)
plt.rc('ytick.minor', size=5, width=1)


def extract_total_energy(ts):
    """Extract the total energy E as a function of time.

    Parameters
    ----------
    ts : yt time series
        the time series to extract energy from

    Returns
    -------
    tuple of numpy.ndarray
        time array and energy deviation array (E - E0)
    """

    E_tot = []
    time = []

    for ds in ts:
        # Get coarse grid data
        cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                               dims=ds.domain_dimensions)

        # Get the sum of the total energy of the whole domain
        E = np.array(cg['boxlib', 'rho_E'].to_ndarray().squeeze())
        vol = np.array(cg['boxlib', 'volume'].to_ndarray().squeeze())
        E_tot.append((vol * E).sum())

        time.append(float(ds.current_time))

    time = np.array(time)
    E_tot = np.array(E_tot)
    dE = np.abs(E_tot - E_tot[0])

    return time, dE


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", nargs="+", metavar="LABEL:PATH",
                        help="runs in the form label:glob, e.g. --runs run1:'/path/to/run1/plt*' run2:'/path/to/run2/plt*'")
    parser.add_argument("-o", "--output", default=None, type=str, metavar="FILENAME",
                        help="Output filename (PNG). If not set, shows interactive plot.")
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(8, 8))

    print(args.runs)

    for entry in args.runs:
        # Get label and the path containing plt files
        label, path = entry.split(":", 1)
        ts = yt.load(path)

        # Get data and do plotting
        time, dE = extract_total_energy(ts)
        ax.plot(time, dE, '-^', linewidth=2, markersize=10, label=label)

    ax.set_ylabel(r"$E - E_0$ [ergs]")
    ax.grid(linestyle=":")
    ax.tick_params(top=True, bottom=True, left=True, right=True)
    ax.legend(frameon=False)
    ax.tick_params(axis="both",direction="in")
    ax.set_xlabel("Time [s]")

    fig.tight_layout()

    if args.output is not None:
        fig.savefig(args.output, format="png", bbox_inches="tight")
    else:
        plt.show()
