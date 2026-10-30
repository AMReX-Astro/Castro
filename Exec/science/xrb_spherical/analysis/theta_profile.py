#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import yt
from yt.frontends.boxlib.api import CastroDataset
from yt.units import (cm, s, g)

def plot_theta_profile(ds, field, r_kms, figsize=(7, 9), outName=None):
    """
    Plot theta profile for different field at multiple radial depths on the same axes.
    """

    ds.force_periodicity()
    dims = ds.domain_dimensions.copy()
    dims[2] = 1

    # Work with level 0 data only for now.
    cg = ds.smoothed_covering_grid(level=0, left_edge=ds.domain_left_edge, dims=dims)

    f = cg["boxlib", field][:, :, 0].to_ndarray()
    r = cg["index", "r"][:, :, 0].to("km").to_ndarray()
    theta = cg["index", "theta"][:, :, 0].to_ndarray()

    r_1d = r[:, 0]
    theta_1d = theta[0, :]

    fig, ax = plt.subplots(figsize=figsize)
    colors  = plt.cm.plasma(np.linspace(0, 1, len(r_kms)))

    for r_km, color in zip(r_kms, colors):
        # find the r data that is closest to the asked r
        r_idx  = np.argmin(np.abs(r_1d - r_km))
        r_actual = r_1d[r_idx]
        f_slice  = f[r_idx, :]

        ax.plot(theta_1d, f_slice, color=color, label=f"r = {r_actual:.4f} km")

    ax.set_xlabel(r"$\theta$ [rad]")
    ax.set_ylabel(field)
    ax.set_yscale("log")
    ax.legend(fontsize=8)

    time = ds.current_time.to("ms").value
    ax.set_title(f"t = {time:.3f} ms")
    ax.grid(linestyle=":")

    # Show ticks on all 4 axes (top, bottom, left, right)
    ax.tick_params(top=True, bottom=True, left=True, right=True)
    plt.tight_layout()
    if outName is not None:
        fig.savefig(outName, format="png", bbox_inches="tight")
    else:
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
        A script to plot theta profile of a quantity over different depth (r).
        """)

    parser.add_argument('fname', nargs=1, type=str,
                        help="dataset file name for plotting")
    parser.add_argument('field', type=str,
                        help="what field quantity to plot")
    parser.add_argument('--rmin', default=None, type=float,
                        help="minimum depth to do the theta profile")
    parser.add_argument('--rmax', default=None, type=float,
                        help="maximum depth to do the theta profile")
    parser.add_argument('--rstep', default=10, type=float,
                        help="Number of theta profiles to plot")
    parser.add_argument("--figsize", nargs=2, type=float, default=[8, 8],
                        metavar=("WIDTH", "HEIGHT"), help="Figure size in inches.")
    parser.add_argument("-o", "--output", default=None, type=str, metavar="FILENAME",
                        help="Output filename (PNG). If not set, shows interactive plot.")

    args = parser.parse_args()

    ds = CastroDataset(args.fname[0])

    rmin = args.rmin
    if rmin is None:
        rel_height= ds.parameters.get("problem.H_star") + 1.5 * ds.parameters.get("problem.atm_delta")
        rmin = ds.domain_left_edge[0].in_units("cm") + rel_height*cm
        rmin = rmin.value * 1e-5 # convert to km

    rmax = args.rmax
    if rmax is None:
        rmax = rmin + 0.05

    r_kms = np.linspace(rmin, rmax, args.rstep)

    plot_theta_profile(ds, args.field, r_kms, figsize=args.figsize, outName=args.output)
