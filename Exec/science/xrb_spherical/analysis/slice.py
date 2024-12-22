#!/usr/bin/env python3

import sys
import os
import yt
import argparse
from typing import List, Optional
import numpy as np
import matplotlib.pyplot as plt
from yt.frontends.boxlib.api import CastroDataset

from yt.units import cm

SMALL_SIZE = 28
MEDIUM_SIZE = 30
BIGGER_SIZE = 32

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('xtick.major', size=7, width=2)
plt.rc('xtick.minor', size=5, width=1)
plt.rc('ytick.major', size=7, width=2)
plt.rc('ytick.minor', size=5, width=1)

def slice(fnames:List[str], fields:List[str],
          loc: str = "top", widthScale: float = 4.0,
          coord: Optional[List[float, float]] = None) -> None:
    """
    A slice plot of the dataset for Spherical2D geometry.

    Parameters
    ==================================================================================
    fnames: A list of file names to plot multiple slice plots between different
            plot files for a given field parameter.
            Note that either fname or field must be single valued.

    fields: A list of field parameters to plot multiple slice plots between different
            field parameters for a given file.
            Note that either fname or field must be single valued.

    loc:    preset center location of the domain. {top, mid, bot}

    widthScale: scaling for the domain width of the slice plot

    coord:  user defined center of the slice in the format [r_center, theta_center]
    """

    ts = CastroDataset(fnames)
    fig = plt.figures(figsize=(16, 9))
    grid = ImageGrid(fig, 111, nrows_ncols=(len(fields)*len(ts), 1),
                     axes_pad=0.25, label_mode="L", cbar_location="right",
                     cbar_mode="each", cbar_size="2.5%", cbar_pad="0%")

    # We will check if all files have the same timestamp later.
    # This matters on how we annotate timestamp for the slice plots.
    if len(fnames) == 1:
        sameTime = True
    else:
        refTimeStamp = ts[0].current_time.in_units("ms")
        sameTime = all(refTimeStamp == ds.current_time.in_units("ms") for ds in ts)

    for j, ds in enumerate(ts):
        # Process information for each dataset

        currentTime = ds.current_time.in_units("ms")
        print(f"Current time of this plot file is {currentTime}")

        # Some geometry properties
        rr = ds.domain_right_edge[0].in_units("km")
        rl = ds.domain_left_edge[0].in_units("km")
        dr = rr - rl
        r_center = 0.5 * (rr + rl)

        thetar = ds.domain_right_edge[1]
        thetal = ds.domain_left_edge[1]
        dtheta = thetar - thetal
        theta_center = 0.5 * (thetar + thetal)

        # Domain width of the slice plot
        width = widthScale * dr
        box_widths = (width, width)

        # Preset centers for the Top, Mid and Bot panels
        centers = {"top":(r_center*np.sin(thetal)+0.5*width, r_center*np.cos(thetal)),
                   "mid":(r_center*np.sin(theta_center), r_center*np.cos(theta_center)),
                   "bot":(r_center*np.sin(thetar)+0.5*width, r_center*np.cos(thetar))}

        if coord is None:
            # Set center vis sp.set_center.
            #This will be physical coordinates in Cartesian
            center = centers[loc]
        else:
            # Set center during SlicePlot call.
            # This will be in Spherical coordinate: [r_center, theta_center, 0]
            # coord will be in [r_center, theta_center], so append 0 here.
            center = coord.append(0)

        for i, field in enumerate(fields):
            # Plot each field parameter

            # Note we can also set center during SlicePlot, however then we would enter in [r_center, theta_center, 0]
            # rather than the physical R-Z coordinate if we do it via sp.set_center

            if coord is None:
                sp = yt.SlicePlot(ds, 'phi', field, width=box_widths)
                sp.set_center(center)
            else:
                sp = yt.SlicePlot(ds, 'phi', field, center=center, width=box_widths)

            sp.set_cmap(field, "viridis")
            if field in ["x_velocity", "y_velocity", "z_velocity"]:
                sp.set_cmap(field, "coolwarm")
            elif field == "Temp":
                sp.set_zlim(field, 5.e7, 2.5e9)
                sp.set_cmap(field, "magma_r")
            elif field == "enuc":
                sp.set_zlim(field, 1.e18, 1.e20)
            elif field == "density":
                sp.set_zlim(field, 1.e-3, 5.e8)

            sp.set_axes_unit("km")
            # sp.annotate_text((0.05, 0.05), f"{currentTime.in_cgs():8.5f} s")

            plot = sp.plots[field]
            plot.figure = fig
            plot.axes = grid[i+j*len(fields)].axes
            plot.cax = grid.cbar_axes[i+j*len(fields)]

            if i+j*len(fields) < len(fnames)+len(fields)-1:
                    grid[i].axes.xaxis.offsetText.set_visible(False)

            sp._setup_plots()

    if sameTime:
        fig.text(0.8, 0.05, f"t = {ds.current_time}",
                 horizontalalignment='right', verticalalignment='center',
                 color="black", transform=fig.transFigure)

    fig.set_size_inches(24, 13.5)
    fig.tight_layout()
    fig.savefig(f"{ds}",format="pdf",bbox_inches="tight")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
        A slice plot script for xrb_spherical problem.
        Given a plot file and field name, it gives slice plots at the top,
        middle, and bottom of the domain (shell).
        """)

    parser.add_argument('--fnames', nargs='+', type=str,
                        help="""dataset file names for plotting. Accepts one or more datasets.
                        If multiple file names are given, a grid of slice plots of different
                        files will be plotted for a given field parameter.
                        Note that either fnames or field must be single valued.""")
    parser.add_argument('--fields', nargs='+', type=str,
                        help="""field parameters for plotting. Accepts one or more datasets.
                        If multiple parameters are given, a grid of slice plots of different
                        field parameters will be plotted for a given fname.
                        Note that either fnames or fields must be single valued.
                        """)
    parser.add_argument('-l', '--loc', default='top', type=str,
                        help="""preset center location of the plot domain.
                        Enter one of the three choices: {top, mid, bot}""")
    parser.add_argument('-c', '--coord', nargs=2, type=float,
                        help="""user defined center location of the plot domain.
                        Enter two numbers in the format of [r_center, theta_center]""")
    parser.add_argument('-w', '--width', default=4.0, type=float,
                        help="scaling for the domain width of the slice plot")


    args = parser.parse_args()

    if len(args.fnames) > 1 and len(args.fields) > 1:
        parser.error("Either fnames or fields must be single valued!")

    loc = args.loc.lower()
    loc_options = ["top", "mid", "bot"]

    if loc not in loc_options:
        parser.error("loc must be one of the three: {top, mid, bot}")

    slice(args.fnames, args.fields,
          loc=loc, coord=args.coord, width_factor=args.width)
