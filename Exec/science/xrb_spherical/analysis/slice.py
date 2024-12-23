#!/usr/bin/env python3

import sys
import os
import yt
import argparse
import math
from typing import List, Optional
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from yt.frontends.boxlib.api import CastroDataset
from yt.units import cm


def slice(fnames:List[str], fields:List[str],
          loc: str = "top", widthScale: float = 4.0,
          theta: Optional[float] = None) -> None:
    """
    A slice plot of the datasets for different field parameters for Spherical2D geometry.

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

    theta:  user defined theta center of the slice plot
    """

    ts = [CastroDataset(fname) for fname in fnames]

    fig = plt.figure(figsize=(16, 9))

    num = len(fields)*len(fnames)
    ny = math.ceil(np.sqrt(num))
    nx = math.ceil(num/ny)

    grid = ImageGrid(fig, 111, nrows_ncols=(nx, ny),
                     axes_pad=1, label_mode="L", cbar_location="right",
                     cbar_mode="each", cbar_size="2.5%", cbar_pad="0%")

    # Output plot file name
    outName = "xrb_spherical_slice.pdf"

    # Determine the appropriate time unit
    refTimeStamp = ts[0].current_time
    timeUnit = "ms"
    if float(refTimeStamp) < 1e-4:
        timeUnit = "us"

    refTimeStamp = ts[0].current_time.in_units(timeUnit)

    # We will check if all files have the same timestamp later.
    # This matters on how we annotate timestamp for the slice plots.
    sameTime = True
    if len(fnames) > 1:
        sameTime = all(refTimeStamp == ds.current_time.in_units(timeUnit) for ds in ts)

    for j, ds in enumerate(ts):
        # Process information for each dataset

        currentTime = ds.current_time.in_units(timeUnit)

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
        # Centers will be physical coordinates in Cylindrical, i.e. R-Z
        centers = {"top":(r_center*np.sin(thetal)+0.5*width, r_center*np.cos(thetal)),
                   "mid":(r_center*np.sin(theta_center), r_center*np.cos(theta_center)),
                   "bot":(r_center*np.sin(thetar)+0.5*width, r_center*np.cos(thetar))}

        if theta is None:
            center = centers[loc]
        else:
            R = r_center*np.sin(theta)
            Z = r_center*np.cos(theta)
            if R < 0.5*width:
                R = 0.5*width
            center = [R, Z]

        for i, field in enumerate(fields):
            # Plot each field parameter

            sp = yt.SlicePlot(ds, 'phi', field, width=box_widths, fontsize=14)
            sp.set_center(center)

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

            sp.set_buff_size((2400,2400))
            sp.set_axes_unit("km")
            # sp.annotate_text((0.05, 0.05), f"{currentTime.in_cgs():8.5f} s")

            plot = sp.plots[field]
            plot.figure = fig
            plot.axes = grid[i+j*len(fields)].axes
            plot.cax = grid.cbar_axes[i+j*len(fields)]

            sp._setup_plots()

    if sameTime:
        # fig.text(0.8, 0.05, f"t = {refTimeStamp:.2f}",
        #          horizontalalignment='right', verticalalignment='center',
        #          color="black", transform=fig.transFigure)

        outName = f"{ts[0]}_slice.pdf"

    fig.set_size_inches(16, 9)
    fig.tight_layout()
    fig.savefig(outName, format="png", bbox_inches="tight")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
        A slice plot script for xrb_spherical problem.
        Given a list of plotfiles or a list of field parameters,
        it plots multiple slice plots.
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
    parser.add_argument('-l', '--loc', default='top', type=str, metavar="{top, mid, bot}",
                        help="""preset center location of the plot domain.
                        Enter one of the three choices: {top, mid, bot}""")
    parser.add_argument('-t', '--theta', type=float,
                        help="""user defined theta center location of the plot domain.
                        Alternative way of defining plotting center""")
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
          loc=loc, theta=args.theta, widthScale=args.width)
