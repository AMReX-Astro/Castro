#!/usr/bin/env python3

import sys
import os
import yt
import argparse
import math
from typing import List, Optional
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1 import ImageGrid
from yt.frontends.boxlib.api import CastroDataset
from yt.units import km
from slice import extract_info, annotate_latitude_lines

def single_slice(fname:str, field:str,
                 loc: str = "top", widthScale: float = 3.0,
                 dr: Optional[float] = None,
                 theta: Optional[float] = None,
                 displace_theta: bool = True,
                 annotate_vline: bool = True,
                 annotate_lat_lines: bool = True,
                 show_full_star: bool = True) -> None:
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

    dr: user defined distance between lower r to upper r boundary. Assumed in unit km.
        This is used to change the center and width of the SlicePlot.

    theta:  user defined theta center of the slice plot
    """

    ds = CastroDataset(fname)

    # Process information for the dataset

    # Some geometry properties
    r, box_widths, center = extract_info(ds,
                                         loc=loc, widthScale=widthScale,
                                         dr=dr,
                                         theta=theta,
                                         displace_theta=displace_theta,
                                         show_full_star=show_full_star)


    # Plot each field parameter
    sp = yt.SlicePlot(ds, 'phi', field, width=box_widths, fontsize=18)
    sp.set_center(center)

    sp.set_cmap(field, "viridis")
    if field in ["x_velocity", "y_velocity", "z_velocity"]:
        sp.set_cmap(field, "coolwarm")
        if field == "z_velocity":
            sp.set_zlim(field, -2.e8, 2.e8)
            sp.set_log(field, False)
    elif field == "Temp":
        sp.set_zlim(field, 5.e7, 2.5e9)
        sp.set_cmap(field, "magma_r")
    elif field == "abar":
        sp.set_zlim(field, 4, 8)
        sp.set_log(field, False)
        sp.set_cmap(field, "plasma_r")
    elif field == "enuc":
        sp.set_zlim(field, 1.e15, 1.e20)
        sp.set_log(field, linthresh=1.e11)
    elif field == "density":
        sp.set_zlim(field, 1.e-3, 5.e8)

    sp.set_buff_size((2400,2400))
    sp.set_axes_unit("km")
    # sp.annotate_text((0.05, 0.05), f"{currentTime.in_cgs():8.5f} s")

    # Plot a vertical to indicate flame front
    if theta is not None and annotate_vline:
        sp.annotate_line([r[0]*np.sin(theta), r[0]*np.cos(theta)],
                         [r[2]*np.sin(theta), r[2]*np.cos(theta)],
                         coord_system="plot",
                         color="k",
                         linewidth=1.5,
                         linestyle="-.")

    ### Annotate Latitude Lines
    if annotate_lat_lines:
        annotate_latitude_lines(sp, center, box_widths, r,
                                show_full_star=show_full_star)

    sp._setup_plots()
    return sp


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
    This script plots the full-star along with a zoom-in plot
    displayed in the center. This works with one field only.
        """)

    parser.add_argument('fname', type=str,
                        help="""A single dataset file name for plotting.""")
    parser.add_argument('-f', '--field', type=str,
                        help="""A single field parameter for plotting slice plot.
                        """)
    parser.add_argument('-l', '--loc', default='top', type=str, metavar="{top, mid, bot}",
                        help="""preset center location of the plot domain.
                        Enter one of the three choices: {top, mid, bot}""")
    parser.add_argument('-t', '--theta', type=float,
                        help="""user defined theta center location of the plot domain.
                        Alternative way of defining plotting center""")
    parser.add_argument('-r', '--dr', type=float,
                        help="""Distance between upper r and lower r shown in the SlicePlot.
                        Assumed in unit km. This is used to control center and width of the SlicePlot""")
    parser.add_argument('-w', '--width', default=2.0, type=float,
                        help="scaling for the domain width of the slice plot")
    parser.add_argument('--displace_theta', action='store_true',
                        help="""whether to displace the theta that defines the center of the frame.
                        This is useful when theta represents the flame front position.""")
    parser.add_argument('--annotate_vline', action='store_true',
                        help="whether to annotate a vertical line along the given theta")
    parser.add_argument('--annotate_lat_lines', action='store_true',
                        help="whether to annotate latitude lines")

    args = parser.parse_args()

    loc = args.loc.lower()
    loc_options = ["top", "mid", "bot"]

    if loc not in loc_options:
        parser.error("loc must be one of the three: {top, mid, bot}")

    # First get the slice plot of the full-star.
    full_star_slice = single_slice(args.fname, args.field, loc=loc,
                            widthScale=args.width, dr=args.dr, theta=args.theta,
                            displace_theta=args.displace_theta, annotate_vline=args.annotate_vline,
                            annotate_lat_lines=args.annotate_lat_lines, show_full_star=True)
    # full_star_slice.render()

    # Extract the figure of the full-star slice and use that as the main figure for plotting.
    fig = full_star_slice.plots[args.field].figure

    # Add an inset ax in the middle of the star
    rect = (0.24, 0.35, 0.4, 0.4)  # Left, Bottom, Width, Height
    inset_ax = fig.add_axes(rect)

    # Get the slice of the zoom-in plot
    zoom_in_slice = single_slice(args.fname, args.field, loc=loc,
                            widthScale=args.width, dr=args.dr, theta=args.theta,
                            displace_theta=args.displace_theta, annotate_vline=args.annotate_vline,
                            annotate_lat_lines=args.annotate_lat_lines, show_full_star=False)
    zoom_in_slice.hide_colorbar()
    # zoom_in_slice.render()

    # Export to mpl figure, otherwise it has problems putting it in the inset ax
    fig_zoom = zoom_in_slice.export_to_mpl_figure((1,1))
    # Need to render so that it displays.
    zoom_in_slice.render()

    # Find the image of the zoom-in slice plot, and replot it using imshow onto the inset_ax
    inset_img = fig_zoom.axes[0].images[0]
    inset_ax.imshow(
        inset_img.get_array(),
        extent=inset_img.get_extent(),
        origin=inset_img.origin,
        cmap=inset_img.get_cmap(),
        norm=inset_img.norm
    )

    # Add back annotation lines, i.e. latitude lines
    for line in fig_zoom.axes[0].lines:
        inset_ax.plot(
            line.get_xdata(),
            line.get_ydata(),
            linestyle=line.get_linestyle(),
            linewidth=line.get_linewidth(),
            color=line.get_color(),
            marker=line.get_marker(),
            markersize=line.get_markersize(),
            alpha=line.get_alpha()
        )

    # Add back annotation text, i.e. latitude labeling:
    for text in fig_zoom.axes[0].texts:
        inset_ax.text(
            text.get_position()[0],
            text.get_position()[1],
            text.get_text(),
            fontsize=text.get_fontsize(),
            color=text.get_color(),
            rotation=text.get_rotation(),
            ha=text.get_ha(),
            va=text.get_va(),
            alpha=text.get_alpha(),
            clip_on=True,
            transform=inset_ax.transData,  # Match data coordinates
            bbox=dict(
                boxstyle="square,pad=0.1",
                facecolor="white",
                edgecolor="white",
                )
        )

    ### Now annotate inset box lines ###
    # loc1 and loc2: corners to connect (1: upper right, 2: upper left, 3: lower left, 4: lower right)
    # Probably need to change loc1 and loc2 depending on where the zoom-in plot is.
    if args.theta < np.pi / 3.0 or (args.theta is None and loc=="top"):
        loc1 = 1
        loc2 = 2
    elif args.theta > 2.0*np.pi/3.0 or (args.theta is None and loc=="bot"):
        loc1 = 3
        loc2 = 4
    else:
        # mid
        loc1 = 1
        loc2 = 4

    mark_inset(fig.axes[0], inset_ax, loc1=loc1, loc2=loc2, fc="none", ec="red",
               linestyle="--", linewidth=1.0)

    # Save the figure
    fig.savefig("full_star_slice.png", format='png', bbox_inches='tight')
