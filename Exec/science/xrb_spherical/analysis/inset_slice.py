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

def single_slice(ds, field:str,
                 loc: str = "top", widthScale: float = 3.0,
                 theta: Optional[float] = None,
                 position: float | None = None,
                 displace_theta: bool = True,
                 show_full_star: bool = True) -> None:
    """
    A slice plot a single dataset for a single field parameters for Spherical2D geometry.
    This is mainly a helper function plot to plot the full-star slice plot with another
    inset zoom-in plot in the center.

    Parameters
    ==================================================================================
    ds:
      A single Castro dataset for plotting the slice plot

    field:
      A single field name for plotting the slice plot.

    loc:
      preset center location of the domain. {top, mid, bot}

    widthScale:
      scaling for the domain width of the slice plot

    theta:
      user defined theta center of the slice plot

    position:
      draw a vertical line on user defined front position in theta (in radian)

    displace_theta:
      When theta is explicitly defined, do we want to displace the center
      of the slice plot by ~0.7. This is helpful when tracking the flame front.

    show_full_star:
      do we want to plot the full star instead of a zoom-in slice plot.
    """

    # Process information for the dataset

    # Some geometry properties
    r, box_widths, center = extract_info(ds,
                                         loc=loc, widthScale=widthScale,
                                         widthRatio=1.0,
                                         theta=theta,
                                         displace_theta=displace_theta,
                                         show_full_star=show_full_star)


    # Plot each field parameter
    sp = yt.SlicePlot(ds, 'phi', field, width=box_widths, fontsize=28)
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
    if position is not None:
        sp.annotate_line([r[0]*np.sin(position), r[0]*np.cos(position)],
                         [r[2]*np.sin(position), r[2]*np.cos(position)],
                         coord_system="plot",
                         color="k",
                         linewidth=1.5,
                         linestyle="-.")

    ### Annotate Latitude Lines
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
    parser.add_argument('-p', '--position', type=float,
                        help="""user defined front position in theta,
                        this will draw a vertical line to annotate the front position""")
    parser.add_argument('-w', '--width', default=2.0, type=float,
                        help="scaling for the domain width of the slice plot")
    parser.add_argument('--displace_theta', action='store_true',
                        help="""whether to displace the theta that defines the center of the frame.
                        This is useful when theta represents the flame front position.""")

    args = parser.parse_args()

    loc = args.loc.lower()
    loc_options = ["top", "mid", "bot"]

    if loc not in loc_options:
        parser.error("loc must be one of the three: {top, mid, bot}")

    # First load the data
    ds = CastroDataset(args.fname)

    # First get the slice plot of the full-star.
    full_star_slice = single_slice(ds, args.field, loc=loc,
                                   widthScale=args.width, theta=args.theta, position=args.position,
                                   displace_theta=args.displace_theta,
                                   show_full_star=True)
    # full_star_slice.render()

    # Extract the figure of the full-star slice and use that as the main figure for plotting.
    fig = full_star_slice.plots[args.field].figure

    # Add an inset ax in the middle of the star
    rect = (0.31, 0.4, 0.31, 0.31)  # Left, Bottom, Width, Height
    inset_ax = fig.add_axes(rect)

    # Get the slice of the zoom-in plot
    zoom_in_slice = single_slice(ds, args.field, loc=loc,
                                 widthScale=args.width, theta=args.theta, position=args.position,
                                 displace_theta=args.displace_theta,
                                 show_full_star=False)
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

    # Increase inset ax tick label size
    inset_ax.tick_params(labelsize=16)

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
                boxstyle="round,pad=0.1",
                facecolor="white",
                edgecolor="white",
                )
        )

    ### Now annotate inset box lines ###
    # loc1 and loc2: corners to connect (1: upper right, 2: upper left, 3: lower left, 4: lower right)
    if args.theta < np.pi/3.0 or (args.theta is None and loc=="top"):
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

    # Add the current simulation time of the data
    time = ds.current_time.in_units("ms")
    if float(time) < 1e-1:
        time = ds.current_time.in_units("us")

    fig.text(0.9, 0.05, f't = {time:.2f}', fontsize=18,
             horizontalalignment='center', verticalalignment='center',
             color='black', transform=fig.transFigure)

    # Save the figure
    fig.savefig(f"{args.fname}_slice.png", format='png', bbox_inches='tight')
