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
    rl = ds.domain_left_edge[0].in_units("km")
    if dr is None:
        rr = ds.domain_right_edge[0].in_units("km")
        dr = rr - rl
    else:
        dr = dr * rl.units
        rr = rl + dr

    if show_full_star:
        r_center = 0.5 * rr
    else:
        r_center = 0.5 * dr + rl

    thetar = ds.domain_right_edge[1]
    thetal = ds.domain_left_edge[1]
    dtheta = thetar - thetal
    theta_center = 0.5 * (thetar + thetal)
        
    # Domain width of the slice plot
    width = widthScale * dr
    box_widths = (width, width)

    # Now determine center of the frame

    if show_full_star:
            # If we want to show the full domain in the background.
            # Change box_widths and center to full star
            center = [r_center*np.sin(theta_center), r_center*np.cos(theta_center)]
            if thetar < 0.5 * np.pi:
                box_widths = (rr*np.sin(thetar), rr*np.cos(thetal))
            else:
                box_widths = ( rr, rr*(np.abs(np.cos(thetar)) + np.cos(thetal)) )

    elif theta is None:
        # Preset centers for the Top, Mid and Bot panels
        # Centers will be physical coordinates in Cylindrical, i.e. R-Z
        centers = {"top":(r_center*np.sin(thetal)+0.5*width, r_center*np.cos(thetal)),
                   "mid":(r_center*np.sin(theta_center), r_center*np.cos(theta_center)),
                   "bot":(r_center*np.sin(thetar)+0.5*width, r_center*np.cos(thetar))}
        
        center = centers[loc]

    else:
        if displace_theta:
            # Optionally displace the center by ~0.7 of the plotting width
            # This is helpful when tracking the flame front.
            # This keeps the front at ~0.7 of the plotting width.
            
            # Determine dtheta that displaces from center to ~0.7 of the plotting domain
            oSevenTheta = np.arcsin(0.7 * width / r_center)
            halfTheta = np.arcsin(0.5 * width / r_center)
            dtheta = oSevenTheta - halfTheta
        else:
            dtheta = 0

        # Determine center using theta but also displace it by dtheta
        R = r_center*np.sin(theta - dtheta)
        Z = r_center*np.cos(theta - dtheta)
        if R < 0.5 * width:
            R = 0.5 * width
        center = [R, Z]


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
        sp.annotate_line([rl*np.sin(theta), rl*np.cos(theta)],
                         [rr*np.sin(theta), rr*np.cos(theta)],
                         coord_system="plot",
                         color="k",
                         linewidth=1.5,
                         linestyle="-.")

    ### Annotate Latitude Lines

    # Start from the theta center [deg] of the slice plot frame.
    thetac = round(math.degrees(np.arcsin(center[0] / r_center)))
    latitude_thetas = [thetac]

    # Determine the upper and lower bound of the frame
    lobnd_r = center[0] - 0.5 * box_widths[0]
    lobnd_z = center[1] - 0.5 * box_widths[1]
    hibnd_r = center[0] + 0.5 * box_widths[0]
    hibnd_z = center[1] + 0.5 * box_widths[1]

    if show_full_star:
        # Plot latitude line every 15 degrees
        start = 15
        end = 195
        step = 15
    else:
        # Plot latitude line every degree
        start = 1
        end = 181
        step = 1

    # Find what theta to latitude lines
    for theta_increment in range(start, end, step):
        # This is actually in theta-coordinate.
        latitude_thetar = thetac + theta_increment
        latitude_thetal = thetac - theta_increment

        # Now find the RZ position of different latitude thetas
        # and see if they're out of the frame.
        lo_r = rl * np.sin(math.radians(latitude_thetal))
        lo_z = rl * np.cos(math.radians(latitude_thetal))

        hi_r = rl * np.sin(math.radians(latitude_thetar))
        hi_z = rl * np.cos(math.radians(latitude_thetar))

        # Check if the point is within the frame and append point
        if (0 <= latitude_thetal < 180 and lo_r >= lobnd_r
            and lo_z >= lobnd_z and lo_z <= hibnd_z):
            latitude_thetas.append(latitude_thetal)
        if (0 <= latitude_thetar < 180 and hi_r < hibnd_r
            and hi_z > lobnd_z and hi_z < hibnd_z):
            latitude_thetas.append(latitude_thetar)

        # If outside the frame, then breakout
        if ((lo_r < lobnd_r or lo_z < lobnd_z or lo_z > hibnd_z) and
            (hi_r >= hibnd_r or hi_z >= hibnd_z or hi_z <= lobnd_z)):
            break

    # Now annotate latitude lines and do the labeling.
    for latitude_theta in latitude_thetas:
        # Break out if we don't want to annotate latitude lines.
        if not annotate_lat_lines:
            break
        
        if latitude_theta == 0:
            continue

        latitude_radian = math.radians(latitude_theta)
        linewidth = 2.0 if not latitude_theta % 5 else 1.0

        # Label latitude
        if show_full_star:
            r_label = rl - 3*dr
        else:
            r_label = rl - 0.15*dr

        sp.annotate_text([r_label*np.sin(latitude_radian),
                          r_label*np.cos(latitude_radian)],
                         f"{int(90 - latitude_theta)}\u00B0",
                         text_args={"color": "silver",
                                    "size": "12",
                                    "family": "monospace",
                                    "horizontalalignment": "center",
                                    "verticalalignment": "center",
                                    },
                         inset_box_args={
                             "boxstyle": "square,pad=0.3",
                             "facecolor": "white",
                             "edgecolor": "white",
                         },
                         coord_system="plot")

        # Find the upper and lower bound of the latitude lines
        if latitude_radian > 0.5*np.pi:
            rll = max(lobnd_r / np.sin(latitude_radian),
                      hibnd_z / np.cos(latitude_radian))
            rrr = min(hibnd_r / np.sin(latitude_radian),
                      lobnd_z / np.cos(latitude_radian))
        else:
            rll = max(lobnd_r / np.sin(latitude_radian),
                      lobnd_z / np.cos(latitude_radian))
            rrr = min(hibnd_r / np.sin(latitude_radian),
                      hibnd_z / np.cos(latitude_radian))

        # First do a line to the lower half of the shell
        sp.annotate_line([rll*np.sin(latitude_radian),
                          rll*np.cos(latitude_radian)],
                         [rl*np.sin(latitude_radian),
                          rl*np.cos(latitude_radian)],
                         coord_system="plot",
                         color="k",
                         alpha=0.2,
                         linewidth=linewidth,
                         linestyle="-")
        
        # Then do a line to the upper half of the shell
        sp.annotate_line([rr*np.sin(latitude_radian),
                          rr*np.cos(latitude_radian)],
                         [rrr*np.sin(latitude_radian),
                          rrr*np.cos(latitude_radian)],
                         coord_system="plot",
                         color="k",
                         alpha=0.2,
                         linewidth=linewidth,
                         linestyle="-")
    
    sp._setup_plots()
    return sp


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
        A slice plot script for xrb_spherical problem.
        Given a list of plotfiles or a list of field parameters,
        it plots multiple slice plots.
        """)

    parser.add_argument('fname', type=str,
                        help="""dataset file names for plotting. Accepts one or more datasets.
                        If multiple file names are given, a grid of slice plots of different
                        files will be plotted for a given field parameter.
                        Note that either fnames or field must be single valued.""")
    parser.add_argument('-f', '--field', type=str,
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
    parser.add_argument('-r', '--dr', type=float,
                        help="""Distance between upper r and lower r shown in the SlicePlot.
                        Assumed in unit km. This is used to control center and width of the SlicePlot""")
    parser.add_argument('-w', '--width', default=3.0, type=float,
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

    full_star_slice = single_slice(args.fname, args.field, loc=loc,
                            widthScale=args.width, dr=args.dr, theta=args.theta,
                            displace_theta=args.displace_theta, annotate_vline=args.annotate_vline,
                            annotate_lat_lines=args.annotate_lat_lines, show_full_star=True)
    # full_star_slice.render()

    fig = full_star_slice.plots[args.field].figure
    rect = (0.24, 0.35, 0.4, 0.4)  # Left, Bottom, Width, Height
    inset_ax = fig.add_axes(rect)
    
    zoom_in_slice = single_slice(args.fname, args.field, loc=loc,
                            widthScale=args.width, dr=args.dr, theta=args.theta,
                            displace_theta=args.displace_theta, annotate_vline=args.annotate_vline,
                            annotate_lat_lines=args.annotate_lat_lines, show_full_star=False)
    zoom_in_slice.hide_colorbar()
    # zoom_in_slice.render()
    fig_zoom = zoom_in_slice.export_to_mpl_figure((1,1))
    zoom_in_slice.render()

    inset_img = fig_zoom.axes[0].images[0]
    inset_ax.imshow(
        inset_img.get_array(),
        extent=inset_img.get_extent(),
        origin=inset_img.origin,
        cmap=inset_img.get_cmap(),
        norm=inset_img.norm
    )

    # Add back annotation lines, i.e. latitude liens
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
    mark_inset(fig.axes[0], inset_ax, loc1=2, loc2=1, fc="none", ec="red",
               linestyle="--", linewidth=1.0)
    
    # fig_inset.savefig("zoom_in_slice.png", format='png', bbox_inches='tight')
    fig.savefig("full_star_slice.png", format='png', bbox_inches='tight')
