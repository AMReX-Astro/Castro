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
from yt.units import km

def extract_info(ds,
                 loc: str = "top", widthScale: float = 3.0,
                 dr: Optional[float] = None,
                 theta: Optional[float] = None,
                 displace_theta: bool = True,
                 show_full_star: bool = True):
    '''
    Extracts relevant infos of a plotting script
    '''
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

    r = [rl, r_center, rr]

    thetar = ds.domain_right_edge[1]
    thetal = ds.domain_left_edge[1]
    theta_center = 0.5 * (thetar + thetal)

    # Domain width of the slice plot
    width = widthScale * dr
    box_widths = (width, width)

    # Now determine center of the frame
    if show_full_star:
        # If we want to show the full domain in the background.
        # Change box_widths and center to full star
        center = [r[1]*np.sin(theta_center), r[1]*np.cos(theta_center)]
        if thetar < 0.5 * np.pi:
            box_widths = (r[2]*np.sin(thetar), r[2]*np.cos(thetal))
        else:
            box_widths = ( r[2], r[2]*(np.abs(np.cos(thetar)) + np.cos(thetal)) )

    elif theta is None:
        # Preset centers for the Top, Mid and Bot panels
        # Centers will be physical coordinates in Cylindrical, i.e. R-Z
        centers = {"top":(r[1]*np.sin(thetal)+0.5*width, r[1]*np.cos(thetal)),
                   "mid":(r[1]*np.sin(theta_center), r[1]*np.cos(theta_center)),
                   "bot":(r[1]*np.sin(thetar)+0.5*width, r[1]*np.cos(thetar))}

        center = centers[loc]
    else:
        if displace_theta:
            # Optionally displace the center by ~0.7 of the plotting width
            # This is helpful when tracking the flame front.
            # This keeps the front at ~0.7 of the plotting width.

            # Determine dtheta that displaces from center to ~0.7 of the plotting domain
            oSevenTheta = np.arcsin(0.7 * width / r[1])
            halfTheta = np.arcsin(0.5 * width / r[1])
            dtheta = oSevenTheta - halfTheta
        else:
            dtheta = 0

        # Determine center using theta but also displace it by dtheta
        R = r[1]*np.sin(theta - dtheta)
        Z = r[1]*np.cos(theta - dtheta)
        if R < 0.5 * width:
            R = 0.5 * width
        center = [R, Z]

    return r, box_widths, center

def annotate_latitude_lines(sp, center, box_widths, r,
                            show_full_star=False):
    """
    Given the slice plot, annotate latitude lines
    Parameters
    ----------
    sp: yt SlicePlot
    center: center of the frame in [R, Z]
    box_widths: width of the slice plot frame.
    r: A list of r in [rl, r_center, rr]
    show_full_star: whether we're dealing with full star slice plot.
    """

    # Start from the theta center [deg] of the slice plot frame.
    thetac = round(math.degrees(np.arccos(center[1] / r[1])))
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

    for theta_increment in range(start, end, step):
        # This is actually in theta-coordinate.
        latitude_thetar = thetac + theta_increment
        latitude_thetal = thetac - theta_increment

        # Now find the RZ position of different latitude thetas
        # and see if they're out of the frame.
        l_r = r[0] * np.sin(math.radians(latitude_thetal))
        l_z = r[0] * np.cos(math.radians(latitude_thetal))

        r_r = r[0] * np.sin(math.radians(latitude_thetar))
        r_z = r[0] * np.cos(math.radians(latitude_thetar))

        # Check if the point is within the frame and append point
        if (0 <= latitude_thetal < 180 and
            lobnd_r <= l_r < hibnd_r and
            lobnd_z <= l_z < hibnd_z):
            latitude_thetas.append(latitude_thetal)

        if (0 <= latitude_thetal < 180 and
            lobnd_r <= r_r < hibnd_r and
            lobnd_z <= r_z < hibnd_z):
            latitude_thetas.append(latitude_thetar)

        # If outside the frame, then breakout
        if ((l_r < lobnd_r or l_r >= hibnd_r or
             l_z < lobnd_z or l_z >= hibnd_z) and
            (r_r < lobnd_r or r_r >= hibnd_r or
             r_z < lobnd_z or r_z >= hibnd_z)):
            break

    # Now annotate latitude lines and do the labeling.
    for latitude_theta in latitude_thetas:
        if latitude_theta == 0 or latitude_theta == 180:
            continue

        latitude_radian = math.radians(latitude_theta)
        linewidth = 2.0 if not latitude_theta % 5 else 1.0

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
        if rll < r[0]:
            sp.annotate_line([rll*np.sin(latitude_radian),
                              rll*np.cos(latitude_radian)],
                             [r[0]*np.sin(latitude_radian),
                              r[0]*np.cos(latitude_radian)],
                             coord_system="plot",
                             color="k",
                             alpha=0.2,
                             linewidth=linewidth,
                             linestyle="-")

        # Then do a line to the upper half of the shell
        if rrr > r[2]:
            sp.annotate_line([r[2]*np.sin(latitude_radian),
                              r[2]*np.cos(latitude_radian)],
                             [rrr*np.sin(latitude_radian),
                              rrr*np.cos(latitude_radian)],
                             coord_system="plot",
                             color="k",
                             alpha=0.2,
                             linewidth=linewidth,
                             linestyle="-")

        # Label latitude if we plotted the latitude line
        if rll < r[0] or rrr > r[2]:
            dr = r[2] - r[0]
            if show_full_star:
                r_label = r[0] - 3.0*dr
            else:
                r_label = r[0] - 0.2*dr

            sp.annotate_text([r_label*np.sin(latitude_radian),
                              r_label*np.cos(latitude_radian)],
                             f"{int(90 - latitude_theta)}\u00B0",
                             text_args={"color": "silver",
                                        "size": "12",
                                        "family": "monospace",
                                        "horizontalalignment": "center",
                                        "verticalalignment": "center",
                                        "clip_on": True,
                                        },
                             inset_box_args={
                                 "boxstyle": "round,pad=0.1",
                                 "facecolor": "white",
                                 "edgecolor": "white",
                             },
                             coord_system="plot")


def slice(fnames:List[str], fields:List[str],
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
    fnames:
      A list of file names to plot multiple slice plots between different
      plot files for a given field parameter.
      Note that either fname or field must be single valued.

    fields:
      A list of field parameters to plot multiple slice plots between different
      field parameters for a given file.
      Note that either fname or field must be single valued.

    loc:
      preset center location of the domain. {top, mid, bot}

    widthScale:
      scaling for the domain width of the slice plot

    dr:
      user defined distance between lower r to upper r boundary. Assumed in unit km.
      This is used to change the center and width of the SlicePlot.

    theta:
      user defined theta center of the slice plot

    displace_theta:
      whether to displace theta so that the vertical lines that represents
      the flame front is offset by some amount

    annotate_vline:
      whether to plot a vertical line to represent the flame front,
      which is represented by what theta is.

    annotate_lat_lines:
      whether to annotate latitude lines.

    show_full_star:
      whether to plot the full star rather than a zoom-in
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
    outName = "xrb_spherical_slice.png"

    for j, ds in enumerate(ts):
        # Process information for each dataset
        r, box_widths, center = extract_info(ds,
                                             loc=loc, widthScale=widthScale,
                                             dr=dr,
                                             theta=theta,
                                             displace_theta=displace_theta,
                                             show_full_star=show_full_star)
        for i, field in enumerate(fields):
            # Plot each field parameter
            sp = yt.SlicePlot(ds, 'phi', field, width=box_widths, fontsize=16)
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

            plot = sp.plots[field]
            plot.figure = fig
            plot.axes = grid[i+j*len(fields)].axes
            plot.cax = grid.cbar_axes[i+j*len(fields)]

            sp._setup_plots()

    if len(fnames) == 1:
        time = ts[0].current_time.in_units("ms")
        if float(time) < 1e-1:
            time = ts[0].current_time.in_units("us")

        # Determine position of the text on grid
        xyPositions = {(1, 1): (0.78, 0.02),
                       (1, 2): (0.95, 0.075),
                       (2, 2): (0.78, 0.02),
                       (2, 3): (0.9, 0.02),
                       (3, 3): (0.78, 0.02)
                       }
        xPosition, yPosition = xyPositions.get((nx, ny), (0.78, 0.02))

        fig.text(xPosition, yPosition, f"t = {time:.3f}", fontsize=16,
                 horizontalalignment='right', verticalalignment='bottom',
                 color="black", transform=fig.transFigure)

        outName = f"{ts[0]}_slice.png"

    fig.set_size_inches(16, 9)
    fig.tight_layout()
    fig.savefig(outName, format="png", bbox_inches="tight")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
        A slice plot script for xrb_spherical problem.
        Given a list of plotfiles or a list of field parameters,
        it plots multiple slice plots.
        """)

    parser.add_argument('fnames', nargs='+', type=str,
                        help="""dataset file names for plotting. Accepts one or more datasets.
                        If multiple file names are given, a grid of slice plots of different
                        files will be plotted for a given field parameter.
                        Note that either fnames or field must be single valued.""")
    parser.add_argument('-f', '--fields', nargs='+', type=str,
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
    parser.add_argument('--show_full_star', action='store_true',
                        help="whether show the full star in the background")

    args = parser.parse_args()

    if len(args.fnames) > 1 and len(args.fields) > 1:
        parser.error("Either fnames or fields must be single valued!")

    loc = args.loc.lower()
    loc_options = ["top", "mid", "bot"]

    if loc not in loc_options:
        parser.error("loc must be one of the three: {top, mid, bot}")

    slice(args.fnames, args.fields, loc=loc,
          widthScale=args.width, dr=args.dr, theta=args.theta,
          displace_theta=args.displace_theta, annotate_vline=args.annotate_vline,
          annotate_lat_lines=args.annotate_lat_lines, show_full_star=args.show_full_star)
