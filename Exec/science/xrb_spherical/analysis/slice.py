#!/usr/bin/env python3

import sys
import os
import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.frontends.boxlib.api import CastroDataset

from yt.units import cm

"""
Given a plot file and field name, it gives slice plots at the top,
middle, and bottom of the domain (shell).
"""

def slice(fname:str, field:str,
          loc: str = "top", width_factor: float = 3.0) -> None:
    """
    A slice plot of the dataset for Spherical2D geometry.

    Parameter
    =======================
    fname: plot file name
    field: field parameter
    loc: location on the domain. {top, mid, bot}
    """

    ds = CastroDataset(fname)
    currentTime = ds.current_time.in_units("s")
    print(f"Current time of this plot file is {currentTime} s")

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
    width = width_factor * dr
    box_widths = (width, width)

    loc = loc.lower()
    loc_options = ["top", "mid", "bot"]

    if loc not in loc_options:
        raise Exception("loc parameter must be top, mid or bot")

    # Centers for the Top, Mid and Bot panels
    centers = {"top":(r_center*np.sin(thetal)+0.5*width, r_center*np.cos(thetal)),
               "mid":(r_center*np.sin(theta_center), r_center*np.cos(theta_center)),
               "bot":(r_center*np.sin(thetar)+0.5*width, r_center*np.cos(thetar))}

    # Note we can also set center during SlicePlot, however then we would enter in [r_center, theta_center, 0]
    # rather than the physical R-Z coordinate if we do it via sp.set_center

    sp = yt.SlicePlot(ds, 'phi', field, width=box_widths)
    sp.set_center(centers[loc])

    sp.set_cmap(field, "viridis")
    if field in ["x_velocity", "y_velocity", "z_velocity"]:
        sp.set_cmap(field, "coolwarm")
    elif field == "Temp":
        sp.set_zlim(field. 5.e7, 2.5e9)
        sp.set_cmap(field, "magma_r")
    elif field == "enuc":
        sp.set_zlim(field, 1.e18, 1.e20)
    elif field == "density":
        sp.set_zlim(field, 1.e-3, 5.e8)

    # sp.annotate_text((0.05, 0.05), f"{currentTime.in_cgs():8.5f} s")
    sp.save(f"{ds}_{loc}")

if __name__ == "__main__":

    if len(sys.argv) < 3:
        raise Exception("Please enter parameters in order of: fname field_name width_factor[optional] loc[optional]")

    fname = sys.argv[1]
    field = sys.argv[2]
    loc = "top"
    width_factor = 3.0

    if len(sys.argv) == 4:
        width_factor = float(sys.argv[3])
    elif len(sys.argv) > 4:
        loc = sys.argv[4]

    slice(fname, field, loc=loc, width_factor=width_factor)
