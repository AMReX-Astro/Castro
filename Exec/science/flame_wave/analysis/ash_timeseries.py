#!/usr/bin/env python

import argparse
import os
import re
import sys
from functools import reduce

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid

import yt
from yt.frontends.boxlib.api import CastroDataset
# assume that our data is in CGS
from yt.units import amu, cm

matplotlib.use('agg')




# define ash
def _ash(field, data):
    """ash is anything beyond O, excluding Fe and Ni"""

    ash_sum = None
    for f in data.ds.field_list:
        field_name = f[-1]
        # matches names like "X(ne21)" or "X(He4)"
        m = re.match(r"^X\(([A-Za-z]+)(\d+)\)$", field_name)
        if m is None:
            continue
        element = m[1].lower()
        aion = int(m[2])
        if element not in {"h", "he", "c", "n", "o", "fe", "ni"}:
            if ash_sum is None:
                ash_sum = data[f]
            else:
                ash_sum += data[f]
    if ash_sum is None:
        raise ValueError("no ash found")
    return ash_sum


def doit(plotfiles):

    ds = CastroDataset(plotfiles[0])

    xmin = ds.domain_left_edge[0]
    xmax = ds.domain_right_edge[0]
    xctr = 0.5*(xmin + xmax)
    L_x = xmax - xmin

    ymin = 0.0*cm
    ymax = 2.5e4*cm

    yctr = 0.5*(ymin + ymax)
    L_y = ymax - ymin

    fig = plt.figure()
    fig.set_size_inches(12.0, 9.0)

    grid = ImageGrid(fig, 111, nrows_ncols=(len(plotfiles), 1),
                     axes_pad=0.25, label_mode="L", cbar_mode="single", cbar_size="0.5%")


    for i, pf in enumerate(plotfiles):

        ds = CastroDataset(pf)

        ds.add_field(("gas", "ash"), function=_ash, display_name="ash", units="(dimensionless)", sampling_type="cell")

        field = "ash"
        sp = yt.SlicePlot(ds, "theta", field, center=[xctr, yctr, 0.0*cm], width=[L_x, L_y, 0.0*cm], fontsize="12")
        sp.set_buff_size((2400,2400))
        sp.set_zlim(field, 1.e-5, 0.1)
        sp.set_log(field, True)
        sp.set_cmap(field, "plasma_r")

        sp.set_axes_unit("cm")

        sp.annotate_text((0.85, 0.8), f"{1000.0*float(ds.current_time.in_cgs()):5.2f} ms",
                         coord_system="axis", text_args={"color": "black", "size": 10})

        plot = sp.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        if i < len(plotfiles)-1:
            grid[i].axes.xaxis.offsetText.set_visible(False)

        sp._setup_plots()

    fig.set_size_inches(10.0, 15.0)
    fig.subplots_adjust(left=0.075, right=0.925, top=0.95, bottom=0.05)
    fig.savefig("time_series_ash.pdf")

if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--skip", type=int, default=1,
                   help="interval between plotfiles")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")

    args = p.parse_args()

    plot_prefix = args.plotfiles[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.plotfiles], key=int)

    plotfiles = []
    for n in range(0, len(plot_nums), args.skip):
        plotfiles.append(f"{plot_prefix}{plot_nums[n]}")

    doit(plotfiles)
