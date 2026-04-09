#!/usr/bin/env python

import os
import argparse
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np
from functools import reduce

from mpl_toolkits.axes_grid1 import ImageGrid

# assume that our data is in CGS
from yt.units import cm, amu

def doit(plotfiles):

    ds = yt.load(plotfiles[0], hint="castro")

    xmin = ds.domain_left_edge[0]
    xmax = 2*ds.domain_right_edge[0]/3
    xctr = 0.5*(xmin + xmax)
    L_x = xmax - xmin

    ymin = 0.0*cm
    ymax = 1.5e4*cm

    yctr = 0.5*(ymin + ymax)
    L_y = ymax - ymin

    fig = plt.figure()

    grid = ImageGrid(fig, 111, nrows_ncols=(len(plotfiles), 1),
                     axes_pad=0.25, label_mode="L",
                     cbar_mode="single", cbar_size="0.5%")

    for i, pf in enumerate(plotfiles):

        ds = yt.load(pf, hint="castro")

        field = "abar"
        sp = yt.SlicePlot(ds, "theta", field, center=[xctr, yctr, 0.0*cm],
                          width=[L_x, L_y, 0.0*cm], fontsize="9")
        sp.set_buff_size((2400,2400))

        sp.set_zlim(field, 4, 8)
        sp.set_log(field, False)
        sp.set_cmap(field, "plasma_r")

        sp.set_axes_unit("cm")

        sp.annotate_text((0.85, 0.8), f"{1000.0*float(ds.current_time.in_cgs()):5.2f} ms",
                         coord_system="axis", text_args={"color": "black", "size": 10})

        plot = sp.plots["abar"]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        if i < len(plotfiles)-1:
            grid[i].axes.xaxis.offsetText.set_visible(False)

        sp._setup_plots()

    fig.set_size_inches(9.5, 10.0)
    #fig.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    #plt.margins(0.0)
    fig.savefig("time_series.pdf")

if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--skip", type=int, default=50,
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
