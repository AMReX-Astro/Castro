#!/usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse
import os
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np
from functools import reduce

from mpl_toolkits.axes_grid1 import ImageGrid

# assume that our data is in CGS
from yt.units import cm, amu
from yt.frontends.boxlib.api import CastroDataset

#times = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35]
times = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
#times = [0.0, 0.15, 0.3, 0.45]

def find_files(plist):

    mask = np.zeros(len(times))
    files_to_plot = []
    for pfile in plist:
        for k, t in enumerate(times):
            if mask[k]:
                continue
            print(pfile)
            ds = CastroDataset(pfile)
            if ds.current_time >= t:
                files_to_plot.append(pfile)
                mask[k] = 1.0

    return files_to_plot

def doit(field, pfiles):

    print("looking to plot: ", pfiles)


    fig = plt.figure()

    if len(pfiles) > 4:
        nrows = 2
        ncols = (len(pfiles) + 1)//2
    else:
        nrows = 1
        ncols = len(pfiles)

    grid = ImageGrid(fig, 111, nrows_ncols=(nrows, ncols),
                     axes_pad=0.75, cbar_pad=0.05, label_mode="L", cbar_mode="single")


    for i in range(nrows * ncols):

        if i < len(pfiles):
            pf = pfiles[i]
        else:
            grid[i].remove()
            continue

        ds = CastroDataset(pf)

        domain_frac = 0.1

        xmin = ds.domain_left_edge[0]
        xmax = domain_frac * ds.domain_right_edge[0]
        xctr = 0.5 * (xmin + xmax)
        L_x = xmax - xmin

        ymin = ds.domain_left_edge[1]
        ymax = ds.domain_right_edge[1]
        yctr = 0.5 * (ymin + ymax)
        L_y = ymax - ymin
        ymin = yctr - 0.5 * domain_frac * L_y
        ymax = yctr + 0.5 * domain_frac * L_y
        L_y = ymax - ymin

        sp = yt.SlicePlot(ds, "theta", field, center=[xctr, yctr, 0.0*cm], width=[L_x, L_y, 0.0*cm], fontsize="12")
        sp.set_buff_size((2400,2400))
        sp.annotate_text((0.05, 0.05), f"time = {float(ds.current_time):8.3f} s", coord_system="axis", text_args={"color": "black"})

        if field == "Temp":
            sp.set_zlim(field, 5.e7, 4e9)
            sp.set_cmap(field, "magma_r")
        elif field == "enuc":
            sp.set_log(field, True, linthresh=1.e18)
            sp.set_zlim(field, -1.e22, 1.e22)
            sp.set_cmap(field, "bwr")
        elif field == "abar":
            sp.set_zlim(field, 4, 28)
            sp.set_log(field, False)
            sp.set_cmap(field, "plasma_r")

        sp.set_axes_unit("km")

        plot = sp.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        if i < len(pfiles)-1:
            grid[i].axes.xaxis.offsetText.set_visible(False)

        sp._setup_plots()

    fig.set_size_inches(19.2, 10.8)
    plt.tight_layout()
    plt.savefig(f"subch_{field}_sequence.png")

if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--var", type=str, default="Temp",
                   help="variable to plot")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")

    args = p.parse_args()
    print(args)
    plist = find_files(args.plotfiles)

    doit(args.var, plist)
