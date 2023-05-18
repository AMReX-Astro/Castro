#!/usr/bin/env python3

import matplotlib
matplotlib.use('agg')

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

def make_plot(plotfile, prefix="plot", size=(19.2, 10.8), cbar_location="right"):

    ds = CastroDataset(plotfile)

    xmin = ds.domain_left_edge[0]
    xmax = ds.domain_right_edge[0]

    ymin = ds.domain_left_edge[1]
    ymax = ds.domain_right_edge[1]

    xctr = 0.0 * xmin
    L_x = xmax - xmin

    yctr = 0.5 * (ymin + ymax)
    L_y = ymax - ymin


    fig = plt.figure()
    fig.set_size_inches(size)

    width_frac = 0.1

    fields = ["MachNumber", "magvort", "abar", "enuc"]

    sp = yt.SlicePlot(ds, "theta", fields,
                      center=[xmin + 0.5*width_frac*L_x, yctr, 0.0*cm],
                      width=[width_frac*L_x, width_frac*L_y, 0.0*cm], fontsize="12")

    sp.set_buff_size((2400,2400))

    for f in fields:
        if f == "Ye":
            sp.set_zlim(f, 0.46, 0.5)
            sp.set_log(f, False)
            sp.set_cmap(f, "magma_r")
        elif f == "abar":
            sp.set_log(f, False)
            sp.set_cmap(f, "viridis")
        elif f == "enuc":
            sp.set_log(f, True, linthresh=1.e12)
            sp.set_zlim(f, -1.e20, 1.e20)
            sp.set_cmap(f, "bwr")
        elif f == "MachNumber":
            sp.set_zlim(f, 1.e-4, 0.3)
            sp.set_cmap(f, "plasma")
        elif f == "magvel":
            sp.set_zlim(f, 100.0, 2.e7)
            sp.set_cmap(f, "viridis")
        elif f == "magvort":
            sp.set_cmap(f, "magma")
            sp.set_zlim(f, 1.e-2, 5)

        if f == "enuc":
            # now do a contour of NSE
            sp.annotate_contour("in_nse", ncont=1, clim=(0.5, 0.5), take_log=False,
                                plot_args={"colors": "k", "linewidths": 2})

    sp.set_axes_unit("cm")

    fig = sp.export_to_mpl_figure((1, len(fields)), cbar_location=cbar_location, cbar_pad="5%")

    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.025, top=0.975)
    fig.text(0.02, 0.02, "time = {:8.5f} s".format(float(ds.current_time)), transform=fig.transFigure)
    fig.set_size_inches(size)

    fig.savefig(f"{prefix}_{os.path.basename(plotfile)}_slice.png", pad_inches=0)


if __name__ == "__main__":

    plotfile = sys.argv[1]

    make_plot(plotfile, "all")
