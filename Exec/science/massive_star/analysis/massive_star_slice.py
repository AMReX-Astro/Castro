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

plotfile = sys.argv[1]
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
fig.set_size_inches(12.0, 9.0)

width_frac = 0.1

fields = ["Ye", "MachNumber", "enuc"]

grid = ImageGrid(fig, 111, nrows_ncols=(1, len(fields)),
                 axes_pad=1.0, label_mode="L", cbar_mode="each", cbar_pad=0)


for i, f in enumerate(fields):

    sp = yt.SlicePlot(ds, "theta", f,
                      center=[xmin + 0.25*width_frac*L_x, yctr, 0.0*cm],
                      width=[0.5*width_frac*L_x, width_frac*L_y, 0.0*cm], fontsize="12")
    sp.set_buff_size((2400,2400))

    if f == "Ye":
        sp.set_zlim(f, 0.46, 0.5)
        sp.set_log(f, False)
        sp.set_cmap(f, "magma_r")
    elif f == "enuc":
        sp.set_log(f, True, linthresh=1.e12)
        sp.set_zlim(f, -1.e20, 1.e20)
        sp.set_cmap(f, "bwr")
    elif f == "MachNumber":
        sp.set_zlim(f, 1.e-4, 0.3)
        sp.set_cmap(f, "plasma_r")

    if f == "enuc":
        # now do a contour of density
        sp.annotate_contour("in_nse", ncont=1, clim=(0.5, 0.5), take_log=False,
                            plot_args={"colors": "k", "linewidths": 2})

    sp.set_axes_unit("cm")

    #sp.annotate_text((0.05, 0.05), "{:8.5f} s".format(float(ds.current_time.in_cgs())),
    #                 coord_system="figure", text_args={"color": "black"})

    plot = sp.plots[f]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    #if i < len(fields)-1:
    #    grid[i].axes.xaxis.offsetText.set_visible(False)
    #print(type(grid[i].axes))

    sp._setup_plots()

    sp.plots[f].axes.xaxis.set_major_locator(plt.MaxNLocator(4))
    sp.plots[f].axes.ticklabel_format(axis="both", style="scientific", scilimits=(0,0))

fig.text(0.02, 0.02, "time = {:8.5f} s".format(float(ds.current_time)), transform=fig.transFigure)

fig.set_size_inches(19.2, 10.8)
fig.tight_layout()
fig.savefig("{}_slice.png".format(os.path.basename(plotfile)))
