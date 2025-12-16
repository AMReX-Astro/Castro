#!/usr/bin/env python
import matplotlib
matplotlib.use('agg')

import os
import re
import sys
from functools import reduce

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import yt
from mpl_toolkits.axes_grid1 import ImageGrid
# assume that our data is in CGS
from yt.units import cm

matplotlib.use('agg')


# define ash derived field

def _ash(field, data):
    """ash is anything beyond O, excluding Fe and Ni"""
    return data["boxlib", "X(ash)"] * data["gas", "density"]


# load data

plotfile = sys.argv[1]
ds = yt.load(plotfile)

ds.add_field(("gas", "ash"), function=_ash,
             display_name=r"\rho X\left(ash\right)",
             units="auto", sampling_type="cell")

# figure out domain size and the buffer (set to be proportional to the
# number of zones)

xmin = ds.domain_left_edge[0]
xmax = ds.domain_right_edge[0]
xctr = 0.5 * (xmin + xmax)
L_x = xmax - xmin

ymin = 0.0 * cm
ymax = 1.5104e4 * cm
yctr = 0.5 * (ymin + ymax)
L_y = ymax - ymin

nx, ny, nz = ds.domain_dimensions
ref = int(np.prod(ds.ref_factors[0:ds.max_level]))

dx = ds.domain_width[0] / nx / ref
dy = ds.domain_width[1] / ny / ref

thinning = 2

pxls = int(L_x / dx / thinning), int(L_y / dy / thinning)


fig = plt.figure()


fields = ["Temp", "ash", "enuc", "z_velocity"] #, "density"]

grid = ImageGrid(fig, 111, nrows_ncols=(len(fields), 1),
                 axes_pad=0.25, label_mode="L", cbar_mode="each", cbar_size="1%")


for i, f in enumerate(fields):

    sp = yt.SlicePlot(ds, "theta", f, center=[xctr, yctr, 0.0*cm], width=[L_x, L_y, 0.0*cm], fontsize="12")
    sp.set_buff_size(pxls)

    if f == "Temp":
        sp.set_zlim(f, 5.e7, 1.5e9)
        sp.set_cmap(f, "magma")
    elif f == "enuc":
        sp.set_zlim(f, 1.e15, 3.e19)
        # set the background color to be white to show no burning
        sp.set_background_color("enuc", "white")
    elif f == "density":
        sp.set_zlim(f, 1.e-3, 5.e8)
    elif f == "z_velocity":
        sp.set_zlim(f, -2.e8, 2.e8)
        sp.set_log(f, False)
        sp.set_cmap(f, "bwr")
    elif f == "abar":
        sp.set_zlim(f, 4, 8)
        sp.set_log(f, False)
        sp.set_cmap(f, "plasma_r")
    elif f == "ash":
        sp.set_zlim(f, 1.e-2, 1e6)
        sp.set_log(f, True)
        sp.set_cmap(f, "plasma_r")

    sp.set_axes_unit("cm")

    plot = sp.plots[f]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]
    if i < len(fields)-1:
        grid[i].axes.xaxis.offsetText.set_visible(False)
    if i > 0:
        grid[i].axes.yaxis.offsetText.set_visible(False)

    sp._setup_plots()

fig.text(0.025, 0.05, f"time = {1000*float(ds.current_time):5.1f} ms",
         transform=fig.transFigure, fontsize="16")
fig.set_size_inches(19.2, 10.8)
fig.tight_layout()

fig.savefig(f"{os.path.basename(plotfile)}_slice.png")
