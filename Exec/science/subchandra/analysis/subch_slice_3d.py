#!/usr/bin/env python

import os
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




plotfile = sys.argv[1]
ds = CastroDataset(plotfile)

domain_frac = 0.2

xmin = ds.domain_left_edge[0]
xmax = ds.domain_right_edge[0]
xctr = 0.5 * (xmin + xmax)
L_x = xmax - xmin
xmin = xctr - 0.5 * domain_frac * L_x
xmax = xctr + 0.5 * domain_frac * L_x
L_x = xmax - xmin

ymin = ds.domain_left_edge[1]
ymax = ds.domain_right_edge[1]
yctr = 0.5 * (ymin + ymax)
L_y = ymax - ymin
ymin = yctr - 0.5 * domain_frac * L_y
ymax = yctr + 0.5 * domain_frac * L_y
L_y = ymax - ymin

zmin = ds.domain_left_edge[2]
zmax = ds.domain_right_edge[2]
zctr = 0.5 * (zmin + zmax)
L_z = zmax - zmin
zmin = zctr - 0.5 * domain_frac * L_z
zmax = zctr + 0.5 * domain_frac * L_z
L_z = zmax - zmin

fig = plt.figure()


fields = ["Temp", "abar", "enuc"]

grid = ImageGrid(fig, 111, nrows_ncols=(1, len(fields)),
                 axes_pad=0.75, cbar_pad=0.05, label_mode="L", cbar_mode="each")


for i, f in enumerate(fields):

    sp = yt.SlicePlot(ds, "y", f,
                      center=[xctr, yctr, zctr], width=[L_x, L_z],
                      fontsize="12")

    sp.set_buff_size((2400,2400))
    sp.swap_axes()

    if f == "Temp":
        sp.set_zlim(f, 5.e7, 4e9)
        sp.set_cmap(f, "magma")
    elif f == "enuc":
        sp.set_log(f, True, linthresh=1.e18)
        sp.set_zlim(f, -1.e22, 1.e22)
        sp.set_cmap(f, "bwr")
    elif f == "density":
        sp.set_zlim(f, 1.e-3, 5.e8)
    elif f == "z_velocity":
        sp.set_zlim(f, -2.e8, 2.e8)
        sp.set_log(f, False)
        sp.set_cmap(f, "bwr")
    elif f == "abar":
        sp.set_zlim(f, 4, 28)
        sp.set_log(f, False)
        sp.set_cmap(f, "plasma_r")

    sp.set_axes_unit("cm")

    plot = sp.plots[f]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]
    if i < len(fields)-1:
        grid[i].axes.xaxis.offsetText.set_visible(False)

    sp._setup_plots()

fig.text(0.02, 0.02, f"time = {float(ds.current_time):8.5f} s", transform=fig.transFigure)

fig.set_size_inches(19.2, 10.8)
plt.tight_layout()
plt.savefig(f"{os.path.basename(plotfile)}_slice.png")
