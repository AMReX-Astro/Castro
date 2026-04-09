#!/usr/bin/env python

import os
import sys

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

import yt
from yt.frontends.boxlib.api import CastroDataset

matplotlib.use('agg')


plotfile = sys.argv[1]
ds = CastroDataset(plotfile)

t_drive = 0.0
if "[*] castro.drive_initial_convection_tmax" in ds.parameters:
    t_drive = ds.parameters["[*] castro.drive_initial_convection_tmax"]
elif "castro.drive_initial_convection_tmax" in ds.parameters:
    t_drive = ds.parameters["castro.drive_initial_convection_tmax"]
print(t_drive)

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


fields = ["MachNumber", "magvort", "abar", "enuc"]

grid = ImageGrid(fig, 111, nrows_ncols=(2, (len(fields)+1)//2),
                 axes_pad=0.75, cbar_pad=0.05, label_mode="L", cbar_mode="each")


for i, f in enumerate(fields):

    sp = yt.SlicePlot(ds, "y", f,
                      center=[xctr, yctr, zctr], width=[L_x, L_z],
                      fontsize="12")

    sp.set_buff_size((2400,2400))
    sp.swap_axes()

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

    sp.set_axes_unit("km")

    plot = sp.plots[f]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]
    if i < len(fields)-1:
        grid[i].axes.xaxis.offsetText.set_visible(False)

    sp._setup_plots()

fig.text(0.02, 0.02, f"$t - \\tau_\\mathrm{{drive}}$ = {float(ds.current_time) - t_drive:6.1f} s",
         transform=fig.transFigure)

fig.set_size_inches(11, 10)
plt.tight_layout()
plt.savefig(f"{os.path.basename(plotfile)}_slice.png")
