#!/usr/bin/env python3

import os
import sys
import yt
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.ticker import ScalarFormatter

# assume that our data is in CGS
from yt.units import cm

plotfile = sys.argv[1]

# slice plot of temperature
ds = yt.load(plotfile, hint="castro")

xctr = 0.5*(ds.domain_left_edge[0] + ds.domain_right_edge[0])
L_x = ds.domain_right_edge[0] - ds.domain_left_edge[0]

yctr = 0.5*(ds.domain_left_edge[1] + ds.domain_right_edge[1])
L_y = ds.domain_right_edge[1] - ds.domain_left_edge[1]


fig = plt.figure()

nrows = 2
ncols = 3

grid = ImageGrid(fig, 111, nrows_ncols=(nrows, ncols), axes_pad=0.75, cbar_pad="2%",
                 label_mode="L", cbar_mode="each")


fields = ["Ye", "abar", "enuc", "Temp", "MachNumber"]

for i in range(nrows * ncols):

    if i >= len(fields):
        grid[i].remove()
        grid.cbar_axes[i].remove()
        continue
    else:
        f = fields[i]

    sp = yt.SlicePlot(ds, "z", f, center=[xctr, yctr, 0.0*cm],
                      width=[L_x, L_y, 0.0*cm], fontsize="12")

    sp.set_buff_size((2000,2000))

    if f == "Ye":
        sp.set_zlim(f, 0.48, 0.5)
        sp.set_log(f, False)
        sp.set_cmap(f, "magma_r")
    elif f == "abar":
        sp.set_log(f, False)
        sp.set_cmap(f, "viridis")
    elif f == "enuc":
        sp.set_log(f, True, linthresh=1.e10)
        sp.set_zlim(f, -1.e18, 1.e18)
        sp.set_cmap(f, "bwr")
    elif f == "MachNumber":
        sp.set_zlim(f, 1.e-5, 0.1)
        sp.set_cmap(f, "plasma")
    elif f == "magvel":
        sp.set_zlim(f, 100.0, 2.e7)
        sp.set_cmap(f, "viridis")
    elif f == "Temp":
        sp.set_cmap(f, "magma")
        sp.set_zlim(f, 1.e8, 6.e9)

    sp.set_axes_unit("cm")

    plot = sp.plots[f]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]
    if i < len(fields)-1:
        grid[i].axes.xaxis.offsetText.set_visible(False)
    #fmt = grid.cbar_axes[i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    sp._setup_plots()



fig.set_size_inches(12.0, 9.0)
plt.tight_layout()
plt.savefig("{}_slice.png".format(os.path.basename(plotfile)))

