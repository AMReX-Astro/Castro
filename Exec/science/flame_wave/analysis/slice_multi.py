#!/usr/bin/env python3

import os
import sys
import yt
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import ImageGrid

# assume that our data is in CGS
from yt.units import cm

plotfile = sys.argv[1]

# slice plot of temperature
ds = yt.load(plotfile)

xctr = 0.5*(ds.domain_left_edge[0] + ds.domain_right_edge[0])
L_x = ds.domain_right_edge[0] - ds.domain_left_edge[0]

ymin = 0.0*cm
ymax = 1.5e4*cm

yctr = 0.5*(ymin + ymax)
L_y = ymax - ymin


fig = plt.figure()
fig.set_size_inches(12.0, 9.0)

grid = ImageGrid(fig, 111, nrows_ncols=(3, 1), axes_pad=0.25, label_mode="L", cbar_mode="each")


fields = ["Temp", "enuc", "z_velocity"] #, "density"]

for i, f in enumerate(fields):

    sp = yt.SlicePlot(ds, "theta", f, center=[xctr, yctr, 0.0], width=[L_x, L_y, 0.0], fontsize="12")
    sp.set_buff_size((2000,2000))

    if f == "Temp":
        sp.set_zlim(f, 5.e7, 1.5e9)
        sp.set_cmap(f, "magma_r")
    elif f == "enuc":
        sp.set_zlim(f, 1.e18, 1.e20)
    elif f == "density":
        sp.set_zlim(f, 1.e-3, 5.e8)
    elif f == "z_velocity":
        sp.set_zlim(f, -2.e8, 2.e8)
        sp.set_log(f, False)
        sp.set_cmap(f, "bwr")

    #if f != "density":
    #    # now do a contour of density
    #    sp.annotate_contour("density", ncont=2, clim=(1.e2, 2.e6),
    #                        plot_args={"colors": "0.5", "linewidths": 1, "linestyle": ":"})

    sp.set_axes_unit("cm")

    #sp.annotate_text((0.05, 0.05), "{:8.5f} s".format(float(ds.current_time.in_cgs())),
    #                 coord_system="figure", text_args={"color": "black"})

    plot = sp.plots[f]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]
    if i < len(fields)-1:
        grid[i].axes.xaxis.offsetText.set_visible(False)

    sp._setup_plots()

fig.set_size_inches(10.0, 9.0)
plt.tight_layout()
plt.savefig(f"{os.path.basename(plotfile)}_slice.png")

