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
ds = yt.load(plotfile)

xctr = 0.5*(ds.domain_left_edge[0] + ds.domain_right_edge[0])
L_x = ds.domain_right_edge[0] - ds.domain_left_edge[0]

yctr = 0.5*(ds.domain_left_edge[1] + ds.domain_right_edge[1])
L_y = ds.domain_right_edge[1] - ds.domain_left_edge[1]


fig = plt.figure()
fig.set_size_inches(12.0, 9.0)

grid = ImageGrid(fig, 111, nrows_ncols=(1, 2), axes_pad=0.75, cbar_pad="2%",
                 label_mode="L", cbar_mode="each")


fields = ["pressure", "magvel"]

for i, f in enumerate(fields):

    sp = yt.SlicePlot(ds, "z", f, center=[xctr, yctr, 0.0], width=[L_x, L_y, 0.0], fontsize="12")
    sp.set_buff_size((2000,2000))

    if f == "pressure":
        sp.set_log(f, False)

    elif f == "magvel":
        sp.set_log(f, False)
        #sp.set_zlim(f, 1.e-3, 2.5e-2)
        sp.set_cmap(f, "magma")

    sp.set_axes_unit("cm")

    plot = sp.plots[f]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]
    if i < len(fields)-1:
        grid[i].axes.xaxis.offsetText.set_visible(False)
    #fmt = grid.cbar_axes[i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    sp._setup_plots()



fig.set_size_inches(10.0, 5.0)
plt.tight_layout()
plt.savefig("{}_slice.pdf".format(os.path.basename(plotfile)))

