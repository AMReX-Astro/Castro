import os
import sys
from functools import reduce

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid

import yt
from yt.units import cm

matplotlib.use('agg')


def doit(plotfile):

    ds = yt.load(plotfile)

    xmin = ds.domain_left_edge[0]
    xmax = ds.domain_right_edge[0]

    ymin = ds.domain_left_edge[1]
    ymax = ds.domain_right_edge[1]

    xctr = 0.0 * xmin
    L_x = xmax - xmin

    yctr = 0.5 * (ymin + ymax)
    L_y = ymax - ymin

    width_frac = 1./3.

    fig = plt.figure()

    fields = ["density", "Temp"]

    grid = ImageGrid(fig, 111, nrows_ncols=(1, len(fields)),
                     axes_pad=0.75, cbar_pad=0.05,
                     label_mode="L", cbar_mode="each")

    for i, field in enumerate(fields):

        slc = yt.SlicePlot(ds, "z", field,
                           center=[xctr, yctr, 0.0*cm],
                           width=[width_frac*L_x, width_frac*L_y, 0.0*cm],
                           fontsize=14)

        if field == "Temp":
            slc.set_zlim(field, 1.e7, 4e9)
            slc.set_cmap(field, "magma")
        elif field == "density":
            slc.set_zlim(field, 1.e-4, 5.e7)

        slc.set_buff_size((3072, 3072))
        slc.set_axes_unit("cm")

        plot = slc.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        if i < len(fields)-1:
            grid[i].axes.xaxis.offsetText.set_visible(False)

        slc._setup_plots()

    fig.text(0.02, 0.02, "time = {:8.5f} s".format(float(ds.current_time)), transform=fig.transFigure)

    fig.set_size_inches(19.2, 10.8)
    plt.tight_layout()
    plt.savefig("{}_slice.png".format(os.path.basename(plotfile)))


    slc.save(f"wdmerger_slice_grid_{plotfile}.png")


if __name__ == "__main__":

    pf = sys.argv[-1]
    doit(pf)
