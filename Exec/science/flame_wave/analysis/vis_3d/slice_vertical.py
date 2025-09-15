import argparse
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid

import yt
from yt.frontends.boxlib.api import CastroDataset
# assume that our data is in CGS
from yt.units import amu, cm

matplotlib.use('agg')




def doit(plotfiles):

    plotfile = sys.argv[1]
    ds = CastroDataset(plotfile)

    xmin = ds.domain_left_edge[0]
    xmax = ds.domain_right_edge[0]
    xctr = 0.5 * (xmin + xmax)
    L_x = (2./3.) * (xmax - xmin)

    ymin = ds.domain_left_edge[1]
    ymax = ds.domain_right_edge[1]
    yctr = 0.5*(ymin + ymax)

    zmin = 0.0*cm
    zmax = 1.e4*cm

    zctr = 0.5*(zmin + zmax)
    L_z = zmax - zmin

    fig = plt.figure()

    grid = ImageGrid(fig, 111, nrows_ncols=(len(plotfiles), 1),
                     axes_pad=0.25, label_mode="L",
                     cbar_mode="single", cbar_size="0.5%")

    for i, pf in enumerate(plotfiles):

        ds = CastroDataset(pf)

        f = "abar"

        sp = yt.SlicePlot(ds, "y", f, origin="native", center=[xctr, yctr, zctr],
                          width=[L_z, L_x], fontsize="9")
        sp.set_buff_size((4800,4800))
        sp.swap_axes()

        sp.set_zlim(f, 4, 5)
        sp.set_log(f, False)
        sp.set_cmap(f, "plasma_r")

        sp.set_axes_unit("cm")

        sp.annotate_text((0.8, 0.9), f"{1000.0*float(ds.current_time.in_cgs()):5.2f} ms",
                         coord_system="axis", text_args={"color": "black", "size": 9})

        plot = sp.plots["abar"]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        if i < len(plotfiles)-1:
            grid[i].axes.xaxis.offsetText.set_visible(False)

        sp._setup_plots()

    fig.set_size_inches(10.0, 3.5)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    fig.savefig("time_series_3D.pdf")

if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--skip", type=int, default=1,
                   help="interval between plotfiles")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")

    args = p.parse_args()

    plot_prefix = args.plotfiles[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.plotfiles], key=int)

    plotfiles = []
    for n in range(0, len(plot_nums), args.skip):
        plotfiles.append(f"{plot_prefix}{plot_nums[n]}")

    doit(plotfiles)
