#!/usr/bin/env python3

import argparse
import os

import matplotlib
import matplotlib.pyplot as plt

import yt
from yt.frontends.boxlib.api import CastroDataset
# assume that our data is in CGS
from yt.units import cm

matplotlib.use('agg')




def make_plot(plotfile, prefix="plot", width_frac=0.1,
              layout=(1, 4),
              size=(19.2, 10.8), cbar_location="right"):

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

    fig = sp.export_to_mpl_figure((layout[0], layout[1]), axes_pad=(1.0, 0.4),
                                  cbar_location=cbar_location, cbar_pad="2%")

    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.025, top=0.975)
    fig.text(0.02, 0.02, f"time = {float(ds.current_time):8.2f} s",
             transform=fig.transFigure)
    fig.set_size_inches(size)
    fig.tight_layout()
    extra = ""
    if layout[0] >= layout[1]:
        extra = "_vertical"

    fig.savefig(f"{prefix}_{os.path.basename(plotfile)}_w{width_frac:04.2f}{extra}.pdf", pad_inches=0)


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("--vertical", action="store_true",
                   help="plot 2x2 or 1x4")
    p.add_argument("--width_fraction", type=float, default=0.1,
                   help="fraction of domain to show")
    p.add_argument("plotfile", type=str, nargs=1,
                   help="plotfile to plot")

    args = p.parse_args()
    plotfile = args.plotfile[0]

    if args.vertical:
        size = (7.5, 11.0)
        layout = (2, 2)
    else:
        size = (19.2, 8.5)
        layout = (1, 4)

    make_plot(plotfile, "all", layout=layout, width_frac=args.width_fraction, size=size)
