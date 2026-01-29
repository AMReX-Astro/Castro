#!/usr/bin/env python

import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid

import yt
from yt.frontends.boxlib.api import CastroDataset
# assume that our data is in CGS
from yt.units import cm

matplotlib.use('agg')

times = [50, 100, 150, 200, 250, 300]

def find_files(plist):

    mask = np.zeros(len(times))
    files_to_plot = []
    for pfile in plist:
        for k, t in enumerate(times):
            if mask[k]:
                continue
            print(pfile)
            ds = CastroDataset(pfile)
            if ds.current_time >= t:
                files_to_plot.append(pfile)
                mask[k] = 1.0

    return files_to_plot


def make_plot(field, pfiles, *,
              width_frac=0.1,
              size=(19.2, 10.8)):

    fig = plt.figure()

    if len(pfiles) > 8:
        nrows = 3
        ncols = (len(pfiles) + 1)//3
    elif len(pfiles) >= 4:
        nrows = 2
        ncols = (len(pfiles) + 1)//2
    else:
        nrows = 1
        ncols = len(pfiles)

    grid = ImageGrid(fig, 111, nrows_ncols=(nrows, ncols),
                     axes_pad=0.2, cbar_pad=0.05,
                     label_mode="L", cbar_mode="single")

    # get some metadata
    ds = CastroDataset(pfiles[0])

    t_drive = 0.0
    if "[*] castro.drive_initial_convection_tmax" in ds.parameters:
        t_drive = ds.parameters["[*] castro.drive_initial_convection_tmax"]
    elif "[*] castro.drive_initial_convection_tmax" in ds.parameters:
        t_drive = ds.parameters["[*] castro.drive_initial_convection_tmax"]
    print(t_drive)

    xmin = ds.domain_left_edge[0]
    xmax = ds.domain_right_edge[0]

    ymin = ds.domain_left_edge[1]
    ymax = ds.domain_right_edge[1]

    xctr = 0.0 * xmin
    L_x = xmax - xmin

    yctr = 0.5 * (ymin + ymax)
    L_y = ymax - ymin

    f = field

    if f == "Temp" or f == "abar":
        text_color = "white"
    else:
        text_color = "black"

    for i in range(nrows * ncols):

        if i < len(pfiles):
            pf = pfiles[i]
        else:
            grid[i].remove()
            continue

        ds = CastroDataset(pf)

        sp = yt.SlicePlot(ds, "theta", f,
                          center=[xmin + 0.5*width_frac*L_x, yctr, 0.0*cm],
                          width=[width_frac*L_x, width_frac*L_y, 0.0*cm],
                          fontsize="12")

        sp.set_buff_size((2400,2400))

        sp.annotate_text((0.05, 0.05),
                         f"$t - \\tau_\\mathrm{{drive}}$ = {float(ds.current_time) - t_drive:8.2f} s",
                         coord_system="axis",
                         text_args={"color": text_color})


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

        # now do a contour of NSE
        if ("boxlib", "in_nse") in ds.derived_field_list:
            sp.annotate_contour("in_nse", levels=1, clim=(0.5, 0.5), take_log=False,
                                plot_args={"colors": "k", "linewidths": 2})

        sp.set_axes_unit("cm")

        plot = sp.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        if i < len(pfiles)-1:
            grid[i].axes.xaxis.offsetText.set_visible(False)
        if i > 0:
            grid[i].axes.yaxis.offsetText.set_visible(False)

        sp._setup_plots()

    fig.set_size_inches(size)
    fig.tight_layout()
    extra = ""

    fig.savefig(f"massive_star_{f}_w{width_frac:04.2f}.pdf", pad_inches=0.1, bbox_inches="tight")


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("--var", type=str, default="Temp",
                   help="variable to plot")
    p.add_argument("--vertical", action="store_true",
                   help="plot 2x2 or 1x4")
    p.add_argument("--width-fraction", type=float, default=0.1,
                   help="fraction of domain to show")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")

    args = p.parse_args()

    if args.vertical:
        size = (7.0, 8.0)
    else:
        size = (19.2, 8.5)

    plist = find_files(args.plotfiles)

    make_plot(args.var, plist,
              width_frac=args.width_fraction, size=size)
