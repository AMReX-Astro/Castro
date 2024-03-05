#!/usr/bin/env python3

import argparse
import os
import sys
from functools import reduce

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid

import yt
from yt.fields.derived_field import ValidateSpatial
from yt.frontends.boxlib.api import CastroDataset
from yt.funcs import just_one
# assume that our data is in CGS
from yt.units import amu, cm

matplotlib.use('agg')


#times = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35]
times = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
#times = [0.0, 0.15, 0.3, 0.45]

clip_val = -35
max_val = -19

# how much to coarsen for the contouring
blocking_factor = 8

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

def _lap_rho(field, data):
    dr = just_one(data["index", "dr"]).d
    r = data["index", "r"].d
    rl = r - 0.5 * dr
    rr = r + 0.5 * dr

    dz = just_one(data["index", "dz"]).d
    dens = data["gas", "density"].d

    _lap = np.zeros_like(dens)

    lapl_field = data.ds.arr(np.zeros(dens.shape, dtype=np.float64), None)

    # r-component
    _lap[1:-1, :] = 1 / (r[1:-1, :] * dr**2) * (
        - 2.0 * r[1:-1, :] * dens[1:-1:, :] +
        rl[1:-1, :] * dens[:-2, :] + rr[1:-1, :] * dens[2:, :])

    _lap[:, 1:-1] += 1 / dz**2 * (dens[:, 2:] + dens[:, :-2] - 2.0 * dens[:, 1:-1])
    lapl_field[1:-1, 1:-1] = np.log(np.abs(_lap[1:-1, 1:-1] / dens[1:-1, 1:-1]))
    lapl_field[lapl_field < clip_val] = clip_val
    return lapl_field

def doit(field, add_contours, pfiles):

    print("looking to plot: ", pfiles)


    fig = plt.figure()

    if len(pfiles) > 8:
        nrows = 3
        ncols = (len(pfiles) + 1)//3
    elif len(pfiles) > 4:
        nrows = 2
        ncols = (len(pfiles) + 1)//2
    else:
        nrows = 1
        ncols = len(pfiles)

    grid = ImageGrid(fig, 111, nrows_ncols=(nrows, ncols),
                     axes_pad=0.3, cbar_pad=0.05, label_mode="L", cbar_mode="single")


    for i in range(nrows * ncols):

        if i < len(pfiles):
            pf = pfiles[i]
        else:
            grid[i].remove()
            continue

        ds = CastroDataset(pf)

        if field == "lap_rho" or add_contours:
            ds.force_periodicity()
            ds.add_field(name=("gas", "lap_rho"), sampling_type="local",
                         display_name=r"$\log_{10}(|\rho^{-1}\nabla^2\rho|)$",
                         function=_lap_rho, units="",
                         validators=[ValidateSpatial(1)])

        domain_frac = 0.2

        xmin = ds.domain_left_edge[0]
        xmax = domain_frac * ds.domain_right_edge[0]
        xctr = 0.5 * (xmin + xmax)
        L_x = xmax - xmin

        ymin = ds.domain_left_edge[1]
        ymax = ds.domain_right_edge[1]
        yctr = 0.5 * (ymin + ymax)
        L_y = ymax - ymin
        ymin = yctr - 0.5 * domain_frac * L_y
        ymax = yctr + 0.5 * domain_frac * L_y
        L_y = ymax - ymin

        sp = yt.SlicePlot(ds, "theta", field, center=[xctr, yctr, 0.0*cm], width=[L_x, L_y, 0.0*cm], fontsize="14")
        sp.set_buff_size((2400,2400))
        if field == "Temp":
            text_color = "white"
        else:
            text_color = "black"

        sp.annotate_text((0.05, 0.05), f"time = {float(ds.current_time):8.3f} s", coord_system="axis", text_args={"color": text_color})

        if field == "Temp":
            sp.set_zlim(field, 5.e7, 4e9)
            sp.set_cmap(field, "magma")
        elif field == "enuc":
            sp.set_log(field, True, linthresh=1.e15)
            sp.set_zlim(field, -1.e22, 1.e22)
            sp.set_cmap(field, "bwr")
        elif field == "abar":
            sp.set_zlim(field, 4, 28)
            sp.set_log(field, False)
            sp.set_cmap(field, "plasma_r")
        elif field == "lap_rho":
            sp.set_zlim(field, clip_val, max_val)
            sp.set_log(field, False)
            sp.set_cmap(field, "bone_r")

        if add_contours:
            sp.annotate_contour(("gas", "lap_rho"), take_log=False, ncont=12, clim=(clip_val, max_val), factor=blocking_factor)

        sp.set_axes_unit("km")

        plot = sp.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        if i < len(pfiles)-1:
            grid[i].axes.xaxis.offsetText.set_visible(False)

        sp._setup_plots()

    fig.set_size_inches(11, 14)
    plt.tight_layout()
    plt.savefig(f"subch_{field}_sequence.pdf")

if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--var", type=str, default="Temp",
                   help="variable to plot")
    p.add_argument("--add_contours", action="store_true", help="add lap{rho}/rho contours")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")

    args = p.parse_args()
    print(args)
    plist = find_files(args.plotfiles)

    doit(args.var, args.add_contours,  plist)
