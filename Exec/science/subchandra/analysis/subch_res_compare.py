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


clip_val = -35
max_val = -19

# how much to coarsen for the contouring
blocking_factor = 8

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

def doit(rows, field):

    nrows = len(rows)
    ncols = len(rows[0][0])


    fig = plt.figure()

    grid = ImageGrid(fig, 111, nrows_ncols=(nrows, ncols),
                     axes_pad=0.25, cbar_pad=0.05, label_mode="L", cbar_mode="single")


    i = 0

    for row in rows:

        plotfiles, label = row

        for irow, pf in enumerate(plotfiles):

            ds = CastroDataset(pf)

            if field == "lap_rho":
                ds.force_periodicity()
                ds.add_field(name=("gas", "lap_rho"), sampling_type="local",
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

            sp.annotate_text((0.05, 0.05), f"time = {float(ds.current_time):8.3f} s", coord_system="axis", text_args={"color": text_color, "fontsize": "12"})
            if (irow == 0):
                sp.annotate_text((0.05, 0.925), f"{label}", coord_system="axis", text_args={"color": text_color, "fontsize": "14"})

            if (irow == 0):
                sp.annotate_grids(max_level=10, cmap="tab10")

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

            sp.set_axes_unit("km")

            plot = sp.plots[field]
            plot.figure = fig
            plot.axes = grid[i].axes
            plot.cax = grid.cbar_axes[i]
            if irow < len(plotfiles)-1:
                grid[i].axes.xaxis.offsetText.set_visible(False)

            sp._setup_plots()

            i += 1

    fig.set_size_inches(10, 14)
    plt.tight_layout()
    plt.savefig(f"subch_{field}_res_compare.pdf")

if __name__ == "__main__":


    subch_40km = [("subch_sdc_40km/subch_plt02123",
                   "subch_sdc_40km/subch_plt04184",
                   "subch_sdc_40km/subch_plt08509",
                   "subch_sdc_40km/subch_plt17272"), "40 km"]

    subch_20km = [("subch_sdc/subch_plt02123",
                   "subch_sdc/subch_plt04196",
                   "subch_sdc/subch_plt08614",
                   "subch_sdc/subch_plt17412"), "20 km"]

    subch_10km = [("subch_sdc_10km_3lev/subch_plt02123",
                   "subch_sdc_10km_3lev/subch_plt04197",
                   "subch_sdc_10km_3lev/subch_plt08582",
                   "subch_sdc_10km_3lev/subch_plt17526"), "10 km"]

    subch_5km = [("subch_sdc_5km_4lev/subch_plt02123",
                   "subch_sdc_5km_4lev/subch_plt04197",
                   "subch_sdc_5km_4lev/subch_plt08581",
                   "subch_sdc_5km_4lev/subch_plt17472"), "5 km"]

    field = "Temp"

    doit([subch_40km, subch_20km, subch_10km, subch_5km], field)
