#!/usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import os
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np
from functools import reduce

from mpl_toolkits.axes_grid1 import ImageGrid

# assume that our data is in CGS
from yt.units import cm, amu
from yt.fields.derived_field import ValidateSpatial
from yt.funcs import just_one
from yt.frontends.boxlib.api import CastroDataset


def _lap_rho(field, data):
    clip_val = -35

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


plotfile = sys.argv[1]
ds = CastroDataset(plotfile)

# Used to plot lap_rho plots
ds.force_periodicity()
ds.add_field(name=("gas", "lap_rho"), sampling_type="local",
             display_name=r"$\log_{10}(|\rho^{-1}\nabla^2\rho|)$",
             function=_lap_rho, units="",
             validators=[ValidateSpatial(1)])

domain_frac = 0.22

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

fig = plt.figure()

fields = ["Temp", "abar", "enuc", "lap_rho"]

grid = ImageGrid(fig, 111, nrows_ncols=(1, len(fields)),
                 axes_pad=0.75, cbar_pad=0.05, label_mode="L", cbar_mode="each")


for i, f in enumerate(fields):

    sp = yt.SlicePlot(ds, "theta", f, center=[xctr, yctr, 0.0*cm], width=[L_x, L_y, 0.0*cm], fontsize="12")
    sp.set_buff_size((2400,2400))

    if f == "Temp":
        sp.set_zlim(f, 5.e7, 6e9)
        sp.set_cmap(f, "magma")
    elif f == "enuc":
        sp.set_log(f, True, linthresh=1.e14)
        sp.set_zlim(f, -1.e22, 1.e22)
        sp.set_cmap(f, "bwr")
    elif f == "density":
        sp.set_zlim(f, 1.e-3, 5.e8)
    elif f == "z_velocity":
        sp.set_zlim(f, -2.e8, 2.e8)
        sp.set_log(f, False)
        sp.set_cmap(f, "bwr")
    elif f == "abar":
        sp.set_zlim(f, 4, 28)
        sp.set_log(f, False)
        sp.set_cmap(f, "plasma_r")
    elif f == "lap_rho":
        sp.set_zlim(f, -35, -19)
        sp.set_log(f, False)
        sp.set_cmap(f, "bone_r")


    #if f != "density":
    #    # now do a contour of density
    #    sp.annotate_contour("density", ncont=2, clim=(1.e2, 2.e6),
    #                        plot_args={"colors": "0.5", "linewidths": 1, "linestyle": ":"})

    if ("boxlib", "in_nse") in ds.derived_field_list:
        sp.annotate_contour("in_nse", levels=1, clim=(0.5, 0.5), take_log=False,
                            plot_args={"colors": "k", "linewidths": 2})

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

fig.text(0.02, 0.02, "time = {:8.5f} s".format(float(ds.current_time)), transform=fig.transFigure)

fig.set_size_inches(19.2, 10.8)
plt.tight_layout()
plt.savefig("{}_slice.png".format(os.path.basename(plotfile)))
