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

def doit(pf, size=(19.2, 10.8)):

    ds = CastroDataset(pf)

    xmin = ds.domain_left_edge[0]
    xmax = ds.domain_right_edge[0]
    xctr = 0.5 * (xmin + xmax)
    L_x = xmax - xmin

    ymin = ds.domain_left_edge[1]
    ymax = ds.domain_right_edge[1]
    yctr = 0.5*(ymin + ymax)

    zmin = ds.domain_left_edge[2]
    zmax = ds.domain_right_edge[2]
    zctr = 0.5*(zmin + zmax)
    L_z = zmax - zmin

    fields = ["magvort", "enuc"]

    sp = yt.SlicePlot(ds, "y", fields,
                      origin="native", center=[xctr, yctr, zctr],
                      width=[L_z, L_x], fontsize="20")
    sp.set_buff_size((4800,4800))
    sp.swap_axes()

    for f in fields:
        if f == "Temp":
            sp.set_zlim(f, 5.e5, 5.e8)
            sp.set_cmap(f, "magma_r")
        elif f == "enuc":
            sp.set_zlim(f, 1.e16, 1.e21)
            sp.set_log(f, True)
            sp.set_cmap(f, "plasma_r")
        elif f == "density":
            sp.set_zlim(f, 1.e-5, 5.e6)
        elif f == "abar":
            sp.set_zlim(f, 1, 8)
            sp.set_log(f, False)
            sp.set_cmap(f, "plasma_r")
        elif f == "magvort":
            sp.set_zlim(f, 2.0, 2.e5)
            sp.set_log(f, True)
        elif f.startswith(r"X("):
            sp.set_zlim(f, 1.e-5, 1.0)
            sp.set_log(f, True)

    sp.set_axes_unit("cm")

    layout = (1, len(fields))

    fig = sp.export_to_mpl_figure((layout[0], layout[1]), axes_pad=(1.0, 0.4),
                                  cbar_location="right", cbar_pad="2%")

    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.025, top=0.975)
    fig.text(0.05, 0.01, f"t = {float(ds.current_time) * 1000:6.2f} ms",
             fontsize="20",
             transform=fig.transFigure)
    fig.set_size_inches(size)
    fig.tight_layout()

    fig.savefig(f"slice_{pf}.png")

if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("plotfile", type=str, nargs=1,
                   help="plotfiles to plot")

    args = p.parse_args()

    doit(args.plotfile[0])
