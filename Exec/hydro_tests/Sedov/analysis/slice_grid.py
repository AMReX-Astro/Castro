#!/usr/bin/env python

import sys
import yt
from matplotlib.colors import Normalize
from matplotlib.patches import (Rectangle, Polygon)
import matplotlib.colors as mcolors
from yt.frontends.boxlib.api import CastroDataset
import matplotlib.pyplot as plt
import numpy as np

"""
Give a pressure slice plot for 2d/3d sedov plot file
"""

fname = sys.argv[1]

outname = None
if len(sys.argv) == 3:
    outname = sys.argv[2]

ds = CastroDataset(fname)

slice_dirs = {"cartesian":"z",
              "cylindrical":"theta",
              "spherical":"phi"}
slice_dir = slice_dirs[ds.geometry]

f = "pressure"
cmap = "plasma_r"

slc = yt.SlicePlot(ds, slice_dir, f)
slc.set_buff_size(1600)
if ds.geometry == "cylindrical" or ds.geometry == "cartesian":
    slc.annotate_grids()
    slc.set_figure_size(12)
    slc.set_font_size(24)
    slc.set_cmap(f, cmap)
    slc.set_log(f, False)
    slc.set_axes_unit("cm")
    if outname is None:
        outname = "slice_grid.png"
    slc.save(outname)

else:
    # Spherical Case.
    # Let's extract data from SlicePlot using FRB and annotate grid
    # yt doesn't support annotating grids in spherical case.

    # Get the frb directly from the slice plot
    # This gives 2D array the shape as the dims
    slc.render()
    frb = slc.frb
    var = np.array(frb[f])

    # Get bounds and compute the coordinates
    R_min, R_max, Z_min, Z_max = [b.value for b in frb.bounds]
    nr, nth = var.shape
    R = np.linspace(R_min, R_max, nr + 1)
    Z = np.linspace(Z_min, Z_max, nth  + 1)

    fig, ax = plt.subplots(figsize=(8, 8))
    pcm = ax.pcolormesh(R, Z, var,
                        norm=mcolors.Normalize(vmin=np.nanmin(var),
                                               vmax=np.nanmax(var)),
                        cmap="plasma_r", shading="auto")

    # Now annotate grids
    level_colors = plt.get_cmap("Set1", ds.max_level + 1)
    for g in ds.index.grids:
        r_lo  = g.LeftEdge[0].value
        r_hi  = g.RightEdge[0].value

        th_lo = g.LeftEdge[1].value
        th_hi = g.RightEdge[1].value
        # sample curved edges
        npts = 100

        theta_outer = np.linspace(th_lo, th_hi, npts)
        theta_inner = np.linspace(th_hi, th_lo, npts)

        # outer arc
        R_outer = r_hi * np.sin(theta_outer)
        Z_outer = r_hi * np.cos(theta_outer)

        # inner arc
        R_inner = r_lo * np.sin(theta_inner)
        Z_inner = r_lo * np.cos(theta_inner)

        # combine into closed polygon
        R_poly = np.concatenate([R_outer, R_inner])
        Z_poly = np.concatenate([Z_outer, Z_inner])

        verts = np.column_stack([R_poly, Z_poly])

        poly = Polygon(verts, closed=True, linewidth=1.5,
                       edgecolor=level_colors(g.Level), facecolor="none")

        ax.add_patch(poly)

    # Deduplicated legend entries
    handles = [Rectangle((0,0), 1, 1, edgecolor=level_colors(l), facecolor="none", label=f"Level {l}")
               for l in range(ds.max_level + 1)]
    ax.legend(handles=handles, loc="upper right", framealpha=0.5)

    cbar = fig.colorbar(pcm, ax=ax)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label(r"$P \ (\frac{\mathrm{dyn}}{\mathrm{cm}^{2}})$", fontsize=16)
    cbar.minorticks_on()

    ax.tick_params(labelsize=14)
    ax.set_ylabel(r"$Z$ (cm)", size=16)
    ax.set_xlabel(r"$R$ (cm)", size=16)
    ax.minorticks_on()

    # Set equal scaling for R-Z coordinate
    ax.set_aspect("equal")

    fig.tight_layout()
    if outname is not None:
        fig.savefig(outname, format="png", bbox_inches="tight")
    else:
        plt.show()
