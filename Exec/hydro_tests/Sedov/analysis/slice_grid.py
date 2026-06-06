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

    # # Get refinement ratio between finest to coarest level
    # ref_ratio = ds.relative_refinement(0, ds.max_level)

    # # Get the max possible number of pixels for the full domain
    # dims = ds.domain_dimensions * ref_ratio
    # pixels = (dims[0], dims[1])

    # Set FRB resolution
    slc.render()

    # Now get the frb directly from the slice plot
    # This gives 2D array the shape as the dims
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

    # fig, ax = plt.subplots(figsize=(8, 8))

    # # Needed for smoothed_covering_grid to use fine level data
    # ds.force_periodicity()

    # # Extract data coarsest to finest level and plot
    # r_min = ds.domain_left_edge[0].value
    # r_max = ds.domain_right_edge[0].value
    # th_min = ds.domain_left_edge[1].value
    # th_max = ds.domain_right_edge[1].value

    # max_level = ds.max_level
    # # First Do the Plotting
    # for level in range(max_level + 1):
    #     grids = [g for g in ds.index.grids if g.Level == level]
    #     if not grids:
    #         continue

    #     # Get left grid edge to get the covering grid
    #     # Make sure that we're within the actual domain and exclude ghost cells
    #     r_left = max(min(g.LeftEdge[0].value for g in grids), r_min)
    #     th_left = max(min(g.LeftEdge[1].value for g in grids), th_min)
    #     r_right = min(max(g.RightEdge[0].value for g in grids), r_max)
    #     th_right = min(max(g.RightEdge[1].value for g in grids), th_max)
    #     left_edge = ds.arr([r_left, th_left, ds.domain_left_edge[2].value], "code_length")
    #     right_edge = ds.arr([r_right, th_right, ds.domain_right_edge[2].value], "code_length")

    #     # Get covering grid corresponding to each individual levels
    #     # Compute the refinement ratio between the coarsest level and current level
    #     # ref_ratio = ds.relative_refinement(0, level)
    #     ref_ratio = ds.refine_by**level
    #     dims = ds.domain_dimensions * ref_ratio
    #     dims[2] = 1
    #     print(ds.refinement_factors)
    #     cg = ds.covering_grid(level=level, left_edge=left_edge, dims=dims)

    #     # Get data and plot
    #     var   = cg["boxlib", f][:, :, 0].to_ndarray()
    #     r     = cg["index", "r"][:, :, 0].to("cm").to_ndarray()
    #     theta = cg["index", "theta"][:, :, 0].to_ndarray()

    #     # Convert to cylindrical R-Z coordinate -- This is cell-centered coordinate
    #     R = r * np.sin(theta)
    #     Z = r * np.cos(theta)

    #     # yt can generate data outside of the simulation domain -- ghost cells?
    #     # Make sure that we're within the actual domain for the field we're plotting
    #     mask = (r >= r_min) & (r <= r_max) & (theta >= th_min) & (theta <= th_max)
    #     masked_var = np.where(mask, var, np.nan)

    #     pcm = ax.pcolormesh(R, Z, masked_var, cmap=cmap,
    #                         norm=Normalize(vmin=var.min(), vmax=var.max()), shading="auto")

    #     if level == 0:
    #         # Compute the edge values so that we get proper x-y lim
    #         dr = r[1, 0] - r[0, 0]
    #         dtheta = theta[0, 1] - theta[0, 0]

    #         Rr = (r + 0.5*dr) * np.sin(theta + 0.5*dtheta)
    #         Rmax = max(Rr[mask])
    #         Rl = (r - 0.5*dr) * np.sin(theta + 0.5*dtheta)
    #         Rmin = min(Rl[mask])

    #         Zr = (r + 0.5*dr) * np.cos(theta + 0.5*dtheta)
    #         Zmax = max(Zr[mask])
    #         Zl = (r - 0.5*dr) * np.cos(theta - 0.5*dtheta)
    #         Zmin = min(Zl[mask])

    #         ax.set_xlim(Rmin, Rmax)
    #         ax.set_ylim(Zmin, Zmax)

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
