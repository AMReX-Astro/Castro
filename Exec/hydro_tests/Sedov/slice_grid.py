#!/usr/bin/env python

import sys
import yt
from yt.frontends.boxlib.api import CastroDataset

"""
Give a pressure slice plot for 2d/3d sedov plot file
"""

fname = sys.argv[1]
if len(sys.argv) == 3:
    outname = sys.argv[2]
else:
    outname = "sedov_slice.png"

ds = CastroDataset(fname)

slice_dirs = {"cartesian":"z",
              "cylindrical":"theta",
              "spherical":"phi"}
slice_dir = slice_dirs[ds.geometry]

f = "pressure" #"rho_E"
slc = yt.SlicePlot(ds, slice_dir, f)
if ds.geometry == "cylindrical" or ds.geometry == "cartesian":
    slc.annotate_grids()
slc.set_figure_size(12)
slc.set_buff_size(1600)
slc.set_font_size(24)
slc.set_cmap(f, "plasma_r")
slc.set_log(f, False)
slc.set_axes_unit("cm")

slc.save(outname)
