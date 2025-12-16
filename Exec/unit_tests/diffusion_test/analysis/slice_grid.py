#!/usr/bin/env python

import sys
import yt
from yt.frontends.boxlib.api import CastroDataset

"""
Give a temperature slice plot for 2d plot file
"""

fname = sys.argv[1]
ds = CastroDataset(fname)

slice_dirs = {"cylindrical":"theta",
             "spherical":"phi"}
slice_dir = slice_dirs[ds.geometry]

slc = yt.SlicePlot(ds, slice_dir, "Temp")
if ds.geometry == "cylindrical":
    slc.annotate_grids()
slc.set_figure_size(12)
slc.set_buff_size(1600)
slc.set_font_size(24)
slc.set_cmap("Temp", "plasma_r")
slc.set_log("Temp", False)
slc.set_axes_unit("cm")

slc.save("diffusion_temp.png")
