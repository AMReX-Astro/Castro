#!/usr/bin/env python

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


plotfile = "subch_plt00000"

fig = plt.figure()

ds = CastroDataset(plotfile)

xmin = 0 * cm
xmax = 1.e8 * cm

xctr = 0.5 * (xmin + xmax)
L_x = xmax - xmin

ymin = 5.42e9 * cm
ymax = 5.58e9 * cm
yctr = 0.5 * (ymin + ymax)
L_y = ymax - ymin

field = "Temp"

sp = yt.SlicePlot(ds, "theta", field, center=[xctr, yctr, 0.0*cm], width=[L_x, L_y, 0.0*cm], fontsize="14")
sp.set_buff_size((2400,2400))

sp.set_zlim(field, 5.e7, 4e9)
sp.set_cmap(field, "magma")

sp.annotate_contour(("gas", "density"), take_log=True, ncont=3, clim=(1.e4, 1.e6), plot_args={"colors": "white", "linestyles": ":"})

sp.set_axes_unit("km")

sp.save(f"subch_{field}_zoom.pdf")
