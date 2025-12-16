#!/usr/bin/env python

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
from yt.frontends.boxlib.api import CastroDataset

def _nrhs(field, data):
    return data["boxlib", "burn_weights_firsthalf"] + data["boxlib", "burn_weights_secondhalf"] + 1

yt.add_field(
    name = ("boxlib", "nrhs"),
    function = _nrhs,
    sampling_type = "local",
    units = None
)


plotfile = sys.argv[1]
ds = CastroDataset(plotfile)

xmin = ds.domain_left_edge[0]
xmax = 0.5*ds.domain_right_edge[0]
xctr = 0.5*(xmin + xmax)
L_x = xmax - xmin

ymin = ds.domain_left_edge[1]
ymax = ds.domain_right_edge[1]
yctr = 0.5*(ymin + ymax)
L_y = ymax - ymin
ymin = yctr - 0.25*L_y
ymax = yctr + 0.25*L_y
L_y = ymax - ymin

fig = plt.figure()

f = ("boxlib", "nrhs")

sp = yt.SlicePlot(ds, "theta", f, center=[xctr, yctr, 0.0*cm], width=[L_x, L_y, 0.0*cm], fontsize="12")
sp.set_buff_size((2400,2400))
sp.set_zlim(f, 1.e1, 1.e4)
sp.set_log(f, True) #False) #True, linthresh=1)
sp.save("strang_nrhs.png")
