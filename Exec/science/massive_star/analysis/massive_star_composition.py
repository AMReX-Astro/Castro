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
from yt.visualization.image_writer import multi_image_composite


def _iron_group(field, data):
    return (
        data["boxlib", "X(Cr48)"] +
        data["boxlib", "X(Fe52)"] +
        data["boxlib", "X(Fe54)"] +
        data["boxlib", "X(Ni56)"] + 1.e-10)

def _si_group(field, data):
    return (
        data["boxlib", "X(Si28)"] +
        data["boxlib", "X(S32)"] +
        data["boxlib", "X(Ar36)"] +
        data["boxlib", "X(Ca40)"] +
        data["boxlib", "X(Ti44)"] + 1.e-10)

def _light_nuclei(field, data):
    return (
        data["boxlib", "X(H1)"] +
        data["boxlib", "X(He3)"] +
        data["boxlib", "X(He4)"] +
        data["boxlib", "X(C12)"] +
        data["boxlib", "X(N14)"] +
        data["boxlib", "X(O16)"] +
        data["boxlib", "X(Ne20)"] +
        data["boxlib", "X(Mg24)"] +
        data["boxlib", "X(n)"] +
        data["boxlib", "X(p)"] + 1.e-10)

yt.add_field(
    name=("gas", "iron_group"),
    function=_iron_group,
    sampling_type="local",
    units="",
)

yt.add_field(
    name=("gas", "si_group"),
    function=_si_group,
    sampling_type="local",
    units="",
)

yt.add_field(
    name=("gas", "light_nuclei"),
    function=_light_nuclei,
    sampling_type="local",
    units="",
)

plotfile = sys.argv[1]
ds = CastroDataset(plotfile)

xmin = ds.domain_left_edge[0]
xmax = ds.domain_right_edge[0]

ymin = ds.domain_left_edge[1]
ymax = ds.domain_right_edge[1]

xctr = 0.0 * xmin
L_x = xmax - xmin

yctr = 0.5 * (ymin + ymax)
L_y = ymax - ymin


fig = plt.figure()
fig.set_size_inches(12.0, 9.0)

width_frac = 0.1

center=[xmin + 0.25*width_frac*L_x, yctr, 0.0*cm]
width=[0.5*width_frac*L_x, width_frac*L_y]

slc = yt.SlicePlot(ds, "theta",
                   fields=[("gas", "iron_group"), ("gas", "si_group"), ("gas", "light_nuclei")],
                   center=[xmin + 0.25*width_frac*L_x, yctr, 0.0*cm],
                   width=[0.5*width_frac*L_x, width_frac*L_y, 0.0*cm], fontsize="12")

res = (1024, 512)
frb = slc.data_source.to_frb(width[0], res, height=width[1]) #width, res)#, center=center)

multi_image_composite("multi_channel1.png",
                      np.transpose(np.log10(frb["iron_group"])),
                      np.transpose(np.log10(frb["si_group"])),
                      np.transpose(np.log10(frb["light_nuclei"])))

