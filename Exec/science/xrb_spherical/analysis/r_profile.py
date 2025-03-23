#!/usr/bin/env python3

# Spherical R profile at different theta

import os
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np
from functools import reduce
import itertools

import matplotlib.ticker as ptick
from yt.frontends.boxlib.api import CastroDataset
from yt.units import cm


plotfile = sys.argv[1]
ds = CastroDataset(plotfile)

rmin = ds.domain_left_edge[0]
rmax = rmin + 5000.0*cm
#rmax = ds.domain_right_edge[0]
print(ds.domain_left_edge[1])
fig, _ax = plt.subplots(2,2)

axes = list(itertools.chain(*_ax))

fig.set_size_inches(7.0, 8.0)

fields = ["Temp", "density", "x_velocity", "y_velocity"]
nice_names = [r"$T$ (K)", r"$\rho$ (g/${cm}^3$)", r"$u$ (cm/s)", r"$v$ (cm/s)"]

# 4 rays at different theta values
thetal = ds.domain_left_edge[1]
thetar = ds.domain_right_edge[1]
thetas = [thetal, 0.25*thetar, 0.5*thetar, 0.75*thetar]

for i, f in enumerate(fields):

    for theta in thetas:
        # simply go from (rmin, theta) -> (rmax, theta). Doesn't need to convert to physical R-Z
        ray = ds.ray((rmin, theta, 0*cm), (rmax, theta, 0*cm))

        # sort by "t", which goes from 0 to 1 representing the spatial order.
        isrt = np.argsort(ray["t"])
        axes[i].plot(ray['r'][isrt], ray[f][isrt], label=r"$\theta$ = {:.4f}".format(float(theta)))

    axes[i].set_xlabel(r"$r$ (cm)")
    axes[i].set_ylabel(nice_names[i])
    axes[i].set_yscale("symlog")

    if i == 0:
        axes[0].legend(frameon=False, loc="lower left")

#fig.set_size_inches(10.0, 9.0)
plt.tight_layout()
plt.savefig("{}_profiles.png".format(os.path.basename(plotfile)))
