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

fig, _ax = plt.subplots(2,2)

axes = list(itertools.chain(*_ax))

fig.set_size_inches(7.0, 8.0)

fields = ["Temp", "density", "x_velocity", "y_velocity"]
nice_names = [r"$T$ (K)", r"$\rho$ (g/${cm}^3$)", r"$u$ (cm/s)", r"$v$ (cm/s)"]

# 4 rays at different theta values
thetas = [0, 0.125*np.pi, 0.25*np.pi, 0.5*np.pi]

for i, f in enumerate(fields):

    for theta in thetas:
        ray = ds.ray((rmin*np.sin(theta), rmin*np.cos(theta), 0*cm),
                     (rmax*np.sin(theta), rmax*np.cos(theta), 0*cm))
        isrt = np.argsort(ray["t"])
        axes[i].plot(ray['r'][isrt], ray[f][isrt], label="theta = {}".format(theta))

    axes[i].set_xlabel(r"$r$ (cm)")
    axes[i].set_ylabel(nice_names[i])
    axes[i].set_yscale("symlog")

    if i == 0:
        axes[0].legend(frameon=False)

#fig.set_size_inches(10.0, 9.0)
plt.tight_layout()
plt.savefig("{}_profiles.png".format(os.path.basename(plotfile)))
