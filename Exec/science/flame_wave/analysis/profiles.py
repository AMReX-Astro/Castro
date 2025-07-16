#!/usr/bin/env python3

# vertical profiles at fixed r

import os
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np
from functools import reduce
import itertools

import matplotlib.ticker as ptick


# assume that our data is in CGS
from yt.units import cm, amu


# Define Abar
def _Abar(field, data):
    """ Mean atomic mass. """

    sum = None

    for i, f in enumerate(mfrac_fields):

       mfracs = data[f]
       A = atomic_masses[i]

       if sum is None: sum = mfracs / A
       else: sum += mfracs / A

    return 1 / sum * amu



plotfile = sys.argv[1]
ds = yt.load(plotfile)

# Get mass fraction fields and atomic masses, assuming an identical field list for each dataset.
xfilt = lambda f: f.startswith("X(") and f.endswith(")")
fields = map(lambda f: f[1], ds.field_list)
mfrac_fields = np.array(list(filter(xfilt, fields)))

def to_atomic_mass(mfrac_field):
    """Conversion function from field names to atomic masses."""

    numeric = filter(lambda char: char.isnumeric(), mfrac_field)
    return int(reduce(lambda a, b: a + b, numeric))

atomic_masses = np.array(list(map(to_atomic_mass, mfrac_fields)))
indices = np.argsort(atomic_masses)

mfrac_fields = mfrac_fields[indices]
atomic_masses = atomic_masses[indices]


# slice plot of temperature
ds.add_field(("gas", "Abar"), function=_Abar, display_name=r"$\bar{A}$", units="amu", sampling_type="cell")

xctr = 0.5*(ds.domain_left_edge[0] + ds.domain_right_edge[0])
L_x = ds.domain_right_edge[0] - ds.domain_left_edge[0]

ymin = 0.0*cm
ymax = 1.5e4*cm

yctr = 0.5*(ymin + ymax)
L_y = ymax - ymin


fig, _ax = plt.subplots(2,2)

axes = list(itertools.chain(*_ax))


fig.set_size_inches(7.0, 8.0)


fields = ["Temp", "Abar", "enuc", "z_velocity"] #, "density"]
nice_names = [r"$T$ (K)", r"$\bar{A}$ (amu)", r"$\dot{e}_\mathrm{nuc}$ (erg/g/s)", r"$w$ (cm/s)"]
log = [1, 0, 1, 0]

# we'll do 3 rays, at r = 0, 1.e4, 2.e4
rvals = [0, 2.e4, 4.e4, 6.e4]

for i, f in enumerate(fields):

    for r in rvals:
        ray = ds.ray((r, 0, 0), (r, 1.5e4, 0))
        isrt = np.argsort(ray["t"])
        axes[i].plot(ray['z'][isrt], ray[f][isrt], label=f"r = {r} cm")

    axes[i].set_xlabel(r"$z$ (cm)")
    axes[i].set_ylabel(nice_names[i])
    if log[i] == 1:
        axes[i].set_yscale("log")
    else:
        # make the offset text look nice and latex-y
        # https://stackoverflow.com/questions/31517156/adjust-exponent-text-after-setting-scientific-limits-on-matplotlib-axis
        axes[i].yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
    if i == 0:
        axes[0].legend(frameon=False)


#fig.set_size_inches(10.0, 9.0)
plt.tight_layout()
plt.savefig(f"{os.path.basename(plotfile)}_profiles.png")

