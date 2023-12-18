#!/usr/bin/env python3
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from numpy.lib import recfunctions as rfn

from diag_parser import deduplicate, read_diag_file

data = deduplicate(read_diag_file(sys.argv[1]))

mass_columns = [name for name in rfn.get_names(data.dtype) if name.startswith("Mass ")]

# compute the total mass in the domain
total_mass = rfn.apply_along_fields(np.sum, data[mass_columns])

# cycle through colors then line styles
plt.rcParams["axes.prop_cycle"] = (
    cycler(linestyle=["-", "--"]) * mpl.rcParamsDefault["axes.prop_cycle"]
)

plt.figure(figsize=(19.20, 10.80))
for col in mass_columns:
    plt.plot(data["TIME"], data[col] / total_mass, label=col.removeprefix("Mass "))
plt.xlabel("Time (s)")
plt.ylabel("Mass fraction")
# plt.xscale("log")
plt.yscale("log")

# put the legend to the right of the plot, not overlapping
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))

plt.tight_layout()
plt.savefig("massfrac_plot.png")
