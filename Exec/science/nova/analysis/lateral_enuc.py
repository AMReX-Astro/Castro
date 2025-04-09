import yt
import os
import re
import numpy as np
from yt.utilities.parallel_tools.parallel_analysis_interface import communication_system
import matplotlib.pyplot as plt
from yt.frontends.boxlib.api import CastroDataset

plt.rcParams['font.size'] = 21
plt.rcParams['font.family'] = 'serif'
plt.rcParams["axes.labelsize"] = 21
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['savefig.bbox']='tight'

curr_dir = os.getcwd()

plt_dir = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt+[0-9]{5,8}$",f.name))]
plt_dir.sort(key=lambda x: int(x[3:]))

dir_paths = []
for dir in plt_dir:
    dir_paths.append(os.path.join(curr_dir, dir))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r'$y\,[km]$')
ax.set_ylabel(r'$\dot{\mathrm{e}}_\mathrm{nuc}\,[\mathrm{erg}/\mathrm{g}/\mathrm{s}]$')

times = [300, 931, 1524, 1532]

for time in times:

    if time == 1532:
        ds =  ds = CastroDataset(dir_paths[-2])
    else:
        id = int(time/0.5)
        ds = CastroDataset(dir_paths[id])

    sp = ds.all_data()
    res = 960
    profile = yt.create_profile(
        sp,
        [("boxlib", "y")], [("boxlib", "enuc")],
        weight_field=('gas', 'mass'),
        accumulation=False,
        n_bins=res,
        )

    y = profile.x.to("km")
    val = profile[("boxlib", "enuc")]
    ax.plot(y, val, label="${:.2f}\,[s]$".format(ds.current_time.value))


plt.axvline(x = 225.8, color='k', linestyle='--')
#The CEI is located at 225.8km from the bottom boundary

ax.legend()
ax.set_yscale('log')
ax.set_ylim(1.0e11, 1.0e18)

fig.set_size_inches(8,6)
ax.tick_params(labelsize=21)
plt.tight_layout()

fig.savefig('lateral_enuc_04_ext.pdf', bbox_inches="tight")
