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
ax.set_xlabel(r'$y\,[\mathrm{km}]$')
ax.set_ylabel(r'$v_y \quad\times 10^6\,[\mathrm{cm}\, \mathrm{s}^{-1}]$')

times = [300, 931, 1524,  1540.5]

field = 'y_velocity'
weight = ('index', 'ones')
res = 960

for time in times:

    id_i = int(time/0.5)

    ds_m = CastroDataset(dir_paths[id_i-1])
    ds_c = CastroDataset(dir_paths[id_i])
    ds_p = CastroDataset(dir_paths[id_i+1])

    sp_m = ds_m.all_data()
    sp_c = ds_c.all_data()
    sp_p = ds_p.all_data()

    profile_m = yt.create_profile(
        sp_m,
        [("index", "y")], [('boxlib', field)],
        weight_field=weight,
        accumulation=False,
        n_bins=res,
        )

    profile_c = yt.create_profile(
        sp_c,
        [("index", "y")], [('boxlib', field)],
        weight_field=weight,
        accumulation=False,
        n_bins=res,
        )

    profile_p = yt.create_profile(
        sp_p,
        [("index", "y")], [('boxlib', field)],
        weight_field=weight,
        accumulation=False,
        n_bins=res,
        )

    x_m = profile_m.x.to("km")
    x_c = profile_c.x.to("km")
    x_p = profile_p.x.to("km")

    val_m = profile_m[('boxlib', field)]
    val_c = profile_c[('boxlib', field)]
    val_p = profile_p[('boxlib', field)]

    x = (x_m + x_c + x_p)/3
    val = (val_m + val_c + val_p)/3

    val /= 1.0e6

    time = (ds_m.current_time.value + ds_c.current_time.value + ds_p.current_time.value)/3
    ax.plot(x, val, label=f"${time:.2f}\,[s]$")

ax.legend(loc="upper left")
plt.axvline(x = 225.8, color='k', linestyle='--')
ax.set_ylim(-6, 16)

fig.set_size_inches(8, 6)
ax.tick_params(labelsize=21)
plt.tight_layout()

fig.savefig('lateral_velocity_04_ext.pdf', bbox_inches="tight")