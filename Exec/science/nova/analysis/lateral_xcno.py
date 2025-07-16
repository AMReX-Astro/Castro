import yt
import os
import re
import numpy as np
from yt.utilities.parallel_tools.parallel_analysis_interface import communication_system
import matplotlib.pyplot as plt
from yt.frontends.boxlib.api import CastroDataset
from yt.units import dimensions

plt.rcParams['font.size'] = 21
plt.rcParams['font.family'] = 'serif'
plt.rcParams["axes.labelsize"] = 21
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['savefig.bbox']='tight'

@yt.derived_field(
    name=("gas","Zcno"),
    sampling_type="local",
    units="auto",
    take_log=False,
    display_name=r"Z_{\mathrm{CNO}}",
    dimensions=dimensions.dimensionless,
)
def x_mixing(field, data):
    cno = (data['X(C12)'] + data['X(C13)'] + data['X(N13)']+ data['X(N14)'] + data['X(N15)']
            + data['X(O14)'] +  data['X(O15)'] + data['X(O16)'] + data['X(O17)'] )
    return cno


curr_dir = os.getcwd()

plt_dir = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt+[0-9]{5,8}$",f.name))]
plt_dir.sort(key=lambda x: int(x[3:]))

dir_paths = []
for dir in plt_dir:
    dir_paths.append(os.path.join(curr_dir, dir))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r'$y\,[km]$')
ax.set_ylabel(r'$\mathrm{Z}_{\mathrm{CNO}}$')

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
        [("boxlib","y")], [("gas", "Zcno")],
        weight_field=('gas', 'mass'),
        accumulation=False,
        n_bins=res,
        )

    y = profile.x.to("km")
    val = profile[("gas", "Zcno")]
    ax.plot(y, val, label=f"${ds.current_time.value:.2f}\,[s]$")

ax.legend()
plt.axvline(x = 225.8, color='k', linestyle='--')

fig.set_size_inches(8, 6)
ax.tick_params(labelsize=21)
plt.tight_layout()

fig.savefig('lateral_xcno_04_ext.pdf', bbox_inches="tight")
