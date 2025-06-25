import yt
import matplotlib.pyplot as plt
import numpy as np
from yt.frontends.boxlib.api import CastroDataset
import os
import re

# yt.enable_parallelism()

# set matplotlib formatting to be the same as yt
plt.rcParams['font.size'] = 14
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'stix'

curr_dir = os.getcwd()
plt_dir = [f.name for f in os.scandir(curr_dir) if (f.is_dir() and re.search("^plt[0-9]{5,8}$",f.name))]
plt_dir.sort(key=lambda x: int(x[3:]))

out_dir = curr_dir+"/slices/"
ext_dir = [f.name[13:-4] for f in os.scandir(out_dir) if (f.is_file() and re.search("plt[0-9]{5,8}.png", f.name))]
# ext_dir.sort(key=lambda x: int(x[3:]))

remaining = [x for x in plt_dir if x not in ext_dir]

ts_04 = yt.DatasetSeries(remaining)
field = 'magvel'
text_color = "white"


for ds in ts_04.piter():

    sp = yt.SlicePlot(ds, 'z', ('boxlib',field), buff_size=(7680,3840),origin="domain", fontsize="14")

    if field == "magvel":
        sp.set_cmap(field, "inferno")
        sp.set_zlim(('boxlib', field), 1.0e5, 1.0e8)
    elif field == "enuc":
        sp.set_log(field, True, linthresh=1.e11)
        sp.set_zlim(field, 0., 1.e19)
        sp.set_cmap(field, "plasma")

    sp.annotate_text((0.05, 0.05), f"time = {float(ds.current_time):8.3f} s", coord_system="axis", text_args={"color": text_color, "fontsize": "14"})
    sp.set_axes_unit("km")
    outdir1 = out_dir + f"drive_{field}_{ds.basename}.png"
    sp.save(outdir1)
    ds.index.clear_all_data()
