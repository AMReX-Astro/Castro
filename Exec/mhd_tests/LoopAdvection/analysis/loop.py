# Visualize the magnetic pressure for the loop-advection test

import sys

import yt
from yt import derived_field

@derived_field(name="magp", sampling_type="cell")
def _magp(field, data):
    return 0.5*(data["B_x"] * data["B_x"] +\
           data["B_y"] * data["B_y"] + data["B_z"] * data["B_z"])

plotfile = sys.argv[1]

ds = yt.load(plotfile)
ds.add_field(("boxlib", "magp"),
             function=_magp, sampling_type="cell")

field = ("boxlib", "magp")

sp = yt.SlicePlot(ds, "z", field)
sp.set_log(field, False)
sp.set_zlim(field, 0, 5e-7) #for no log scale
sp.set_cmap(field, cmap='RdBu')
sp.save("loop.png")
