import matplotlib
matplotlib.use('agg')

import os
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.axes_grid1 import ImageGrid

# assume that our data is in CGS
from yt.units import cm, amu
from yt.frontends.boxlib.api import CastroDataset

plotfile = sys.argv[1]
ds = CastroDataset(plotfile)

xmin = ds.domain_left_edge[0]
xmax = ds.domain_right_edge[0]
xctr = 0.5 * (xmin + xmax)
L_x = (2./3.) * (xmax - xmin)


ymin = ds.domain_left_edge[1]
ymax = ds.domain_right_edge[1]
yctr = 0.5*(ymin + ymax)


zmin = 0.0*cm
zmax = 1.e4*cm

zctr = 0.5*(zmin + zmax)
L_z = zmax - zmin


fig = plt.figure()
fig.set_size_inches(16.0, 9.0)


f = "abar"

sp = yt.SlicePlot(ds, "y", f, origin="native", center=[xctr, yctr, zctr],
                  width=[L_z, L_x]) #0.0*cm, L_z]) #, fontsize="12", origin="native")
sp.set_buff_size((4800,4800))
sp.swap_axes()

sp.set_zlim(f, 4, 5)
sp.set_log(f, False)
sp.set_cmap(f, "plasma_r")
sp.set_axes_unit("cm")
sp.set_figure_size(20)
sp.save(f"{os.path.basename(plotfile)}_slice.pdf")
