#!/usr/bin/env python3

import os
import sys
import yt

# assume that our data is in CGS
from yt.units import cm

plotfile = sys.argv[1]

# slice plot of temperature
ds = yt.load(plotfile)

xctr = 0.5*(ds.domain_left_edge[0] + ds.domain_right_edge[0])
L_x = ds.domain_right_edge[0] - ds.domain_left_edge[0]

ymin = 0.0*cm
ymax = 1.5e4*cm

yctr = 0.5*(ymin + ymax)
L_y = ymax - ymin


fields = ["Temp", "enuc", "density"]

for f in fields:
    sp = yt.SlicePlot(ds, "theta", f, center=[xctr, yctr, 0.0], width=[L_x, L_y, 0.0])
    if f == "Temp":
        sp.set_zlim(f, 1.e7, 1.e9)
    elif f == "enuc":
        sp.set_zlim(f, 1.e15, 1.e20)
    elif f == "density":
        sp.set_zlim(f, 1.e-3, 5.e8)

    if f != "density":
        # now do a contour of density
        sp.annotate_contour("density", ncont=2, clim=(1.e4, 2.e6),
                            plot_args={"colors": "red", "linewidths": 2})

    sp.annotate_text((0.05, 0.05), "{}".format(ds.current_time.in_cgs()),
                     coord_system="figure", text_args={"color": "black"})

    sp.save("{}_{}_slice.png".format(os.path.basename(plotfile), f))
