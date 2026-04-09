#!/usr/bin/env python

import os
import sys
import yt

plotfile = sys.argv[1]

# slice plot of temperature
ds = yt.load(plotfile)
sp = yt.SlicePlot(ds, "theta", "Temp")
sp.set_zlim("Temp", 1.e-3, 100)

# now do a contour of density
sp.annotate_contour("density", ncont=2, clim=(0.01, 1.0),
                    plot_args={"colors": "red", "linewidths": 2})

sp.save(f"{os.path.basename(plotfile)}_slice.png")

