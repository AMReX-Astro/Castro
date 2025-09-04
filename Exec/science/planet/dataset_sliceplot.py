#!/usr/bin/env python

import yt
from yt.mods import *
import glob
import matplotlib
import os
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from yt.units import kpc

delta_i=10000
i=0
max_data=90000

every_nthfile=3
profilefilelist = sorted(glob.glob('plt*'), key=lambda name: int(name[3:]))
print profilefilelist


old_time=0
time_difference=0.000001
old_time=-time_difference
j=0
for i in profilefilelist:

  pf=load(i)

  time=float(pf.current_time)
  filename_time=str(int(time))
  filename=filename_time+'s_'+i
  print(filename)
  my_ray = pf.ortho_ray(1,(0,0))

  if time>=old_time+time_difference:

    slicevector=[0,1,0]

    ps= SlicePlot(pf,2,'density',slicevector)
    temp_max=1.0e-1
    temp_min=1.0e-21
    ps.annotate_timestamp()
    ps.set_zlim('density', temp_min, temp_max)
    ps.save('plot_density_'+filename)

    pd= SlicePlot(pf,2,'pressure',slicevector)
    temp_max=1.0e12
    temp_min=1.0e-23
    pd.annotate_timestamp()
    pd.set_zlim('pressure', temp_min, temp_max)
    pd.save('plot_pressure_'+filename)

    pg= SlicePlot(pf,2,'Temp',slicevector)
    temp_max=1.0e4
    temp_min=1.0e-2
    pg.annotate_timestamp()
    pg.set_zlim('Temp', temp_min, temp_max)
    pg.save('plot_Temp_'+filename)



    pe= SlicePlot(pf,2,'entropy',slicevector)
    temp_max=5.0e9
    temp_min=1.0e8
    pe.annotate_timestamp()
    pe.set_zlim('entropy', temp_min, temp_max)
    pe.save('plot_entropy_'+filename)

    ph= SlicePlot(pf,2,'rad',slicevector)
    temp_max=1.0e-1
    temp_min=1.0e-15
    ph.annotate_timestamp()
    ph.set_zlim('rad', temp_min, temp_max)
    ph.save('plot_radiation_'+filename)






