#!/usr/bin/env python

import yt
from yt.mods import *
import glob
import matplotlib
import os
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
np.set_printoptions(threshold=np.inf)

delta_i=10000
i=0
max_data=90000
j=0
every_nthfile=3
profilefilelist = sorted(glob.glob('plt*'), key=lambda name: int(name[3:]))
number_plt=0
for i in profilefilelist:
  number_plt=number_plt+1


average_T=0
average_P=0
Dimension=2


old_time=0
time_difference_x_den=1000.0
time_difference_x_T=1000.0
time_difference_x_P=1000.0
time_difference_x_entropy=1000.0
time_difference_x_rad=1000.0
time_difference_T_P=1000.0
time_difference_x_opacity=1000.0
time_difference_x_Mach=1000.0
j=0
old_time=-time_difference_x_den
max_time=10000.0


for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if time>=old_time+time_difference_x_den and time<max_time :
    old_time=time
    print('Adding lines at t='+filename_time+'s for $\rho$-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')

    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5

    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5



    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])

    my_ray = ds.ortho_ray(1,(0,0))
    plt.semilogy(x_coord,dens,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$\rho\/[\mathrm{g}/\mathrm{cm}^{3}]$')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
  #  plt.xscale('linear')
  #  plt.yscale('linear')
    plt.title('Density vs height')
    plt.savefig("figure_x_den.png")
  j=j+1

print("density plot made")
plt.close()



j=0
old_time=-time_difference_T_P
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))


  kappa=0.18*1e-9
  SB_constant=5.6704*1e-5


  if Dimension==1:
    my_ray = ds.ortho_ray(0,(0,0))
    srt=np.argsort(my_ray['x'])
    x_coord=np.array(my_ray['x'][srt])/1.e5
  elif Dimension==2:
    my_ray = ds.ortho_ray(1,(0,0))
    srt=np.argsort(my_ray['y'])
    x_coord=np.array(my_ray['y'][srt])/1.e5
  temp=np.array(my_ray['Temp'][srt])
  ent=np.array(my_ray['entropy'][srt])
  dens=np.array(my_ray['density'][srt])
  press=np.array(my_ray['pressure'][srt])/1.e6
  average_T=(average_T*float(j)+temp)/float(j+1)
  average_P=(average_P*float(j)+press)/float(j+1)


  if j==0 :
    temp0=np.array(my_ray['Temp'][srt])
    ent0=np.array(my_ray['entropy'][srt])
    pressure0=np.array(my_ray['pressure'][srt])
    density0=np.array(my_ray['density'][srt])
    location1=np.amin(pressure0)/1e10
    location2=np.amax(pressure0)

  if time>=old_time+time_difference_T_P and time<max_time :
    old_time=time
    plt.plot(average_P,average_T,label='t='+filename_time+'s')
    plt.legend()
    print('Adding lines at t='+filename_time+'s')
    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$T\/[\mathrm{K}]$')
    plt.xlim(0.001,2000)
    plt.ylim(1000,2000)
    plt.xscale('log')
    #plt.yscale('log')
    plt.title('<T> vs <P> (average up to t)')
    plt.savefig("figure_T_P_average.png")

  j=j+1
print("T-P plot average made")
plt.close()

j=0
old_time=-time_difference_T_P
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))



  kappa=0.18*1e-9
  SB_constant=5.6704*1e-5


  if Dimension==1:
    my_ray = ds.ortho_ray(0,(0,0))
    srt=np.argsort(my_ray['x'])
    x_coord=np.array(my_ray['x'][srt])/1.e5
  elif Dimension==2:
    my_ray = ds.ortho_ray(1,(0,0))
    srt=np.argsort(my_ray['y'])
    x_coord=np.array(my_ray['y'][srt])/1.e5
  temp=np.array(my_ray['Temp'][srt])
  ent=np.array(my_ray['entropy'][srt])
  dens=np.array(my_ray['density'][srt])
  press=np.array(my_ray['pressure'][srt])/1.e6
  average_T=(average_T*float(j)+temp)/float(j+1)
  average_P=(average_P*float(j)+press)/float(j+1)


  if j==0 :
    temp0=np.array(my_ray['Temp'][srt])
    ent0=np.array(my_ray['entropy'][srt])
    pressure0=np.array(my_ray['pressure'][srt])/1.e6
    density0=np.array(my_ray['density'][srt])
    location1=np.amin(pressure0)/1e10
    location2=np.amax(pressure0)/1e10

  diff_temp=(average_T-temp0)
  diff_pressure=average_P-pressure0

  if time>=old_time+time_difference_T_P and time<max_time :
    old_time=time
    plt.plot(average_P,diff_temp,label='t='+filename_time+'s')
    plt.legend()
    print('Adding lines at t='+filename_time+'s')
    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$\Delta T\/[\mathrm{K}]$')
    plt.xlim(0.001,2000)
    plt.ylim(-50,50)
    plt.figtext(0.5,0.3,r'$\Delta t_{\mathrm{damp}}=7000$s', fontsize=20)
    plt.xscale('log')
    #plt.yscale('log')
    plt.title(r'<$\Delta$ T> vs <P> (average up to t)')
    plt.savefig("figure_T_P_diff_average.png")

  j=j+1
print("T-P difference average plot made")
plt.close()


j=0
old_time=-time_difference_x_P
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if time>=old_time+time_difference_x_P and time<max_time:
    old_time=time
    print('Adding lines at t='+filename_time+'s for P-height plot')
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5


    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6


    if j==0 :
      temp0=np.array(my_ray['Temp'][srt])
      ent0=np.array(my_ray['entropy'][srt])
      pressure0=np.array(my_ray['pressure'][srt])
      density0=np.array(my_ray['density'][srt])
      location1=np.amin(pressure0)/1e10
      location2=np.amax(pressure0)

    plt.plot(x_coord,press,label='t='+filename_time+'s')
    plt.legend()

    plt.ylabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
  # plt.xscale('linear')
    plt.title('P vs height')
    plt.yscale('log')

    plt.savefig("figure_x_P.png")


  j=j+1


print("P plot made")
plt.close()

j=0
old_time=-time_difference_x_T
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if time>=old_time+time_difference_x_T and time<max_time:
    old_time=time
    print('Adding lines at t='+filename_time+'s for T-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')

    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5


    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6

    plt.semilogy(x_coord,temp,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$T$ [K]')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
    plt.title('T vs height')
  #plt.legend(loc=3)
  #plt.yscale('linear')
    plt.savefig("figure_x_T.png")


  j=j+1

print("T plot made")
plt.close()

j=0
old_time=-time_difference_x_opacity
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if time>=old_time+time_difference_x_opacity and time<max_time:
    old_time=time
    print('Adding lines at t='+filename_time+r's for $\kappa$-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')

    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5


    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6
    opacity=np.array(my_ray['pressure'][srt])*0.18/1.e9
    plt.semilogy(x_coord,opacity,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$\kappa$ [$\mathrm{cm}^{2}\mathrm{g}^{-1}$]')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
    plt.title(r'Opacity $\kappa$ vs height')
    #plt.legend(loc=3)
    #plt.yscale('linear')
    plt.savefig("figure_x_opacity.png")


  j=j+1

print("opacity plot made")
plt.close()



j=0
old_time=-time_difference_x_Mach
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if time>=old_time+time_difference_x_Mach and time<max_time:
    old_time=time
    print('Adding lines at t='+filename_time+r's for $\mathca{M}$-height plot')
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5


    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6
    Mach_number=np.array(my_ray['MachNumber'][srt])

    if j==0 :
      temp0=np.array(my_ray['Temp'][srt])
      ent0=np.array(my_ray['entropy'][srt])
      pressure0=np.array(my_ray['pressure'][srt])
      density0=np.array(my_ray['density'][srt])
      location1=np.amin(pressure0)/1e10
      location2=np.amax(pressure0)

    plt.plot(press,Mach_number,label='t='+filename_time+'s')
    plt.legend()

    plt.ylabel(r'Mach number')
    plt.xlabel(r'$P\/[\mathrm{bar}]$')
    #plt.ylim(location1,location2)
    #  plt.xscale('linear')
    plt.xlim(0.001,2000)
    plt.ylim(0,5)
    plt.yscale('linear')
    plt.xscale('log')
    plt.title(r'Mach number \mathcal{M} vs height')
    plt.savefig("figure_P_Mach.png")


  j=j+1


print("Mach plot made")
plt.close()



j=0
old_time=-time_difference_x_opacity
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if time>=old_time+time_difference_x_opacity and time<max_time:
    old_time=time
    print('Adding lines at t='+filename_time+r's for $\kappa$-T plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')

    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5


    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6
    opacity=np.array(my_ray['pressure'][srt])*0.18/1.e9
    plt.semilogy(temp,opacity,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$\kappa$ [$\mathrm{cm}^{2}\mathrm{g}^{-1}$]')
    plt.xlabel(r'$T\/[\mathrm{K}]$')
    #plt.legend(loc=3)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(400,10000)
    plt.ylim(-100,1000)
    plt.title(r'Opacity $\kappa $vs T')
    plt.savefig("figure_T_opacity.png")


  j=j+1

print("opacity-T plot made")
plt.close()



j=0
old_time=-time_difference_x_entropy
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if time>=old_time+time_difference_x_entropy and time<max_time:
    old_time=time
    print('Adding lines at t='+filename_time+'s for entropy-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')

    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5


    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6

    plt.semilogy(x_coord,ent,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$\mathrm{Entropy}$')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
    #plt.legend(loc=3)
    #plt.yscale('linear')
    plt.title('Entropy $vs height')
    plt.savefig("figure_x_entropy.png")


  j=j+1

print("entropy plot made")
plt.close()


j=0
old_time=-time_difference_x_rad
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if time>=old_time+time_difference_x_rad and time<max_time:
    old_time=time
    print('Adding lines at t='+filename_time+'s for $U$-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')

    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5


    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6
    radiation=np.array(my_ray['rad'][srt])
    plt.semilogy(x_coord,radiation,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$U [{\rm erg} {\rm cm}^{-3}]$')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
    #plt.legend(loc=3)
    #plt.yscale('linear')
    plt.figtext(0.5,0.3,r'$\Delta t_{\mathrm{damp}}=7000$s', fontsize=20)
    plt.title('Radiation energy density $U$ vs height')
    plt.savefig("figure_x_radiation.png")


  j=j+1

print("radiation plot made")
plt.close()






j=0
old_time=-time_difference_T_P
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if time>=old_time+time_difference_T_P and time<max_time:
    old_time=time
    print('Adding lines at t='+filename_time+'s for $T-P$ plot')


    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5


    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6


    if j==0 :
      temp0=np.array(my_ray['Temp'][srt])
      ent0=np.array(my_ray['entropy'][srt])
      pressure0=np.array(my_ray['pressure'][srt])
      density0=np.array(my_ray['density'][srt])
      location1=np.amin(pressure0)/1e10
      location2=np.amax(pressure0)

    plt.plot(press,temp,label='t='+filename_time+'s')
    plt.legend()

    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$T\/[\mathrm{K}]$')
    plt.xlim(0.001,2000)
    plt.ylim(1000,2000)
    plt.xscale('log')
    #plt.yscale('log')
    plt.title('T vs P')
    plt.savefig("figure_T_P.png")


  j=j+1
print("T-P plot made")
plt.close()



j=0
old_time=-time_difference_T_P
for i in profilefilelist:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if time>=old_time+time_difference_T_P and time<max_time:
    old_time=time
    print('Adding lines at t='+filename_time+r's for $\Delta T - P$ plot')


    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5


    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6

    if j==0 :
      temp0=np.array(my_ray['Temp'][srt])
      ent0=np.array(my_ray['entropy'][srt])
      pressure0=np.array(my_ray['pressure'][srt])
      density0=np.array(my_ray['density'][srt])
      location1=np.amin(pressure0)/1e10
      location2=np.amax(pressure0)

    diff_temp=(temp-temp0)
    plt.plot(press,diff_temp,label='t='+filename_time+'s')
    plt.legend()

    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$\Delta T\/[\mathrm{K}]$')
    plt.xlim(0.001,2000)
    plt.ylim(-200,200)
    plt.figtext(0.5,0.3,r'$\Delta t_{\mathrm{damp}}=7000$s', fontsize=20)
    plt.title(r'$\Delta$ T vs P')
    plt.xscale('log')
    #plt.yscale('log')

    plt.savefig("figure_T_P_diff.png")
  j=j+1
print("T-P_diff plot made")
