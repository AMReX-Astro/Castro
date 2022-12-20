"""
Another plotting script to use after running trajectories_plot.py
If you've identified several trajectories of interest, this script plots all of them in one figure.
"""

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax
from matplotlib import pylab
import matplotlib.ticker as ticker
from matplotlib import colors

plt.rcParams["font.family"] = ["font.serif"]
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rc('axes', labelsize=25)
plt.rc('legend',fontsize=20)

#Several trajectories selected after looking through the individual trajectory plots
files = ['1_113_info.txt','2_185_info.txt','1_179_info.txt','1_130_info.txt','2_313_info.txt','2_382_info.txt']

fig = plt.figure(figsize=(15, 6))
colors = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown"]

#Marks the hottest region and cooler region (rough estimate, mainly to distinguish where particles are undergoing vortical motion)
plt.axhspan(2e3,5e3,0,1.5e5,color='darkorange',alpha=0.1,zorder=0)
plt.axhspan(5e3,1e4,0,1.5e5,color='indigo',alpha=0.1,zorder=1)

for i,f in enumerate(files):
    print(f)
    # read file, only getting x and y positions for now and storing them as separate arrays
    arr = np.loadtxt(f, delimiter = " ", usecols = (0,1), unpack=False)
    if arr.ndim == 1:
        print("only has one point, skip")
        continue
    else:
        Xs = arr[:,0]
        Ys = arr[:,1]
    plt.plot(Xs, Ys,linewidth=2,color=colors[i],zorder=2,label=f.rsplit("_", 1)[0])
    plt.plot(Xs[0],Ys[0],'k*',markersize=8,zorder=5)
    plt.plot(Xs[-1],Ys[-1],'ro',zorder=10)
    plt.xlabel('r (cm)')
    plt.ylabel('z (cm)')
    plt.xlim([0,1.8e5])
    plt.ylim([2e3,1e4])
    plt.legend(loc='lower right')

fig.savefig('./plots/stacked_trajectory.png',bbox_inches="tight")