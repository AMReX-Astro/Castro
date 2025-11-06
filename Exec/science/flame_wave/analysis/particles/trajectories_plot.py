"""
Creates a trajectory plot for each unique particle from a Castro simulation with tracer particles enabled.
The start and end points of each trajectory are also marked.

Set the xlim and ylim as needed to match the region where particles are distributed (likely a subset of problo and probhi)
"""

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax
from matplotlib import pylab
import matplotlib.ticker as ticker
from matplotlib import colors

files = sorted(glob.glob("./*.txt"))
for f in files:
    print(f)
    # read file, only getting x and y positions for now and storing them as separate arrays
    arr = np.loadtxt(f, delimiter = " ", usecols = (0,1), unpack=False)
    if arr.ndim == 1: #skip file if the particle immediately left the domain
        print("only has one point, skip")
        continue
    else:
        Xs = arr[:,0]
        Ys = arr[:,1]
    fig = plt.figure(figsize=(15, 6))
    plt.plot(Xs, Ys,linewidth=0.5,zorder=0,label="trajectory")
    plt.plot(Xs[0],Ys[0],'k*',zorder=5,label="start")
    plt.plot(Xs[-1],Ys[-1],'ro',zorder=10,label="end")
    plt.xlabel('r (cm)')
    plt.ylabel('z (cm)')
    plt.xlim([0,1.8e5])
    plt.ylim([2e3,1e4])
    plt.legend()
    fig.savefig("plots/"+f.rsplit("_", 1)[0]+'_trajectory.png')
    plt.clf()
    plt.close(fig)