#!/usr/bin/env python3

import sys
import pandas as pd
from scipy.stats import linregress

import matplotlib as mpl
import matplotlib.pyplot as plt

def measure_and_plot(rad, t, stab_ind):
    
    slopes = dict()
    
    rad.plot()

    for col in rad:
        
        radarr = rad[col][stab_ind:]
        tarr = t["time"][stab_ind:]
        
        m, b, r, _, _ = linregress(tarr, radarr)
        print("{:>20}\t{:.2e}\t{:.3f}".format(col, m, r))
        
        if r > 0.8:
            
            plt.plot(t["time"].index.values, m * t["time"] + b, "k--")
            
        slopes["col"] = m
        
    return slopes

if __name__ == "__main__":
    
    try:
        
        rad = pd.read_csv(sys.argv[1], sep="\s+")
        t = pd.read_csv(sys.argv[2], sep="\s+")
        stab_ind = int(sys.argv[3])
        
    except:
        
        print("Usage: ./flame_speed.py radii.dat times.dat starting_index")
    
    plt.style.use("seaborn")
    measure_and_plot(rad, t, stab_ind)
    plt.show()
