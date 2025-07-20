#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import pandas as pd

# Set some fontsize
SMALL_SIZE = 18
MEDIUM_SIZE = 20
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('xtick.major', size=7, width=2)
plt.rc('xtick.minor', size=5, width=1)
plt.rc('ytick.major', size=7, width=2)
plt.rc('ytick.minor', size=5, width=1)

parser = argparse.ArgumentParser(description='''This script uses output
from front_tracker.py to plot the time evolution of flame front θ.
''')

parser.add_argument('tracking_fname', type=str,
                    help='cvs file generated from front_tracker.py to track flame front position')

args = parser.parse_args()

# Read in data
# data has columns: fname, time, front_theta, theta_max_avg, max_avg, theta_max, max_val.
tracking_data = pd.read_csv(args.tracking_fname)

# Get time and theta
times = tracking_data['time[ms]']
front_thetas = tracking_data['front_theta']

# Do plotting
fig, ax = plt.subplots()
ax.plot(times, front_thetas)
ax.set_xlabel("time [ms]")
ax.set_ylabel("θ")
fig.tight_layout()
fig.set_size_inches(8, 8)
fig.savefig("flame_position.png", bbox_inches="tight")
