#!/usr/bin/env python


import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt
import sys
import numpy as np
from scipy.optimize import brentq

yt.funcs.mylog.setLevel(50)

class Profile:
    """read a plotfile using yt and store the 1d profile for T and enuc"""

    def __init__(self, plotfile):
        ds = yt.load(plotfile, hint="castro")

        time = float(ds.current_time)
        ad = ds.all_data()

        # Sort the ray values by 'x' so there are no discontinuities
        # in the line plot
        srt = np.argsort(ad['x'])
        x_coord = np.array(ad['x'][srt])
        temp = np.array(ad['Temp'][srt])
        enuc = np.array(ad['enuc'][srt])

        self.time = time
        self.x = x_coord
        self.T = temp
        self.enuc = enuc

    def get_dx(self):
        """ return the grid dx -- assumes uniform spacing"""
        return self.x[1] - self.x[0]

    def find_x_for_T(self, T_0=1.e9):
        """ given a profile x(T), find the x_0 that corresponds to T_0 """

        # our strategy here assumes that the hot ash is in the early
        # part of the profile.  We then fit a quartic conservative interpolant
        # through the 4 zones surrounding the desired temperature and root find to
        # find the x that corresponds to the T

        # find the index of the first point where T drops below T_0
        idx = np.where(self.T < T_0)[0][0]
        
        # the 4 points for our quartic will then be: idx-2, idx-1, idx, and idx+1
        T1 = self.T[idx-2]
        x1 = self.x[idx-2]

        T2 = self.T[idx-1]
        x2 = self.x[idx-1]

        T3 = self.T[idx]
        x3 = self.x[idx]

        T4 = self.T[idx+1]
        x4 = self.x[idx+1]

        dx = x2 - x1

        s1 = (T4 - 3*T3 + 3*T2 - T1)/(6*dx**3)
        s2 = (T3 - 2*T2 + T1)/(2*dx**2)
        s3 = (-5*T4 + 27*T3 - 15*T2 - 7*T1)/(24*dx)
        s4 = (-T3 + 26*T2 - T1)/24.0

        x_flame = brentq(lambda x: s1*(x-x2)**3 + s2*(x-x2)**2 + s3*(x-x2) + s4 - T_0,
                         x1, x4)

        return x_flame

    def find_flame_width(self):
        """ given a profile T(x), find the flame width """

        gradT = np.max(np.abs( (self.T[1:] - self.T[:-1])/(self.x[1:] - self.x[:-1]) ))
        dT = self.T.max() - self.T.min()

        return dT/gradT

    def find_flame_width_reaction_zone(self, fraction_low=0.2, fraction_high=0.8):
        """ 
    Calculate flame width as ΔT/max(gradT) within the reaction zone.
    The reaction zone is defined dynamically as the region between
    fraction_low and fraction_high of the total temperature range.
    
    For example: fraction_low=0.2, fraction_high=0.8 means the reaction
    zone spans from 20% to 80% of the temperature range.
    """
    
        T_min = self.T.min()
        T_max = self.T.max()
        T_range = T_max - T_min
    
    # Define reaction zone boundaries based on fractions of temperature range
        T_low = T_min + fraction_low * T_range
        T_high = T_min + fraction_high * T_range
    
    # Find indices where temperature is between T_low and T_high
        in_reaction_zone = (self.T >= T_low) & (self.T <= T_high)
    
        if not np.any(in_reaction_zone):
            raise ValueError(f"No points found in reaction zone between {T_low:.2e} and {T_high:.2e} K")
    
    # Calculate temperature gradient at all points
        gradT = np.abs((self.T[1:] - self.T[:-1]) / (self.x[1:] - self.x[:-1]))
    
    # Find max gradient within the reaction zone (note: gradT is one element shorter)
        gradT_in_zone = gradT[in_reaction_zone[:-1]]
    
        if len(gradT_in_zone) == 0:
            raise ValueError("No gradient points in reaction zone")
    
        max_gradT = np.max(gradT_in_zone)
    
    # ΔT is the temperature range in the reaction zone
        dT = T_high - T_low
    
        return dT / max_gradT
  

if __name__ == "__main__":
    p = argparse.ArgumentParser()

    p.add_argument("--skip", type=int, default=1,
                   help="interval between plotfiles")
    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")
    p.add_argument("--end-plotfile", type=int, default=None,
               help="the last plotfile number to process(only number)")
    p.add_argument("--start-plotfile", type=int, default=None,
               help="the first plotfile number to process(only number)")
    args = p.parse_args()

    prefix = args.plotfiles[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.plotfiles], key=int)
    if args.start_plotfile is not None:
        plot_nums = [p for p in plot_nums if int(p) >= args.start_plotfile]
    if args.end_plotfile is not None:
        plot_nums = [p for p in plot_nums if int(p) <= args.end_plotfile]
    t = []
    v = []
    vs = []
    w = []
    w_rz = [] 
    pltfiles = []  # Add this
    ptimes = []    # Add this
    
    x1_old = None
    x2_old = None
    x3_old = None
    width_rz_old = None  # Add this
    for n in range(0, len(plot_nums), args.skip):

        pfile = f"{prefix}{plot_nums[n]}"
        print(f"Processing {pfile}...")  # Add this line
        try:
            p = Profile(pfile)
            x1 = p.find_x_for_T(T_0=1.5e9)
            x2 = p.find_x_for_T(T_0=1.e9)
            x3 = p.find_x_for_T(T_0=2.e9)
            width = p.find_flame_width()
            width_rz = p.find_flame_width_reaction_zone(fraction_low=0.2, fraction_high=0.8)
        except (ValueError, IndexError) as e:
            print(f"Error in {pfile}: {e}")
            print(f"Stopping analysis. Processed {len(t)} time intervals.")
            break

        if x1_old is not None:

            # difference with the previous file to find the det speed
            # note: the corresponding time is centered in the interval
            v1 = (x1 - x1_old)/(p.time - time_old)
            v2 = (x2 - x2_old)/(p.time - time_old)
            v3 = (x3 - x3_old)/(p.time - time_old)
            v_avg = (v1 + v2 + v3)/3.0
            v_sigma = np.std([v1, v2, v3])

            v.append(v_avg)
            vs.append(v_sigma)

            w.append(0.5*(width + width_old))
            t.append(0.5*(time_old + p.time))
            pltfiles.append(pfile)  # Add this
            ptimes.append(p.time)    # Add this
            w_rz.append(0.5*(width_rz + width_rz_old))  # Add this
        time_old = p.time
        x1_old = x1
        x2_old = x2
        x3_old = x3
        width_old = width
        width_rz_old = width_rz  # Add this

    f = plt.figure()
    f.set_size_inches(7.0, 9.0)

    ax_s = f.add_subplot(211)
    ax_w = f.add_subplot(212)

    ax_s.plot(t, v)
    ax_w.plot(t, w)

    ax_s.set_yscale("log")  #symlog->log
    ax_w.set_yscale("log")

    ax_w.set_xlabel("time (s)")
    ax_s.set_ylabel("flame speed (cm/s)")
    ax_w.set_ylabel("flame width (cm)")

    plt.tight_layout()
    plt.savefig("speed.png")

    print(f"{'Plotfile':15s} {'plt_time':>10s} {'avg_time':>10s} : {'flame_speed':>15s} +/- {'std_dev':>15s} | {'width_full':>15s} {'width_rz':>15s}")
    print("-" * 130)
    for pfile, ptime, to, vo, vso, wo, wrz in zip(pltfiles, ptimes, t, v, vs, w, w_rz):
        print(f"{pfile:15s} {ptime:10.3g} {to:10.3g} : {vo:15.8g} +/- {vso:15.8g} | {wo:15.8g} {wrz:15.8g}")
