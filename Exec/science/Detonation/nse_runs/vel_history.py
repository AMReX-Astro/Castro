import glob
import os
import operator
import sys

import numpy as np
import matplotlib.pyplot as plt

import yt

yt.funcs.mylog.setLevel(50)

class Profile:
    """read a plotfile using yt and store the 1d profile for T and enuc"""

    def __init__(self, plotfile):
        ds = yt.load(plotfile)

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

    def find_x_for_T(self, T_0=1.e9):
        """ given a profile x(T), find the x_0 that corresponds to T_0 """

        # our strategy here assumes that the hot ash is in the early
        # part of the profile.  We then find the index of the first
        # point where T drops below T_0
        idx = np.where(self.T < T_0)[0][0]

        T1 = self.T[idx-1]
        x1 = self.x[idx-1]

        T2 = self.T[idx]
        x2 = self.x[idx]

        slope = (x2 - x1)/(T2 - T1)

        return x1 + slope*(T_0 - T1)

def get_velocity(p1, p2):
        """look at the last 2 plotfiles and estimate the velocity by
        finite-differencing"""

        # we'll do this by looking at 3 different temperature
        # thresholds and averaging
        T_ref = [2.e9, 3.e9, 4.e9]

        vs = []
        for T0 in T_ref:
            x1 = p1.find_x_for_T(T0)
            x2 = p2.find_x_for_T(T0)
            vs.append((x1 - x2)/(p1.time - p2.time))

        vs = np.array(vs)
        v = np.mean(vs)
        v_sigma = np.std(vs)
        return v, v_sigma


if __name__ == "__main__":

    run_dir = sys.argv[1]

    print(run_dir)

    # get all the plotfiles
    pfiles = glob.glob("{}/*plt?????".format(run_dir))
    pfiles.sort()

    profiles = []
    for pf in pfiles:
        profiles.append(Profile(pf))

    for n in range(1, len(profiles)):
        p1 = profiles[n-1]
        p2 = profiles[n]
        print(0.5*(p1.time + p2.time), get_velocity(p1, p2))

