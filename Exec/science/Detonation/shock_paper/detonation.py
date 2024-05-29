import glob
import os

import numpy as np
import matplotlib.pyplot as plt

import yt
from yt.frontends.boxlib.api import CastroDataset


yt.funcs.mylog.setLevel(50)

class Profile:
    """read a plotfile using yt and store the 1d profile for T and enuc"""

    def __init__(self, plotfile):
        ds = CastroDataset(plotfile)

        time = float(ds.current_time)
        ad = ds.all_data()

        # Sort the ray values by 'x' so there are no discontinuities
        # in the line plot
        srt = np.argsort(ad['x'])
        x_coord = np.array(ad['x'][srt])
        temp = np.array(ad['Temp'][srt])
        enuc = np.array(ad['enuc'][srt])
        shock = np.array(ad['Shock'][srt])

        self.time = time
        self.x = x_coord
        self.T = temp
        self.enuc = enuc
        self.shock = shock

    def find_x_for_T(self, T_0=1.e9):
        """ given a profile x(T), find the x_0 that corresponds to T_0 """

        # our strategy here assumes that the hot ash is in the early
        # part of the profile.  We then find the index of the first
        # point where T drops below T_0
        try:
            idx = np.where(self.T < T_0)[0][0]
        except IndexError:
            idx = len(self.T)-1

        T1 = self.T[idx-1]
        x1 = self.x[idx-1]

        T2 = self.T[idx]
        x2 = self.x[idx]

        slope = (x2 - x1)/(T2 - T1)

        return x1 + slope*(T_0 - T1)


class Detonation:
    def __init__(self, dirname):
        self.name = dirname

        self.v = None
        self.v_sigma = None

        # find all the output (plot) files
        self.files = glob.glob(f"{dirname}/*plt?????")
        self.files.sort()

        # precompute the velocity and the data profiles
        if len(self.files) >= 3:
            self.v, self.v_sigma = self.get_velocity()
        else:
            self.v, self.v_sigma = 0.0, 0.0

    def get_velocity(self):
        """look at the last 2 plotfiles and estimate the velocity by
        finite-differencing"""

        vs = []
        pairs = [(-2, -1), (-3, -1), (-3, -2)]

        for i1, i2 in pairs:
            p1 = self.get_data(i1)
            p2 = self.get_data(i2)

            # we'll do this by looking at 3 different temperature
            # thresholds and averaging
            T_ref = [1.e9, 2.e9, 3.e9]

            for T0 in T_ref:
                x1 = p1.find_x_for_T(T0)
                x2 = p2.find_x_for_T(T0)
                vs.append((x1 - x2)/(p1.time - p2.time))

        vs = np.array(vs)
        v = np.mean(vs)
        v_sigma = np.std(vs)
        return v, v_sigma

    def get_data(self, num=-1):
        """get the temperature and energy generation rate from the
        num'th plotfile (defaults to the last)"""

        return Profile(os.path.join(self.files[num]))

