import glob
import os
import operator

import numpy as np
import matplotlib.pyplot as plt

import yt

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


class Detonation:
    def __init__(self, name):
        self.name = name

        self.cfl = None
        self.nzones = None
        self.integrator = None
        self.niters = None

        # read the meta data
        with open(os.path.join(name, "run.meta")) as mf:
            for line in mf.readlines():
                k, v = [q.strip() for q in line.split("=")]
                if k == "cfl":
                    self.cfl = float(v)
                elif k == "nzones":
                    self.nzones = int(v)
                elif k == "integrator":
                    self.integrator = v
                elif k == "niters":
                    self.niters = int(v)

        # find all the output (plot) files
        cwd = os.getcwd()
        os.chdir(name)
        self.files = glob.glob("*plt?????")
        self.files.sort()
        os.chdir(cwd)

        # precompute the velocity and the data profiles
        self.v = self.get_velocity()
        self.data = self.get_data()

    def get_velocity(self):
        """look at the last 2 plotfiles and estimate the velocity by
        finite-differencing"""

        f1 = self.files[-2]
        p1 = Profile(os.path.join(self.name, f1))

        f2 = self.files[-1]
        p2 = Profile(os.path.join(self.name, f2))

        # we'll do this by looking at 3 different temperature
        # thresholds and averaging
        T_ref = [2.e9, 3.e9, 4.e9]
        v = 0.0
        for T0 in T_ref:
            x1 = p1.find_x_for_T(T0)
            x2 = p2.find_x_for_T(T0)
            v += (x1 - x2)/(p1.time - p2.time)

        v /= len(T_ref)
        return v

    def get_data(self):
        """get the temperature and energy generation rate from the last
        plotfile"""

        return Profile(os.path.join(self.name, self.files[-1]))


if __name__ == "__main__":

    # get all the data
    run_dirs = glob.glob("det_s*")
    runs = []
    for run in run_dirs:
        try:
            runs.append(Detonation(run))
        except IndexError:
            # the run didn't produce output -- it might still be running?
            print("run {} didn't produce output".format(run))

    print(len(runs))

    # make a plot of speed vs. CFL, grouped by Strang, SDC2, SDC3,
    # SDC4 for 1024 zones
    strang = [q for q in runs if q.integrator == "Strang" and q.nzones == 1024]
    strang.sort(key=operator.attrgetter("cfl"))

    sdc2 = [q for q in runs if q.integrator == "SDC" and q.niters == 2 and q.nzones == 1024]
    sdc2.sort(key=operator.attrgetter("cfl"))

    sdc3 = [q for q in runs if q.integrator == "SDC" and q.niters == 3 and q.nzones == 1024]
    sdc3.sort(key=operator.attrgetter("cfl"))

    sdc4 = [q for q in runs if q.integrator == "SDC" and q.niters == 4 and q.nzones == 1024]
    sdc4.sort(key=operator.attrgetter("cfl"))

    fig = plt.figure(1)
    fig.clear()

    ax = fig.add_subplot(111)
    ax.plot([q.cfl for q in strang], [q.v for q in strang],
            marker="x", label="Strang")
    ax.plot([q.cfl for q in sdc2], [q.v for q in sdc2],
            marker="o", label="SDC (2 iters)")
    ax.plot([q.cfl for q in sdc2], [q.v for q in sdc3],
            marker="*", label="SDC (3 iters)")
    ax.plot([q.cfl for q in sdc2], [q.v for q in sdc4],
            marker="^", label="SDC (4 iters)")

    ax.legend(frameon=False)
    ax.set_xlabel("CFL")
    ax.set_ylabel("velocity (cm/s)")

    fig.savefig("speed_vs_cfl.png")

    # make a plot of speed vs. resolution, grouped by Strang, SDC2,
    # SDC3, SDC4 for CFL = 0.8
    strang = [q for q in runs if q.integrator == "Strang" and q.cfl == 0.8]
    sdc2 = [q for q in runs if q.integrator == "SDC" and q.niters == 2 and q.cfl == 0.8]
    sdc3 = [q for q in runs if q.integrator == "SDC" and q.niters == 3 and q.cfl == 0.8]
    sdc4 = [q for q in runs if q.integrator == "SDC" and q.niters == 4 and q.cfl == 0.8]

    # make a plot of T, enuc vs. x for different Strang CFL

    # make a plot of T, enuc vs. x for different SDC CFL

    # make a plot of T, enuc vs. x for different SDC iters

    # make a plot of T, enuc vs. x for Strang and SDC
