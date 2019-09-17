#!/bin/env python3

import glob
import os
import operator

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
        rho_he4 = np.array(ad['rho_he4'][srt])

        self.time = time
        self.x = x_coord
        self.T = temp
        self.enuc = enuc
        self.rho_he4 = rho_he4

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

        # dissect the name to find the parameters
        # it is of the form det_zXXX_cX.XX_integrator

        if "strang_ctu" in name:
            self.integrator = "strang"
        elif "lobatto_sdc4" in name:
            self.integrator = "lobatto"
        elif "radau_sdc4" in name:
            self.integrator = "radau"

        fields = name.split("_")

        # resolution is in the form z???
        self.nzones = int(fields[1][1:])

        # CFL is in the form c?.??
        self.cfl = float(fields[2][1:])

        # find all the output (plot) files
        cwd = os.getcwd()
        os.chdir(name)
        self.files = glob.glob("*plt?????")
        self.files.sort()
        os.chdir(cwd)

    def __str__(self):
        return "{} CFL={} N={}".format(self.integrator, self.cfl, self.nzones)

    def __repr__(self):
        return self.name

    def __lt__(self, other):
        """sort by CFL number and resolution"""
        if self.cfl == other.cfl:
            return self.nzones < other.nzones
        else:
            return self.cfl < other.cfl

    def get_velocity(self):
        """look at the last 2 plotfiles and estimate the velocity by
        finite-differencing"""

        vs = []
        pairs = [(-2, -1), (-3, -1), (-3, -1)]

        for i1, i2 in pairs:
            f1 = self.files[i1]
            p1 = Profile(os.path.join(self.name, f1))

            f2 = self.files[i2]
            p2 = Profile(os.path.join(self.name, f2))

            # we'll do this by looking at 3 different temperature
            # thresholds and averaging
            T_ref = [2.e9, 3.e9, 4.e9]

            for T0 in T_ref:
                x1 = p1.find_x_for_T(T0)
                x2 = p2.find_x_for_T(T0)
                vs.append((x1 - x2)/(p1.time - p2.time))

        vs = np.array(vs)
        v = np.mean(vs)
        v_sigma = np.std(vs)
        return v, v_sigma

    def get_data(self, n=-1):
        """get the temperature and energy generation rate from the nth
        plotfile (starting to count at 0).""" 

        return Profile(os.path.join(self.name, self.files[n]))


if __name__ == "__main__":

    # get all the data
    run_dirs = glob.glob("det_z*")
    runs = []
    for run in run_dirs:
        try:
            if not os.path.isfile(os.path.join(run, "Backtrace.0")):
                runs.append(Detonation(run))
        except IndexError:
            # the run didn't produce output -- it might still be running?
            print("run {} didn't produce output".format(run))

    runs.sort()
    print(runs)

    # make a plot of the profiles for different CFLs and different resolutions
    nzones = set([q.nzones for q in runs])
    idx = 5
    for nz in nzones:
        for intg in ["strang", "lobatto", "radau"]:
            dset = [q for q in runs if q.integrator == intg and q.nzones == nz]

            if len(dset) == 0:
                continue

            fig = plt.figure(1)
            fig.clear()

            fig.set_size_inches(7.0, 9.0)

            ax_T = fig.add_subplot(311)
            ax_e = fig.add_subplot(312)
            ax_he = fig.add_subplot(313)

            for q in dset:
                pf = q.get_data(idx)
                print("working on {}, time = {}".format(q.name, pf.time))

                ax_T.plot(pf.x, pf.T, label="CFL = {}".format(q.cfl))
                ax_e.plot(pf.x, pf.enuc)
                ax_he.plot(pf.x, pf.rho_he4)

            ax_T.legend(frameon=False)

            ax_T.set_ylabel(r"$T$ (K)")
            ax_e.set_ylabel(r"$H_\mathrm{nuc}$ (erg/g/s)")
            ax_he.set_ylabel(r"$\rho X({}^4\mathrm{He})$ (g/cm${}^3$)")
            ax_he.set_xlabel("x (cm)")

            ax_T.set_xlim(5000, 15000)
            ax_e.set_xlim(5000, 15000)
            ax_he.set_xlim(5000, 15000)

            ax_e.set_yscale("log")
            ax_he.set_yscale("log")

            fig.suptitle("integrator = {}; number of zones = {}".format(intg, nz))

            fig.tight_layout()
            fig.savefig("det_cfl_compare_{}_nz{}.png".format(intg, nz))
