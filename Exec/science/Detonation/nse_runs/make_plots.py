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
        self.dtnuce = None

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
                elif k == "dtnuc_e":
                    # if the dtnuc_e is < 1, then we are limiting
                    self.dtnuce = float(v) < 1


        # find all the output (plot) files
        cwd = os.getcwd()
        os.chdir(name)
        self.files = glob.glob("*plt?????")
        self.files.sort()
        os.chdir(cwd)

        # precompute the velocity and the data profiles
        self.v, self.v_sigma = self.get_velocity()
        self.data = self.get_data()

    def __repr__(self):
        return self.name

    def __lt__(self, other):
        """sort by CFL number and resolution and then # of SDC
        iterations"""
        if self.cfl == other.cfl:
            if self.nzones == other.nzones:
                if self.niters is None:
                    return True
                elif other.niters is None:
                    return False
                else:
                    return self.niters < other.niters
            else:
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
            if not os.path.isfile(os.path.join(run, "Backtrace.0")):
                runs.append(Detonation(run))
        except IndexError:
            # the run didn't produce output -- it might still be running?
            print(f"run {run} didn't produce output")

    runs.sort()
    print(runs)

    # make a plot of speed vs. CFL, grouped by Strang, SDC2, SDC3,
    # SDC4 for the same resolution
    nzones = set([q.nzones for q in runs])
    for nz in nzones:
        strang = [q for q in runs if q.integrator == "Strang" and q.nzones == nz and q.dtnuce == False]
        sdc2 = [q for q in runs if q.integrator == "SDC" and q.niters == 2 and q.nzones == nz]
        sdc3 = [q for q in runs if q.integrator == "SDC" and q.niters == 3 and q.nzones == nz]
        sdc4 = [q for q in runs if q.integrator == "SDC" and q.niters == 4 and q.nzones == nz]

        fig = plt.figure(1)
        fig.clear()

        ax = fig.add_subplot(111)

        dsets = [(strang, "Strang"),
                 (sdc2, "SDC (2 iters)"),
                 (sdc3, "SDC (3 iters)"),
                 (sdc4, "SDC (4 iters)")]

        for dset, label in dsets:
            if len(dset) == 0:
                continue
            ax.errorbar([q.cfl for q in dset], [q.v for q in dset],
                        yerr=[q.v_sigma for q in dset],
                        marker="x", label=label)

        ax.legend(frameon=False)
        ax.set_xlabel("CFL")
        ax.set_ylabel("velocity (cm/s)")
        ax.set_title(f"number of zones: {nz}")

        fig.savefig(f"speed_vs_cfl_{nz}.png")

    # make a plot of speed vs. resolution, grouped by Strang, SDC2,
    # SDC3, SDC4 grouped by CFL
    cfls = set([q.cfl for q in runs])
    for cfl in cfls:
        strang = [q for q in runs if q.integrator == "Strang" and q.cfl == cfl and q.dtnuce == False]
        # CFL doesn't matter for dtnuc_e
        strang_limit = [q for q in runs if q.integrator == "Strang" and q.dtnuce == True]
        sdc2 = [q for q in runs if q.integrator == "SDC" and q.niters == 2 and q.cfl == cfl]
        sdc3 = [q for q in runs if q.integrator == "SDC" and q.niters == 3 and q.cfl == cfl]
        sdc4 = [q for q in runs if q.integrator == "SDC" and q.niters == 4 and q.cfl == cfl]

        fig = plt.figure(1)
        fig.clear()

        ax = fig.add_subplot(111)

        dsets = [(strang, "Strang"),
                 (strang_limit, r"Strang (with energy $\Delta t$ limit)"),
                 (sdc2, "SDC (2 iters)"),
                 (sdc3, "SDC (3 iters)"),
                 (sdc4, "SDC (4 iters)")]

        for dset, label in dsets:
            if len(dset) == 0:
                continue
            ax.errorbar([q.nzones for q in dset], [q.v for q in dset],
                        yerr=[q.v_sigma for q in dset],
                        marker="x", label=label)

        ax.legend(frameon=False)
        ax.set_xlabel("# of zones")
        ax.set_ylabel("velocity (cm/s)")
        ax.set_title(f"CFL = {cfl}")

        fig.savefig(f"speed_vs_nzones_cfl{cfl}.png")

    # make a plot of T, enuc vs. x for different Strang / SDC CFL
    for nz in nzones:
        strang = [q for q in runs if q.integrator == "Strang" and q.nzones == nz]
        sdc = [q for q in runs if q.integrator == "SDC" and q.nzones == nz and q.niters == 2]

        for dset, title, fname in [(strang, "Strang", "strang"), (sdc, "SDC (niters = 2)", "sdc_niter2")]:

            fig, axs = plt.subplots(2, 1, figsize=(7, 10), constrained_layout=True)

            ax1 = axs.flatten()[0]
            ax2 = axs.flatten()[1]

            for p in dset:
                ax1.plot(p.data.x, p.data.T, label=f"CFL = {p.cfl}; dtnuc_e = {p.dtnuce}")
                ax2.plot(p.data.x, np.abs(p.data.enuc), label=f"CFL = {p.cfl}; dtnuc_e = {p.dtnuce}")

            ax1.legend(frameon=False)
            ax1.set_ylabel("T [K]")
            ax1.set_yscale("log")

            ax2.legend(frameon=False)
            ax2.set_xlabel("x [cm]")
            ax2.set_ylabel("enuc [erg/g/s]")
            ax2.set_yscale("log")
            ax2.set_ylim(1.e14)

            fig.suptitle(f"{title}, nzones = {nz}")
            fig.savefig(f"profile_{fname}_{nz}.png")

    # make a plot of T, enuc vs. x for different SDC CFL

    # make a plot of T, enuc vs. x for different SDC iters

    # make a plot of T, enuc vs. x for Strang and SDC
