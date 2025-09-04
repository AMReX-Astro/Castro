#!/bin/env python3

import glob
import os
import operator

import numpy as np
import matplotlib.pyplot as plt

import detonation as dt

import yt
yt.funcs.mylog.setLevel(50)

if __name__ == "__main__":

    # get all the data
    run_dirs = glob.glob("det_[st]*")
    runs = []
    for run in run_dirs:
        try:
            if not os.path.isfile(os.path.join(run, "Backtrace.0")):
                print(run)
                runs.append(dt.Detonation(run))
        except IndexError:
            # the run didn't produce output -- it might still be running?
            print(f"run {run} didn't produce output")

    runs.sort()
    print(runs)

    # make a plot of speed vs. CFL, grouped by Strang, SDC2, SDC3,
    # SDC4 for the same resolution
    nzones = {q.nzones for q in runs}
    for nz in nzones:
        strang = [q for q in runs if q.integrator == "Strang" and q.nzones == nz and q.dtnuce == False]
        sdc2 = [q for q in runs if q.integrator == "simplified-SDC" and q.niters == 2 and q.nzones == nz]
        print(sdc2)
        sdc3 = [q for q in runs if q.integrator == "simplified-SDC" and q.niters == 3 and q.nzones == nz]
        sdc4 = [q for q in runs if q.integrator == "simplified-SDC" and q.niters == 4 and q.nzones == nz]
        truesdc_o2 = [q for q in runs if q.integrator == "true-SDC" and q.order == 2 and q.nzones == nz]
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
    cfls = {q.cfl for q in runs}
    for cfl in cfls:
        strang = [q for q in runs if q.integrator == "Strang" and q.cfl == cfl and q.dtnuce == False]
        # CFL doesn't matter for dtnuc_e
        strang_limit = [q for q in runs if q.integrator == "Strang" and q.dtnuce == True]
        sdc2 = [q for q in runs if q.integrator == "simplified-SDC" and q.niters == 2 and q.cfl == cfl]
        sdc3 = [q for q in runs if q.integrator == "simplified-SDC" and q.niters == 3 and q.cfl == cfl]
        sdc4 = [q for q in runs if q.integrator == "simplified-SDC" and q.niters == 4 and q.cfl == cfl]

        fig = plt.figure(1)
        fig.clear()

        ax = fig.add_subplot(111)

        dsets = [(strang, "Strang"),
                 (strang_limit, r"Strang (with energy $\Delta t$ limit)"),
                 (sdc2, "simplified-SDC (2 iters)"),
                 (sdc3, "simplified-SDC (3 iters)"),
                 (sdc4, "simplified-SDC (4 iters)")]

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
        sdc2 = [q for q in runs if q.integrator == "simplified-SDC" and q.nzones == nz and q.niters == 2]
        sdc3 = [q for q in runs if q.integrator == "simplified-SDC" and q.nzones == nz and q.niters == 3]

        for dset, title, fname in [(strang, "Strang", "strang"),
                                   (sdc2, "simplified-SDC (niters = 2)", "sdc_niter2"),
                                   (sdc3, "simplified-SDC (niters = 3)", "sdc_niter3")]:

            if not dset:
                continue

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
