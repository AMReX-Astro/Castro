#!/usr/bin/env python3

# run as: ./test1-helm.py castro_exec_dir plotfle

# note: this relies on fextract.XXXX.ex being in your path somewhere

import sys
import os
import shutil
import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def process(castro_exec_dir, plotfile):

    run_dir = os.getcwd()

    # 1. find the fextract tool by looking through the User's path
    path = os.environ["PATH"].split(":")

    for d in path:
        full_dir = os.path.expanduser(d)
        if not os.path.isdir(full_dir):
            continue
        for f in os.listdir(full_dir):
            if (os.path.isfile(full_dir+"/"+f) and
                f.startswith("fextract") and f.endswith(".ex")):
                analysis_routine = full_dir+"/"+f
                break

    print("analysis_routine = {}".format(analysis_routine))

    shutil.copy(analysis_routine, run_dir)


    # 2. analyze the data


    # output the average profile
    os.system("./{} -s {} {}".format(os.path.basename(analysis_routine), "sod_slice.out", plotfile))


    analytic = os.path.join(castro_exec_dir, "Exec/hydro_tests/Sod/Verification/sod-exact.out")
    analytic_data = np.loadtxt(analytic)

    # need to be more flexible with the data from the simulations, as the
    # columns shift depending on the dimensionality.  This gets the column
    # names from the header
    data = np.genfromtxt("sod_slice.out", skip_header=2, names=True)

    # 3. make the plot
    plt.subplot(411)

    plt.plot(analytic_data[:,0], analytic_data[:,1])
    plt.scatter(data["x"], data["density"], marker="+", color="r")

    plt.xlabel("x")
    plt.ylabel("density")

    plt.subplot(412)

    # figure out which dimensions are present
    d = 1
    for n in data.dtype.names:
        if n == "ymom" or n == "zmom": d += 1

    dim = "xmom"
    if d >= 2 and data["ymom"].ptp() > data["xmom"].ptp(): dim = "ymom"
    if d == 3 and data["zmom"].ptp() > data[dim].ptp(): dim = "zmom"

    plt.plot(analytic_data[:,0], analytic_data[:,2])
    plt.scatter(data["x"], data[dim]/data["density"], marker="+", color="r")

    plt.xlabel("x")
    plt.ylabel("velocity")

    plt.subplot(413)

    plt.plot(analytic_data[:,0], analytic_data[:,3])
    plt.scatter(data["x"], data["pressure"], marker="+", color="r")

    plt.xlabel("x")
    plt.ylabel("pressure")

    plt.subplot(414)

    plt.plot(analytic_data[:,0], analytic_data[:,4])
    plt.scatter(data["x"], data["rho_e"]/data["density"], marker="+", color="r")

    plt.xlabel("x")
    plt.ylabel("internal energy")


    ax = plt.gca()
    ax.set_yscale("log")

    f = plt.gcf()
    f.set_size_inches(6.0, 9.0)

    plt.tight_layout()

    index = plotfile.rfind("_plt")
    if (index > 0):
        plt.savefig(plotfile[:index] + ".png")
    else:
        plt.savefig("sod_test.png")



if __name__ == "__main__":

    castro_exec_dir = str(sys.argv[1])
    plotfile = str(sys.argv[2])

    process(castro_exec_dir, plotfile)
