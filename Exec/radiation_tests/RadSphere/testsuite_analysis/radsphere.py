#!/usr/bin/env python3
#
# this script is meant to be run from the test suite.  It will read in
# the plotfile data from the RadSphere problem (using the fradsphere
# tool) and compare to the analytic solution, and output a plot.
#
# run as: ./radsphere.py castroradiation_dir plotfle

import sys
import os
import shutil
import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def process(castro_dir, plotfile):

    run_dir = os.getcwd()

    # 1. make sure that the analysis tool is built
    build_dir = castro_dir + "/Diagnostics/Radiation/"
    os.chdir(build_dir)
    os.system("make DIM=1 rad_sphere.ex >& /dev/null")

    # find the executable
    analysis_routine = None
    for file in os.listdir(build_dir):
        if (os.path.isfile(file) and
            file.startswith("rad_sphere") and
            file.endswith(".ex")):
            analysis_routine = file
            break

    if analysis_routine is None:
        sys.exit("Error: failed to build the analysis routine")

    print("analysis_routine = ", analysis_routine)

    shutil.copy(analysis_routine, run_dir)

    os.chdir(run_dir)


    # 2. analyze the data

    # output the average profile -- note: this expects the
    # group_structure.dat file in the same directory as the plotfile
    group_file = "group_structure.dat"
    if not os.path.isfile(group_file):
        group_file = os.path.join(os.path.dirname(os.path.abspath(plotfile)), "group_structure.dat")
        if not os.path.isfile(group_file):
            sys.exit("error: cannot find group_structure.dat")

    variable = "rad"

    cmd = f"./{analysis_routine} -p {plotfile} -r 0.06 -g {group_file} -v {variable}"
    print(cmd)
    os.system(cmd)

    group = []
    rad_energy = []
    with open("rad_sphere.out") as of:
        lines = of.readlines()
        for line in lines:
            if line.strip().startswith("#"):
                continue
            _, group_energy, e_rad, _ = line.strip().split()
            group.append(float(group_energy))
            rad_energy.append(float(e_rad))

    # Now do the same for the analytic data

    variable = "rad_analytic_"

    cmd = f"./{analysis_routine} -p {plotfile} -r 0.06 -g {group_file} -v {variable}"
    print(cmd)
    os.system(cmd)

    rad_energy_analytic = []
    with open("rad_sphere.out") as of:
        lines = of.readlines()
        for line in lines:
            if line.strip().startswith("#"):
                continue
            _, _, e_rad, _ = line.strip().split()
            rad_energy_analytic.append(float(e_rad))

    # 3. make the plot
    plt.plot(np.array(group), np.array(rad_energy_analytic))
    plt.scatter(np.array(group), np.array(rad_energy), marker="+", color="r")

    plt.xlabel(r"group center (Hz)")
    plt.ylabel(r"group radiation energy density (erg/cm$^3$/Hz)")

    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log")

    #plt.tight_layout()

    plt.savefig("radsphere.png")



if __name__ == "__main__":

    if len(sys.argv) < 3:
        sys.exit("usage: ./sedov_3d_sph.py castro_dir plotfile")

    castro_dir = str(sys.argv[1])
    plotfile = str(sys.argv[2])

    process(castro_dir, plotfile)
