#!/usr/bin/env python
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

def process(castrorad_dir, plotfile):

    run_dir = os.getcwd()

    # 1. make sure that the analysis tool is built
    build_dir = castrorad_dir + "/Exec/RadSphere/Tools/"
    os.chdir(build_dir)
    os.system("make programs=fradsphere >& /dev/null")

    # find the executable
    for file in os.listdir(build_dir):
        if (os.path.isfile(file) and 
            file.startswith("fradsphere") and
            file.endswith(".exe")):
            analysis_routine = file
            break

    print "analysis_routine = ", analysis_routine

    shutil.copy(analysis_routine, run_dir)

    os.chdir(run_dir)


    # 2. analyze the data
    
    
    # output the average profile
    os.system("./{} {} > {}".format(analysis_routine, plotfile, "radsphere_testsuite.out"))


    analytic = castrorad_dir + "/Exec/RadSphere/testsuite_analysis/radsphere_analytic.out"
    analytic_data = np.loadtxt(analytic)

    data = np.loadtxt("radsphere_testsuite.out")


    # 3. make the plot
    plt.plot(analytic_data[:,1], analytic_data[:,2])
    plt.scatter(data[:,1], data[:,2], marker="+", color="r")

    plt.xlabel(r"group center (Hz)")
    plt.ylabel(r"group radiation energy density (erg/cm$^3$/Hz)")

    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log")

    #plt.tight_layout()

    plt.savefig("radsphere.png")



if __name__ == "__main__":

    castrorad_dir = str(sys.argv[1])
    plotfile = str(sys.argv[2])

    process(castrorad_dir, plotfile)

