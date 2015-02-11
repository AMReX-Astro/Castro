#!/usr/bin/env python

# run as: ./test2-helm.py castro_exec_dir plotfle

# note: this relies on fextract.XXXX.exe being in your path somewhere

import sys
import os
import shutil
import numpy as np

import matplotlib
matplotlib.use('agg')
import pylab

def process(castro_exec_dir, plotfile):

    run_dir = os.getcwd()

    # 1. find the fextract tool by looking through the User's path
    path = os.environ["PATH"].split(":")
    
    for d in path:
        if not os.path.isdir(d): continue
        for f in os.listdir(d):
            if (os.path.isfile(d+"/"+f) and 
                f.startswith("fextract") and f.endswith(".exe")):
                analysis_routine = d+"/"+f
                break

    print "analysis_routine = ", analysis_routine

    shutil.copy(analysis_routine, run_dir)


    # 2. analyze the data
    

    # output the average profile
    os.system("./{} -s {} {}".format(os.path.basename(analysis_routine), "test2-helm.out", plotfile))


    analytic = castro_exec_dir + "/Sod_stellar/Verification/test2.exact.128.out"
    analytic_data = np.loadtxt(analytic)

    # need to be more flexible with the data from the simulations, as the 
    # columns shift depending on the dimensionality.  This gets the column
    # names from the header
    data = np.genfromtxt("test2-helm.out", skip_header=2, names=True)

    # 3. make the plot
    pylab.subplot(411)

    pylab.plot(analytic_data[:,1], analytic_data[:,2])
    pylab.scatter(data["x"], data["density"], marker="+", color="r")

    pylab.xlabel("x")
    pylab.ylabel("density")
    pylab.xlim(0,1.e5)


    pylab.subplot(412)

    # figure out which dimensions are present
    d = 1
    for n in data.dtype.names:
        if n == "ymom" or n == "zmom": d += 1

    dim = "xmom"
    if d >= 2 and data["ymom"].ptp() > data["xmom"].ptp(): dim = "ymom"
    if d == 3 and data["zmom"].ptp() > data[dim].ptp(): dim = "zmom"
    
    pylab.plot(analytic_data[:,1], analytic_data[:,3])
    pylab.scatter(data["x"], data[dim]/data["density"], marker="+", color="r")

    pylab.xlabel("x")
    pylab.ylabel("velocity")
    pylab.xlim(0,1.e5)


    pylab.subplot(413)

    pylab.plot(analytic_data[:,1], analytic_data[:,4])
    pylab.scatter(data["x"], data["pressure"], marker="+", color="r")

    pylab.xlabel("x")
    pylab.ylabel("pressure")
    pylab.xlim(0,1.e5)


    pylab.subplot(414)

    pylab.plot(analytic_data[:,1], analytic_data[:,5])
    pylab.scatter(data["x"], data["Temp"], marker="+", color="r")

    pylab.xlabel("x")
    pylab.ylabel("temperature")
    pylab.xlim(0,1.e5)


    #ax = pylab.gca()
    #ax.set_yscale("log")

    f = pylab.gcf()
    f.set_size_inches(6.0, 9.0)

    pylab.tight_layout()

    index = plotfile.rfind("_plt")
    if (index > 0):
        pylab.savefig(plotfile[:index] + ".png")
    else:
        pylab.savefig("test2-helm.png")



if __name__ == "__main__":

    castro_exec_dir = str(sys.argv[1])
    plotfile = str(sys.argv[2])

    process(castro_exec_dir, plotfile)

