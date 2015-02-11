#!/usr/bin/env python                                                           

# make a plot of the exact model

import sys
import math
import numpy
import numpy.linalg
import pylab
import string
import dataRead
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib
import fnmatch
import os

class dataObj:
    # just a simply container

    def __init__(self):
        self.x = None
        self.rho = None
        self.u = None
        self.p = None
        self.T = None


def model():

    problems = ['test1', 'test2', 'test3']

    runs = ['exact', 'CW', 'CW-ev', 'flash', 'flash-nochar', 'MC', 'MC-ev', 'MC-ppmT-I-ev', 'MC-ppmT-II-ev', 'MC-ppmT-III-ev']


    for p in problems:

        print "working on problem: ", p

        # read in all the data, store in a dictionary
        data = {}


        for r in runs:
            if r == "exact":
                modelData = dataRead.getData("exact/%s.exact.128.out" % (p))

                vars = dataObj()

                vars.x = modelData[:,1]
                vars.rho = modelData[:,2]
                vars.u = modelData[:,3]
                vars.p = modelData[:,4]
                vars.T = modelData[:,5]

                data[r] = vars

            elif r == "flash":

                modelData = dataRead.getData("flash/%s/flash.par.%s.data" % (p, p))

                vars = dataObj()

                vars.x = modelData[:,0]
                vars.rho = modelData[:,1]
                vars.T = modelData[:,3]
                vars.u = modelData[:,4]
                vars.p = modelData[:,2]

                data[r] = vars
                

            elif r == "flash-nochar":

                modelData = dataRead.getData("flash/%s/flash.par.%s.nocharlimit.data" % (p, p))

                vars = dataObj()

                vars.x = modelData[:,0]
                vars.rho = modelData[:,1]
                vars.T = modelData[:,3]
                vars.u = modelData[:,4]
                vars.p = modelData[:,2]

                data[r] = vars
                

            else:

                # find the last slice file with our pattern
                files = []
                
                for f in os.listdir(r):
                    if fnmatch.fnmatch(f, "%s*plt?????.slice" % (p)):
                        files.append(f)

                files.sort()
                dataFile = files[len(files)-1]
                modelData = dataRead.getData("%s/%s" % (r, dataFile))

                vars = dataObj()

                vars.x = modelData[:,0]
                vars.rho = modelData[:,1]
                vars.u = modelData[:,2]/modelData[:,1]
                vars.p = modelData[:,10]
                vars.T = modelData[:,6]

                data[r] = vars


        # done reading

        print " "
        print "problem: {}".format(p)

        print "L2 norm"
        print "{:^20s} {:^20s} {:^20s} {:^20s} {:^20s}".format("run", "rho", "u", "p", "T")
        for r in runs:
            e_rho = numpy.linalg.norm((data[r].rho - data["exact"].rho)/data["exact"].rho)
            e_u   = numpy.linalg.norm((data[r].u   - data["exact"].u)/data["exact"].u)
            e_p   = numpy.linalg.norm((data[r].p   - data["exact"].p)/data["exact"].p)
            e_T   = numpy.linalg.norm((data[r].T   - data["exact"].T)/data["exact"].T)

            print "{:20s} {:20g} {:20g} {:20g} {:20g}".format(r, e_rho, e_u, e_p, e_T)

        print " "
        print "L-infty norm"
        print "{:^20s} {:^20s} {:^20s} {:^20s} {:^20s}".format("run", "rho", "u", "p", "T")
        for r in runs:
            e_rho = numpy.max(numpy.abs(data[r].rho - data["exact"].rho)/data["exact"].rho)
            e_u   = numpy.max(numpy.abs(data[r].u   - data["exact"].u)/data["exact"].u)
            e_p   = numpy.max(numpy.abs(data[r].p   - data["exact"].p)/data["exact"].p)
            e_T   = numpy.max(numpy.abs(data[r].T   - data["exact"].T)/data["exact"].T)

            print "{:20s} {:20g} {:20g} {:20g} {:20g}".format(r, e_rho, e_u, e_p, e_T)




if __name__== "__main__":

    model()

