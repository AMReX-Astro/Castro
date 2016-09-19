#!/usr/bin/env python                                                           

# make a plot of the exact model

import sys
import math
import numpy
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

    runs = ['exact', 'MC', 'CW', 'ppmT-I', 'ppmT-II', 'ppmT-III']

    markers = ["o", "x", "+", "*", "D", "h"]
    colors = ["r", "b", "c", "g", "m", "0.5"]

    xmax = {"test1":1.e6, "test2":1.e5, "test3":2.e5}

    for p in problems:

        print "working on problem: ", p

        # read in all the data, store in a dictionary
        data = {}

        print "  ..reading in data"

        for r in runs:
            if r == "exact":
                modelData = dataRead.getData("exact/%s.exact.out" % (p))

                vars = dataObj()

                vars.x = modelData[:,1]
                vars.rho = modelData[:,2]
                vars.u = modelData[:,3]
                vars.p = modelData[:,4]
                vars.T = modelData[:,5]

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

        print "  ..making plots"

        pylab.rc("font", size=9)
        pylab.rc("legend", loc="best")

        #pylab.ticklabel_format(style='sci', axis='x', scilimits=(-3,3), useMathText=True )
        fmt = pylab.ScalarFormatter(useMathText=True, useOffset=False)
        fmt.set_powerlimits((-3,3))

        pylab.clf()

        #-------------------------------------------------------------------------
        # density plot
        #-------------------------------------------------------------------------
        
        pylab.subplot(221)

        isym = 0
        for r in runs:
            if (r == "exact"):
                pylab.plot(data[r].x, data[r].rho, label=r, c="k")
            else:
                pylab.scatter(data[r].x, data[r].rho, label=r, 
                              marker=markers[isym], c=colors[isym], s=7, edgecolor=colors[isym])
                isym += 1


        pylab.xlabel("x")
        pylab.ylabel("density (g/cc)")
    
        pylab.legend(frameon=False, fontsize=9)

        ax = pylab.gca()

        pylab.xlim(0, xmax[p])

        ax.xaxis.set_major_formatter(fmt)
        ax.yaxis.set_major_formatter(fmt)


        #-------------------------------------------------------------------------
        # velocity plot
        #-------------------------------------------------------------------------

        pylab.subplot(222)

        isym = 0
        for r in runs:
            if (r == "exact"):
                pylab.plot(data[r].x, data[r].u, label=r, c="k")
            else:
                pylab.scatter(data[r].x, data[r].u, label=r,
                              marker=markers[isym], c=colors[isym], s=7, edgecolor=colors[isym])
                isym += 1
                
        pylab.xlabel("x")
        pylab.ylabel("velocity (cm/s)")
            
        pylab.legend(frameon=False, fontsize=9)
            
        ax = pylab.gca()

        pylab.xlim(0, xmax[p])
            
        ax.xaxis.set_major_formatter(fmt)
        ax.yaxis.set_major_formatter(fmt)
            
                
        #-------------------------------------------------------------------------
        # pressure plot
        #-------------------------------------------------------------------------
        
        pylab.subplot(223)
            
        isym = 0
        for r in runs:
            if (r == "exact"):
                pylab.plot(data[r].x, data[r].p, label=r, c="k")
            else:
                pylab.scatter(data[r].x, data[r].p, label=r,
                              marker=markers[isym], c=colors[isym], s=7, edgecolor=colors[isym])
                isym += 1
                    
        pylab.xlabel("x")
        pylab.ylabel("pressure (erg/cc)")
            
        pylab.legend(frameon=False, fontsize=9)
                    
        ax = pylab.gca()

        pylab.xlim(0, xmax[p])

        ax.xaxis.set_major_formatter(fmt)
        ax.yaxis.set_major_formatter(fmt)


        #-------------------------------------------------------------------------
        # temperature plot
        #-------------------------------------------------------------------------

        pylab.subplot(224)
        
        isym = 0
        for r in runs:
            if (r == "exact"):
                pylab.plot(data[r].x, data[r].T, label=r, c="k")
            else:
                pylab.scatter(data[r].x, data[r].T, label=r,
                              marker=markers[isym], c=colors[isym], s=7, edgecolor=colors[isym])
                isym += 1


        pylab.xlabel("x")
        pylab.ylabel("temperature (K)")
            
        pylab.legend(frameon=False, fontsize=9)
            
        ax = pylab.gca()

        ax.set_yscale('log')

        pylab.xlim(0, xmax[p])
        
        ax.xaxis.set_major_formatter(fmt)
        #ax.yaxis.set_major_formatter(fmt)
            

        f = pylab.gcf()
        f.set_size_inches(7.0,9.0)

        pylab.tight_layout()

        print "saving figure: %s.png" % (p)
        pylab.savefig("%s.png" % (p))
        pylab.savefig("%s.eps" % (p))
        


if __name__== "__main__":

    model()

