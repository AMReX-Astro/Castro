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

    problems = ['test1', 'test2', 'test3', 'test4']

    runs = ['exact', 'MC-ev', 'flash', 'flash-nochar'] #, 'flash-nosteep']
    labels = ['exact', 'Castro + CGF', 'Flash', 'Flash (no char limiting)']
    
    markers = ["o", "x", "+", "*", "D", "h"]
    colors = ["r", "b", "c", "g", "m", "0.5"]
    symsize = [12, 12, 25, 15, 10, 10]

    xmax = {"test1":1.e6, "test2":1.e5, "test3":2.e5, "test4":1.e5}

    for p in problems:

        print "working on problem: ", p

        # read in all the data, store in a dictionary
        data = {}

        print "  ..reading in data"

        for r in runs:
            if r == "exact":
                modelData = dataRead.getData("exact/%s.exact.128.out" % (p))

                vars = dataObj()

                vars.x = modelData[:,1]
                vars.rho = modelData[:,2]
                vars.u = modelData[:,3]
                vars.p = modelData[:,4]
                vars.T = modelData[:,5]
                vars.e = modelData[:,6]
                
                data[r] = vars

            elif r == "flash":

                modelData = dataRead.getData("flash/%s/flash.par.%s.data" % (p, p))

                vars = dataObj()

                vars.x = modelData[:,0]
                vars.rho = modelData[:,1]
                vars.T = modelData[:,3]
                vars.u = modelData[:,4]
                vars.p = modelData[:,2]
                vars.e = modelData[:,5]
                
                data[r] = vars
                

            elif r == "flash-nochar":

                modelData = dataRead.getData("flash/%s/flash.par.%s.nocharlimit.data" % (p, p))

                vars = dataObj()

                vars.x = modelData[:,0]
                vars.rho = modelData[:,1]
                vars.T = modelData[:,3]
                vars.u = modelData[:,4]
                vars.p = modelData[:,2]
                vars.e = modelData[:,5]
                
                data[r] = vars
                

            elif r == "flash-nosteep":

                modelData = dataRead.getData("flash/%s/flash.par.%s.nosteep.data" % (p, p))

                vars = dataObj()

                vars.x = modelData[:,0]
                vars.rho = modelData[:,1]
                vars.T = modelData[:,3]
                vars.u = modelData[:,4]
                vars.p = modelData[:,2]
                vars.e = modelData[:,5]
                
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
                vars.e = modelData[:,5]/vars.rho
                
                data[r] = vars


        # done reading

        print "  ..making plots"

        pylab.rc("font", size=9)
        pylab.rc("legend", loc="best")

        #pylab.ticklabel_format(style='sci', axis='x', scilimits=(-3,3), useMathText=True )
        fmt = pylab.ScalarFormatter(useMathText=True, useOffset=False)
        fmt.set_powerlimits((-3,3))


        #----------------------------------------------------------------------
        # plots 
        #----------------------------------------------------------------------

        pylab.clf()
        
        vars = ["density", "velocity", "pressure", "temperature", "gammae"]
        units = ["(g/cc)", "(cm/s)", "(erg/cc)", "(K)", ""]

        for v in vars:

            if v == "density":
                pylab.subplot(321)
            elif v == "velocity":
                pylab.subplot(322)
            elif v == "pressure":
                pylab.subplot(323)
            elif v == "temperature":
                pylab.subplot(324)
            elif v == "gammae":
                pylab.subplot(325)
                
            isym = 0
            for r, l in zip(runs, labels):

                if v == "density":
                    varData = data[r].rho
                elif v == "velocity":
                    varData = data[r].u
                elif v == "pressure":
                    varData = data[r].p
                elif v == "temperature":
                    varData = data[r].T
                elif v == "gammae":
                    varData = data[r].p/(data[r].rho * data[r].e) + 1.0                    

                if (r == "exact"):
                    pylab.plot(data[r].x, varData, label=l, c="k")
                else:
                    pylab.plot(data[r].x, varData, c=colors[isym], ls=":", zorder=-100, alpha=0.75)
                    pylab.scatter(data[r].x, varData, label=l, 
                                  marker=markers[isym], c=colors[isym], s=7, edgecolor=colors[isym])
                    isym += 1


            pylab.xlabel("x")
            if v == "gammae":
                vl = r"$\gamma_e$"
            else:
                vl = v
            
            pylab.ylabel(vl + " " + units[vars.index(v)])

            if v == "density":
                pylab.legend(frameon=False, fontsize=9)

            ax = pylab.gca()

            pylab.xlim(0, xmax[p])

            ax.xaxis.set_major_formatter(fmt)
            if v == "temperature":
                ax.set_yscale('log')
            else:
                ax.yaxis.set_major_formatter(fmt)


            if p == "test4" and v in ["density", "pressure", "temperature"]:
                ax.set_yscale('log')


        f = pylab.gcf()
        f.set_size_inches(7.0,9.0)

        pylab.tight_layout()

        print "saving figure: %s-MC-final.png" % (p)
        pylab.savefig("%s-MC-final.png" % (p))
        pylab.savefig("%s-MC-final.eps" % (p))
        

        #----------------------------------------------------------------------    
        # residual plots
        #----------------------------------------------------------------------
        pylab.clf()

        for v in vars:

            if v == "density":
                pylab.subplot(321)
            elif v == "velocity":
                pylab.subplot(322)
            elif v == "pressure":
                pylab.subplot(323)
            elif v == "temperature":
                pylab.subplot(324)
            elif v == "gammae":
                pylab.subplot(325)
                
            isym = 0
            for r, l in zip(runs, labels):

                if v == "density":
                    varData = data[r].rho
                    refData = data["exact"].rho
                elif v == "velocity":
                    varData = data[r].u
                    refData = data["exact"].u
                elif v == "pressure":
                    varData = data[r].p
                    refData = data["exact"].p
                elif v == "temperature":
                    varData = data[r].T
                    refData = data["exact"].T
                elif v == "gammae":
                    varData = data[r].p/(data[r].rho * data[r].e) + 1.0
                    refData = data["exact"].p/(data["exact"].rho * data["exact"].e) + 1.0
                    
                if (r == "exact"):
                    pass
                else:
                    # sanity check
                    if not numpy.max(data[r].x - data["exact"].x) == 0.0:
                        print "grid differences with {}: max error = {}".format(r, numpy.max(data[r].x - data["exact"].x))

                    pylab.plot(data[r].x, varData-refData, c=colors[isym], ls=":", zorder=-100, alpha=0.75)
                    pylab.scatter(data[r].x, varData-refData, label=l, 
                                  marker=markers[isym], c=colors[isym], s=7, edgecolor=colors[isym])
                    isym += 1


            pylab.xlabel("x")
            if v == "gammae":
                vl = r"$\gamma_e$"
            else:
                vl = v

            pylab.ylabel("{} error {}".format(vl, units[vars.index(v)]))

            if v == "density":            
                pylab.legend(frameon=False, fontsize=9)

            ax = pylab.gca()

            pylab.xlim(0, xmax[p])

            ax.xaxis.set_major_formatter(fmt)
            ax.yaxis.set_major_formatter(fmt)



        f = pylab.gcf()
        f.set_size_inches(7.0,9.0)

        pylab.tight_layout()

        print "saving figure: %s-MC-final-resid.png" % (p)
        pylab.savefig("%s-MC-final-resid.png" % (p))
        pylab.savefig("%s-MC-final-resid.eps" % (p))





if __name__== "__main__":

    model()

