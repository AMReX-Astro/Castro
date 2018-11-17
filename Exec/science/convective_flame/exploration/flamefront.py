# When copied into a directory with plotfiles, this code makes slices of those
# plotfiles at specified times and plots them on top of one another.
# Ideally this shows the movement of a flame front.
#
# Written by Blaire Ness
#
# Last updated May 2018

from __future__ import print_function

import os
import shlex
import subprocess
import matplotlib
import numpy
import matplotlib.pyplot as plt

class PlotFile:
    def __init__(self, run, value, plotfile, full_path=""):
        self.run = run                 # variable being altered
        self.value = value             # low, high
        self.plotfile = plotfile
        self.full_path = full_path     # path to file
        self.find_time = -1            # time of file

def run(string):

    # shlex.split will preserve inner quotes
    prog = shlex.split(string)
    p0 = subprocess.Popen(prog, stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE)

    stdout0, stderr0 = p0.communicate()
    rc = p0.returncode
    p0.stdout.close()
    p0.stderr.close()

    return stdout0, stderr0, rc


def find_time(pf):
    stdout, _, _ = run("ftime.Linux.gfortran.exe {}".format(pf))
    return float(stdout.split()[-1])


def process(top_dir=".", prefix="plt"):

    all_plot_files = []

    current_run = None
    current_value = None

    top_dir = os.getcwd()

    # loop through directories looking for plotfiles
    for root, dirs, files in os.walk(top_dir):
        for d in dirs:
            if prefix in d:
                path = root.split("/")
                if path[0] == ".":
                    current_run = "reference"
                    current_value = None

                if path[-1] == "lowres":
                    current_value = "low"
                    current_run = path[-1]
                elif path[-1] == "highres":
                    current_value = "high"
                    current_run = path[-1]

                all_plot_files.append(PlotFile(current_run, current_value, d, full_path=os.path.join(root, d)))

    # find time for plot files    
    for pf in all_plot_files:
        pf.time = find_time(pf.full_path)
    

    #process reference here
    for pf in all_plot_files:
        command = "fextract.Linux.gfortran.exe -y 2.5 -v Temp %s" % pf.full_path
        #print(command)
        if (pf.time == 0.0):
            run(command)
        elif (0.000099 < pf.time < 0.000105):
            run(command)
        elif (0.000199 < pf.time < 0.000205):
            run(command)
        elif (0.000299 < pf.time < 0.000305):
            run(command)

if __name__ == "__main__":
    process()

for pfile in os.listdir("."):
    if pfile.endswith(".slice"):

        with open(pfile) as f:
            for i, line in enumerate(f):
                if i == 1:
                    time = line
                elif i > 1:
                    break
                    
                    print(time)

        x, y = numpy.loadtxt(pfile, unpack="True")
        plt.plot(x,y, label=time)
plt.xlim([0,30])
plt.show()
plt.savefig("flamefronts.png")
