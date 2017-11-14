from __future__ import print_function

import os
import shlex
import subprocess

class PlotFile:
    def __init__(self, run, value, plotfile, full_path=""):
        self.run = run
        self.value = value
        self.plotfile = plotfile
        self.full_path = full_path
        self.find_time = -1

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

    # find directories with runs

    current_run = None
    current_value = None

    all_plot_files = []

    for root, dirs, files in os.walk(top_dir):
        for d in dirs:
            if prefix in d:
                path = root.split("/")
                if path[-1] == "reference":
                    current_run = "reference"
                    current_value = None
                else:
                    if path[-1] == "low":
                        current_value = "low"
                        current_run = path[-2]
                    elif path[-1] == "high":
                        current_value = "high"
                        current_run = path[-2]                        
                
                all_plot_files.append(PlotFile(current_run, current_value, d, full_path=os.path.join(root, d)))

    # find all the unique runs
    runs = set([pf.run for pf in all_plot_files])    

    # file files and their output times
    for pf in all_plot_files:
        pf.time = find_time(pf.full_path)

    ref = [pf for pf in all_plotfiles if pf.run == "reference"]

    # loop over runs and write html
    for r in runs:

        # find "low" and "high" files for this run
        low = [pf for pf in all_plotfiles if pf.run == r and pf.value == "low"]
        high = [pf for pf in all_plotfiles if pf.run == r and pf.value == "high"]

        # create images (using run())

        # write to html



if __name__ == "__main__":
    process()
