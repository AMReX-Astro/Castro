#!/usr/bin/env python

import os
import shutil
import subprocess
import shlex
from multiprocessing import Pool
import time

CFL = [0.8, 0.4, 0.2]
SDC_ITERS = [2, 3]
DTNUC_E = [1.e200, 0.5, 0.25]

NZONES = [256, 512, 1024] #, 2048, 4096, 8192]

job_list = []

NEEDED_SDC_FILES = ["helm_table.dat", "probin-det-x.nse_disabled", "Castro1d.gnu.SDC.ex"]
NEEDED_STRANG_FILES = ["helm_table.dat", "probin-det-x.nse_disabled", "Castro1d.gnu.ex"]

def setup_runs():

    with open("inputs.template.sdc", "r") as tf:
        sdc_template = tf.readlines()

    with open("inputs.template.strang", "r") as tf:
        strang_template = tf.readlines()

    # setup the Strang runs
    for nz in NZONES:
        for dtn in DTNUC_E:
            for c in CFL:

                # don't do all the CFLs when using the nuc energy dt limiter
                if dtn < 1.0 and c != CFL[0]:
                    continue

                # make the output directory
                odir = "det_strang_cfl{}_dtnuce{}_nzones{}".format(c, dtn, nz)
                os.mkdir(odir)

                job_list.append(odir)

                # dump the metdata file
                with open("{}/run.meta".format(odir), "w") as meta:
                    meta.write("cfl = {}\n".format(c))
                    meta.write("nzones = {}\n".format(nz))
                    meta.write("integrator = Strang\n")
                    meta.write("dtnuc_e = {}\n".format(dtn))

                # write the inputs file
                with open("{}/inputs".format(odir), "w") as inpf:
                    for line in strang_template:
                        inpf.write(line.replace("@@CFL@@", str(c)).replace("@@NZONES@@", str(nz)).replace("@@method@@", "0").replace("@@DTNUC_E@@", str(dtn)))

                # copy the remaining files
                for f in NEEDED_STRANG_FILES:
                    shutil.copy(f, odir)


    # setup the simplified SDC runs
    for nz in NZONES:
        for c in CFL:
            for niter in SDC_ITERS:

                # make the output directory
                odir = "det_sdc_iter{}_cfl{}_nzones{}".format(niter, c, nz)
                os.mkdir(odir)

                job_list.append(odir)

                # dump the metdata file
                with open("{}/run.meta".format(odir), "w") as meta:
                    meta.write("cfl = {}\n".format(c))
                    meta.write("nzones = {}\n".format(nz))
                    meta.write("integrator = SDC\n")
                    meta.write("niters = {}\n".format(niter))

                # write the inputs file
                with open("{}/inputs".format(odir), "w") as inpf:
                    for line in sdc_template:
                        inpf.write(line.replace("@@CFL@@", str(c)).replace("@@NZONES@@", str(nz)).replace("@@SDC_ITERS@@", str(niter)).replace("@@method@@", "3"))

                # copy the remaining files
                for f in NEEDED_SDC_FILES:
                    shutil.copy(f, odir)

    return job_list

def run(string):
    """run a command and capture the output and return code"""

    # shlex.split will preserve inner quotes
    prog = shlex.split(string)
    p0 = subprocess.Popen(prog, stdin=None, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)
    stdout0, stderr0 = p0.communicate()
    rc = p0.returncode
    stdout = stdout0.decode('utf-8')

    return stdout, rc


def do_run(name):
    print("running {}...".format(name))

    if "sdc" in name:
        command = "./Castro1d.gnu.MPI.SMPLSDC.ex inputs"
    else:
        command = "./Castro1d.gnu.MPI.ex inputs"

    cwd = os.getcwd()
    os.chdir(name)

    stdout, rc = run(command)

    # add a file indicated the job is complete
    with open("job_status", "w") as jf:
        jf.write("completed\n")

    # output stdout
    with open("stdout", "w") as sout:
        for line in stdout:
            sout.write(line)

    os.chdir(cwd)


if __name__ == "__main__":
    job_list = setup_runs()

    # do the runs
    p = Pool(8)
    p.map(do_run, job_list)
