import os
import shutil
import subprocess
import shlex
from multiprocessing import Pool
import time

CFL = [0.8, 0.4] #, 0.2, 0.1]
SDC_ITERS = [2, 3] #, 4]
NZONES = [256, 512] #, 1024, 2048]

job_list = []

NEEDED_SDC_FILES = ["helm_table.dat", "probin.nse_test.sdc", "Castro1d.gnu.SDC.ex"]
NEEDED_STRANG_FILES = ["helm_table.dat", "probin.nse_test.strang", "Castro1d.gnu.ex"]

def setup_runs():

    with open("inputs.big.template.sdc", "r") as tf:
        sdc_template = tf.readlines()

    with open("inputs.big.template.strang", "r") as tf:
        strang_template = tf.readlines()

    # setup the Strang runs
    for c in CFL:
        for nz in NZONES:
            # make the output directory
            odir = "det_strang_cfl{}_nzones{}".format(c, nz)
            os.mkdir(odir)

            job_list.append(odir)

            # dump the metdata file
            with open("{}/run.meta".format(odir), "w") as meta:
                meta.write("cfl = {}\n".format(c))
                meta.write("nzones = {}\n".format(nz))
                meta.write("integrator = Strang\n")

            # write the inputs file
            with open("{}/inputs".format(odir), "w") as inpf:
                for line in strang_template:
                    inpf.write(line.replace("@@CFL@@", str(c)).replace("@@NZONES@@", str(nz)).replace("@@SDC_ITERS@@", "1").replace("@@method@@", "0"))

            # copy the remaining files
            for f in NEEDED_STRANG_FILES:
                shutil.copy(f, odir)

    # setup the simplified SDC runs
    for c in CFL:
        for nz in NZONES:
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
                          stderr=subprocess.PIPE)
    stdout0, stderr0 = p0.communicate()
    rc = p0.returncode
    stdout = stdout0.decode('utf-8')
    stderr = stderr0.decode('utf-8')

    return stdout, stderr, rc


def do_run(name):
    print("running {}...".format(name))

    if "sdc" in name:
        command = "./Castro1d.gnu.SDC.ex inputs"
    else:
        command = "./Castro1d.gnu.ex inputs"

    cwd = os.getcwd()
    os.chdir(name)

    stdout, stderr, rc = run(command)

    # add a file indicated the job is complete
    with open("job_status", "w") as jf:
        jf.write("completed\n")

    os.chdir(cwd)


if __name__ == "__main__":
    job_list = setup_runs()

    # do the runs
    p = Pool(16)
    p.map(do_run, job_list)
