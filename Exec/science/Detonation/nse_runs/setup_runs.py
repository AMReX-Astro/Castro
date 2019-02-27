import os
import shutil

CFL = [0.8, 0.4, 0.2, 0.1]
SDC_ITERS = [2, 3, 4]
NZONES = [512, 1024, 2048]

job_list = []

NEEDED_SDC_FILES = ["helm_table.dat", "probin-det-x.SDC", "Castro1d.gnu.SDC.ex"]
NEEDED_STRANG_FILES = ["helm_table.dat", "probin-det-x.SDC", "Castro1d.gnu.ex"]


with open("inputs.template", "r") as tf:
    template = tf.readlines()

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
            for line in template:
                inpf.write(line.replace("@@CFL@@", str(c)).replace("@@NZONES@@", str(nz)).replace("@@SDC_ITERS@@", "1"))

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
                for line in template:
                    inpf.write(line.replace("@@CFL@@", str(c)).replace("@@NZONES@@", str(nz)).replace("@@SDC_ITERS@@", str(niter)))

            # copy the remaining files
            for f in NEEDED_SDC_FILES:
                shutil.copy(f, odir)


# do the runs
print("jobs: ", job_list)
