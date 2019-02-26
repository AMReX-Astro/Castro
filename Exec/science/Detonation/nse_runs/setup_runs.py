import os

CFL = [0.8, 0.4, 0.2, 0.1]
SDC_ITERS = [2, 3, 4]
NZONES = [512, 1024, 2048]

with open("inputs.template", "r") as tf:
    template = tf.readlines()

# setup the Strang runs
for c in CFL:
    for nz in NZONES:
        # dump the metdata file

        # make the output directory
        odir = "det_strang_cfl{}_nzones{}".format(c, nz)
        os.mkdir(odir)

        # write the inputs file
        with open("{}/inputs".format(odir), "w") as inpf:
            for line in template:
                inpf.write(line.replace("@@CFL@@", str(c)).replace("@@NZONES@@", str(nz)).replace("@@SDC_ITERS@@", "1"))

        # copy the remaining files


# setup the simplified SDC runs
for c in CFL:
    for nz in NZONES:
        for niter in SDC_ITERS:
            # dump the metdata file

            # make the output directory
            odir = "det_sdc_iter{}_cfl{}_nzones{}".format(niter, c, nz)
            os.mkdir(odir)

            # write the inputs file
            with open("{}/inputs".format(odir), "w") as inpf:
                for line in template:
                    inpf.write(line.replace("@@CFL@@", str(c)).replace("@@NZONES@@", str(nz)).replace("@@SDC_ITERS@@", str(niter)))

            # copy the remaining files

