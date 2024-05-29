#!/usr/bin/env python

import os
import shutil

# tuples of label, nx, and max_level
RES_INFO = [("12.288km", "48", "0"),
            ("1.536km", "384", "0"),
            ("0.192km", "3072", "0"),
            ("0.024km", "6144", "1"),
            ("0.003km", "12288", "2")]

INPUTS_TEMPLATE = "inputs-shock-burn.template"

COMMON_FILES = ["helm_table.dat",
                "Castro1d.gnu.MPI.SMPLSDC.ex"]

shock_flag = "0"
shock_thresh = "0.6666"

def doit():

    # read in the template
    with open(INPUTS_TEMPLATE) as tf:
        template = tf.readlines()

    # loop over resolutions
    for label, nx, max_level in RES_INFO:

        # create output direct
        odir = f"res{label}"
        if shock_flag == 1:
            odir += "_noshockburn"

        os.mkdir(odir)

        # copy files
        for f in COMMON_FILES:
            shutil.copy(f, odir)

        # modify inputs
        with open(f"{odir}/inputs", "w") as fin:
            for line in template:
                fin.write(line.replace("@nx@", nx).replace("@nlevel@", max_level).replace("@shock_flag@", shock_flag).replace("@shock_thresh@", shock_thresh))


if __name__ == "__main__":
    doit()
