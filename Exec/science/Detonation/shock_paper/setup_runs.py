#!/usr/bin/env python

import os
import shutil

# tuples of label, nx, and max_level
RES_INFO = [("24.0km", "48", "0"),
            ("12.0km", "96", "0"),
            ("6.0km", "192", "0"),
            ("3.0km", "384", "0"),
            ("1.5km", "768", "0"),
            ("0.1875km", "6144", "0"),
            ("0.0234375km", "12288", "1")]

INPUTS_TEMPLATE = "inputs-shock-burn.template"

COMMON_FILES = ["helm_table.dat",
                "Castro1d.gnu.MPI.SMPLSDC.ex"]

shock_flag = "1"
shock_thresh = "0.666"

def doit():

    # read in the template
    with open(INPUTS_TEMPLATE) as tf:
        template = tf.readlines()

    # loop over resolutions
    for label, nx, max_level in RES_INFO:

        # create output direct
        odir = f"res{label}"
        if shock_flag == "1":
            odir += f"_noshockburn_{shock_thresh}"

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
