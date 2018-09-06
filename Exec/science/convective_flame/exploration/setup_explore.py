# This is a code to set up a test suite for the convective flame problem.
#
# The code takes template inputs and probin files and replaces placeholder values
# with the values specified in 'params'.
#
# A "reference" directory is created, along with directories for each of the variables
# to be modified. There is a "high" and a "low" subdirectory in each, but these names
# names are arbitrary and the values of the variables can be whatever you wish.
#
# NOTE: After compiling the problem code it needs to be copied into each directory by hand
#
# Written by Blaire Ness and Michael Zingale
#
# Last Updated 9/5/18

import os
import shutil

# specify the names of the template files and the model file
inputs_file = "inputs.2d.template"
probin_file = "probin_temp"
model_file = "convective_flame.hse.tanh.delta_0.040cm.dx_0.050cm"
exe_file = "Castro2d.gnu.MPI.ex"

# function for finding and replacing placeholder values within the template files
def write_files(odir, param_dict):

    # Look for existing directories for the different variables
    if os.path.exists(odir) == False:
        os.mkdir(odir)

    # Replace placeholders in inputs template file
    with open("{}/inputs.2d".format(odir), "w") as in_file:
        for line in inputs_lines:
            for key in param_dict:
                line = line.replace("@@{}@@".format(key), "{}".format(param_dict[key]))
            in_file.write(line)

    # Replace placeholders in probin template file
    with open("{}/probin".format(odir), "w") as prob_file:
        for line in probin_lines:
            for key in param_dict:
                line = line.replace("@@{}@@".format(key), "{}".format(param_dict[key]))
            prob_file.write(line)

    with open("{}/convective_flame.hse.tanh.delta_0.040cm.dx_0.050cm".format(odir), "w") as mod_file:
        for line in mod_lines:
            mod_file.write(line)

    shutil.copy(exe_file, odir)
    print("exe copied ",odir)


# read in the templates and model file
with open(inputs_file, "r") as in_file:
    inputs_lines = in_file.readlines()

with open(probin_file, "r") as prob_file:
    probin_lines = prob_file.readlines()

with open(model_file, "r") as mod_file:
    mod_lines = mod_file.readlines()


# values to be replaced in the template files.
# the format is (low, reference, high)
params = {"ROT_PERIOD": (0.001, 0.01, 0.1),
          "X_PERT_LOC": (0.5, 1.0, 10.0),
          "NU": (0.5, 1, 2),
          "Q_BURN": (1.24e9, 1.24e9, 1.24e9),
          "COND": (5.e8, 5.e9, 5.e10),
          "PERT_WIDTH": (0.0025, 0.025, 0.25),
          "T_BURN_REF": (0.1, 1.0, 10.0),
          "PERT_FACTOR": (0.1, 1.0, 10.0),
          "RTILDE": (4.0, 40.0, 20.0),
          "T_BASE": (0.01, 1.0, 10.0)}

# by definition, the "reference" run is the middle set of parameter
reference = {}
for key in params:
    reference[key] = params[key][1]

# setup reference directory
odir = "reference"
write_files(odir, reference)

# setup high and low directories
for key in params:

    out_dict = {}
    for k in reference:
        out_dict[k] = reference[k]

    if os.path.exists(key) == False:
        os.mkdir(key)

    low, _, high = params[key]

    # work on low
    odir = "{}/{}".format(key, "low")
    out_dict[key] = low
    write_files(odir, out_dict)

    # work on high
    odir = "{}/{}".format(key, "high")
    out_dict[key] = high
    write_files(odir, out_dict)
