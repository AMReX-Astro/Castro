import os

inputs_file = "inputs.2d.template"
probin_file = "probin.template"
model_file = "convective_flame.hse.tanh.delta_0.040cm.dx_0.050cm"

def write_files(odir, param_dict):

    os.mkdir(odir)
    with open("{}/inputs.2d".format(odir), "w") as in_file:
        for line in inputs_lines:
            for key in param_dict:
                line = line.replace("@@{}@@".format(key), "{}".format(param_dict[key]))
            in_file.write(line)

    with open("{}/probin".format(odir), "w") as prob_file:
        for line in probin_lines:
            for key in param_dict:
                line = line.replace("@@{}@@".format(key), "{}".format(param_dict[key]))
            prob_file.write(line)

# read in the templates
with open(inputs_file, "r") as in_file:
    inputs_lines = in_file.readlines()

with open(probin_file, "r") as prob_file:
    probin_lines = prob_file.readlines()

# parameters to loop over
params = {"ROT_PERIOD": (0.001, 0.01, 0.1),
          "X_PERT_LOC": (0.5, 1.0, 10.0),
          "NU": (1, 4, 16),
          "Q_BURN": (1.24e7, 1.24e8, 1.24e9),
          "COND": (5.e8, 5.e9, 5.e10)}

# by definition, the "reference" run is the middle set of parameter
reference = {}
for key in params:
    reference[key] = params[key][1]

# setup reference
odir = "reference"
write_files(odir, reference)

# create the remainder
for key in params:

    out_dict = {}
    for k in reference:
        out_dict[k] = reference[k]

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



