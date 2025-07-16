"""Helper functions for working with Castro diagnostic files (*_diag.out)

To use these in a standalone script, you can do one of the following:

* append $CASTRO_HOME/Util/scripts to sys.path at the top of your script:
    sys.path.append("<path to Castro>/Util/scripts")

* add a symlink to this file in the same directory as your script:
    $ ln -s "$CASTRO_HOME/Util/scripts/diag_parser.py" .

* copy this file into the same directory as your script

Then you can do `from diag_parser import deduplicate, read_diag_file`.
"""

from pathlib import Path

import numpy as np

""" Format notes
files are opened in Castro.cpp, data is written in sum_integrated_quantities.cpp

data_logs[0]: grid_diag.out
intwidth, fixwidth, datwidth*

data_logs[1]: gravity_diag.out
- this was previously missing the last column number (8), which we handle for
  backwards compatibility
intwidth, fixwidth, datwidth, datwidth, datwidth, datwidth, datwidth, datwidth

data_logs[2]: species_diag.out
intwidth, fixwidth, datwidth*

data_logs[3]: amr_diag.out
- if compiled with GPU support, this will have two additional integer fields at
  the end with size `datwidth` for the GPU memory usage
- column 5 (max number of subcycles) is an integer
intwidth, fixwidth, fixwidth, intwidth, fixwidth, datwidth
"""

datwidth = 25  # Floating point data in scientific notation
fixwidth = 25  # Floating point data not in scientific notation
intwidth = 12  # Integer data

# Any additional columns after these are assumed to be floating point values in
# scientific notation (amr_diag.out gets special handling)
FIELD_WIDTHS = {
    "grid_diag.out": [intwidth, fixwidth],
    "gravity_diag.out": [intwidth, fixwidth] + [datwidth] * 6,
    "species_diag.out": [intwidth, fixwidth],
    "amr_diag.out": [intwidth, fixwidth, fixwidth, intwidth, fixwidth, datwidth],
}


def read_diag_file(file_path):
    """Reads a Castro diagnostic file into a numpy structured array.

    Currently only supports the default files that Castro generates.
    """
    if not isinstance(file_path, Path):
        file_path = Path(file_path)
    filetype = file_path.name
    if filetype not in FIELD_WIDTHS:
        raise ValueError("Unsupported file name")
    widths = FIELD_WIDTHS[filetype]
    with open(file_path) as f:
        # try getting the number of columns from the first line
        first_line = f.readline().rstrip("\n")
        if filetype == "gravity_diag.out":
            # gravity_diag.out is missing the last column number, but it
            # fortunately has a fixed number of columns
            num_columns = 8
        else:
            num_columns = int(first_line.split()[-1])
        # pad out the widths list on the right if necessary
        widths.extend([datwidth] * (num_columns - len(widths)))
        # infer datatypes from the widths
        dtypes = [int if w == intwidth else float for w in widths]
        # amr_diag.out has several integer columns with long names
        if filetype == "amr_diag.out":
            dtypes[4] = int  # max number of subcycles
            if num_columns >= 8:
                dtypes[6] = int  # maximum gpu memory used
                dtypes[7] = int  # minimum gpu memory free
        # already read the first header line, so we don't need to skip any rows
        data = np.genfromtxt(
            f, delimiter=widths, comments="#", dtype=dtypes, names=True
        )
    return data


def deduplicate(data):
    """Deduplicate based on the timestep, keeping the only last occurrence."""
    # get the unique indices into the reversed timestep array, so we find the
    # final occurrence of each timestep
    _, rev_indices = np.unique(data["TIMESTEP"][::-1], return_index=True)
    # np.unique() sorts by value, so we don't need to un-reverse rev_indices
    unique_indices = data.shape[0] - rev_indices - 1
    return data[unique_indices]
