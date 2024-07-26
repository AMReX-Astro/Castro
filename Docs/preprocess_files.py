"""
This script strips out preprocessor directives from C++ headers
and saves the results in source/preprocessed_files
"""

import os
import re

# directory of the source files
rootdir = "../Source"

outdir = "source/preprocessed_files"


def strip_directives(filename, filepath, outpath):
    """
    Read in file, remove all preprocessor directives and output.

    This is also going to switch square brackets initializing arrays to
    parentheses and remove the new-line characters
    """

    with open(os.path.join(filepath, filename)) as infile:
        txt = infile.read()

        outtxt = re.sub(r"(^#.*$\n)", '', txt, flags=re.M)
        outtxt = re.sub(r"(&\n)\s*", '', outtxt)
        outtxt = re.sub(r"\[", r"(\\", outtxt)
        outtxt = re.sub(r"\]", r'\\)', outtxt)

        with open(os.path.join(outpath, filename), 'w') as outfile:
            outfile.write(outtxt)

if __name__ == "__main__":

    # make the output directory if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # loop over source dir
    for subdir in sorted(os.listdir(rootdir)):
        if not os.path.isdir(os.path.join(rootdir, subdir)):
            continue

        # loop over files in subdirectories and run strip_directives on all
        # C++ header files
        for f in sorted(os.listdir(os.path.join(rootdir, subdir))):
            if (f[-2:] == ".H"):
                strip_directives(f, os.path.join(rootdir, subdir), outdir)
