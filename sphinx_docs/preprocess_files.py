"""
This script strips out preprocessor directives from C++ headers and Fortran files
and saves the results in source/preprocessed_files
"""

import os
import re

# directory of the source files
rootdir = "../Source"

outdir = "source/preprocessed_files"


def strip_directives(filename, filepath, outpath):
    """
    Read in file, remove all preprocessor directives and output
    """

    # r = re.compile(r"(^#.*$\n)")

    with open(os.path.join(filepath, filename)) as infile:
        txt = infile.read()

        outtxt = re.sub(r"(^#.*$\n)", '', txt, flags=re.M)

        with open(os.path.join(outpath, filename), 'w') as outfile:
            outfile.write(outtxt)


def delete_lines(filename, filepath):
    """
    For some reason sphinx-fortran does not like these lines, so we're going to
    delete them out by hand in the preprocessed files
    """

    if filename == 'sdc_util.F90':
        txt = ""
        with open(os.path.join(filepath, filename)) as infile:
            txt = infile.read()

        txt = txt.replace("integer, parameter :: lrw = 22 + 9*(nspec_evolve+2) + 2*(nspec_evolve+2)**2",
                          "")

        with open(os.path.join(filepath, filename), 'w') as outfile:
            outfile.write(txt)


if __name__ == "__main__":

    excl_files = ['HABEC_1D.F90', 'HABEC_2D.F90', 'HABEC_3D.F90', 'RAD_1D.F90', 'RAD_2D.F90', 'RAD_3D.F90', 'bc_fill_3d.F90',
                  'problem_tagging_3d.f90',
                  'bc_ext_fill_3d.F90',
                  'Prob_1d.f90',
                  'Prob_3d.f90',
                  'bc_ext_fill_1d.F90',
                  'bc_fill_1d.F90',
                  'bc_ext_fill_2d.F90',
                  'problem_tagging_1d.f90',
                  'bc_fill_2d.F90',
                  'Prob_2d.f90',
                  'problem_tagging_2d.f90']

    # make the output directory if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # loop over source dir
    for subdir in sorted(os.listdir(rootdir)):
        if not os.path.isdir(os.path.join(rootdir, subdir)):
            continue

        # loop over files in subdirectories and run strip_directives on all
        # C++ header files and Fortran files
        for f in sorted(os.listdir(os.path.join(rootdir, subdir))):
            if ((f[-2:] == ".H" and f[-4:] != "_F.H") or f[-4:] == ".F90" or f[-4:] == ".f90") and f not in excl_files:
                strip_directives(f, os.path.join(rootdir, subdir), outdir)
                delete_lines(f, outdir)
