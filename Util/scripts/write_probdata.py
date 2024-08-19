#!/usr/bin/env python3

"""This routine parses plain-text parameter files that list runtime
parameters for use in our codes.  The general format of a parameter
is:

max_step                            integer            1
small_dt                            real               1.d-10
xlo_boundary_type                   character          ""
octant                              logical            .false.

This specifies the parameter name, datatype, and default
value.

The optional 4th column indicates whether the parameter appears
in the fortin namelist ("y" or "Y").

If the parameter is an array, the optional 5th column indicates the size of
the array.

Note: there are two types of parameters here, the ones that are in the
namelist are true runtime parameters.  The non-namelist parameters
should be avoided if at all possible.

"""

import argparse
import os
import re
import sys

import runtime_parameters as rp

CXX_HEADER = """
#ifndef PROBLEM_PARAMETERS_H
#define PROBLEM_PARAMETERS_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <network_properties.H>

"""

CXX_FOOTER = """
#endif
"""


def get_next_line(fin):
    """return the next, non-blank line, with comments stripped"""
    line = fin.readline()

    pos = line.find("#")

    while (pos == 0 or line.strip() == "") and line:
        line = fin.readline()
        pos = line.find("#")

    if pos == -1:
        return line.strip()
    return line[:pos]


def parse_param_file(params_list, param_file):
    """read all the parameters in the prob_param_files and add valid
    parameters to the params list.  This returns the parameter list.

    """

    namespace = "problem"

    try:
        f = open(param_file)
    except FileNotFoundError:
        sys.exit(f"write_probdata.py: ERROR: file {param_file} does not exist")

    line = get_next_line(f)

    err = 0

    while line and not err:

        # this splits the line into separate fields.  A field is a
        # single word or a pair in parentheses like "(a, b)"
        fields = re.findall(r'[\w\"\+\./\-]+|\([\w+\./\-]+\s*,\s*[\w\+\.\-]+\)', line)

        if len(fields) < 3:
            print("write_probdata.py: ERROR: missing one or more fields in parameter definition.")
            err = 1
            continue

        name = fields[0]
        dtype = fields[1]
        default = fields[2]

        current_param = rp.Param(name, dtype, default,
                                 namespace=namespace)

        # optional field: in namelist
        try:
            in_namelist_in = fields[3]
            if in_namelist_in in ["y", "Y"]:
                in_namelist = True
            else:
                in_namelist = False

        except IndexError:
            in_namelist = False


        # optional field: size
        try:
            size = fields[4]
        except IndexError:
            size = 1

        current_param.in_namelist = in_namelist
        current_param.size = size

        # check to see if this parameter is defined in the current
        # list if we delete the old one and take the new one (we
        # assume that later files automatically have higher
        # priority)
        p_names = [p.name for p in params_list]
        try:
            idx = p_names.index(current_param.name)
        except ValueError:
            pass
        else:
            params_list.pop(idx)

        if not err == 1:
            params_list.append(current_param)

        line = get_next_line(f)

    return err


def abort(outfile):
    """ abort exits when there is an error.  A dummy stub file is written
    out, which will cause a compilation failure """

    fout = open(outfile, "w")
    fout.write("There was an error parsing the parameter files")
    fout.close()
    sys.exit(1)


def write_probin(prob_param_files, cxx_prefix):
    """write_probin will read through the list of parameter files and
    C++ header and source files to manage the runtime parameters"""

    params = []

    print(" ")
    print(f"write_probdata.py: creating prob_param C++ files")

    # read the parameters defined in the parameter files

    for f in prob_param_files:
        err = parse_param_file(params, f)
        if err:
            abort(f"Error parsing {f}")


    # now handle the C++ -- we need to write a header and a .cpp file
    # for the parameters

    cxx_base = os.path.basename(cxx_prefix)

    ofile = f"{cxx_prefix}_parameters.H"
    with open(ofile, "w") as fout:
        fout.write(CXX_HEADER)

        fout.write(f"  void init_{cxx_base}_parameters();\n\n")

        fout.write("  namespace problem {\n\n")

        for p in params:
            if p.dtype == "string":
                fout.write(f"  extern std::string {p.name};\n\n")
            else:
                if p.is_array():
                    if p.size == "nspec":
                        fout.write(f"  extern AMREX_GPU_MANAGED {p.get_cxx_decl()} {p.name}[NumSpec];\n\n")
                    else:
                        fout.write(f"  extern AMREX_GPU_MANAGED {p.get_cxx_decl()} {p.name}[{p.size}];\n\n")
                else:
                    fout.write(f"  extern AMREX_GPU_MANAGED {p.get_cxx_decl()} {p.name};\n\n")

        fout.write("  }\n\n")

        fout.write(CXX_FOOTER)

    # now the C++ job_info tests
    ofile = f"{cxx_prefix}_job_info_tests.H"
    with open(ofile, "w") as fout:
        for p in params:
            if not p.is_array():
                if p.in_namelist:
                    fout.write(p.get_job_info_test())

    # now the C++ initialization routines
    ofile = f"{cxx_prefix}_parameters.cpp"
    with open(ofile, "w") as fout:
        fout.write(f"#include <{cxx_base}_parameters.H>\n")
        fout.write("#include <AMReX_ParmParse.H>\n")
        fout.write("#include <AMReX_REAL.H>\n\n")
        for p in params:
            if p.dtype == "string":
                fout.write(f"  std::string problem::{p.name};\n\n")
            else:
                if p.is_array():
                    if p.size == "nspec":
                        fout.write(f"  AMREX_GPU_MANAGED {p.get_cxx_decl()} problem::{p.name}[NumSpec];\n\n")
                    else:
                        fout.write(f"  AMREX_GPU_MANAGED {p.get_cxx_decl()} problem::{p.name}[{p.size}];\n\n")
                else:
                    fout.write(f"  AMREX_GPU_MANAGED {p.get_cxx_decl()} problem::{p.name};\n\n")

        fout.write("\n")
        fout.write(f"  void init_{cxx_base}_parameters() {{\n")


        # now write the parmparse code to get the value from the C++
        # inputs.

        # open namespace
        fout.write("    {\n")

        # we need access to _rt
        fout.write("        using namespace amrex;\n")

        fout.write(f"        amrex::ParmParse pp(\"problem\");\n\n")
        for p in params:
            if p.is_array():
                size = p.size
                if (size == "nspec"):
                    size = "NumSpec"
                fout.write(f"        for (auto & e : problem::{p.name}) {{\n")
                fout.write(f"            e = {p.default_format(lang='C++')};\n")
                fout.write(f"        }}\n")
            else:
                fout.write(f"        {p.get_default_string()}")

            if p.in_namelist:
                fout.write(f"        {p.get_query_string()}")
            fout.write("\n")
        fout.write("    }\n")

        fout.write("  }\n")


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', type=str, help='problem parameter file names (space separated in quotes)')
    parser.add_argument('--cxx_prefix', type=str, default="prob",
                        help="a name to use in the C++ file names")

    args = parser.parse_args()

    prob_params = args.p.split()

    write_probin(prob_params, args.cxx_prefix)

if __name__ == "__main__":
    main()
