#!/usr/bin/env python3

"""
This script parses the list of C++ runtime parameters and writes the
necessary header files and Fortran routines to make them available
in Castro's C++ routines and (optionally) the Fortran routines
through meth_params_module.

parameters have the format:

  name  type  default  need-in-fortran?  ifdef

the first three (name, type, default) are mandatory:

  name: the name of the parameter.  This will be the same name as the
    variable in C++ unless a pair is specified as (name, cpp_name)

  type: the C++ data type (int, bool, Real, string)

  default: the default value.  If specified as a pair, (a, b), then
    the first value is the normal default and the second is for
    debug mode (#ifdef AMREX_DEBUG)

the next are optional:

   need-in-fortran: if "y" then we do a pp.query() in meth_params_nd.F90

   ifdef: only define this parameter if the name provided is #ifdef-ed

Any line beginning with a "#" is ignored

Commands begin with a "@":

   @namespace: sets the namespace that these will be under (see below)

     e.g. @namespace castro

Note: categories listed in the input file aren't used for code generation
but are used for the documentation generation


For a namespace, name, we write out:

  -- name_params.H  (for castro, included in Castro.H):
     sets up the namespace and extern parameters

  -- name_declares.H  (for castro, included in Castro.cpp):
     declares the runtime parameters

  -- name_queries.H  (for castro, included in Castro.cpp):
     does the parmparse query to override the default in C++

  -- name_job_info_tests.H
     this tests the current value against the default and outputs
     into a file

we write out a single copy of:

  -- meth_params_nd.F90
     does the parmparse query to override the default in Fortran,
     and sets a number of other parameters specific to the F90 routines

"""

import argparse
import re
import sys

import runtime_parameters as rp

FWARNING = """
! This file is automatically created by parse_castro_params.py at build time.
! To update or add runtime parameters, please edit _cpp_parameters and rebuild.\n
"""

CWARNING = """
// This file is automatically created by parse_castro_params.py at build time.
// To update or add runtime parameters, please edit _cpp_parameters and rebuild.\n
"""

def write_meth_module(plist, meth_template, out_directory):
    """this writes the meth_params_module, starting with the meth_template
       and inserting the runtime parameter declaration in the correct
       place
    """

    try:
        mt = open(meth_template, "r")
    except IOError:
        sys.exit("invalid template file")

    try:
        mo = open(f"{out_directory}/meth_params_nd.F90", "w")
    except IOError:
        sys.exit("unable to open meth_params_nd.F90 for writing")

    mo.write(FWARNING)

    param_decls = [p.get_f90_decl_string() for p in plist if p.in_fortran == 1]
    params = [p for p in plist if p.in_fortran == 1]

    decls = ""

    for p in param_decls:
        decls += f"  {p}"

    for line in mt:
        if line.find("@@f90_declarations@@") > 0:
            mo.write(decls)

        elif line.find("@@set_castro_params@@") >= 0:

            namespaces = {q.namespace for q in params}
            print("Fortran namespaces: ", namespaces)
            for nm in namespaces:
                params_nm = [q for q in params if q.namespace == nm]
                ifdefs = {q.ifdef for q in params_nm}

                for ifdef in ifdefs:
                    if ifdef is None:
                        for p in [q for q in params_nm if q.ifdef is None]:
                            mo.write(p.get_f90_default_string())
                    else:
                        mo.write(f"#ifdef {ifdef}\n")
                        for p in [q for q in params_nm if q.ifdef == ifdef]:
                            mo.write(p.get_f90_default_string())
                        mo.write("#endif\n")

                mo.write("\n")

                mo.write(f'    call amrex_parmparse_build(pp, "{nm}")\n')

                for ifdef in ifdefs:
                    if ifdef is None:
                        for p in [q for q in params_nm if q.ifdef is None]:
                            mo.write(p.get_query_string("F90"))
                    else:
                        mo.write(f"#ifdef {ifdef}\n")
                        for p in [q for q in params_nm if q.ifdef == ifdef]:
                            mo.write(p.get_query_string("F90"))
                        mo.write("#endif\n")

                mo.write('    call amrex_parmparse_destroy(pp)\n')

                mo.write("\n\n")

        elif line.find("@@free_castro_params@@") >= 0:

            params_free = [q for q in params if q.in_fortran == 1]

            for p in params_free:
                if p.dtype != "string":
                    mo.write(f"    if (allocated({p.name})) then\n")
                    mo.write(f"        deallocate({p.name})\n")
                    mo.write("    end if\n")

            mo.write("\n\n")


        else:
            mo.write(line)

    mo.close()
    mt.close()


def parse_params(infile, meth_template, out_directory):

    params = []

    namespace = None

    try:
        f = open(infile)
    except IOError:
        sys.exit("error opening the input file")

    for line in f:
        if line[0] == "#":
            continue

        if line.strip() == "":
            continue

        if line[0] == "@":
            # this is a command
            cmd, value = line.split(":")
            if cmd == "@namespace":
                fields = value.split()
                namespace = fields[0]

            else:
                sys.exit("invalid command")

            continue

        # this splits the line into separate fields.  A field is a
        # single word or a pair in parentheses like "(a, b)"
        fields = re.findall(r'[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)', line)

        name = fields[0]
        if name[0] == "(":
            name, cpp_var_name = re.findall(r"\w+", name)
        else:
            cpp_var_name = name

        dtype = fields[1].lower()

        default = fields[2]
        if default[0] == "(":
            default, debug_default = re.findall(r"\w+", default)
        else:
            debug_default = None

        try:
            in_fortran_string = fields[3]
        except IndexError:
            in_fortran = 0
        else:
            if in_fortran_string.lower().strip() == "y":
                in_fortran = 1
            else:
                in_fortran = 0

        try:
            ifdef = fields[4]
        except IndexError:
            ifdef = None

        if namespace is None:
            sys.exit("namespace not set")

        params.append(rp.Param(name, dtype, default,
                               cpp_var_name=cpp_var_name,
                               namespace=namespace,
                               debug_default=debug_default,
                               in_fortran=in_fortran,
                               ifdef=ifdef))



    # output

    # find all the namespaces
    namespaces = {q.namespace for q in params}

    for nm in namespaces:

        params_nm = [q for q in params if q.namespace == nm]
        ifdefs = {q.ifdef for q in params_nm}

        # write name_declares.H
        try:
            cd = open(f"{out_directory}/{nm}_declares.H", "w")
        except IOError:
            sys.exit(f"unable to open {nm}_declares.H for writing")

        cd.write(CWARNING)
        cd.write(f"#ifndef _{nm.upper()}_DECLARES_H_\n")
        cd.write(f"#define _{nm.upper()}_DECLARES_H_\n")

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    cd.write(p.get_declare_string())
            else:
                cd.write(f"#ifdef {ifdef}\n")
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    cd.write(p.get_declare_string())
                cd.write("#endif\n")

        cd.write("#endif\n")
        cd.close()

        # write name_params.H
        try:
            cp = open(f"{out_directory}/{nm}_params.H", "w")
        except IOError:
            sys.exit(f"unable to open {nm}_params.H for writing")

        cp.write(CWARNING)
        cp.write(f"#ifndef _{nm.upper()}_PARAMS_H_\n")
        cp.write(f"#define _{nm.upper()}_PARAMS_H_\n")

        cp.write("\n")
        cp.write(f"namespace {nm} {{\n")

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    cp.write(p.get_decl_string())
            else:
                cp.write(f"#ifdef {ifdef}\n")
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    cp.write(p.get_decl_string())
                cp.write("#endif\n")
        cp.write("}\n\n")
        cp.write("#endif\n")
        cp.close()

        # write castro_queries.H
        try:
            cq = open(f"{out_directory}/{nm}_queries.H", "w")
        except IOError:
            sys.exit(f"unable to open {nm}_queries.H for writing")

        cq.write(CWARNING)

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    cq.write(p.get_default_string())
                    cq.write(p.get_query_string("C++"))
                    cq.write("\n")
            else:
                cq.write(f"#ifdef {ifdef}\n")
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    cq.write(p.get_default_string())
                    cq.write(p.get_query_string("C++"))
                    cq.write("\n")
                cq.write("#endif\n")
            cq.write("\n")
        cq.close()

        # write the job info tests
        try:
            jo = open(f"{out_directory}/{nm}_job_info_tests.H", "w")
        except IOError:
            sys.exit(f"unable to open {nm}_job_info_tests.H")

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    jo.write(p.get_job_info_test())
            else:
                jo.write(f"#ifdef {ifdef}\n")
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    jo.write(p.get_job_info_test())
                jo.write("#endif\n")

        jo.close()

    # write the Fortran module
    write_meth_module(params, meth_template, out_directory)


def main():
    """the main driver"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", type=str, default=None,
                        help="template for the meth_params module")
    parser.add_argument("-o", type=str, default=None,
                        help="output directory for the generated files")
    parser.add_argument("input_file", type=str, nargs=1,
                        help="input file containing the list of parameters we will define")

    args = parser.parse_args()

    parse_params(args.input_file[0], args.m, args.o)

if __name__ == "__main__":
    main()
