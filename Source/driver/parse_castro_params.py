#!/usr/bin/env python

from __future__ import print_function

# This script parses the list of C++ runtime parameters and writes the
# necessary header files and Fortran routines to make them available
# in Castro's C++ routines and (optionally) the Fortran routines
# through meth_params_module.
#
# parameters have the format:
#
#   name  type  default  need-in-fortran?  ifdef fortran-name  fortran-type
#
# the first three (name, type, default) are mandatory:
#
#   name: the name of the parameter.  This will be the same name as the
#     variable in C++ unless a pair is specified as (name, cpp_name)
#
#   type: the C++ data type (int, Real, string)
#
#   default: the default value.  If specified as a pair, (a, b), then
#     the first value is the normal default and the second is for
#     debug mode (#ifdef DEBUG)
#
# the next are optional:
#
#    need-in-fortran: if "y" then we do a pp.query() in meth_params.F90
#
#    ifdef: only define this parameter if the name provided is #ifdef-ed
#
#    fortran-name: if a different variable name in Fortran, specify here
#
#    fortran-type: if a different data type in Fortran, specify here
#
# Any line beginning with a "#" is ignored
#
# Commands begin with a "@":
#
#    @namespace: sets the namespace that these will be under (see below)
#      it also gives the C++ class name.
#      if we include the keyword "static" after the name, then the parameters
#      will be defined as static member variables in C++
#
#      e.g. @namespace castro Castro static
#
# Note: categories listed in the input file aren't used for code generation
# but are used for the documentation generation
#
#
# For a namespace, name, we write out:
#
#   -- name_params.H  (for castro, included in Castro.H):
#      declares the static variables of the Castro class
#
#   -- name_defaults.H  (for castro, included in Castro.cpp):
#      sets the defaults of the runtime parameters
#
#   -- name_queries.H  (for castro, included in Castro.cpp):
#      does the parmparse query to override the default in C++
#
# we write out a single copy of:
#
#   -- meth_params.F90
#      does the parmparse query to override the default in Fortran,
#      and sets a number of other parameters specific to the F90 routinse
#

import argparse
import re
import sys

FWARNING = """
! This file is automatically created by parse_castro_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh\n
"""

CWARNING = """
// This file is automatically created by parse_castro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh\n
"""

param_include_dir = "param_includes/"


class Param(object):
    """ the basic parameter class.  For each parameter, we hold the name,
        type, and default.  For some parameters, we also take a second
        value of the default, for use in debug mode (delimited via
        #ifdef DEBUG)

    """

    def __init__(self, name, dtype, default,
                 cpp_var_name=None,
                 namespace=None, cpp_class=None, static=None,
                 debug_default=None,
                 in_fortran=0, f90_name=None, f90_dtype=None,
                 ifdef=None):

        self.name = name
        self.dtype = dtype
        self.default = default
        self.cpp_var_name = cpp_var_name

        self.namespace = namespace
        self.cpp_class = cpp_class

        if static is None:
            self.static = 0
        else:
            self.static = static

        self.debug_default = debug_default
        self.in_fortran = in_fortran

        if ifdef == "None":
            self.ifdef = None
        else:
            self.ifdef = ifdef

        if f90_name is None:
            self.f90_name = name
        else:
            self.f90_name = f90_name

        if f90_dtype is None:
            self.f90_dtype = dtype
        else:
            self.f90_dtype = f90_dtype

    def get_default_string(self):
        # this is the line that goes into castro_defaults.H included
        # into Castro.cpp

        if self.dtype == "int":
            tstr = "int         {}::{}".format(self.cpp_class, self.cpp_var_name)
        elif self.dtype == "Real":
            tstr = "amrex::Real {}::{}".format(self.cpp_class, self.cpp_var_name)
        elif self.dtype == "string":
            tstr = "std::string {}::{}".format(self.cpp_class, self.cpp_var_name)
        else:
            sys.exit("invalid data type for parameter {}".format(self.name))

        ostr = ""

        if not self.ifdef is None:
            ostr = "#ifdef {}\n".format(self.ifdef)

        if not self.debug_default is None:
            ostr += "#ifdef DEBUG\n"
            ostr += "{} = {};\n".format(tstr, self.debug_default)
            ostr += "#else\n"
            ostr += "{} = {};\n".format(tstr, self.default)
            ostr += "#endif\n"
        else:
            ostr += "{} = {};\n".format(tstr, self.default)

        if not self.ifdef is None:
            ostr += "#endif\n"

        return ostr

    def get_f90_default_string(self):
        # this is the line that goes into set_castro_method_params()
        # to set the default value of the variable

        ostr = ""

        # convert to the double precision notation Fortran knows
        # if the parameter is already of the form "#.e###" then
        # it is easy as swapping out "e" for "d"; if it is a number
        # like 0.1 without a format specifier, then add a d0 to it
        # because the C++ will read it in that way and we want to
        # give identical results (at least to within roundoff)

        if self.debug_default is not None:
            debug_default = self.debug_default
            if self.dtype == "Real":
                if "e" in debug_default:
                    debug_default = debug_default.replace("e", "d")
                else:
                    debug_default += "d0"

        default = self.default
        if self.dtype == "Real":
            if "e" in default:
                default = default.replace("e", "d")
            else:
                default += "d0"

        name = self.f90_name

        # for a character, we need to allocate its length.  We allocate
        # to 1, and the Fortran parmparse will resize
        if self.dtype == "string":
            ostr += "    allocate(character(len=1)::{})\n".format(name)

        if not self.debug_default is None:
            ostr += "#ifdef DEBUG\n"
            ostr += "    {} = {};\n".format(name, debug_default)
            ostr += "#else\n"
            ostr += "    {} = {};\n".format(name, default)
            ostr += "#endif\n"
        else:
            ostr += "    {} = {};\n".format(name, default)

        return ostr

    def get_query_string(self, language):
        # this is the line that queries the ParmParse object to get
        # the value of the runtime parameter from the inputs file.
        # This goes into castro_queries.H included into Castro.cpp

        ostr = ""
        if not self.ifdef is None:
            ostr += "#ifdef {}\n".format(self.ifdef)

        if language == "C++":
            ostr += "pp.query(\"{}\", {});\n".format(self.name, self.cpp_var_name)
        elif language == "F90":
            ostr += "    call pp%query(\"{}\", {})\n".format(self.name, self.f90_name)
        else:
            sys.exit("invalid language choice in get_query_string")

        if not self.ifdef is None:
            ostr += "#endif\n".format(self.ifdef)

        return ostr

    def get_decl_string(self):
        # this is the line that goes into castro_params.H included
        # into Castro.H

        static = ""
        if self.static:
            static = "static"

        if self.dtype == "int":
            tstr = "{} int {};\n".format(static, self.cpp_var_name)
        elif self.dtype == "Real":
            tstr = "{} amrex::Real {};\n".format(static, self.cpp_var_name)
        elif self.dtype == "string":
            tstr = "{} std::string {};\n".format(static, self.cpp_var_name)
        else:
            sys.exit("invalid data type for parameter {}".format(self.name))

        ostr = ""

        if not self.ifdef is None:
            ostr = "#ifdef {}\n".format(self.ifdef)

        ostr += tstr

        if not self.ifdef is None:
            ostr += "#endif\n"

        return ostr

    def get_f90_decl_string(self):
        # this is the line that goes into meth_params.f90

        if not self.in_fortran:
            return None

        if self.f90_dtype == "int":
            tstr = "integer         , save :: {}\n".format(self.f90_name)
        elif self.f90_dtype == "Real":
            tstr = "real(rt), save :: {}\n".format(self.f90_name)
        elif self.f90_dtype == "logical":
            tstr = "logical         , save :: {}\n".format(self.f90_name)
        elif self.f90_dtype == "string":
            tstr = "character (len=:), allocatable, save :: {}\n".format(self.f90_name)
        else:
            sys.exit("unsupported datatype for Fortran: {}".format(self.name))

        return tstr


def write_meth_module(plist, meth_template):
    """this writes the meth_params_module, starting with the meth_template
       and inserting the runtime parameter declaration in the correct
       place
    """

    try: mt = open(meth_template, "r")
    except:
        sys.exit("invalid template file")

    try: mo = open("meth_params.F90", "w")
    except:
        sys.exit("unable to open meth_params.F90 for writing")


    mo.write(FWARNING)

    param_decls = [p.get_f90_decl_string() for p in plist if p.in_fortran == 1]
    params = [p for p in plist if p.in_fortran == 1]

    decls = ""

    for p in param_decls:
        decls += "  {}".format(p)

    for line in mt:
        if line.find("@@f90_declarations@@") > 0:
            mo.write(decls)

            # Now do the OpenACC declarations

            mo.write("\n")
            mo.write("  !$acc declare &\n")
            mo.write("  !$acc create(")

            for n, p in enumerate(params):
                if p.f90_dtype == "string": 
                    print("warning: string parameter {} will not be available on the GPU".format(p.name),
                          file=sys.stderr)
                    continue

                mo.write("{}".format(p.f90_name))

                if n == len(params)-1:
                    mo.write(")\n")
                else:
                    if n % 3 == 2:
                        mo.write(") &\n  !$acc create(")
                    else:
                        mo.write(", ")

        elif line.find("@@set_castro_params@@") >= 0:

            namespaces = list(set([q.namespace for q in params]))
            print("namespaces: ", namespaces)
            for nm in namespaces:
                params_nm = [q for q in params if q.namespace == nm]

                for p in params_nm:
                    mo.write(p.get_f90_default_string())

                mo.write("\n")

                mo.write('    call amrex_parmparse_build(pp, "{}")\n'.format(nm))

                for p in params_nm:
                    mo.write(p.get_query_string("F90"))

                mo.write('    call amrex_parmparse_destroy(pp)\n')
                
                mo.write("\n\n")

            # Now do the OpenACC device updates

            mo.write("\n")
            mo.write("    !$acc update &\n")
            mo.write("    !$acc device(")

            for n, p in enumerate(params):
                if p.f90_dtype == "string": continue
                mo.write("{}".format(p.f90_name))

                if n == len(params)-1:
                    mo.write(")\n")
                else:
                    if n % 3 == 2:
                        mo.write(") &\n    !$acc device(")
                    else:
                        mo.write(", ")

        elif line.find("@@free_castro_params@@") >= 0:

            params_free = [q for q in params if q.in_fortran == 1 and q.f90_dtype == "string"]

            for p in params_free:
                mo.write("    deallocate({})\n".format(p.f90_name))
                
            mo.write("\n\n")


        else:
            mo.write(line)

    mo.close()
    mt.close()


def parse_params(infile, meth_template):

    params = []

    namespace = None
    cpp_class = None
    static = None

    try: f = open(infile)
    except:
        sys.exit("error openning the input file")


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
                cpp_class = fields[1]

                try: static = fields[2]
                except: static = ""

                # do we have the static keyword?
                if "static" in static:
                    static = 1
                else:
                    static = 0

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

        dtype = fields[1]

        default = fields[2]
        if default[0] == "(":
            default, debug_default = re.findall(r"\w+", default)
        else:
            debug_default = None

        try: in_fortran_string = fields[3]
        except: in_fortran = 0
        else:
            if in_fortran_string.lower().strip() == "y":
                in_fortran = 1
            else:
                in_fortran = 0

        try: ifdef = fields[4]
        except: ifdef = None

        try: f90_name = fields[5]
        except: f90_name = None

        try: f90_dtype = fields[6]
        except: f90_dtype = None

        if namespace is None:
            sys.exit("namespace not set")

        params.append(Param(name, dtype, default,
                            cpp_var_name=cpp_var_name,
                            namespace=namespace,
                            cpp_class=cpp_class,
                            static=static,
                            debug_default=debug_default,
                            in_fortran=in_fortran, f90_name=f90_name, f90_dtype=f90_dtype,
                            ifdef=ifdef))



    # output

    # find all the namespaces
    namespaces = list(set([q.namespace for q in params]))

    for nm in namespaces:

        params_nm = [q for q in params if q.namespace == nm]

        # write name_defaults.H
        try: cd = open("{}/{}_defaults.H".format(param_include_dir, nm), "w")
        except:
            sys.exit("unable to open {}_defaults.H for writing".format(nm))

        cd.write(CWARNING)

        for p in params_nm:
            cd.write(p.get_default_string())

        cd.close()

        # write name_params.H
        try: cp = open("{}/{}_params.H".format(param_include_dir, nm), "w")
        except:
            sys.exit("unable to open {}_params.H for writing".format(nm))

        cp.write(CWARNING)

        for p in params_nm:
            cp.write(p.get_decl_string())

        cp.close()

        # write castro_queries.H
        try: cq = open("{}/{}_queries.H".format(param_include_dir, nm), "w")
        except:
            sys.exit("unable to open {}_queries.H for writing".format(nm))

        cq.write(CWARNING)

        for p in params_nm:
            cq.write(p.get_query_string("C++"))

        cq.close()


    # write the Fortran module
    write_meth_module(params, meth_template)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", type=str, default=None,
                        help="template for the meth_params module")
    parser.add_argument("input_file", type=str, nargs=1,
                        help="input file containing the list of parameters we will define")

    args = parser.parse_args()

    parse_params(args.input_file[0], args.m)
