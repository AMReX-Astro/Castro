#!/usr/bin/env python


# This script parses the list of C++ runtime parameters and writes the
# necessary header files and Fortran routines to make them available
# in Castro's C++ routines and (optionally) the Fortran routines
# through meth_params_module.
#
# specifically, we write out:
#
#   -- castro_params.H  (included in Castro.H):
#      declares the static variables of the Castro class
#
#   -- castro_defaults.H  (included in Castro.cpp):
#      sets the defaults of the runtime parameters
#
#   -- castro_queries.H  (included in Castro.cpp):
#      does the parmparse query to override the default in C++
#
#   -- meth_params.F90
#      does the parmparse query to override the default in Fortran,
#      and sets a number of other parameters specific to the F90 routinse

import argparse
import re
import sys

class Param(object):
    def __init__(self, name, dtype, default,
                 debug_default=None,
                 in_fortran=0, f90_name=None, f90_dtype=None,
                 ifdef=None):

        self.name = name
        self.dtype = dtype
        self.default = default
        self.debug_default=debug_default
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
            tstr = "int         Castro::{}".format(self.name)
        elif self.dtype == "Real":
            tstr = "Real        Castro::{}".format(self.name)
        elif self.dtype == "string":
            tstr = "std::string Castro::{}".format(self.name)
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

        if not self.debug_default is None:
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
            ostr += "pp.query(\"{}\", {});\n".format(self.name, self.name)
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

        if self.dtype == "int":
            tstr = "static int {};\n".format(self.name)
        elif self.dtype == "Real":
            tstr = "static Real {};\n".format(self.name)
        elif self.dtype == "string":
            tstr = "static std::string {};\n".format(self.name)
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
            tstr = "double precision, save :: {}\n".format(self.f90_name)
        elif self.f90_dtype == "logical":
            tstr = "logical         , save :: {}\n".format(self.f90_name)
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

    try: mo = open("Src_nd/meth_params.F90", "w")
    except:
        sys.exit("unable to open meth_params.F90 for writing")

    param_decls = [p.get_f90_decl_string() for p in plist if p.in_fortran == 1]
    params = [p for p in plist if p.in_fortran == 1]

    decls = ""

    for p in param_decls:
        decls += "  {}".format(p)

    for line in mt:
        if line.find("@@f90_declarations@@") > 0:
            mo.write(decls)

        elif line.find("@@set_castro_params@@") >= 0:
            for p in params:
                mo.write(p.get_f90_default_string())

            mo.write("\n")

            for p in params:
                mo.write(p.get_query_string("F90"))

        else:
            mo.write(line)

    mo.close()
    mt.close()


def parse_params(infile, meth_template):

    params = []

    try: f = open(infile)
    except:
        sys.exit("error openning the input file")


    for line in f:
        if line[0] == "#":
            continue

        if line.strip() == "":
            continue

        fields = re.findall(r'[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)', line)

        name = fields[0]
        dtype = fields[1]

        default = fields[2]
        if default[0] == "(":
            default, debug_default = re.findall(r'\w+', default)
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

        params.append(Param(name, dtype, default, debug_default=debug_default,
                            in_fortran=in_fortran, f90_name=f90_name, f90_dtype=f90_dtype,
                            ifdef=ifdef))



    # output

    # write castro_defaults.H
    try: cd = open("castro_defaults.H", "w")
    except:
        sys.exit("unable to open castro_defaults.H for writing")

    for p in params:
        cd.write(p.get_default_string())

    cd.close()

    # write castro_params.H
    try: cp = open("castro_params.H", "w")
    except:
        sys.exit("unable to open castro_params.H for writing")

    for p in params:
        cp.write(p.get_decl_string())

    cp.close()

    # write castro_queries.H
    try: cq = open("castro_queries.H", "w")
    except:
        sys.exit("unable to open castro_queries.H for writing")

    for p in params:
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
