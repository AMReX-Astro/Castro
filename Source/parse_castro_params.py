#!/usr/bin/env python


# This script parses the list of C++ runtime parameters and writes the
# necessary header files and Fortran routines to make them available
# in Castro's C++ routines and (optionally) the Fortran routines
# through meth_params_module.
#
# specifically, we write out:
#
#   -- castro_set_meth.H  (included in Castro_F.H):
#      defines the prototype for set_castro_method_params
#
#   -- castro_call_set_meth.H  (included in Castro_setup.cpp):
#      implements the actual call to set_castro_method_params
#
#   -- castro_params.H  (included in Castro.H):
#      declares the static variables of the Castro class
#
#   -- castro_defaults.H  (included in Castro.cpp):
#      sets the defaults of the runtime parameters
#
#   -- castro_queries.H  (included in Castro.cpp):
#      does the parmparse query to override the default
#

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

    def get_query_string(self):
        # this is the line that queries the ParmParse object to get
        # the value of the runtime parameter from the inputs file.
        # This goes into castro_queries.H included into Castro.cpp
        ostr = ""
        if not self.ifdef is None:
            ostr += "#ifdef {}\n".format(self.ifdef)

        ostr += "pp.query(\"{}\", {});\n".format(self.name, self.name)

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

    def get_static_default_string(self):
        # this is used before the call to set_castro_meth_param to 
        # give defaults when the ifdef is not defined

        if self.dtype == "int":
            tstr = "static int {} = {};\n".format(self.name, self.default)
        elif self.dtype == "Real":
            tstr = "static Real {} = {};\n".format(self.name, self.default)
        elif self.dtype == "string":
            tstr = "static std::string {} = {};\n".format(self.name, self.default)
        else:
            sys.exit("invalid data type for parameter {}".format(self.name))

        return tstr

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

    def get_prototype_decl(self):
        # this give the line in the set_castro_method_params prototype
        # for this runtime parameter that will be written into
        # castro_set_meth.H and defined in Castro_F.H

        if not self.in_fortran:
            return None

        if self.dtype == "int":
            tstr = "const int& {}".format(self.name)
        elif self.dtype == "Real":
            tstr = "const Real& {}".format(self.name)
        else:
            sys.exit("unsupported datatype for Fortran: {}".format(self.name))

        return tstr


def get_prototype(plist):

    decls = [p.get_prototype_decl() for p in plist if not p.get_prototype_decl() is None]

    indent = 4

    prototype = "void set_castro_method_params\n"
    prototype += (indent-1)*" " + "("

    for n, d in enumerate(decls):
        prototype += "{}".format(d)
        if n == len(decls)-1:
            prototype += ");\n"
            break

        else:
            prototype += ", "

        if n % 2 == 1:
            prototype += "\n"
            prototype += indent*" "

    return prototype


def get_cpp_call(plist):

    args = [p.name for p in plist if p.in_fortran == 1]

    indent = 4

    call = ""

    # any special ifdefs to deal with?
    for p in plist:
        if p.ifdef is None or not p.in_fortran: continue
        
        call += "#ifndef {}\n".format(p.ifdef)
        call += p.get_static_default_string()
        call += "#endif\n"
        

    call += "set_castro_method_params\n"

    call += (indent-1)*" " + "("

    for n, a in enumerate(args):
        call += "{}".format(a)
        if n == len(args)-1:
            call += ");\n"
            break

        else:
            call += ", "

        if n % 2 == 1:
            call += "\n"
            call += indent*" "

    return call


def write_meth_module(plist, meth_template):
    """this writes the meth_params_module, starting with the meth_template
       and inserting the runtime parameter declaration in the correct
       place
    """

    try: mt = open(meth_template, "r")
    except:
        sys.exit("invalid template file")

    try: mo = open("Src_nd/meth_params.f90", "w")
    except:
        sys.exit("unable to open meth_params.f90 for writing")

    param_decls = [p.get_f90_decl_string() for p in plist if p.in_fortran == 1]

    decls = ""

    for p in param_decls:
        decls += "  {}".format(p)

    for line in mt:
        if line.find("@@f90_declaractions@@") > 0:
            mo.write(decls)
        else:
            mo.write(line)

    mo.close()
    mt.close()


def write_set_meth_sub(plist, set_template):
    """this writes the set_castro_meth_params routine, starting with the
       set_template and inserting the runtime parameter info where needed
    """

    try: st = open(set_template, "r")
    except:
        sys.exit("invalid template file")

    try: so = open("Src_nd/set_castro_params.f90", "w")
    except:
        sys.exit("unable to open set_castro_params.f90 for writing")

    params = [p for p in plist if p.in_fortran == 1]

    for line in st:
        if line.find("@@start_sub@@") >= 0:
            so.write("subroutine set_castro_method_params( &\n  ")
            for n, p in enumerate(params):
                so.write("{}_in".format(p.f90_name))

                if n == len(params)-1:
                    so.write(") bind(C)\n")
                else:
                    so.write(", ")

                if n % 3 == 2:
                    so.write(" &\n  ")

        elif line.find("@@set_params@@") >= 0:

            # first the declarations
            for p in params:
                if p.dtype == "int":
                    so.write("  integer,          intent(in) :: {}_in\n".format(p.f90_name))
                elif p.dtype == "Real":
                    so.write("  double precision, intent(in) :: {}_in\n".format(p.f90_name))
                else:
                    sys.exit("unsupported datatype for Fortran: {}".format(p.name))

            so.write("\n")

            # now store it in the module -- be aware of changing data types
            for p in params:
                if p.dtype == p.f90_dtype:
                    so.write("  {} = {}_in\n".format(p.f90_name, p.f90_name))
                else:
                    if p.dtype == "int" and p.f90_dtype == "logical":
                        so.write("  {} = {}_in .ne. 0\n".format(p.f90_name, p.f90_name))
                    else:
                        sys.exit("unsupported combination of data types")


        elif line.find("@@end_sub@@") >= 0:
            so.write("end subroutine set_castro_method_params\n")

        else:
            so.write(line)

    so.close()
    st.close()


def parse_params(infile, meth_template, set_template):

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

    # write castro_set_meth.H
    prototype = get_prototype(params)
    try: csm = open("castro_set_meth.H", "w")
    except:
        sys.exit("unable to open castro_set_meth.H for writing")
        
    csm.write(prototype)
    csm.close()

    # write castro_call_set_meth.H
    call = get_cpp_call(params)
    try: ccsm = open("castro_call_set_meth.H", "w")
    except:
        sys.exit("unable to open castro_call_set_meth.H for writing")
        
    ccsm.write(call)
    ccsm.close()

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
        cq.write(p.get_query_string())

    cq.close()


    # write the Fortran routines
    write_meth_module(params, meth_template)
    write_set_meth_sub(params, set_template)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", type=str, default=None,
                        help="template for the meth_params module")
    parser.add_argument("-s", type=str, default=None,
                        help="template for the castro_set_meth_params subroutine")
    parser.add_argument("input_file", type=str, nargs=1,
                        help="input file containing the list of parameters we will define")

    args = parser.parse_args()

    parse_params(args.input_file[0], args.m, args.s)
