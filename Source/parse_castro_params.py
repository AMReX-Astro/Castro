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
        return "pp.query(\"{}\", {});\n".format(self.name, self.name)

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


def write_prototype(plist):

    decls = [p.get_prototype_decl() for p in plist if not p.get_prototype_decl() is None]

    indent = 4

    prototype = "BL_FORT_PROC_DECL(SET_CASTRO_METHOD_PARAMS, set_castro_methd_params)\n"
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


def write_cpp_call(plist):

    args = [p.name for p in plist if p.in_fortran == 1]

    indent = 4

    call = "BL_FORT_PROC_CALL(SET_CASTRO_METHOD_PARAMS, set_castro_method_params)\n"

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


def write_fortran_decl(plist):

    params = [p.get_f90_decl_string() for p in plist if p.in_fortran == 1]

    decls = ""
    
    for p in params:
        decls += p

    return decls


def parser(infile):

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


    for p in params:
        print p.get_default_string()

    print write_prototype(params)
    print write_cpp_call(params)
    print write_fortran_decl(params)

if __name__ == "__main__":

    try: infile = sys.argv[1]
    except:
        sys.exit("need to specify an input file")

    parser(infile)
