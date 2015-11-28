import re
import sys

class Param(object):
    def __init__(self, name, dtype, default, 
                 debug_default=None, in_fortran=0, ifdef=None):
        self.name = name
        self.dtype = dtype
        self.default = default
        self.debug_default=debug_default
        self.in_fortran = in_fortran
        self.ifdef = ifdef

    def get_default_string(self):
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
            ostr += "{} = {}\n".format(tstr, self.debug_default)
            ostr += "#else\n"
            ostr += "{} = {}\n".format(tstr, self.default)
            ostr += "#endif\n"
        else:
            ostr += "{} = {}\n".format(tstr, self.default)
            
        if not self.ifdef is None:
            ostr += "#endif\n"

        return ostr


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

        try: in_fortran = fields[3]
        except: in_fortran = 0

        try: ifdef = fields[4]
        except: ifdef = None

        params.append(Param(name, dtype, default, debug_default=debug_default, ifdef=ifdef))


    for p in params:
        print p.get_default_string()



if __name__ == "__main__":
    
    try: infile = sys.argv[1]
    except: 
        sys.error("need to specify an input file")

    parser(infile)

