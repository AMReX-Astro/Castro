#!/usr/bin/env python3

# walk through all of the Fortran files and convert any double
# precision declarations to
#
#   real(c_real) :: x

import os
import re
import sys


class Unit(object):
    def __init__(self, name="", ftype="subroutine", has_decls=False, in_module=False):
        self.name = name
        self.ftype = ftype
        self.has_decls = has_decls
        self.in_module = in_module
        

def find_files(top_dir, extension):
    """ find files with a given extension -- return a list """

    sfiles = []

    for dir_name, _, files in os.walk(top_dir):
        sfiles += [os.path.join(dir_name, f) for f in files
                   if f.endswith(extension) and os.path.isfile(os.path.join(dir_name, f))]

    return sfiles


def main():

    top_dir = os.getcwd()

    extensions = [".f90", ".F90"]

    # define the regular expressions for double precision declaration -- these are
    # what we seek to replace
    
    # match declarations like "real (kind=dp_t)"
    r1 = re.compile(r"(real)\s*(\(\s*kind\s*=\s*dp_t\s*\))", re.IGNORECASE|re.DOTALL)

    # match declarations like "real (dp_t)"
    r2 = re.compile(r"(real)\s*(\(\s*dp_t\s*\))", re.IGNORECASE|re.DOTALL)

    # match declarations like "double precision"
    r3 = re.compile(r"(double)\s*(precision)", re.IGNORECASE|re.DOTALL)

    regexs = [r1, r2, r3]

    # find the source files
    sfiles = []
    for e in extensions:
        sfiles += find_files(top_dir, e)


    # main loop over the source files
    for sf in sfiles:

        # the tricky part of the conversion is that we need to add the
        # "use bl_fort_module" to the source in any program unit (or
        # scope) that has double precision declarations.
        
        # keep track of all the program units in file file (modules,
        # subroutines, etc.)
        units = []
        
        # read the file
        try:
            f = open(sf, "r")
        except IOError:
            sys.exit("error reading file {}".format(sf))
        else:
            lines = f.readlines()
            f.close()

        # move the source file to a backup
        try:
            os.rename(sf, "{}_orig".format(sf))
        except:
            sys.exit("error renaming {}".format(sf))
            
        # parse it first looking for subroutine and function
        # definitions, marking which have any double precision
        # definitions
        current_unit = None


        # now write out the file, line by line.  When we encounter a
        # subroutine or function, add the necessary module line.  When
        # we encounter a declaration, convert it to the new form
        new_lines = []
        for line in lines:
            for r in regexs:
                decl = r.search(line)
                if decl:
                    line = line.replace(decl.group(0), "real(c_real)")
            new_lines.append(line)

        # skip comments
        try:
            f = open(sf, "w")
        except IOError:
            sys.exit("error openning {} for writing".format(sf))
        else:
            f.writelines(new_lines)
            f.close()
            
if __name__ == "__main__":
    main()
