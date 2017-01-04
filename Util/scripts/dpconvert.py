#!/usr/bin/env python3

# walk through all of the Fortran files and convert any double
# precision declarations to
#
#   real(c_real) :: x

import os
import re
import sys


# routine signatures
routine_re = re.compile(r"^(subroutine|function|module)\s+(\w*)", re.IGNORECASE|re.DOTALL)

# declaration signatures
decl_re = re.compile(r"^(real|integer|double\s*precision|logical|character|type)", re.IGNORECASE|re.DOTALL)

# implicit none signature
implno_re = re.compile(r"^(implicit\s+none)", re.IGNORECASE|re.DOTALL)

# use module signature
use_re = re.compile(r"^(use)\s+(\w*)", re.IGNORECASE|re.DOTALL)


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


    # also we want to convert numeric constants that are of the form X.Y_dp_t, etc.
    c_re = re.compile(r"([0-9edED.\+\-]+)(_dp_t)", re.IGNORECASE|re.DOTALL)

    # and... we want to replace and "d" scientific notation with the new style
    # this matches stuff like -1.25d-10, and gives us separate groups for the
    # prefix and exponent
    d_re = re.compile(r"([\+\-0-9.]+)[dD]([\+\-0-9]+)", re.IGNORECASE|re.DOTALL)

    # find the source files
    sfiles = []
    for e in extensions:
        sfiles += find_files(top_dir, e)


    # main loop over the source files
    for sf in sfiles:

        # the tricky part of the conversion is that we need to add the
        # "use bl_fort_module" to the source in any program unit (or
        # scope) that has double precision declarations.

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

        # keep track of all the program units in file file (modules,
        # subroutines, etc.)
        units = {}

        current_unit = None
        has_decl = False

        for line in lines:
            tline = line.strip()

            # are we a routine?
            rout = routine_re.search(tline)
            if rout:
                # save the old info
                if current_unit is not None:
                    units[current_unit] = has_decl

                # start the next unit
                print("here: {}".format(rout.group(0)))
                print(tline)
                current_unit = rout.group(2)
                has_decl = False

            # are there any declarations that we need to replace?
            for r in regexs:
                decl = r.search(line)
                if decl:
                    has_decl = True

        # we never stored the last found unit
        if current_unit is not None and current_unit not in units:
            units[current_unit] = has_decl


        # now write out the file, line by line.  When we encounter a
        # subroutine or function, add the necessary module line.  When
        # we encounter a declaration, convert it to the new form
        current_unit = None
        has_decl = False

        print("keys:", units.keys())

        new_lines = []
        while lines:

            line = lines.pop(0)

            # is this the start of a routine?
            rout = routine_re.search(line.strip())
            if rout:
                current_unit = rout.group(2)
                has_decl = units[current_unit]

                new_lines.append(line)

                while has_decl:
                    # we need to add the use line before any declarations or 
                    # and implicit none.  We don't simply add it after first 
                    # detecting the start of a routine, since the function
                    # definition might span multiple lines
                    
                    line = lines.pop(0)
                    
                    is_decl = decl_re.search(line.strip())
                    is_implno = implno_re.search(line.strip())

                    if is_decl or is_implno:
                        # we need to add the module use statement
                        new_lines.append("use bl_fort_module\n")
                        has_decl = False

                        # this line will be added at the end of the while

                    else:
                        new_lines.append(line)


            # replace declarations
            for r in regexs:
                decl = r.search(line)
                if decl:
                    old_decl = decl.group(0)
                    new_decl = "real(c_real)"
                    # preserve spacing, if possible
                    lo = len(old_decl) - len(new_decl)
                    if lo > 0:
                        new_decl = new_decl + lo*" "
                    line = line.replace(old_decl, new_decl)

            # replace constants
            const = c_re.search(line)
            if const:
                num_value = const.group(1)
                new_const = "{}_rt".format(num_value)
                old_const = const.group(0)
                line = line.replace(old_const, new_const)

            # update "d" scientific notation -- allow for multiple constants on a single line
            for dd in d_re.finditer(line):
                prefix = dd.group(1)
                exponent = dd.group(2)
                new_num = "{}e{}_rt".format(prefix, exponent)
                old_num = dd.group(0)
                line = line.replace(old_num, new_num)

            new_lines.append(line)


        try:
            f = open(sf, "w")
        except IOError:
            sys.exit("error openning {} for writing".format(sf))
        else:
            f.writelines(new_lines)
            f.close()

if __name__ == "__main__":
    main()
