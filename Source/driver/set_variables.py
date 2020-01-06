#!/usr/bin/env python3

# parse the _variables file and write the set of functions that will
# define the indices.  We write two files with the following functions:
#
# 1. set_indices.F90:
#
#    * ca_set_auxiliary_indices: the auxiliary state information
#
#    * ca_set_conserved_indices: the conserved state
#
#    * ca_set_godunov_indices: the interface state
#
#    * ca_set_primitive_indices: the primitive variable state
#
# 2. set_conserved.H, set_primitive.H, set_godunov.H
#
#    This simply sets the C++ indices
#

import argparse
import os
import re
import sys

HEADER = """
! DO NOT EDIT!!!

! This file is automatically created by set_variables.py.  To update
! or add variable indices, please edit _variables and then rerun the
! script.

"""

CHECK_EQUAL = """
subroutine check_equal(index1, index2)

  use castro_error_module

  implicit none

  integer, intent(in) :: index1, index2

#ifndef AMREX_USE_CUDA
  if (index1 /= index2) then
    call castro_error("ERROR: mismatch of indices")
  endif
#endif

end subroutine check_equal


"""

LINE_SPLIT = r"[\w\"\+\.-]+|\[[\w\.\-\(\)\s,]+\]|\([\w\.\-\s,]+\)"
ADDS_TO_PAIR = r"[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)"

def split_pair(pair_string):
    """given an option of the form "(val1, val2)", split it into val1 and
    val2"""
    return pair_string.replace("(", "").replace(")", "").replace(" ", "").split(",")

def split_adds_to_list(pair_string):
    """given an option of the form "[val1, val2]", split it into val1 and
    val2.  Note, each of val1 and val2 could itself be a pair like "(a, b)". """
    return re.findall(ADDS_TO_PAIR, pair_string)

class Index:
    """an index that we want to set"""

    def __init__(self, name, f90_var, default_group=None, iset=None,
                 also_adds_to=None, count=1, cxx_var=None, ifdef=None):
        """ parameters:
               name: a descriptive name for the quantity
               f90_var: name of the variable in Fortran
               default_group: the name of a counter that we increment (e.g., NVAR)
               iset: a descriptive name for the set of the variables this belongs to
                     (e.g., conserved)
               also_adds_to: any other counters that we increment
               count: the number of variables in this group
               cxx_var: the name of the variable in C++
               ifdef: any ifdef that wraps this variable
        """
        self.name = name
        self.cxx_var = cxx_var
        self.f90_var = f90_var
        self.iset = iset
        self.default_group = default_group
        self.adds_to = also_adds_to

        # count may have different names in Fortran and C++
        if count.startswith("("):
            self.count, self.count_cxx = split_pair(count)
        else:
            self.count = count
            self.count_cxx = count

        self.ifdef = ifdef

    def __str__(self):
        return self.f90_var

    def get_set_string(self, val, set_default=None):
        """return the Fortran code that sets this variable index (to val).
        Here set_default is a value to set the key to in the case that
        a string value (like nspec) is 0

        """
        sstr = ""
        if self.ifdef is not None:
            sstr += "#ifdef {}\n".format(self.ifdef)

        if set_default is not None and self.count != "1":
            # this adds a test that the count is greater than 0
            sstr += "  if ({} > 0) then\n".format(self.count)
            sstr += "    {} = {}\n".format(self.f90_var, val)
            sstr += "  else\n"
            sstr += "    {} = {}\n".format(self.f90_var, set_default)
            sstr += "  endif\n"
        else:
            sstr += "  {} = {}\n".format(self.f90_var, val)

        if self.cxx_var is not None:
            sstr += "  call check_equal({},{}_in+1)\n".format(self.f90_var, self.cxx_var)
        if self.ifdef is not None:
            sstr += "#endif\n"
        sstr += "\n"
        return sstr

    def get_cxx_set_string(self):
        """get the C++ code that sets the variable index and increments the
        counters"""

        if self.iset == "primitive":
            counter = "qcnt"
        elif self.iset == "godunov":
            counter = "gcnt"
        else:
            counter = "cnt"

        sstr = ""
        if self.ifdef is not None:
            sstr += "#ifdef {}\n".format(self.ifdef)

        if self.count != "1":
            sstr += "  if ({} > 0) {{\n".format(self.count_cxx)
            sstr += "    {} = {};\n".format(self.cxx_var, counter)
            sstr += "    {} += {};\n".format(counter, self.count_cxx)
            sstr += "  }\n"
        else:
            sstr += "  {} = {};\n".format(self.cxx_var, counter)
            sstr += "  {} += {};\n".format(counter, self.count_cxx)

        if self.ifdef is not None:
            sstr += "#endif\n"
        sstr += "\n"
        return sstr


class Counter:
    """a simple object to keep track of how many variables there are in a
    set"""

    def __init__(self, name, starting_val=1):
        """name: the name of that counter (this will be used in Fortran)"""

        self.name = name
        self.numeric = starting_val
        self.strings = []

        self.starting_val = starting_val

    def increment(self, value):
        """increment the pointer to the variable start by value"""
        try:
            i = int(value)
        except ValueError:
            self.strings.append(value.strip())
        else:
            self.numeric += i

    def get_value(self, offset=0):
        """return the current value of the counter"""
        if self.strings:
            val = "{} + {}".format(self.numeric-offset, " + ".join(self.strings))
        else:
            val = "{}".format(self.numeric-offset)

        return val

    def get_set_string(self):
        """return the Fortran needed to set this as a parameter"""
        return "integer, parameter :: {} = {}".format(
            self.name, self.get_value(offset=self.starting_val))


def doit(variables_file, odir, defines, nadv,
<<<<<<< HEAD
         ngroups,
         n_neutrino_species, neutrino_groups, debug=0):
    """the main driver for reading _variables and writing the include files"""
=======
         ngroups):
>>>>>>> development

    # are we doing radiation?
    if not "RADIATION" in defines:
        ngroups = None

    # read the file and create a list of indices
    indices = []
    default_set = {}
    with open(variables_file, "r") as f:
        current_set = None
        default_group = None
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            elif line.startswith("@"):
                _, current_set, default_group = line.split()
                default_set[current_set] = default_group
            else:

                # this splits the line into separate fields.  A field is a
                # single word or a pair in parentheses like "(a, b)"
                fields = re.findall(LINE_SPLIT, line)

                name = fields[0]
                cxx_var = fields[1]
                f90_var = fields[2]
                adds_to_raw = fields[3]
                count = fields[4]
                ifdef = fields[5]

                # the adds_to column can use [] to hold a list of
                # variables that we need to add to.  Split this into
                # a list now
                if adds_to_raw == "None":
                    adds_to = []
                else:
                    adds_to = split_adds_to_list(adds_to_raw)

                # each item in the adds_to list may be of the form
                # (SET, DEFINE), in which case we only add to SET if
                # we define DEFINE
                all_adds = []
                for a in adds_to:
                    if a.startswith("("):
                        add_set, define = split_pair(a)
                        if define in defines:
                            all_adds.append(add_set)
                    else:
                        all_adds.append(a)

                if cxx_var == "None":
                    cxx_var = None
                if ifdef == "None":
                    ifdef = None

                indices.append(Index(name, f90_var, default_group=default_group,
                                     iset=current_set, also_adds_to=all_adds,
                                     count=count, cxx_var=cxx_var, ifdef=ifdef))


    # find the set of set names
    unique_sets = {q.iset for q in indices}

    # we'll keep track of all the counters across all the sets.  This
    # will be used later to write a module that sets parameters with
    # the size of each set
    all_counters = []

    # all these routines will live in a single file
    with open(os.path.join(odir, "set_indices.F90"), "w") as f:

        f.write(HEADER)
        f.write(CHECK_EQUAL)

        # loop over sets and create the functions
        for s in sorted(unique_sets):
            subname = "ca_set_{}_indices".format(s)

            if debug: print("working on set: {}".format(s))

            set_indices = [q for q in indices if q.iset == s]

            # cxx names in this set
            cxx_names = [q.cxx_var for q in set_indices if q.cxx_var is not None]

            # add to
            adds_to = []
            for idx in set_indices:
                adds_to += [q for q in idx.adds_to]
            adds_to = set(adds_to)
            if debug: print("adds to = ", adds_to)

            # write the function heading
            sub = ""

            # arg list will be C++ names to compare to (if any)
            if not cxx_names:
                sub += "subroutine {}()\n".format(subname)
            else:
                sub += "subroutine {}( &\n".format(subname)
                # we need to put all of the arguments that are ifdef-ed up front
                # so we can properly close the argument list (no hanging commas)
                # note: the argument list will always include the ifdefs, so we have
                # a consistent interface we can call.  Below we will only set those
                # variables that have an ifdef in defines
                cxx_with_ifdef = [q for q in set_indices if q.cxx_var is not None and q.ifdef is not None]
                cxx_wo_ifdef = [q for q in set_indices if q.cxx_var is not None and q.ifdef is None]

                cxx_all = cxx_with_ifdef + cxx_wo_ifdef
                for n, i in enumerate(cxx_all):
                    if i.cxx_var is not None:
                        if i.ifdef is not None:
                            sub += "#ifdef {}\n".format(i.ifdef)
                        if n == len(cxx_all)-1:
                            sub += "           {} {}_in &\n".format(" "*len(subname), i.cxx_var)
                        else:
                            sub += "           {} {}_in, &\n".format(" "*len(subname), i.cxx_var)
                        if i.ifdef is not None:
                            sub += "#endif\n"

                sub += "           {})\n".format(" "*len(subname))

            # done with the subroutine interface, now include the modules we need
            sub += "\n\n"
            sub += "  use meth_params_module\n"
            sub += "  use network, only: naux, nspec\n"
            sub += "#ifdef RADIATION\n  use rad_params_module, only : ngroups\n#endif\n"
            sub += "  implicit none\n"

            # declare the arguments
            for i in set_indices:
                if i.cxx_var is None:
                    continue
                if i.ifdef is not None:
                    sub += "#ifdef {}\n".format(i.ifdef)
                sub += "  integer, intent(in) :: {}_in\n".format(i.cxx_var)
                if i.ifdef is not None:
                    sub += "#endif\n"
            sub += "\n"

            # initialize the counters
            counter_main = Counter(default_set[s])
            counter_adds = []
            for a in adds_to:
                counter_adds.append(Counter(a))

            # write the lines to set the indices
            for i in set_indices:

                if debug:
                    print("considering variable {}".format(i.f90_var))
                    print("   should add to: ", i.adds_to)

                # if this variable has an ifdef, make sure it is in
                # defines, otherwise skip
                if i.ifdef is not None:
                    if i.ifdef not in defines:
                        print(" ")
                        continue

                # get the index value for the main counter
                val = counter_main.get_value()

                # increment the counters
                counter_main.increment(i.count)
                for add_var in i.adds_to:
                    for ca in counter_adds:
                        if ca.name == add_var:
                            if debug: print("   incrementing counter {} by {}".format(ca.name, i.count))
                            ca.increment(i.count)


                # for variables in the "conserved", primitive, or godunov, sets,
                # it may be the case that the variable that defines
                # the count is 0 (e.g. for nadv).  We need to
                # initialize it specially then.
                if s in ["conserved", "primitive", "godunov"]:
                    sub += i.get_set_string(val, set_default=0)
                else:
                    sub += i.get_set_string(val)

                if debug: print(" ")

            # end the function
            sub += "end subroutine {}\n\n".format(subname)

            # store the counters for later writing
            all_counters += [counter_main]
            all_counters += counter_adds

            f.write(sub)


    # write the module containing the size of the sets
    with open(os.path.join(odir, "state_sizes.f90"), "w") as ss:
        ss.write("module state_sizes_module\n")
        ss.write("   use network, only : nspec, naux\n")
        ss.write("   implicit none\n")
        ss.write("   integer, parameter :: nadv = {}\n".format(nadv))
        if ngroups is not None:
            ss.write("   integer, parameter :: ngroups = {}\n".format(ngroups))
        for ac in all_counters:
            ss.write("   {}\n".format(ac.get_set_string()))
        ss.write("end module state_sizes_module\n")


    # write the C++ includes
    conserved_indices = [q for q in indices if q.iset == "conserved" and q.cxx_var is not None]

    with open(os.path.join(odir, "set_conserved.H"), "w") as f:
        f.write("  int cnt = 0;\n")
        for c in conserved_indices:
            f.write(c.get_cxx_set_string())

    primitive_indices = [q for q in indices if q.iset == "primitive" and q.cxx_var is not None]

    with open(os.path.join(odir, "set_primitive.H"), "w") as f:
        f.write("  int qcnt = 0;\n")
        for p in primitive_indices:
            f.write(p.get_cxx_set_string())

    godunov_indices = [q for q in indices if q.iset == "godunov" and q.cxx_var is not None]

    with open(os.path.join(odir, "set_godunov.H"), "w") as f:
        f.write("  int gcnt = 0;\n")
        for g in godunov_indices:
            f.write(g.get_cxx_set_string())

def main():
    """read the commandline arguments and pass them to the driver"""

    # note: you need to put a space at the start of the string
    # that gives defines so that the '-' is not interpreted as
    # an option itself
    # https://stackoverflow.com/questions/16174992/cant-get-argparse-to-read-quoted-string-with-dashes-in-it

    parser = argparse.ArgumentParser()
    parser.add_argument("--odir", type=str, default="",
                        help="output directory")
    parser.add_argument("--defines", type=str, default="",
                        help="preprocessor defines to interpret")
    parser.add_argument("--nadv", type=int, default=0,
                        help="the number of pure advected quantities")
    parser.add_argument("--ngroups", type=int, default=1,
                        help="the number of radiation groups")
    parser.add_argument("variables_file", type=str, nargs=1,
                        help="input variable definition file")
    args = parser.parse_args()

    if args.odir != "" and not os.path.isdir(args.odir):
        os.makedirs(args.odir)

    doit(args.variables_file[0], args.odir, args.defines, args.nadv,
         args.ngroups)


if __name__ == "__main__":
    main()
