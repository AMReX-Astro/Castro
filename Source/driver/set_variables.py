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
# 2. state_indices.H
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

def split_pair(pair_string):
    """given an option of the form "(val1, val2)", split it into val1 and
    val2"""
    return pair_string.replace("(", "").replace(")", "").replace(" ", "").split(",")

class Index:
    """an index that we want to set"""

    def __init__(self, name, var, default_group=None, iset=None,
                 also_adds_to=None, count=1, ifdef=None, exists=True):
        """ parameters:
               name: a descriptive name for the quantity
               var: name of the variable
               default_group: the name of a counter that we increment (e.g., NVAR)
               iset: a descriptive name for the set of the variables this belongs to
                     (e.g., conserved)
               also_adds_to: any other counters that we increment
               count: the number of variables in this group
               ifdef: any ifdef that wraps this variable
               exists: true if the ifdef says it is defined
        """
        self.name = name
        self.var = var
        self.iset = iset
        self.default_group = default_group
        self.adds_to = also_adds_to

        # this will be the integer value of the index -- this can only be
        # set once we add the index to the counter
        self.value = None
        self.cxx_value = None

        # count may have different names in Fortran and C++
        if count.startswith("("):
            self.count, self.count_cxx = split_pair(count)
        else:
            self.count = count
            self.count_cxx = count

        self.ifdef = ifdef

        self.exists = exists

    def __str__(self):
        return self.var

    def set_value(self, val, cxx_val):
        self.value = val
        self.cxx_value = cxx_val

    def get_f90_set_string(self, set_default=None):
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
            sstr += "    {} = {}\n".format(self.var, self.value)
            sstr += "  else\n"
            sstr += "    {} = {}\n".format(self.var, set_default)
            sstr += "  endif\n"
        else:
            sstr += "  {} = {}\n".format(self.var, self.value)

        if self.ifdef is not None:
            sstr += "#endif\n"
        sstr += "\n"
        return sstr

    def get_cxx_set_string(self, set_default=None):
        """return the C++ code that sets this variable index.  Note: since C++
        is 0-based, we subtract 1, so we sync with the Fortran
        value
        """
        sstr = "  constexpr int {} = {};\n".format(self.var, self.cxx_value)
        return sstr


class Counter:
    """a simple object to keep track of how many variables there are in a
    set"""

    def __init__(self, name, cxx_name=None, starting_val=1):
        """name: the name of that counter (this will be used in Fortran)"""

        self.name = name
        if cxx_name is None:
            self.cxx_name = name
        else:
            self.cxx_name = cxx_name
        self.numeric = starting_val

        self.strings = []
        self.cxx_strings = []

        self.starting_val = starting_val
        self.cxx_starting_val = starting_val-1

    def add_index(self, index):
        """increment the counter"""
        try:
            i = int(index.count)
        except ValueError:
            self.strings.append(index.count.strip())
            self.cxx_strings.append(index.count_cxx.strip())
        else:
            self.numeric += i

    def get_value(self, offset=0):
        """return the current value of the counter"""
        if self.strings:
            val = "{} + {}".format(self.numeric-offset, " + ".join(self.strings))
        else:
            val = "{}".format(self.numeric-offset)

        return val

    def get_cxx_value(self, offset=0):
        """return the current value of the counter for C++ (0-based)"""
        if self.strings:
            val = "{} + {}".format(self.numeric - offset - 1, " + ".join(self.cxx_strings))
        else:
            val = "{}".format(self.numeric - offset - 1)

        return val

    def get_f90_set_string(self):
        """return the Fortran needed to set this as a parameter"""
        return "integer, parameter :: {} = {}".format(
            self.name, self.get_value(offset=self.starting_val))

    def get_cxx_set_string(self):
        """return the C++ needed to set this as a parameter"""
        return "constexpr int {} = {};".format(
            self.cxx_name, self.get_cxx_value(offset=self.cxx_starting_val))


def doit(variables_file, odir, defines, nadv,
         ngroups):

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
                # this stores the total number of state variables in each default set,
                # as a tuple, with the Fortran variable and C++ variable,
                # e.g. default_set[conserved] = (NVAR, NUM_STATE)
                _, current_set, default_group, cxx_group = line.split()
                default_set[current_set] = (default_group, cxx_group)
            else:

                # this splits the line into separate fields.  A field is a
                # single word or a pair in parentheses like "(a, b)"
                fields = re.findall(r'[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)', line)

                name = fields[0]
                var = fields[1]
                adds_to = fields[2]
                count = fields[3]
                ifdef = fields[4]

                # we may be fed a pair of the form (SET, DEFINE),
                # in which case we only add to SET if we define
                # DEFINE
                if adds_to.startswith("("):
                    add_set, define = split_pair(adds_to)
                    if not define in defines:
                        adds_to = None
                    else:
                        adds_to = add_set

                if adds_to == "None":
                    adds_to = None

                if ifdef == "None" or ifdef in defines:
                    exists = True
                else:
                    exists = False

                if ifdef == "None":
                    ifdef = None


                indices.append(Index(name, var, default_group=default_group,
                                     iset=current_set, also_adds_to=adds_to,
                                     count=count, ifdef=ifdef, exists=exists))


    # find the set of set names
    unique_sets = {q.iset for q in indices}

    # we'll keep track of all the counters across all the sets.  This
    # will be used later to write a module that sets parameters with
    # the size of each set
    all_counters = []

    # loop over sets, create the counters, and store any indicies that belong to those

    for s in sorted(unique_sets):

        # these are the indices that belong to the default set s.
        # s will be one of our counters
        set_indices = [q for q in indices if q.iset == s]

        # these indices may also add to other counters
        adds_to = {q.adds_to for q in set_indices if q.adds_to is not None}

        # initialize the counters
        counter_main = Counter(default_set[s][0], cxx_name=default_set[s][1])
        counter_adds = []
        for a in adds_to:
            counter_adds.append(Counter(a))

        # add the indices to the respective counters
        for i in set_indices:

            if not i.exists:
                continue

            # set the integer value for this index to the current
            # counter value
            i.set_value(counter_main.get_value(), counter_main.get_cxx_value())

            # increment the counters
            counter_main.add_index(i)
            if i.adds_to is not None:
                for ca in counter_adds:
                    if ca.name == i.adds_to:
                        ca.add_index(i)

        # store the counters for later writing
        all_counters += [counter_main]
        all_counters += counter_adds

    # now loop over counters and write out code to set the indices
    # all these routines will live in a single file

    # first the Fortran
    with open(os.path.join(odir, "set_indices.F90"), "w") as f:

        f.write(HEADER)

        # loop over sets and create the functions
        for s in sorted(unique_sets):
            subname = "ca_set_{}_indices".format(s)

            set_indices = [q for q in indices if q.iset == s]


            # write the function heading
            sub = ""
            sub += "subroutine {}()\n".format(subname)

            # done with the subroutine interface, now include the modules we need
            sub += "\n\n"
            sub += "  use meth_params_module\n"
            sub += "  use network, only: naux, nspec\n"
            sub += "#ifdef RADIATION\n  use rad_params_module, only : ngroups\n#endif\n"
            sub += "  implicit none\n"

            sub += "\n"

            # write the lines to set the indices
            for i in set_indices:

                # if this variable has an ifdef, make sure it is in
                # defines, otherwise skip
                if not i.exists:
                    continue



                # for variables in the "conserved", primitive,
                # auxillary, or godunov, sets, it may be the case that
                # the variable that defines the count is 0 (e.g. for
                # nadv).  We need to initialize it specially then.
                if s in ["conserved", "primitive", "godunov", "auxiliary"]:
                    sub += i.get_f90_set_string(set_default=0)
                else:
                    sub += i.get_f90_set_string()

            # end the function
            sub += "end subroutine {}\n\n".format(subname)

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
            ss.write("   {}\n".format(ac.get_f90_set_string()))
        ss.write("end module state_sizes_module\n")


    # now the C++
    with open(os.path.join(odir, "state_indices.H"), "w") as f:

        # first write out the counter sizes
        f.write("#ifndef _state_indices_H_\n")
        f.write("#define _state_indices_H_\n")

        f.write("#include <actual_network.H>\n\n")

        f.write("  constexpr int NumAdv = {};\n".format(nadv))
        for ac in all_counters:
            f.write("  {}\n".format(ac.get_cxx_set_string()))
        f.write("  constexpr int npassive = NumSpec + NumAux + NumAdv;\n")

        # we only loop over the default sets for setting indices, not the
        # "adds to", so we don't set the same index twice
        for s in unique_sets:
            set_indices = [q for q in indices if q.iset == s]
            f.write("\n   // {}\n".format(s))
            for i in set_indices:
                if not i.exists:
                    continue
                f.write(i.get_cxx_set_string(set_default=0))

        f.write("\n#endif\n")

def main():

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
        try:
            os.makedirs(args.odir)
        except FileExistsError:
            # this exception is needed in case of a race condition
            # to create the directory by another make target
            pass

    doit(args.variables_file[0], args.odir, args.defines, args.nadv,
         args.ngroups)


if __name__ == "__main__":
    main()
