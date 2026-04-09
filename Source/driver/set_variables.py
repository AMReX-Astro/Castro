#!/usr/bin/env python

"""parse the _variables file and write a C++ header that defines the indices
 and size of state arrays.  We write:

  * state_indices.H

The indices are all 0-based.

For the adds-to column, we can take a list of counters or tuples, e.g.,
  [N1, (N2, DEFINE)]
where we only add to N2 is DEFINE is defined as a preprocessor variable.
"""

import argparse
import os
import re


class Index:
    """an index that we want to set"""

    def __init__(self, name, var, default_group=None, iset=None,
                 also_adds_to=None, count=1):
        """ parameters:
               name: a descriptive name for the quantity
               var: name of the variable
               default_group: the name of a counter that we increment (e.g., NVAR)
               iset: a descriptive name for the set of the variables this belongs to
                     (e.g., conserved)
               also_adds_to: any other counters that we increment
               count: the number of variables in this group
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

        self.count_cxx = count

    def __str__(self):
        return self.var

    def set_value(self, cxx_val):
        self.cxx_value = cxx_val

    def get_cxx_set_string(self):
        """return the C++ code that sets this variable index.  Note: since C++
        is 0-based, we subtract 1, so we sync with the Fortran
        value
        """
        sstr = f"  constexpr int {self.var} = {self.cxx_value};\n"
        return sstr


class Counter:
    """a simple object to keep track of how many variables there are in a
    set"""

    def __init__(self, name, starting_val=1):
        """name: the name of that counter (this will be used in Fortran)"""

        self.name = name

        self.numeric = starting_val

        self.cxx_strings = []

        self.cxx_starting_val = starting_val-1

    def add_index(self, index):
        """increment the counter"""

        try:
            i = int(index.count_cxx)
        except ValueError:
            self.cxx_strings.append(index.count_cxx.strip())
        else:
            self.numeric += i

    def get_cxx_value(self, offset=0):
        """return the current value of the counter for C++ (0-based)"""
        if self.cxx_strings:
            val = "{} + {}".format(self.numeric - offset - 1, " + ".join(self.cxx_strings))
        else:
            val = f"{self.numeric - offset - 1}"

        return val

    def get_cxx_set_string(self):
        """return the C++ needed to set this as a parameter"""
        return "constexpr int {} = {};".format(
            self.name, self.get_cxx_value(offset=self.cxx_starting_val))


def doit(variables_file, odir, defines, nadv):

    # read the file and create a list of indices
    indices = []

    # default_set is the main category to which this index belongs
    # (e.g., conserved, primitive, ...)
    default_set = {}

    with open(variables_file) as f:
        current_set = None
        default_group = None
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            if line.startswith("@"):
                # this stores the total number of state variables in each default set,
                # as a tuple, with the Fortran variable and C++ variable,
                # e.g. default_set[conserved] = (NVAR, NUM_STATE)
                _, current_set, cxx_group = line.split()
                default_set[current_set] = cxx_group
            else:

                # this splits the line into separate fields.  A field is a
                # single word or, for the "also adds to" section, a group in []
                fields = re.findall(r'[\w-]+|\[.*?\]', line)

                name = fields[0]
                var = fields[1]
                adds_to_temp = fields[2]
                count = fields[3]
                ifdef = fields[4]

                # we may be fed a pair of the form (SET, DEFINE),
                # in which case we only add to SET if we define
                # DEFINE
                if adds_to_temp.startswith("["):
                    # now split it into entities that are either single
                    # (counter, ifdef) or just (counter)
                    subfields = re.findall(r'[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)', adds_to_temp)

                    # loop over the subfields and build a list of tuples
                    adds_to = []
                    for s in subfields:
                        if s.startswith("("):
                            counter, define = re.findall(r'[\w]+', s)
                            if define in defines:
                                adds_to.append(counter)
                        else:
                            adds_to.append(s)
                else:
                    if adds_to_temp == "None":
                        adds_to = None
                    else:
                        adds_to = [adds_to_temp]

                # only recognize the index if we defined any required preprocessor variable
                if ifdef == "None" or ifdef in defines:
                    indices.append(Index(name, var, default_group=default_group,
                                         iset=current_set, also_adds_to=adds_to,
                                         count=count))

    # find the set of set names
    unique_sets = {q.iset for q in indices}

    # we'll keep track of all the counters across all the sets.  This
    # will be used later to write a module that sets parameters with
    # the size of each set
    all_counters = []

    # loop over sets, create the counters, and store any indices that belong to those
    for s in sorted(unique_sets):

        # these are the indices that belong to the default set s.
        # s will be one of our counters
        set_indices = [q for q in indices if q.iset == s]

        # these indices may also add to other counters
        adds_to = []
        for q in set_indices:
            if q.adds_to is not None:
                for v in q.adds_to:
                    adds_to.append(v)

        adds_to = sorted(set(adds_to))

        # initialize the counters
        counter_main = Counter(default_set[s])
        counter_adds = []
        for a in adds_to:
            counter_adds.append(Counter(a))

        # add the indices to the respective counters
        for i in set_indices:

            # set the integer value for this index to the current
            # counter value
            i.set_value(counter_main.get_cxx_value())

            # increment the counters
            counter_main.add_index(i)
            if i.adds_to:
                for ca in counter_adds:
                    if ca.name in i.adds_to:
                        ca.add_index(i)

        # store the counters for later writing
        all_counters += [counter_main]
        all_counters += counter_adds

    # now loop over counters and write out code to set the indices
    # all these routines will live in a single file

    with open(os.path.join(odir, "state_indices.H"), "w") as f:

        # first write out the counter sizes
        f.write("#ifndef STATE_INDICES_H\n")
        f.write("#define STATE_INDICES_H\n")

        f.write("#include <network_properties.H>\n\n")

        f.write(f"  constexpr int NumAdv = {nadv};\n")
        for ac in all_counters:
            f.write(f"  {ac.get_cxx_set_string()}\n")
        f.write("  constexpr int npassive = NumSpec + NumAux + NumAdv;\n")

        # we only loop over the default sets for setting indices, not the
        # "adds to", so we don't set the same index twice
        for s in sorted(unique_sets):
            set_indices = [q for q in indices if q.iset == s]
            f.write(f"\n   // {s}\n")
            for i in set_indices:
                f.write(i.get_cxx_set_string())

        f.write("\n#endif\n")


def main():
    """the main driver"""

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

    doit(args.variables_file[0], args.odir, args.defines, args.nadv)


if __name__ == "__main__":
    main()
