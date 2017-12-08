# parse the _variables file and write the set of functions that will
# define the indices
#
# ca_set_conserved_indices: the conserved state
#
# ca_set_primitive_indices: the primitive variable state
#
# ca_set_godunov_indices: the interface state
#
# ca_set_auxillary_indices: the auxillary state information

class Index(object):
    """an index that we want to set"""
    def __init__(self, name, f90_var, default_group=None, iset=None,
                 also_adds_to=None, count=1, cxx_var=None, ifdef=None):
        self.name = name
        self.cxx_var = cxx_var
        self.f90_var = f90_var
        self.iset = iset
        self.default_group = default_group
        self.adds_to = also_adds_to
        self.count = count
        self.ifdef = ifdef

    def __str__(self):
        return self.f90_var

    def get_set_string(self):
        sstr = ""
        sstr += "  {} = {}\n".format(self.f90_var, self.default_group)
        sstr += "  {} = {} + {}\n".format(self.default_group, self.default_group, self.count)
        if self.adds_to is not None:
            sstr += "  {} = {} + {}\n".format(self.adds_to, self.adds_to, self.count)
        if self.cxx_var is not None:
            sstr += "  call check_equal({},{}+1)\n".format(self.f90_var, self.cxx_var)
        sstr += "\n"
        return sstr

def doit():

    # read the file and create a list of indices
    indices = []
    default_set = {}
    with open("_variables", "r") as f:
        current_set = None
        default_group = None
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            elif line.startswith("@"):
                print(line)
                _, current_set, default_group = line.split()
                default_set[current_set] = default_group
            else:
                name, cxx_var, f90_var, adds_to, count, ifdef = line.split()
                if adds_to == "None":
                    adds_to = None
                if cxx_var == "None":
                    cxx_var = None
                if ifdef == "None":
                    ifdef = None

                indices.append(Index(name, f90_var, default_group=default_group,
                                     iset=current_set, also_adds_to=adds_to,
                                     count=count, cxx_var=cxx_var, ifdef=ifdef))


    # find the set of set names
    unique_sets = set([q.iset for q in indices])

    # all these routines will live in a single file
    with open("set_indices.F90", "w") as f:

        # loop over sets and create the function
        # arg list will be C++ names to compare to
        for s in unique_sets:
            subname = "ca_set_{}_indices".format(s)

            set_indices = [q for q in indices if q.iset == s]

            # within this set, find the unique ifdef names
            ifdefs = set([q.ifdef for q in set_indices])

            # cxx names in this set
            cxx_names = set([q.cxx_var for q in set_indices])

            # add to
            adds_to = set([q.adds_to for q in set_indices])

            # write the function heading
            sub = ""
            sub += "subroutine {}( &\n".format(subname)
            for n, cxx in enumerate(cxx_names):
                if cxx is None:
                    continue
                sub += "          {} {}, &\n".format(" "*len(subname), cxx)

            sub += "          {})\n\n".format(" "*len(subname), cxx)
            sub += "  use meth_params_module\n"

            for cxx in cxx_names:
                if cxx is None:
                    continue
                sub += "  integer, intent(in) :: {}\n".format(cxx)
            sub += "\n"

            # initialize the counters
            sub += "  {} = 1\n\n".format(default_set[s])
            for a in adds_to:
                if a is None:
                    continue
                sub += "  {} = 1\n\n".format(a)

            # write the lines to set the indices
            # first do those without an ifex
            for i in [q for q in set_indices if q.ifdef is None]:
                sub += i.get_set_string()

            for ifd in ifdefs:
                if ifd is None:
                    continue
                sub += "#ifdef {}\n".format(ifd)
                for i in [q for q in set_indices if q.ifdef == ifd]:
                    sub += i.get_set_string()
                sub += "#endif\n"

            # finalize the counters
            sub += "  {} = {} - 1\n".format(default_set[s], default_set[s])
            for a in adds_to:
                if a is None:
                    continue
                sub += "  {} = {} - 1\n".format(a, a)

            # end the function
            sub += "end subroutine {}\n\n".format(subname)

            f.write(sub)

    # write the C++ function signatures


if __name__ == "__main__":
    doit()
