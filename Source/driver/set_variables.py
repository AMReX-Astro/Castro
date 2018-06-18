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
# 2. set_conserved.H
#
#    This simply sets the C++ indices
#

import re

HEADER = """
! DO NOT EDIT!!!

! This file is automatically created by set_variables.py.  To update
! or add variable indices, please edit _variables and then rerun the
! script.

"""

CHECK_EQUAL = """
subroutine check_equal(index1, index2)

  use amrex_error_module

  implicit none

  integer, intent(in) :: index1, index2

#ifndef AMREX_USE_CUDA
  if (index1 /= index2) then
    call amrex_error("ERROR: mismatch of indices")
  endif
#endif

end subroutine check_equal


"""


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

        # count may have different names in Fortran and C++
        if count.startswith("("):
            self.count, self.count_cxx = count.replace("(", "").replace(")", "").split(",")
        else:
            self.count = count
            self.count_cxx = count

        self.ifdef = ifdef

    def __str__(self):
        return self.f90_var

    def get_set_string(self, set_default=None):
        sstr = ""
        if self.ifdef is not None:
            sstr += "#ifdef {}\n".format(self.ifdef)

        if set_default is not None and self.count != "1":
            sstr += "  if ({} > 0) then\n".format(self.count)
            sstr += "    {} = {}\n".format(self.f90_var, self.default_group)
            sstr += "    {} = {} + {}\n".format(self.default_group, self.default_group, self.count)
            sstr += "  else\n"
            sstr += "    {} = {}\n".format(self.f90_var, set_default)
            sstr += "  endif\n"
        else:
            sstr += "  {} = {}\n".format(self.f90_var, self.default_group)
            sstr += "  {} = {} + {}\n".format(self.default_group, self.default_group, self.count)

        if self.adds_to is not None:
            sstr += "  {} = {} + {}\n".format(self.adds_to, self.adds_to, self.count)
        if self.cxx_var is not None:
            sstr += "  call check_equal({},{}+1)\n".format(self.f90_var, self.cxx_var)
        if self.ifdef is not None:
            sstr += "#endif\n"
        sstr += "\n"
        return sstr

    def get_cxx_set_string(self):
        sstr = ""
        if self.ifdef is not None:
            sstr += "#ifdef {}\n".format(self.ifdef)

        if self.count != "1":
            sstr += "  if ({} > 0) {{\n".format(self.count_cxx)
            sstr += "    {} = cnt;\n".format(self.cxx_var)
            sstr += "    cnt += {};\n".format(self.count_cxx)
            sstr += "  }\n"
        else:
            sstr += "  {} = cnt;\n".format(self.cxx_var)
            sstr += "  cnt += {};\n".format(self.count_cxx)

        if self.ifdef is not None:
            sstr += "#endif\n"
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
                _, current_set, default_group = line.split()
                default_set[current_set] = default_group
            else:

                # this splits the line into separate fields.  A field is a
                # single word or a pair in parentheses like "(a, b)"
                fields = re.findall(r'[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)', line)

                name = fields[0]
                cxx_var = fields[1]
                f90_var = fields[2]
                adds_to = fields[3]
                count = fields[4].replace(" ","").strip()
                ifdef = fields[5]

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

        f.write(HEADER)

        # loop over sets and create the function
        # arg list will be C++ names to compare to
        f.write(CHECK_EQUAL)

        for s in sorted(unique_sets):
            subname = "ca_set_{}_indices".format(s)

            set_indices = [q for q in indices if q.iset == s]

            # cxx names in this set
            cxx_names = [q.cxx_var for q in set_indices if q.cxx_var is not None]

            # add to
            adds_to = set([q.adds_to for q in set_indices if q.adds_to is not None])

            # write the function heading
            sub = ""

            if not cxx_names:
                sub += "subroutine {}()\n".format(subname)
            else:
                sub += "subroutine {}( &\n".format(subname)
                # we need to put all of the arguments that are ifdef-ed up front
                # so we can properly close the argument list (no hanging commas)
                cxx_with_ifdef = [q for q in set_indices if q.cxx_var is not None and q.ifdef is not None]
                cxx_wo_ifdef = [q for q in set_indices if q.cxx_var is not None and q.ifdef is None]

                cxx_all = cxx_with_ifdef + cxx_wo_ifdef
                for n, i in enumerate(cxx_all):
                    if i.cxx_var is not None:
                        if i.ifdef is not None:
                            sub += "#ifdef {}\n".format(i.ifdef)
                        if n == len(cxx_all)-1:
                            sub += "           {} {} &\n".format(" "*len(subname), i.cxx_var)
                        else:
                            sub += "           {} {}, &\n".format(" "*len(subname), i.cxx_var)
                        if i.ifdef is not None:
                            sub += "#endif\n"

                sub += "           {})\n".format(" "*len(subname))

            sub += "\n\n"
            sub += "  use meth_params_module\n"
            sub += "  use network, only: naux, nspec\n"
            sub += "#ifdef RADIATION\n  use rad_params_module, only : ngroups\n#endif\n"
            sub += "  implicit none\n"

            for i in set_indices:
                if i.cxx_var is None:
                    continue
                if i.ifdef is not None:
                    sub += "#ifdef {}\n".format(i.ifdef)
                sub += "  integer, intent(in) :: {}\n".format(i.cxx_var)
                if i.ifdef is not None:
                    sub += "#endif\n"
            sub += "\n"

            # initialize the counters
            sub += "  {} = 1\n\n".format(default_set[s])
            for a in adds_to:
                sub += "  {} = 1\n\n".format(a)

            # write the lines to set the indices
            # first do those without an ifex
            for i in set_indices:
                if s == "conserved":
                    sub += i.get_set_string(set_default=0)
                else:
                    sub += i.get_set_string()

            # finalize the counters
            sub += "  {} = {} - 1\n".format(default_set[s], default_set[s])
            for a in adds_to:
                if a is None:
                    continue
                sub += "  {} = {} - 1\n".format(a, a)

            # end the function
            sub += "end subroutine {}\n\n".format(subname)

            f.write(sub)

    # write the C++ include
    conserved_indices = [q for q in indices if q.iset == "conserved" and q.cxx_var is not None]

    with open("set_conserved.H", "w") as f:
        f.write("  int cnt = 0;\n")
        for c in conserved_indices:
            f.write(c.get_cxx_set_string())


if __name__ == "__main__":
    doit()
