import argparse
import importlib
import os
import re
from pathlib import Path
import sys

# some parameters that are not defined in _cpp_parameters

whitelist = ["castro.lo_bc",
             "castro.hi_bc",
             "gravity.abs_tol",
             "gravity.rel_tol"]

# we don't have all of the radiation parametrs in the _cpp_parameters
# yet, so we won't check these namespaces

namespace_ignore = ["radiation", "radsolve"]

def doit(castro_dir):

    castro = Path(os.path.abspath(castro_dir))

    # import the module that defines the Castro runtime params
    sys.path.append(str(castro / "Source" / "driver/"))
    import parse_castro_params

    # read in the parameters defined in _cpp_parameters
    param_file = castro / "Source" / "driver" / "_cpp_parameters"
    params = parse_castro_params.read_param_file(str(param_file))

    namespaces = set(p.namespace for p in params)
    runtime_parameters = [f"{p.namespace}.{p.name}" for p in params]

    pattern = re.compile(r"[A-Za-z0-9_]+\.[A-Za-z0-9_]+", re.IGNORECASE)

    # loop over all the inputs files
    exec_path = castro / "Exec"
    for f in exec_path.glob("**/inputs*"):

        if os.path.isdir(f):
            continue

        # find all the params in each namespace
        with open(f) as infile:
            print(f"working on {f}")
            for line in infile:
                # remove comments
                idx = line.find("#")
                if idx > 0:
                    line = line[:idx]

                found_param = pattern.match(line)
                if not found_param:
                    continue

                p = found_param.group(0)
                nm = p.split(".")[0]
                if nm in namespaces and nm not in namespace_ignore:
                    if not (p in runtime_parameters or p in whitelist):
                        sys.exit(f"Error: {p} not valid")


if __name__ == "__main__":

    # we need the top-level Castro directory

    p = argparse.ArgumentParser()
    p.add_argument("castro_dir", type=str, nargs=1,
                   help="top level Castro directory")

    args = p.parse_args()

    doit(args.castro_dir[0])


