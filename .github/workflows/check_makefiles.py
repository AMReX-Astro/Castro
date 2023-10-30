import re
import sys
from pathlib import Path

correct_params = {
    "DEBUG": "FALSE",
    "USE_MPI": "TRUE",
    "USE_OMP": "FALSE",
    "COMP": "gnu",
    "USE_CUDA": "FALSE",
    "USE_HIP": "FALSE",
    "PRECISION": "DOUBLE",
    "PROFILE": "FALSE"}


def find_source_files():
    p = Path("./Exec")
    files = list(p.glob(r"**/GNUmakefile"))
    return files

def check_makefile(makefile):

    with open(makefile) as mf:
        for _line in mf:
            if idx := _line.find("#") >= 0:
                line = _line[:idx]
            else:
                line = _line

            for key in correct_params:
                if key in line:
                    try:
                        k, v = re.split(":=|\?=|=", line)
                    except ValueError:
                        sys.exit(f"invalid line: {line}")

                    if not v.strip() == correct_params[key]:
                        sys.exit(f"invalid param {key} in {makefile}")

if __name__ == "__main__":

    for f in find_source_files():
        check_makefile(f)



