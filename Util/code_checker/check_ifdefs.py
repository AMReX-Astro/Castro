#!/bin/env python

import re

from pathlib import Path


def find_source_files():
    p = Path("./")
    files = list(p.glob(r"**/*.cpp"))
    files += list(p.glob(r"**/*.H"))
    return files

def check_file(filename):

    # this is a general check to see if we should further examine a line
    if_general_re = re.compile(r"^(?:#if|#elif)", re.IGNORECASE|re.DOTALL)

    # this checks something of the form "#ifdef NAME"
    ifdef_re = re.compile(r"^#if[n]*def\s+([a-z_0-9]+)", re.IGNORECASE|re.DOTALL)

    # this checks something of the form
    # #if (NAME == X)
    if_re = re.compile(r"^(?:#if|#elif)\s+[\(]?([a-z_0-9]+)", re.IGNORECASE|re.DOTALL)

    # together these check something of the form
    # #if defined(NAME1) || !defined(NAME2)
    if_defined_re1 = re.compile(r"^(?:#if|#elif)\s+[!]?(?:defined)", re.IGNORECASE|re.DOTALL)
    if_defined_re2 = re.compile(r"[!]?(?:defined)\s*[\(]?([a-z_0-9]+)[\)]?", re.IGNORECASE|re.DOTALL)

    ierr = 0
    defines = []

    with open(filename) as cf:
        for line in cf:
            if if_general_re.search(line):
                # check each of the patterns
                if g := ifdef_re.search(line):
                    defines.append(g.group(1))
                    continue

                if g := if_re.search(line):
                    defines.append(g.group(1))
                    continue

                if if_defined_re1.search(line):
                    g = if_defined_re2.findall(line)
                    defines += g
                    continue

                # if we made it here, then we didn't handle things
                ierr = 1
                print(f"unhandled, file: {filename} | {line}")

    return ierr, set(defines)

if __name__ == "__main__":
    all_defines = []
    for f in find_source_files():
        ierr, defines = check_file(f)
        all_defines += defines

    # remove and header guards
    defines = []
    for d in all_defines:
        print(f"checking {d}")
        if not (d.endswith("_H") or d.endswith("_H_")):
            defines.append(d)

    for d in sorted(set(defines)):
        print(d)

