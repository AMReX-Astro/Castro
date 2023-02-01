#!/bin/env python

import re
import sys

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
                if if_defined_re1.search(line):
                    g = if_defined_re2.findall(line)
                    defines += g
                    continue

                if g := ifdef_re.search(line):
                    defines.append(g.group(1))
                    continue

                if g := if_re.search(line):
                    defines.append(g.group(1))
                    continue


                # if we made it here, then we didn't handle things
                ierr = 1
                print(f"unhandled, file: {filename} | {line}")

    return ierr, set(defines)

if __name__ == "__main__":

    good_defines_file = sys.argv[1]

    # read in the list of good defines

    good_defines = []
    with open(good_defines_file) as gd:
        for line in gd:
            good_defines.append(line.strip())

    all_defines = []
    total_errors = 0
    for f in find_source_files():
        if "tmp_build_dir" in f.parts:
            # skip generated files
            continue
        ierr, defines = check_file(f)
        all_defines += defines
        total_errors += ierr

    # remove any header guards

    defines = []
    for d in all_defines:
        if d.endswith("_H") or d.endswith("_H_"):
            continue
        if len(d) == 1 and d.isdigit():
            continue
        defines.append(d)

    defines = sorted(set(defines))

    print("found defines:")
    for d in defines:
        print(d)

    # now check to make sure that all the defines we found are okay

    invalid = []
    for d in defines:
        if d not in good_defines:
            invalid.append(d)

    if invalid or total_errors > 0:
        if invalid:
            print("\ninvalid defines:")
        for bad in invalid:
            print(bad)
        sys.exit(1)
