#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import textwrap

main_header = """
+----------------------------------+---------------------------------------------------------+---------------+
| parameter                        | description                                             | default value |
+==================================+=========================================================+===============+
"""

separator = """
+----------------------------------+---------------------------------------------------------+---------------+
"""

entry = """
| {:32} | {:55} | {:13} |
"""

WRAP_LEN = 55

class Parameter(object):
    # container class for the parameters

    def __init__(self):
        self.var = ""
        self.default = ""
        self.description = []
        self.category = ""

    def value(self):
        """ the value is what we sort based on """
        return self.category + "." + self.var

    def __cmp__(self, other):
        return cmp(self.value(), other.value())


def make_rest_table(param_files):

    params_list=[]

    for pf in param_files:

        # each file is a category
        category = os.path.basename(os.path.dirname(pf)).replace("_", "\_")

        # open the file
        try: f = open(pf, "r")
        except IOError:
            sys.exit("ERROR: {} does not exist".format(pf))

        descr = r""

        # read in the file
        line = f.readline()
        while line:

            # we assume that parameters have an optional descriptive
            # heading before them without any blank line between the
            # description and the parameter definition.  Therefore,
            # if we encounter a blank line, zero out the description.
            if line.strip() == "":
                descr = r""
                line = f.readline()
                continue

            if line.startswith("#------"):
                line = f.readline()
                continue

            # find the description
            if line.startswith("#"):
                # handle descriptions here
                descr += line[1:].rstrip().replace("@@",r"\newline")
                line = f.readline()
                continue

            else:
                current_param = Parameter()
                line_list = line.split()

                current_param.var = line_list[0]
                current_param.default = line_list[2].replace("_", "\_")
                current_param.description = descr
                current_param.category = category

                descr = r""

                # store the current parameter in the list
                params_list.append(current_param)

            line = f.readline()


    categories = sorted (set([q.category for q in params_list]))

    for c in categories:

        # print the heading

        params = [q for q in params_list if q.category == c]

        clen = len(c)
        print(c)
        print(clen*"=" + "\n")

        print(main_header.strip())

        for p in params:
            desc = list(textwrap.wrap(p.description.strip(), WRAP_LEN))
            if len(desc) == 0:
                desc = [""]

            for n, d in enumerate(desc):
                if n == 0:
                    print(entry.format("``"+p.var+"``", d, p.default).strip())
                else:
                    print(entry.format(" ", d, " ").strip())

            print(separator.strip())

        print("\n\n")

if __name__ == "__main__":

    # find all of the _parameter files
    top_dir = "../"

    param_files = []
    for root, dirs, files in os.walk(top_dir):
        for f in files:
            if f == "_parameters":
                param_files.append(os.path.normpath("/".join([root, f])))

    make_rest_table(param_files)
