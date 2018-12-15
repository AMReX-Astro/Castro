#!/usr/bin/env python3

from __future__ import print_function

import os
import re
import sys
import textwrap
from more_itertools import unique_everseen

main_header = """
+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
"""

separator = """
+----------------------------------------+---------------------------------------------------------+---------------+
"""

entry = """
| {:38} | {:55} | {:13} |
"""

WRAP_LEN = 55

class Parameter(object):
    # container class for the parameters

    def __init__(self):
        self.var = ""
        self.default = ""
        self.description = []
        self.category = ""
        self.namespace=""

    def value(self):
        """ the value is what we sort based on """
        return self.category + "." + self.var

    def __lt__(self, other):
        self.value() < other.value()


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
        category = ""

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

            if line.startswith("@"):
                # this is a command -- we only know namespace
                fields = line.split()
                if fields[0].startswith("@namespace"):
                    namespace = fields[1].strip()
                    category = ""
                    line = f.readline()
                    continue


            # look for category definition
            if line.startswith("#------"):

                # the next line should be the category definition
                line = f.readline()
                index = line.find(":")
                category = line[index+1:]

                # following this is another #---------
                line = f.readline()
                if not line.startswith("#------"):
                    sys.exit("ERROR: category block not formatted correctly")

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

                # this splits the line into separate fields.  A field
                # is a single word or a pair in parentheses like "(a,
                # b)"
                fields = re.findall(r'[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)', line) 

                current_param.var = fields[0]
                if current_param.var.startswith("("):
                    current_param.var = re.findall(r"\w+", fields[0])[0]

                current_param.default = fields[2].replace("_", "\_")
                current_param.description = descr
                current_param.category = category.strip()
                current_param.namespace = namespace.strip()
                descr = r""

                # store the current parameter in the list
                params_list.append(current_param)

            line = f.readline()


    namespaces = list(unique_everseen([q.namespace for q in params_list]))

    for nm in sorted(namespaces):

        # print the heading
        heading_name = r"namespace: ``{}``".format(nm)
        nmlen = len(heading_name)
        print(heading_name)
        print(nmlen*"-" + "\n")

        # now group by category
        categories = list(unique_everseen([q.category for q in params_list if q.namespace == nm]))

        for c in categories:

            # print the subheading
            if c != "":
                print("**{}**\n".format(c))

            params = [q for q in params_list if q.namespace == nm and q.category == c]

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
    param_files = ["../Source/driver/_cpp_parameters"]

    make_rest_table(param_files)
