#!/usr/bin/env python

from __future__ import print_function

import re
import sys

# tex format stuff
Mheader = r"""
\label{ch:parameters}


%%%%%%%%%%%%%%%%
% symbol table
%%%%%%%%%%%%%%%%

\begin{landscape}
"""

header = r"""
{\small

\renewcommand{\arraystretch}{1.5}
%
\begin{center}
\begin{longtable}{|l|p{5.25in}|l|}
\caption[@@catname@@]{@@catname@@} \label{table: @@catname@@ runtime} \\
%
\hline \multicolumn{1}{|c|}{\textbf{parameter}} & 
       \multicolumn{1}{ c|}{\textbf{description}} & 
       \multicolumn{1}{ c|}{\textbf{default value}} \\ \hline 
\endfirsthead

\multicolumn{3}{c}%
{{\tablename\ \thetable{}---continued}} \\
\hline \multicolumn{1}{|c|}{\textbf{parameter}} & 
       \multicolumn{1}{ c|}{\textbf{description}} & 
       \multicolumn{1}{ c|}{\textbf{default value}} \\ \hline 
\endhead

\multicolumn{3}{|r|}{{\em continued on next page}} \\ \hline
\endfoot

\hline 
\endlastfoot

"""

footer = r"""

\end{longtable}
\end{center}

} % ends \small
"""

Mfooter = r"""
\end{landscape}

%

"""

param_file = "../../Source/driver/_cpp_parameters"


class Parameter(object):
    # container class for the parameters

    def __init__(self):
        self.var=""
        self.default=""
        self.description=[]
        self.category=""
        self.namespace=""

    def value(self):
        """ the value is what we sort based on """
        return self.category + "." + self.var

    def __cmp__(self, other):
        return cmp(self.value(), other.value())


def make_tex_table():

    # open the file
    try: f = open(param_file, "r")
    except IOError:
        sys.exit("ERROR: {} does not exist".format(param_file))

    # local storage for the parameters
    params_list = []
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

            else:
                print(fields[0])
                sys.exit("invalid command in parameter file")
                

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

            # this splits the line into separate fields.  A field is a                              
            # single word or a pair in parentheses like "(a, b)"                                    
            fields = re.findall(r'[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)', line) 

            current_param.var = fields[0]
            if current_param.var.startswith("("):
                current_param.var = re.findall(r"\w+", fields[0])[0]

            current_param.default = fields[2].replace("_", "\_")
            current_param.description = descr
            current_param.category = category
            current_param.namespace = namespace
            descr = r""

        # store the current parameter in the list
        params_list.append(current_param)
                
        line = f.readline()


    namespaces = list(set([q.namespace for q in params_list]))

    for nm in sorted(namespaces):
  
        print(r"\section{{ {{\tt {} }} Namespace}}".format(nm))

        # dump the main header
        print(Mheader)

        # sort the parameters and dump them in latex-fashion.  Group
        # things by category
        current_category = -1
        start = 1

        nm_params = [q for q in params_list if q.namespace == nm]

        for param in sorted(nm_params):

            if not param.category == current_category:
                if not start == 1:
                    print(footer)

                current_category = param.category
                odd = 1
                if param.category != "":
                    cat_header = header.replace("@@catname@@", "{} : {} parameters".format(
                        param.namespace, param.category))
                else:
                    cat_header = header.replace("@@catname@@", "{} parameters".format(
                        param.namespace))

                print(cat_header)
                start = 0

            if odd == 1:
                print(r"\rowcolor{tableShade}")
                odd = 0
            else:
                odd = 1

            print(r"\runparamNS{{{}}}{{{}}} & {} & {} \\".format(
                param.var.replace("_", r"\_"), param.namespace,
                param.description, param.default))

        # dump the footer
        print(footer)
        print(Mfooter)

if __name__ == "__main__":
    make_tex_table()
