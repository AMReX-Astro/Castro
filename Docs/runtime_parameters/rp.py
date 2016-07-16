#!/usr/bin/env python

from __future__ import print_function

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

param_file = "../../Source/_cpp_parameters"


class Parameter(object):
    # container class for the parameters

    def __init__(self):
        self.var=""
        self.default=""
        self.description=[]
        self.category=""

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
            line_list = line.split()

            current_param.var = line_list[0]
            current_param.default = line_list[2].replace("_", "\_")
            current_param.description = descr
            current_param.category = category

            descr = r""

        # store the current parameter in the list
        params_list.append(current_param)
                
        line = f.readline()

    
    # dump the main header
    print(Mheader)

    # sort the parameters and dump them in latex-fashion.  Group
    # things by category
    current_category = ""
    start = 1

    for param in sorted(params_list):

        if not param.category == current_category:
            if not start == 1:
                print(footer)

            current_category = param.category
            odd = 1
            cat_header = header.replace("@@catname@@", param.category + " parameters.")
            print(cat_header)
            start = 0

        if odd == 1:
            print(r"\rowcolor{tableShade}")
            odd = 0
        else:
            odd = 1

        print(r"\verb= {} = & {} & {} \\".format(
            param.var, param.description, param.default))

    # dump the footer
    print(footer)
    print(Mfooter)

if __name__ == "__main__":
    make_tex_table()
