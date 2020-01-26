#!/usr/bin/env python3

import sys
import argparse

class Species:
    """the species class holds the properties of a single species"""
    def __init__(self):
        self.name = ""
        self.short_name = ""
        self.A = -1
        self.Z = -1

    def __str__(self):
        return "species {}, (A,Z) = {},{}".format(self.name, self.A, self.Z)

class AuxVar:
    """convenience class for an auxilliary variable"""
    def __init__(self):
        self.name = ""

    def __str__(self):
        return "auxillary variable {}".format(self.name)


def get_next_line(fin):
    """get_next_line returns the next, non-blank line, with comments
    stripped"""
    line = fin.readline()

    pos = str.find(line, "#")

    while (pos == 0 or str.strip(line) == "") and line:
        line = fin.readline()
        pos = str.find(line, "#")

    line = line[:pos]

    return line


def get_object_index(objs, name):
    """look through the list and returns the index corresponding to the
    network object (species or auxvar) specified by name

    """

    index = -1

    for n, o in enumerate(objs):
        if o.name == name:
            index = n
            break

    return index


def parse_net_file(species, aux_vars, net_file):
    """parse_net_file read all the species listed in a given network
    inputs file and adds the valid species to the species list

    """

    err = 0

    try:
        f = open(net_file, "r")
    except IOError:
        sys.exit("write_network.py: ERROR: file "+str(net_file)+" does not exist")


    print("write_network.py: working on network file "+str(net_file)+"...")

    line = get_next_line(f)

    while line and not err:

        fields = line.split()

        # read the species or auxiliary variable from the line
        net_obj, err = parse_network_object(fields)
        if net_obj is None:
            return err

        objs = species
        if isinstance(net_obj, AuxVar):
            objs = aux_vars

        # check to see if this species/auxvar is defined in the current list
        index = get_object_index(objs, net_obj.name)

        if index >= 0:
            print("write_network.py: ERROR: {} already defined.".format(net_obj))
            err = 1
        # add the species or auxvar to the appropriate list
        objs.append(net_obj)

        line = get_next_line(f)

    return err


def parse_network_object(fields):
    """parse the fields in a line of the network file for either species
    or auxiliary variables.  Aux variables are prefixed by '__aux_' in
    the network file

    """

    err = 0

    # check for aux variables first
    if fields[0].startswith("__aux_"):
        ret = AuxVar()
        ret.name = fields[0][6:]

    # check for missing fields in species definition
    elif not len(fields) == 4:
        print(" ".join(fields))
        print("write_network.py: ERROR: missing one or more fields in species definition.")
        ret = None
        err = 1
    else:
        ret = Species()

        ret.name = fields[0]
        ret.short_name = fields[1]
        ret.A = fields[2]
        ret.Z = fields[3]

    return ret, err


def abort(outfile):
    """exit when there is an error.  A dummy stub file is written out,
    which will cause a compilation failure

    """

    fout = open(outfile, "w")
    fout.write("There was an error parsing the network files")
    fout.close()
    sys.exit(1)



def write_network(network_template, header_template,
                  net_file,
                  network_file, header_file):
    """read through the list of species and output the new out_file

    """

    species = []
    aux_vars = []


    #-------------------------------------------------------------------------
    # read the species defined in the net_file
    #-------------------------------------------------------------------------
    err = parse_net_file(species, aux_vars, net_file)

    if err:
        abort(out_file)


    #-------------------------------------------------------------------------
    # write out the Fortran and C++ files based on the templates
    #-------------------------------------------------------------------------
    templates = [(network_template, network_file), (header_template, header_file)]

    for tmp, out_file in templates:

        print("writing {}".format(out_file))

        # read the template
        try:
            template = open(tmp, "r")
        except IOError:
            sys.exit("write_network.py: ERROR: file {} does not exist".format(tmp))
        else:
            template_lines = template.readlines()
            template.close()

        # output the new file, inserting the species info in between the @@...@@
        fout = open(out_file, "w")

        for line in template_lines:

            index = line.find("@@")

            if index >= 0:
                index2 = line.rfind("@@")

                keyword = line[index+len("@@"):index2]
                indent = index*" "

                if keyword == "NSPEC":
                    fout.write(line.replace("@@NSPEC@@", str(len(species))))

                elif keyword == "NAUX":
                    fout.write(line.replace("@@NAUX@@", str(len(aux_vars))))

                elif keyword == "SPEC_NAMES":
                    for n, spec in enumerate(species):
                        fout.write("{}spec_names({}) = \"{}\"\n".format(indent, n+1, spec.name))

                elif keyword == "SHORT_SPEC_NAMES":
                    for n, spec in enumerate(species):
                        fout.write("{}short_spec_names({}) = \"{}\"\n".format(indent, n+1, spec.short_name))

                elif keyword == "AION":
                    for n, spec in enumerate(species):
                        fout.write("{}aion({}) = {}\n".format(indent, n+1, spec.A))

                elif keyword == "ZION":
                    for n, spec in enumerate(species):
                        fout.write("{}zion({}) = {}\n".format(indent, n+1, spec.Z))

                elif keyword == "AUX_NAMES":
                    for n, aux in enumerate(aux_vars):
                        fout.write("{}aux_names({}) = \"{}\"\n".format(indent, n+1, aux.name))
                        fout.write("{}short_aux_names({}) = \"{}\"\n".format(indent, n+1, aux.name))

            else:
                fout.write(line)

        print(" ")
        fout.close()


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", type=str, default="",
                        help="fortran template for the network")
    parser.add_argument("-o", type=str, default="",
                        help="fortran module output file name")
    parser.add_argument("--header_template", type=str, default="",
                        help="C++ header template file name")
    parser.add_argument("--header_output", type=str, default="",
                        help="C++ header output file name")
    parser.add_argument("-s", type=str, default="",
                        help="network file name")

    args = parser.parse_args()

    if args.t == "" or args.o == "":
        sys.exit("write_probin.py: ERROR: invalid calling sequence")

    write_network(args.t, args.header_template,
                  args.s,
                  args.o, args.header_output)

if __name__ == "__main__":
    main()
