#!/usr/bin/python

from numpy import *

def read_group_struct(filenm):
    f = open(filenm, 'r')

    line = f.readline()
    ngroups = int(line.split('=')[1])

    line = f.readline()

    nu = []
    dnu = []
    for i in range(ngroups):
        line = f.readline()
        words = line.split()
        nu.append(float(words[0]))
        dnu.append(float(words[1]))

    line = f.readline()
    line = f.readline()

    xnu = []
    for i in range(ngroups+1):
        line = f.readline()
        words = line.split('=')
        xnu.append(float(words[1]))

    f.close()
    return ngroups, array(nu), array(dnu), array(xnu)

def read_neut_groups(filenm):
    f = open(filenm, 'r')

    line = f.readline()
    ngroups = int(line.split('=')[1])

    line = f.readline()
    words = line.split(':')[1].split(',')
    ng0 = int(words[0])
    ng1 = int(words[1])
    ng2 = int(words[2])

    line = f.readline()

    nu = []
    dnu = []
    for i in range(ngroups):
        line = f.readline()
        words = line.split()
        nu.append(float(words[0]))
        dnu.append(float(words[1]))

    f.close()
    return ng0, ng1, ng2, array(nu), array(dnu)


if __name__ == "__main__":
    ngroups, nu, dnu, xnu = read_group_struct('../run-mg/group_structure.dat')
    print ngroups
    print nu
    print dnu
    print xnu
