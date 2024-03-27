#!/usr/bin/python

import sys, getopt
import numpy as np
from phys import c, a
import cPickle as pickle

def SOunits(eps, printmessg=0):

    kap = 1.0
    T0 = 1.0e6

    alpha = 4.*a/eps

    ul = 1.0/kap
    utime = 1./(eps*c*kap)
    utemp = T0
    urhoe = alpha * T0**4 / 4.0
    uEr = a * T0**4

    units = {'L':ul, 'time':utime, 'Temp': utemp, 'Eg':uEr, 'rhoe':urhoe, 'alpha': alpha}
    if printmessg:
        print 'units:'
        print units

    fid = open("SOunits.p", 'w')
    pickle.dump(units,fid)
    fid.close()

    return units


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["eps="])
    except getopt.GetoptError:
        print '???'
        sys.exit(1)

    eps = 1.0

    for o, a in opts:
        if o == "--eps":
            eps = float(a)

    SOunits(eps, 1)

