#!/usr/bin/python

import sys, getopt
import numpy as np
from phys import c, h, k, eV
import cPickle as pickle

def SBunits(T0, rho0, kap0, printmessg=0):

    B0 = 8.*np.pi*h/c**3
    nu0 = k*T0/h
    l0 = nu0**3/kap0
    x0 = l0/np.sqrt(3.0)
    t0 = l0/c
    u0 = B0*nu0**3
    E0 = u0*nu0
    cv = k*u0/(h*rho0)

    units = {'L':x0, 'time':t0, 'Temp': T0, 'Enu':u0, 'nu':nu0, 'Eg':E0, 'cv': cv}
    if printmessg:
        print 'units:'
        print units

    fid = open("SBunits.p", 'w')
    pickle.dump(units,fid)
    fid.close()

    return units


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["rho0=", "kap0=", "T0="])
    except getopt.GetoptError:
        print '???'
        sys.exit(1)

    rho0 = 1.8212111e-5
    kap0 = 4.0628337e43
    T0 = 0.1  # keV

    for o, a in opts:
        if o == "--T0":
            T0 = float(a)
        elif o == "--rho0":
            rho0 = float(a)
        elif o == "--kap0":
            kap0 = float(a)

    T0 = T0 * 1.e3 * eV / k  # now in erg

    SBunits(T0, rho0, kap0, 1)

