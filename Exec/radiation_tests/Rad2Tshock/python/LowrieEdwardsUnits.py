#!/usr/bin/python

import sys, getopt
import numpy as np
from phys import c, R
from phys import a as aRad
import cPickle as pickle

def LEunits(P0=1.0e-4, gamma=5./3., uL=1.0e5, uT=100.0, mu=1.0, printmessg=0):

    L0 = uL
    T0 = uT
    a0 = np.sqrt(gamma*R*T0/mu)
    rho0 = aRad * T0**4 / (P0 * a0**2)

    ut = L0/a0
    urho = rho0
    uv = a0
    ue = a0**2
    up = rho0*a0**2
    usigma = a0/(L0*c)
    ukappa = a0*L0

    units = {'L':uL, 't':ut, 'rho':urho, 'v':uv, 'e':ue, 'p':up, 'T':uT, 'sigma':usigma, 'kappa':ukappa}
    if printmessg:
        print 'units for'
        print 'uL=', uL
        print 'uT=', uT
        print 'gamma=', gamma
        print 'P0=', P0
        print 'mu=', mu
        print units

    fid = open("units-shock.p", 'w')
    pickle.dump(units,fid)
    fid.close()

    return units


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "L:T:", ["gamma=", "P0=", "mu="])
    except getopt.GetoptError:
        print '???'
        sys.exit(1)

    mu = 1.0
    P0 = 1.0e-4
    gamma = 5./3.
    uL = 1.0e5
    uT = 100.0

    for o, a in opts:
        if o == "-L":
            uL = float(a)
        elif o == "-T":
            uT = float(a)
        elif o == "--gamma":
            gamma = float(a)
        elif o == "--P0":
            P0 = float(a)
        elif o == "--mu":
            mu = float(a)

    LEunits(P0=P0, gamma=gamma, uL=uL, uT=uT, mu=mu, printmessg=1)

