#!/usr/bin/python

import sys, getopt
from LowrieEdwardsUnits import LEunits
from phys import R, c
from RadShock import OneTShock

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "L:T:", ["P0","gamma=","sigma=","kappa=", "M0=", "mu="])
    except getopt.GetoptError:
        print '???'
        sys.exit(1)

    P0 = 1.0e-4
    gamma = 5./3.
    sigma = 1.0e6
    kappa = 1.0
    M0 = 2.0
    uL = 1.0e5
    uT = 100.0
    mu = 1.0
    for o, a in opts:
        if o == "--P0":
            P0 = float(a)
        elif o == "--gamma":
            gamma = float(a)
        elif o == "--sigma":
            sigma = float(a)
        elif o == "--kappa":
            kappa = float(a)
        elif o == "--M0":
            M0 = float(a)
        elif o == "-L":
            uL = float(a)
        elif o == "-T":
            uT = float(a)
        elif o == "-mu":
            mu = float(a)

    print "In Lowrie-Edwards units and shock frame: (Here sigma is \kappa_P and kappa is c\lambda/\Chi)"
    print 'P0=', P0, 'gamma=', gamma, 'mu=', mu, 'sigma=', sigma, 'kappa=', kappa, 'M0=', M0
    print ''

    units = LEunits(P0=P0, gamma=gamma, uL=uL, uT=uT, mu=mu)
    L = units['L']
    T = units['T']
    rho = units['rho']
    v = M0 * units['v']
    kappa_p = sigma * units['sigma']
    kappa_r = c * (1./3.) / (kappa * units['kappa'])
    cv = R / (gamma-1.0)/mu

    rho1, v1, T1 = OneTShock(P0=P0, gamma=gamma, M0=M0)
    rho1 = rho1 * units['rho']
    v1 = v1 * units['v']
    T1 = T1 * units['T']

    print "In cgs units:"
    print 'Length unit: ', L
    print 'Pre-shock: ', 'T=', T, 'rho=', rho, 'v=', v, 'v/c=', v/c
    print 'Post-shock: ', 'T=', T1, 'rho=', rho1, 'v=', v1, 'v/c=', v1/c
    print 'kappa_p=', kappa_p, 'kappa_r=', kappa_r
    print 'cv=', cv, 'gamma=', gamma, 'mu=', mu


if __name__ == "__main__":
    main()
