#!/usr/bin/python

import sys, getopt, phys
from numpy import *
from pylab import *

def main():
    kappa = 1.0
    eps = 0.1
    Finc = 1.0

    xu001 = [0.1, 0.25, 0.5, 0.75, 1.0]
    u001 = [0.17979, 0.11006, 0.04104, 0.01214, 0.00268]
    xv001 = [0.1, 0.25, 0.5, 0.75]
    v001 = [0.00110, 0.00055, 0.00012, 0.00003]

    x03 = [0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5., 7.5]
    u03 = [0.44289, 0.38544, 0.30500, 0.24062, 0.18922, 0.04167, 0.00238, 0.00008]
    v03 = [0.10124, 0.08551, 0.06437, 0.04830, 0.03612, 0.00584, 0.00020, 0.00001]

    Erfilea = '../run-paper/Er_0035.gnu'
    Tfilea = '../run-paper/Temp_0035.gnu'
    Erfileb = '../run-paper/Er_1002.gnu'
    Tfileb = '../run-paper/Temp_1002.gnu'

    def read_gnu_file(filenm):
        x = []
        y = []
        f = open(filenm, 'r')
        line = f.readline()
        t = float(line.split('"')[1].split('=')[2])
        for line in f.readlines():
            if not line[0] == ";":
                words = line.split()
                x.append(float(words[0]))
                y.append(float(words[1]))
        f.close()
        return array(y), array(x), t

    Ta, za, ta = read_gnu_file(Tfilea)
    Era, za, ta = read_gnu_file(Erfilea)
    ua = Era * phys.c / (4.*Finc)
    va = phys.c/4. * (phys.a * Ta**4) / Finc
    xa = sqrt(3.) * kappa * za
    taua = ta * phys.c * kappa * eps
    print 'taua=', taua

    Tb, zb, tb = read_gnu_file(Tfileb)
    Erb, zb, tb = read_gnu_file(Erfileb)
    ub = Erb * phys.c / (4.*Finc)
    vb = phys.c/4. * (phys.a * Tb**4) / Finc
    xb = sqrt(3.) * kappa * zb
    taub = tb * phys.c * kappa * eps
    print 'taub=', taub

    loglog(xa, ua, '-k', label=r'$u(x,0.01)$')
    loglog(xa, va, '--k', label=r'$v(x,0.01)$')
    loglog(xb, ub, '-.k', label=r'$u(x,0.3)$')
    loglog(xb, vb, ':k', label=r'$v(x,0.3)$')

    legend(loc='best', prop=matplotlib.font_manager.FontProperties(size=18),
           labelspacing=0.01)

    loglog(xu001, u001, 'o', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
    loglog(xv001, v001, 'o', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
    loglog(x03, u03, 'o', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
    loglog(x03, v03, 'o', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)

    xlim(0.03, 10.0)
    ylim(0.3e-5,1.0)
    xlabel('x')
    draw()
    show()


if __name__ == "__main__":
    main()
