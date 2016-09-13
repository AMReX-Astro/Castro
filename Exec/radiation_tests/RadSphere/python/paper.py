#!/usr/bin/python

from pylab import *
from RadGroup import read_group_struct
from phys import h, eV

d='../run-paper/'
f = open(d+'castro.out', 'r')
for i in range(4):
    line = f.readline()
nu_c = []
Enu_c = []
for line in f.readlines():
    words = line.split()
    nu_c.append(float(words[1]))
    Enu_c.append(float(words[2]))
f.close()

f = open(d+'analytic.out', 'r')
line = f.readline()
nu_e = []
Enu_e = []
for line in f.readlines():
    words = line.split()
    nu_e.append(float(words[1]))
    Enu_e.append(float(words[2]))
f.close()

ngroups, nu, dnu, xnu = read_group_struct(d+'group_structure.dat')
Enu_e = array(Enu_e) / dnu  
Enu_c = array(Enu_c) / dnu  
nu_e = array(nu_e) 
nu_c = array(nu_c) 

figure(1, figsize=(7,6.5))
subplots_adjust(left=0.15, bottom=0.12, right=0.97, top=0.97, wspace=0,hspace=0)

loglog(nu_c, Enu_c, 'o', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
loglog(nu_e, Enu_e, 'k')

xlim(1.0/(h/eV), 10**5.2/(h/eV) )
ylim(10.**-26.9*(h/eV), 10.**12.9 *(h/eV) )

xlabel(r'$\nu\ (\mathrm{Hz})$')
ylabel(r'$E_\nu\ (\mathrm{erg}\,\mathrm{cm}^{-3}\,\mathrm{Hz}^{-1})$')

legend(('numerical', 'analytic'), loc='lower left')

draw()
show()
