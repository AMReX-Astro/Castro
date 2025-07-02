#!/usr/bin/python

import cPickle as pickle
from pylab import *
from read_gnu import *
import phys

n = '2000'
Er, x, t = read_gnu_file('../run-M2-mg/Er_'+n+'.gnu')
rho, x, t = read_gnu_file('../run-M2-mg/dens_'+n+'.gnu')
v, x, t = read_gnu_file('../run-M2-mg/vx_'+n+'.gnu')
Temp, x, t = read_gnu_file('../run-M2-mg/Temp_'+n+'.gnu')
print t

x += 10

fid = open('M2shock.p', 'rb')
xa = pickle.load(fid)
Ta = pickle.load(fid)
thetaa = pickle.load(fid)
rhoa = pickle.load(fid)
va = pickle.load(fid)
parsa = pickle.load(fid)
fid.close()

fid = open('units-M2shock.p', 'rb')
units =pickle.load(fid)
fid.close()

xa = xa * units['L']
Ta = Ta * units['T']
thetaa = thetaa * units['T']
rhoa = rhoa * units['rho']

figure(1, figsize=(7,8))
subplots_adjust(left=0.13, bottom=0.08, right=0.96, top=0.97,
                hspace=0)

subplot(211)
plot(xa,rhoa/1.e-12, 'k',lw=1.8)
plot(x,rho/1.e-12, 'ko',markerfacecolor='none')
xlim(-800,400)
ylim(0.5001, 1.3)
xticks(visible=False)
ylabel(r'$\rho\ (10^{-12} \mathrm{g\,cm}^{-3})$')
minorticks_on()

subplot(212)
plot(xa,Ta,'k',lw=1.5)
plot(xa,thetaa, 'k',lw=1.8)
plot(x,Temp, 'k+', markersize=9.)
plot(x,(Er/phys.a)**0.25, 'ks', markerfacecolor='none')
xlim(-800,400)
ylim(90,230)
xlabel(r'$x\ (\mathrm{cm})$')
ylabel(r'$T\ (\mathrm{K})$')
minorticks_on()


draw()
show()


