#!/usr/bin/python

import cPickle as pickle
from pylab import *
from read_gnu import *
import phys

Er, x, t = read_gnu_file('../run-M5/Er_1735.gnu')
rho, x, t = read_gnu_file('../run-M5/dens_1735.gnu')
v, x, t = read_gnu_file('../run-M5/vx_1735.gnu')
Temp, x, t = read_gnu_file('../run-M5/Temp_1735.gnu')
print t

x -= 208.

fid = open('M5shock.p', 'rb')
xa = pickle.load(fid)
Ta = pickle.load(fid)
thetaa = pickle.load(fid)
rhoa = pickle.load(fid)
va = pickle.load(fid)
parsa = pickle.load(fid)
fid.close()

fid = open('M5shock_hi.p', 'rb')
xhi = pickle.load(fid)
Thi = pickle.load(fid)
thetahi = pickle.load(fid)
rhohi = pickle.load(fid)
vhi = pickle.load(fid)
parshi = pickle.load(fid)
fid.close()

fid = open('units-M5shock.p', 'rb')
units =pickle.load(fid)
fid.close()

xa = xa * units['L']
xa = append(xa, [2000.])
Ta = Ta * units['T']
Ta = append(Ta, Ta[-2:-1])
thetaa = thetaa * units['T']
thetaa = append(thetaa, thetaa[-2:-1])
rhoa = rhoa * units['rho']
rhoa = append(rhoa, rhoa[-2:-1])

xhi = xhi * units['L']
Thi = Thi * units['T']
thetahi = thetahi * units['T']
rhohi = rhohi * units['rho']

figure(1, figsize=(7,8))
subplots_adjust(left=0.14, bottom=0.08, right=0.97, top=0.97,
                hspace=0)

subplot(211)
plot(x,rho/1.e-12, 'ko',markerfacecolor='none')
#plot(xa,rhoa/1.e-12, 'r',lw=1.8)
plot(xa,rhoa/1.e-12, 'k',lw=2.4)
xlim(-3000,900)
xticks(visible=False)
ylabel(r'$\rho\ (10^{-12} \mathrm{g\,cm}^{-3})$')
minorticks_on()

subplot(212)
plot(x,Temp, 'k+', markersize=10.)
plot(x,(Er/phys.a)**0.25, 'ks', markerfacecolor='none')
#plot(xa,Ta,'b',lw=1.8)
plot(xa,Ta,'k-',lw=2.4)
#plot(xa,thetaa, 'g',lw=1.8)
plot(xa,thetaa, 'k--',lw=2.4, dashes=(20,10))
xlim(-3000,900)
ylim(0,1180)
xlabel(r'$x\ (\mathrm{cm})$')
ylabel(r'$T\ (\mathrm{K})$')
minorticks_on()

ax = axes([0.6, 0.13, 0.345, 0.19])
plot(x,Temp, 'k+', markersize=10.,markeredgewidth=2.)
plot(x,(Er/phys.a)**0.25, 'ks', markersize=8., markerfacecolor='none',markeredgewidth=1.5)
#plot(xhi,Thi,'b',lw=1.8)
plot(xhi,Thi,'k-',lw=2.4)
#plot(xhi,thetahi, 'g',lw=1.8)
plot(xhi,thetahi, 'k--',lw=2.4, dashes=(20,10))
xlim(-10,15)
ylim(810,1140)

draw()
show()

