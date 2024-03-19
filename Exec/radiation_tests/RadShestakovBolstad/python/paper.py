#!/usr/bin/env python

from pylab import *
from read_gnu import *
import cPickle as pickle

n = '0200'
ew, x, t = read_gnu_file('../run-paper/eint_'+n+'.gnu')
Ew, x, t = read_gnu_file('../run-paper/Er_'+n+'.gnu')

fid = open('SBunits.p', 'rb')
units = pickle.load(fid)
fid.close()

Tw = (ew/units['cv']) / units['Temp']
Ew = Ew / units['Eg']
x = x / units['L']
t = t / units['time']

print 't = ', t

ex = array([0.0, 0.2, 0.4, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.6, 0.8, 1.0])
eT = array([9.9373253e-1, 9.9339523e-1, 9.8969664e-1, 9.8060848e-1, 9.7609654e-1, 9.6819424e-1,
            9.5044751e-1, 4.9704000e-1, 4.3632445e-2, 2.5885608e-2, 1.7983134e-2, 1.3470947e-2,
            4.3797848e-3, 6.4654865e-4, 1.9181546e-4])
eE = array([5.6401674e-3, 5.5646351e-3, 5.1047352e-3, 4.5542134e-3, 4.3744933e-3, 4.1294850e-3,
            3.7570008e-3, 2.9096931e-3, 2.0623647e-3, 1.6898183e-3, 1.4447063e-3, 1.2648409e-3,
            7.1255738e-4, 2.3412650e-4, 1.0934921e-4])

figure(1, figsize=(7,6.5))
subplots_adjust(left=0.15, bottom=0.08, right=0.97, top=0.97, wspace=0,hspace=0)

ax1 = subplot(211)
plot(ex, eE, 'o', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
plot(x,Ew,'k')
xticks(visible=False)
xlim(0,1.05)
#ylim(0.0,0.14)
ylabel('$E_r$')
minorticks_on()
text(0.9,0.86,'(a)', transform = ax1.transAxes)

ax2 = subplot(212)
plot(ex,eT,'o', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
plot(x,Tw,'k')
legend(['Exact', 'Numerical'], loc='lower left')
xlabel('x')
xlim(0,1.05)
ylim(-0.04,1.14)
ylabel('$T$')
minorticks_on()
text(0.9,0.86,'(b)', transform = ax2.transAxes)

draw()
show()
