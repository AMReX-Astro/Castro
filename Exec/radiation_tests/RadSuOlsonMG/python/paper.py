#!/usr/bin/env python

from pylab import *
from read_gnu import *
import cPickle as pickle

Er0_3, x, t_3, = read_gnu_file('../run-paper/Er0_0031.gnu')
Er1_3, x, t_3, = read_gnu_file('../run-paper/Er1_0031.gnu')
eint_3, x, t_3, = read_gnu_file('../run-paper/eint_0031.gnu')

Er0_30, x, t_30, = read_gnu_file('../run-paper/Er0_0301.gnu')
Er1_30, x, t_30, = read_gnu_file('../run-paper/Er1_0301.gnu')
eint_30, x, t_30, = read_gnu_file('../run-paper/eint_0301.gnu')

fid = open('SOunits.p', 'rb')
units =pickle.load(fid)
fid.close()

U1_3 = Er0_3 / units['Eg']
U2_3 = Er1_3 / units['Eg']
V_3  = eint_3 / units['rhoe']

U1_30 = Er0_30 / units['Eg']
U2_30 = Er1_30 / units['Eg']
V_30  = eint_30 / units['rhoe']

x = x / units['L']
t_3 = t_3 / units['time']
t_30 = t_30 / units['time']

print 't = ', t_3, t_30

ex_3 = array([0.0, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.50, 5.00, 7.50, 10.0, 15.0, 20.0, 25.0])
eU1_3 = array([0.11433, 0.11343, 0.11073, 0.10717, 0.10367, 0.09692, 0.09047, 0.08433, 0.07295, 0.05802, 0.03838, 0.02435, 0.00858, 0.00250, 0.00060])
eU2_3 = array([0.44368, 0.40645, 0.28679, 0.15591, 0.08076, 0.01931, 0.00442, 0.00137, 0.00064, 0.00045, 0.00025, 0.00013, 0.00003, 0.00001, 0.00000])
eV_3  = array([0.67399, 0.61576, 0.42201, 0.21349, 0.10346, 0.02258, 0.00566, 0.00256, 0.00170, 0.00124, 0.00072, 0.00040, 0.00011, 0.00003, 0.00000])

ex_30 = array([0.00, 0.50, 1.00, 1.50, 2.50, 5.00, 6.00, 7.50, 10.0, 17.5, 25.0, 35.0, 50.0, 65.0, 80.0])
eU1_30 = array([0.06128, 0.06123, 0.06108, 0.06083, 0.06013, 0.05785, 0.05684, 0.05529, 0.05257, 0.04372, 0.03447, 0.02297, 0.01025, 0.00359, 0.00099])
eU2_30 = array([0.35004, 0.33536, 0.29535, 0.24013, 0.12910, 0.02108, 0.01550, 0.01319, 0.01141, 0.00735, 0.00455, 0.00225, 0.00067, 0.00017, 0.00003])
eV_30  = array([0.71204, 0.68082, 0.59613, 0.48038, 0.25242, 0.04111, 0.03091, 0.02663, 0.02309, 0.01493, 0.00929, 0.00463, 0.00141, 0.00036, 0.00007])

figure(1, figsize=(7,8))

subplots_adjust(left=0.13, bottom=0.08, right=0.97, top=0.97, wspace=0,hspace=0)

subplot(311)
line1=loglog(ex_3, eU1_3, 'ko', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
line2=loglog(ex_30, eU1_30, 'kD', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
line3 = loglog(x,U1_3, 'k')
line4 = loglog(x,U1_30, 'k--')
xticks(visible=False)
ylabel('$U_1$')
minorticks_on()
xlim(0.2,90)
ylim(4.e-4, 0.2)
legend((line3,line1,line4,line2),
       (r'$\tau = 3$'+' numerical',
        r'$\tau = 3$'+' analytic',
        r'$\tau = 30$'+' numerical',
        r'$\tau = 30$'+' analytic'),
       loc='lower left',
       prop=matplotlib.font_manager.FontProperties(size=16),
       handlelength=2.5,
       labelspacing=0.01)

subplot(312)
loglog(ex_3, eU2_3, 'ko', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
loglog(ex_30, eU2_30, 'kD', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
loglog(x,U2_3, 'k')
loglog(x,U2_30, 'k--')
xticks(visible=False)
ylabel('$U_2$')
minorticks_on()
xlim(0.2,90)
ylim(3.e-6, 0.9)

subplot(313)
loglog(ex_3, eV_3, 'ko', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
loglog(ex_30, eV_30, 'kD', markeredgecolor='k', markerfacecolor='w', markeredgewidth=2, markersize=8)
loglog(x,V_3, 'k')
loglog(x,V_30, 'k--')
xlabel('$x$')
ylabel('$V$')
minorticks_on()
xlim(0.2,90)
ylim(5.e-5, 3.)

draw()
show()

