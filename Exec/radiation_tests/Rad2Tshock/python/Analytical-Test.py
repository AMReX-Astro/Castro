#!/usr/bin/python
"""
This is a script to run the analytical solution of Rad2TShock. This will take that solution and
with a user inputted plotfile, will produce graphs of the solutions to density and temperature
vs the plot files.
"""
import argparse
import getopt
import os
import sys

#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import yt
from scipy import integrate, optimize

plt.switch_backend('agg')

parser = argparse.ArgumentParser()
parser.add_argument('plotfile', type=str, help='Path to plotfile for which we will compare')
args = parser.parse_args()

def get_x_var(pf, var):
    ds = yt.load(pf)
    time = float(ds.current_time)
    ad = ds.all_data()
    srt = np.argsort(ad['x'])
    x_coord = np.array(ad['x'][srt])
    var = np.array(ad[var][srt])
    return time, x_coord, var


def OneTShock(P0, gamma, M0):  # Lowrie & Rauenzahn 07
    def f1_f(T1):
        return 3.*(gamma+1.)*(T1-1.) - P0*gamma*(gamma-1.)*(7.+T1**4)
    def f2_f(T1):
        return 12.*(gamma-1.)**2*T1*(3.+gamma*P0*(1.+7.*T1**4))
    def rho1_f(T1):
        f1 = f1_f(T1)
        return (f1+np.sqrt(f1**2+f2_f(T1))) / (6.*(gamma-1.)*T1)
    def Eq13(T1):
        rho1 = rho1_f(T1)
        return 3.0*rho1*(rho1*T1-1.) + gamma*P0*rho1*(T1**4-1.) - 3.*gamma*(rho1-1.)*M0**2
    def T1_smallP0():
        return (1.-gamma+2.*gamma*M0**2) * (2.+(gamma-1.)*M0**2) / ((gamma+1.)*M0)**2
    def T1_bigP0():
        M0s = M0 / np.sqrt(4.*P0/9.)
        return (8.*M0s**2-1.)/7.
    T1min = T1_smallP0()
    T1max = T1_bigP0()
    if (T1max < T1min):
        print("T1min=", T1min, "T1max=", T1max)
        sys.exit(1)
    while Eq13(T1min) > 0.0:
        T1min *= 0.99
    while Eq13(T1max) < 0.0:
        T1max *= 1.01
    T1 = optimize.brentq(Eq13, T1min, T1max)
    rho1 = rho1_f(T1)
    v1 = M0 / rho1
    return rho1, v1, T1

# references: Lowrie & Edwards 2008,  Lowrie & Rauenzahn 2007
def TwoTShock(P0=1.0e-4, gamma=5./3., sigma=1.0e6, kappa=1., M0=5):

    def peos(rho, T):
        return rho*T/gamma
    def eeos(T):
        return T/(gamma*(gamma-1.))

    Cp = 1./(gamma-1.)
    M_ISP = 1.0/np.sqrt(gamma)  # Eq. 28
    Km = 3.*(gamma*M0**2+1.) + gamma*P0 # Eq. 16

    def rho_f(T, theta, s): # Eq. 17
        b = Km - gamma*P0*theta**4
        ac4 = 36.*gamma*M0**2*T
        return (b + s*np.sqrt(b**2-ac4))/(6.*T)

    def dTdtheta(T, theta, rho, M, v, s):  # Eq. 54
        b = Km - gamma*P0*theta**4
        d = np.sqrt(b**2 - 36.*gamma*M0**2*T)
        drhodT = - (rho + s*(3.*gamma*M0**2)/d) / T
        drhodtheta = -(2./3.)*P0*gamma*theta**3/T*(1.+s*b/d)
        c1 = M0/(24.*P0*kappa*rho**2*theta**3)
        c2 = P0/(3.*Cp*M0*(M**2-1.))
        dGdT = c1 * (6.*Cp*rho*(2.*drhodT*(T-1.)+rho) - 6.*M0**2*rho*drhodT + 8.*P0*drhodT*(theta**4-2.*rho))
        dGdtheta = c1 * (12.*Cp*drhodtheta*rho*(T-1.) - 6.*M0**2*rho*drhodtheta + 8.*P0*(drhodtheta*(theta**4-2.*rho)+4.*rho*theta**3))
        dFdT = c2 * (4.*v*theta**3*dGdT - 12.*sigma*(gamma*M**2-1.)*T**3)
        dFdtheta = c2 * (4.*v*theta**3*dGdtheta + 12.*sigma*(gamma*M**2-1.)*theta**3)
        S1 = dFdT - dGdtheta
        S2 = np.sqrt((dFdT-dGdtheta)**2 + 4.*dGdT*dFdtheta)
        dTdtheta_tmp = (S1 + S2) / (2.*dGdT)
        drhodtheta_full = drhodtheta + drhodT * dTdtheta_tmp
        if (drhodtheta > 0.0):
            return dTdtheta_tmp
        else:
            return (S1 - S2) / (2.*dGdT)

    def thetaprime_f(T, theta, rho, v): # Eq.18
        return v*(6.*Cp*rho*(T-1.)+3.*rho*(v**2-M0**2)+8.*P0*(theta**4-rho))/(24.*P0*kappa*theta**3)

    def dxTdM(xT, M): # Eq. 37
        T = xT[1]
        rho = M0 / (M * np.sqrt(T))  # Eq. 40
        v = M0/rho
        th4 = (Km-3.*gamma*M0**2/rho-3.*T*rho) / (gamma*P0)  # Eq. 41
        r = 3.*rho*sigma*(th4-T**4) # Eq. 25
        theta = th4**0.25
        th3 = th4/theta
        thetaprime = thetaprime_f(T,theta,rho,v)
        foo = 4.*M0*th3*thetaprime
        ZD = foo + (gamma-1.)/(gamma+1.)*(gamma*M**2+1.)*r  # Eq. 39
        ZN = foo + (gamma*M**2-1.)*r  # Eq. 24
        dxdM =  -6.*M0*rho*T/((gamma+1.)*P0*M) * (M**2-1.)/ZD
        dTdM = -2.*(gamma-1.)/(gamma+1.) * T/M * ZN/ZD
        return [dxdM, dTdM]

    # Step 1 of LE08
    T0 = 1.0
    theta0 = 1.0
    v0 = M0
    rho0 = 1.0

    rho1, v1, T1 = OneTShock(P0=P0, gamma=gamma, M0=M0)
    M1 = v1 / np.sqrt(T1)
    theta1 = T1
    print('THIS IS A TEST rho1, v1, theta1, T1=', rho1, v1, theta1, T1)
    # which branch?
    s0 = -1.0
    if M1 < M_ISP:
        s1 = 1.
    else:
        s1 = -1.

    # Step 2a of LE08
    # i= 0
    eps = 1.0e-10
    theta0_eps = theta0 + eps
    T0_eps = T0 + eps * dTdtheta(T0, theta0, rho0, M0, v0, s0)
    rho0_eps = rho_f(T0_eps, theta0_eps, s0)
    v0_eps = M0 / rho0_eps
    M0_eps = v0_eps / np.sqrt(T0_eps)
    print('       %10s %10s %10s %10s %10s' % ('rho', 'v', 'M', 'T', 'theta'))
    print('0    : %.8f %.8f %.8f %.8f %.8f' % (rho0, v0, M0, T0, theta0))
    print('0,eps: %.8f %.8f %.8f %.8f %.8f' % (rho0_eps, v0_eps, M0_eps, T0_eps, theta0_eps))

    # Step 2b
    x0_eps = 0.0
    x0 = x0_eps - eps/thetaprime_f(T0_eps, theta0_eps, rho0_eps, v0_eps)  # Eq. 44
    print('x0=', x0)

    # Step 2c
    eps_ASP = 1.0e-6
    Mpre = np.linspace(M0_eps, 1.+eps_ASP, num=50000)
    xT0 = [x0_eps, T0_eps]
    xTpre = integrate.odeint(dxTdM, xT0, Mpre, mxstep=100000)
    xpre = xTpre[:,0]
    Tpre = xTpre[:,1]
    rhopre = M0 / (Mpre * np.sqrt(Tpre))  # Eq. 40
    vpre = M0/rhopre
    thetapre = ((Km-3.*gamma*M0**2/rhopre-3.*Tpre*rhopre) / (gamma*P0))**0.25  # Eq. 41


    # Step 3a
    # i= 1
    eps = -1.0e-10
    theta1_eps = theta1 + eps
    T1_eps = T1 + eps * dTdtheta(T1, theta1, rho1, M1, v1, s1)
    rho1_eps = rho_f(T1_eps, theta1_eps, s1)
    v1_eps = M0 / rho1_eps
    M1_eps = v1_eps / np.sqrt(T1_eps)
    print('')
    print('       %10s %10s %10s %10s %10s' % ('rho', 'v', 'M', 'T', 'theta'))
    print('1    : %.8f %.8f %.8f %.8f %.8f' % (rho1, v1, M1, T1, theta1))
    print('1,eps: %.8f %.8f %.8f %.8f %.8f' % (rho1_eps, v1_eps, M1_eps, T1_eps, theta1_eps))

    # Step 3b
    x1_eps = 0.0
    x1 = x1_eps - eps/thetaprime_f(T1_eps, theta1_eps, rho1_eps, v1_eps)
    print('x1=', x1)

    # Step 3c
    # resolution
    Mrel = np.linspace(M1_eps, 1.-eps_ASP, num=50000)
    xT0 = [x1_eps, T1_eps]
    xTrel = integrate.odeint(dxTdM, xT0, Mrel, mxstep=100000)
    xrel = xTrel[:,0]
    Trel = xTrel[:,1]
    rhorel = M0 / (Mrel * np.sqrt(Trel))  # Eq. 40
    vrel = M0/rhorel
    thetarel = ((Km-3.*gamma*M0**2/rhorel-3.*Trel*rhorel) / (gamma*P0))**0.25  # Eq. 41
    # Reverse the solution
    xrel = xrel[::-1]
    Trel = Trel[::-1]
    thetarel = thetarel[::-1]
    Mrel = Mrel[::-1]
    vrel = vrel[::-1]
    rhorel = rhorel[::-1]


    # Step 4
    print('')
    print('thetapre[-1]=', thetapre[-1], '  thetarel[0]=', thetarel[0])
    if thetapre[-1] < thetarel[0] :
        print('Continuous case')
        print('todo')
    else:
        print('Embedded shock case')
        irel=0
        Lrel = len(thetarel)
        Lpre = len(thetapre)
        ilasts=[]
        irels=[]
        errs=[]
        while irel<Lrel and thetarel[irel] <= thetapre[-1]:
            index, = np.where(thetapre >= thetarel[irel])
            ilast = index[0]-1

            rhop = rhopre[ilast]
            vp = vpre[ilast]
            Tp = Tpre[ilast]
            thetap = thetapre[ilast]
            pp = peos(rhop, Tp)
            ep = eeos(Tp)
            Ep = ep + 0.5*vp**2

            rhos = rhorel[irel]
            vs = vrel[irel]
            Ts = Trel[irel]
            thetas = thetarel[irel]
            ps = peos(rhos, Ts)
            es = eeos(Ts)
            Es = es + 0.5*vs**2

            foo = rhop*vp**2 + pp  # Eq. 12b
            bar = rhos*vs**2 + ps
            err1 = abs(foo-bar)/min(foo,bar)

            foo = vp*(rhop*Ep+pp)  # Eq. 12c
            bar = vs*(rhos*Es+ps)
            err2 = abs(foo-bar)/min(foo,bar)

            errs.append(np.sqrt(err1**2+err2**2))
            ilasts.append(ilast)
            irels.append(irel)
            irel += 1

        ierr = np.array(errs).argmin()
        ilast = ilasts[ierr]
        irel = irels[ierr]
        print('min relative error in jump conditions =', min(errs))
        print('theta near the shock:', thetapre[ilast], thetarel[irel])
        print('Jump in T:', Tpre[ilast], Trel[irel])
        print('Jump in rho:', rhopre[ilast], rhorel[irel])


        #This shifts it so x=0 is the shock.

        x0   =  x0 - xpre[ilast+1]
        xpre =  xpre - xpre[ilast+1]
        x1   =  x1 - xrel[irel]
        xrel =  xrel - xrel[irel]
        x    =  np.append(xpre[0:ilast+1], xrel[irel:])



        T = np.append(Tpre[0:ilast+1], Trel[irel:])
        theta = np.append(thetapre[0:ilast+1], thetarel[irel:])
        M = np.append(Mpre[0:ilast+1], Mrel[irel:])
        v = np.append(vpre[0:ilast+1], vrel[irel:])
        rho = np.append(rhopre[0:ilast+1], rhorel[irel:])

        print('test:')
        print('testx0',x0)
        print('testx1', x1)

        x = np.append(x0,x)
        T = np.append(T0,T)
        theta = np.append(theta0,theta)
        rho = np.append(rho0,rho)
        v = np.append(v0,v)

        x = np.append(x,x1)
        T = np.append(T,T1)
        theta = np.append(theta,theta1)
        rho = np.append(rho,rho1)
        v = np.append(v,v1)



        return x, T, theta, rho




if __name__ == "__main__":

        P0 = 1.0e-4
        gamma = 5./3.
        sigma = 1.0e6
        kappa = 1.
        M0 = 5.

        x_sol, T_sol, theta_sol, rho_sol  = TwoTShock(P0, gamma, sigma, kappa, M0)


        t1, x1, temp1 = get_x_var(args.plotfile, "temperature")
        max_value = np.max(temp1)
        max_index = np.where(temp1 == max_value)

        x1 = np.array(x1) - x1[max_index]
        plt.figure(1)
        plt.plot(x1,temp1, '.', label="Simulation")
        plt.plot(x_sol*1e5, T_sol*100, '-', label="Solution")
        plt.legend()
        plt.savefig("{}_Temp_Test.png".format(os.path.basename(args.plotfile)))


        t1, x1, rho1 = get_x_var(args.plotfile, "density")
        x1 = np.array(x1) - x1[max_index]
        plt.figure(2)
        plt.plot(x1, rho1, '.', label='Simulation')
        plt.plot(x_sol*1e5, rho_sol*5.459690277750901e-13,'-',label='Solution')
        plt.legend()
        plt.savefig("{}_Density_Test.png".format(os.path.basename(args.plotfile)))
