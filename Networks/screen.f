      subroutine screenz (t,d,z1,z2,a1,a2,ymass,aion,zion,nion,
     1                    scfac, dscfacdt)


      implicit none

      integer nion
      double precision t, d, z1, z2, a1, a2
      double precision, dimension(nion) :: ymass, aion, zion
      double precision scfac
      double precision dscfacdt

!.... this subroutine calculates screening factors for nuclear reaction
!.... rates in the weak, intermediate , and strong regimes
!.... given the temperature (t--degk), the density (d--g/cc), the
!.... atomic numbers and weights of the elements in the reaction
!.... channel with the largest coulomb barrier (z1,z2,a1,a2),
!.... and the mean plasma parameters
!.... calculated in main and passed over in common aver:
!.... (mean atomic number--zbar, mean square of the atomic number--z2bar,
!.... mean atomic weight--abar, and total number of moles of
!.... nuclei per gram--ytot1).
!.... the unscreened rate is to be multiplied by the dimensionless
!.... the treatment is based on graboske, dewit, grossman, and cooper
!.... ap j. 181,457 (1973) for weak screening and on
!.... alastuey and jancovici, ap.j. 226, 1034, 1978, with plasma
!.... parameters from itoh, totsuji, setsuo, and dewitt, ap.j. 234,
!.... 1079,1979,  for strong screening (rkw modification).

!.... last revision 15 nov 1982

      double precision x13, x14, x53, x512
      parameter (x13=1.d+0/3.d+0,
     &           x14=1.d+0/4.d+0,
     &           x53=5.d+0/3.d+0,
     &           x512=5.d+0/12.d+0)

      double precision abar, zbar, ytot1, z2bar, theta

      integer iy

      double precision qlam0, ztilda, qlam0z, gamp, taufac
      double precision dqlam0dt, dqlam0zdt, dgampdt, dtaufacdt
      double precision zhat, zhat2, gamef, tau12, alph12
      double precision dgamefdt, dtau12dt, dalph12dt 
      double precision h12w, h12, c, h12fac
      double precision dh12wdt, dh12dt, dcdt

!... calculate averages for screening routine
!... nb  y = x/a with x the mass fraction
!... zi and ai are the nuclear chage and atomic mass number
!... respectively

! this part came in through a common block in Kepler -- do it
! directly here
      abar=0.d0
      zbar=0.d0
      ytot1=0.d0
      z2bar=0.d0

      do iy = 1, nion
         ytot1 = ytot1 + ymass(iy)
         z2bar = z2bar + zion(iy)**2*ymass(iy)
         abar = abar + aion(iy)*ymass(iy)
         zbar = zbar + zion(iy)*ymass(iy)
      enddo

      z2bar=z2bar/ytot1
      abar=abar/ytot1
      zbar=zbar/ytot1

! resume original Kepler screen...
      theta=1.d0
      ytot1=1.d0/abar
      
!.... calculate average plasma parameters
!....
      if((z1*z2).le.0.d0) go to 50
      qlam0=1.88d+8*sqrt(d/(abar*t**3))
      dqlam0dt=-1.5d+0*qlam0 / t

      ztilda=sqrt(z2bar+zbar*theta)

      qlam0z=qlam0*ztilda
      dqlam0zdt=ztilda*dqlam0dt

      gamp=2.27493d+5*(d*zbar*ytot1)**x13/t
      dgampdt=-gamp/t

      taufac=4.248719d+3/t**x13
      dtaufacdt=-x13*taufac/t

!.... calculate screening factor
!.... approx. for strong screening only good for alpha .lt. 1.6

      zhat=(z1+z2)**x53-z1**x53-z2**x53
      zhat2=(z1+z2)**x512-z1**x512-z2**x512

      gamef=2.d0**x13*gamp*z1*z2/(z1+z2)**x13
      dgamefdt=gamef*dgampdt/gamp

      tau12=taufac*(z1**2*z2**2*a1*a2/(a1+a2))**x13
      dtau12dt=tau12*dtaufacdt/taufac

      alph12=3.d0*gamef/tau12
      dalph12dt=alph12*(dgamefdt/gamef - dtau12dt/tau12)

!....
!.... limit alph12 to 1.6 to prevent unphysical behavior
!.... (h dec. as rho inc.) at high rho.  this should really
!.... be replaced by a pycnonuclear reaction rate formula.
!....
      if(alph12.le.1.6d0) go to 60
      alph12=1.6d0
      dalph12dt=0.0d0

      gamef=1.6d0*tau12/3.d0
      dgamefdt=gamef*dtau12dt/tau12

      gamp=gamef*(z1+z2)**x13/(2.d0**x13*z1*z2)
      dgampdt=gamp*dgamefdt/gamef

 60   h12w=z1*z2*qlam0z
      dh12wdt=h12w*dqlam0zdt/qlam0z

      h12=h12w
      dh12dt=dh12wdt
      if(gamef.gt.0.3d0) go to 30
      go to 40
      
 30   c=0.896434d0*gamp*zhat-3.44740d0*gamp**x14*zhat2- 
     1     0.5551d0*(log(gamp)+x53*log(z1*z2/(z1+z2)))-2.996d0
      dcdt=0.896434d0*dgampdt*zhat-
     1     3.44740d0*x14*gamp**(x14-1.0d0)*zhat2*dgampdt-
     2     0.5551d0*dgampdt/gamp

      h12=c-(tau12/3.d0)*(5.d0*alph12**3/32.d0-0.014d0*alph12**4 
     1     -0.0128d0*alph12**5)-gamef*(0.0055d0*alph12**4 
     2     -0.0098d0*alph12**5+0.0048d0*alph12**6)
      dh12dt=dcdt - ((dtau12dt*alph12**3 + 3.0d0*tau12*alph12**2*
     1       dalph12dt)*(5.d0/32.d0 - 0.014d0*alph12 - 
     2       0.0128d0*alph12**2) + tau12*alph12**3*dalph12dt*(-0.014d0
     3       - 2.d0*0.0128d0*alph12))/3.d0 -(dgamefdt*alph12**4 + 4.d0
     4       *gamef*alph12**3*dalph12dt)*(0.0055d0 - 0.0098d0*alph12 -
     5       0.0048d0*alph12**2) - gamef*alph12**4*dalph12dt*(-0.0098d0
     6       + 2.d0*0.0048d0*alph12)

      h12fac=0.77d0

      h12=log(max(1.d+0-0.0562d+0*alph12**3,h12fac))+h12
      if (1.d+0-0.0562d+0*alph12**3 .gt. h12fac) then
         dh12dt=(-3.d0*0.0562d0*alph12**2*dalph12dt)/
     1          (1.d0-0.0562d0*alph12**3) + dh12dt
      endif

      if(gamef.gt.0.8d0) go to 40

      h12=h12w*((0.8d0-gamef)/0.5d0)+h12*((gamef-0.3d0)/0.5d0)
      dh12dt=((dh12wdt*(0.8d0-gamef) - h12w*dgamefdt + dh12dt*
     1       (gamef-0.3d0) + h12*dgamefdt)/0.5d0)

 40   if(h12.gt.300.d0) then
         h12=300.d0
         dh12dt=0.d0
      endif
      if(h12.lt.0.d0) then
         h12=0.d0
         dh12dt=0.d0
      endif

      scfac=exp(h12)
      dscfacdt=scfac*dh12dt
      return

 50   scfac=1.d0
      dscfacdt=0.d0
      return

      end subroutine screenz
