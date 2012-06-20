! Frank Timmes Sedov solution program, from:
!
! http://cococubed.asu.edu/code_pages/sedov.shtml

      program sedov3
      implicit none

c..exercises the sedov solver

c..declare
      character*80     outfile,string
      integer          i,nstep,iargc,nmax
      parameter        (nmax = 1000)
      real*16          time,zpos(nmax),
     1                 eblast,rho0,omega,vel0,ener0,pres0,cs0,gamma,
     2                 xgeom,
     3                 den(nmax),ener(nmax),pres(nmax),vel(nmax),
     4                 cs(nmax),
     5                 zlo,zhi,zstep,value


c..popular formats
 01   format(1x,t4,a,t8,a,t22,a,t36,a,t50,a,t64,a,t78,a,t92,a)
 02   format(1x,i4,1p8e12.4)
 03   format(1x,i4,1p8e14.6)



c..if your compiler/linker handles command line arguments
c..get the number of spatial points, blast energy, geometry type,
c..density exponent,  and output file name

c      i = iargc()
c      if (i. lt. 2) stop 'too few arguments'

c      call getarg(1,string)
c      nstep = int(value(string))

c      call getarg(2,string)
c      eblast = value(string)

c      call getarg(3,string)
c      xgeom = value(string)

c      call getarg(4,string)
c      omega = value(string)

c      call getarg(5,outfile)


c..otherwise explicitly set stuff
c..standard cases
c..spherical constant density should reach r=1 at t=1

      nstep = 120
      eblast = 1.0
      xgeom  = 3.0q0
      omega  = 0.0q0
      outfile = 'spherical_sedov.dat'


c..input parameters in cgs
      time   = 0.01q0
      rho0   = 1.0q0
      vel0   = 0.0q0
      ener0  = 0.0q0
      pres0  = 0.0q0
      cs0    = 0.0q0
      gamma  = 1.4q0



c..number of grid points, spatial domain, spatial step size.
c..to match hydrocode output, use the mid-cell points.
      zlo   = 0.0q0
      zhi   = 1.2q0
      zstep = (zhi - zlo)/float(nstep)
      do i=1,nstep
       zpos(i)   = zlo + 0.5q0*zstep + float(i-1)*zstep
      enddo


c..get the solution for all spatial points at once

       call sed_1d(time,nstep,zpos,
     1             eblast,omega,xgeom,
     2             rho0,vel0,ener0,pres0,cs0,gamma,
     3             den,ener,pres,vel,cs)


c..output file 
      open(unit=2,file=outfile,status='unknown')
      write(2,02) nstep,time
      write(2,01) 'i','x','den','ener','pres','vel','cs'
      do i=1,nstep
       write(2,03) i,zpos(i),den(i),ener(i),pres(i),vel(i),cs(i)
      enddo
      close(unit=2)

      write(6,*)

c..close up shop
      end







      subroutine sed_1d(time,nstep,xpos,
     1                  eblast,omega_in,xgeom_in,
     2                  rho0,vel0,ener0,pres0,cs0,gam0,
     3                  den,ener,pres,vel,cs)
      implicit none


c..this routine produces 1d solutions for a sedov blast wave propagating
c..through a density gradient rho = rho**(-omega)
c..in planar, cylindrical or spherical geometry 
c..for the standard, singular and vaccum cases.

c..standard case: a nonzero solution extends from the shock to the origin, 
c..               where the pressure is finite.
c..singular case: a nonzero solution extends from the shock to the origin, 
c..               where the pressure vanishes.
c..vacuum case  : a nonzero solution extends from the shock to a boundary point,
c..               where the density vanishes making the pressure meaningless.


c..input: 
c..time     = temporal point where solution is desired seconds
c..xpos(i)  = spatial points where solution is desired cm
c..eblast   = energy of blast erg
c..rho0     = ambient density g/cm**3    rho = rho0 * r**(-omega_in)
c..omegain  = density power law exponent rho = rho0 * r**(-omega_in)
c..vel0     = ambient material speed cm/s
c..pres0    = ambient pressure erg/cm**3
c..cs0      = ambient sound speed cm/s
c..gam0   = gamma law equation of state
c..xgeom_in = geometry factor, 3=spherical, 2=cylindircal, 1=planar


c..for efficiency reasons (doing the energy integrals only once),
c..this routine returns the solution for an array of spatial points 
c..at the desired time point.

c..output:
c..den(i)  = density  g/cm**3
c..ener(i) = specific internal energy erg/g
c..pres(i) = presssure erg/cm**3
c..vel(i)  = velocity cm/s
c..cs(i)   = sound speed cm/s


c..this routine is based upon two papers:
c.."evaluation of the sedov-von neumann-taylor blast wave solution"
c..jim kamm, la-ur-00-6055
c.."the sedov self-similiar point blast solutions in nonuniform media"
c..david book, shock waves, 4, 1, 1994

c..although the ordinary differential equations are analytic, 
c..the sedov expressions appear to become singular for various 
c..combinations of parameters and at the lower limits of the integration
c..range. all these singularies are removable and done so by this routine.

c..these routines are written in real*16 precision because the 
c..real*8 implementations simply run out of precision "near" the origin 
c..in the standard case or the transition region in the vacuum case.


c..declare the pass
      integer          nstep
      real*16          time,xpos(*),
     1                 eblast,rho0,omega_in,vel0,ener0,pres0,cs0,
     1                 gam0,xgeom_in,den(*),ener(*),pres(*),
     3                 vel(*),cs(*)


c..local variables
      external         midpnt,midpowl,midpowl2,sed_v_find,sed_r_find,
     1                 efun01,efun02
      integer          i
      real*16          efun01,efun02,eval1,eval2
      real*16          v0,v2,vstar,vmin,
     1                 alpha,vstep,us,u2,rho2,p2,e2,cs2,
     2                 zeroin,sed_v_find,sed_r_find,
     3                 vat,l_fun,dlamdv,f_fun,g_fun,h_fun,
     4                 denom2,denom3,rho1


c..eps controls the integration accuracy, don't get too greedy or the number
c..of function evaluations required kills.
c..eps2 controls the root find accuracy
c..osmall controls the size of transition regions

      real*16          iprint,eps,eps2,osmall,pi
      parameter        (iprint = 1,
     1                  eps    = 1.0q-10, 
     2                  eps2   = 1.0q-30, 
     3                  osmall = 1.0q-4,
     4                  pi     = 3.1415926535897932384626433832795029q0)



c..common block communication
      logical          lsingular,lstandard,lvacuum,lomega2,lomega3
      real*16          gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv
      common /slap/    gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv,
     3                 lsingular,lstandard,lvacuum,lomega2,lomega3




c..common block communication with the integration stepper
      real*16          gam_int
      common /cmidp/   gam_int


c..popular formats
 87   format(1x,1p10e14.6)
 88   format(1x,8(a7,1pe14.6,' '))


c..initialize the solution
      do i=1,nstep
       den(i)  = 0.0q0
       vel(i)  = 0.0q0
       pres(i) = 0.0q0
       ener(i) = 0.0q0
       cs(i)   = 0.0q0
      end do


c..return on unphysical cases
c..infinite mass
      if (omega_in .ge. xgeom_in) return



c..transfer the pass to common block and create some frequent combinations
      gamma  = gam0
      gamm1  = gamma - 1.0q0
      gamp1  = gamma + 1.0q0
      gpogm  = gamp1 / gamm1
      xgeom  = xgeom_in
      omega  = omega_in
      xg2    = xgeom + 2.0q0 - omega
      denom2 = 2.0q0*gamm1 + xgeom - gamma*omega
      denom3 = xgeom * (2.0q0 - gamma) - omega


c..post shock location v2 and location of singular point vstar
c..kamm equation 18 and 19

      v2    = 4.0q0 / (xg2 * gamp1)
      vstar = 2.0q0 / (gamm1*xgeom + 2.0q0)


c..set two logicals that determines the type of solution

      lstandard = .false.
      lsingular = .false.
      lvacuum   = .false.
      if (abs(v2 - vstar) .le. osmall) then
       lsingular = .true.
       if (iprint .eq. 1) write(6,*) 'singular'
      else if (v2 .lt. vstar - osmall) then
       lstandard = .true.
       if (iprint .eq. 1) write(6,*) 'standard'
      else if (v2 .gt. vstar + osmall) then
       lvacuum = .true.
       if (iprint .eq. 1) write(6,*) 'vacuum'
      end if

c..two apparent singularies, book's notation for omega2 and omega3
      lomega2 = .false.
      lomega3 = .false.
      if (abs(denom2) .le. osmall) then
       lomega2 = .true.
       denom2  = 1.0q-8
       if (iprint .eq. 1) write(6,*) 'omega2 case'
      else if (abs(denom3) .le. osmall) then
       lomega3 = .true.
       denom3  = 1.0q-8
       if (iprint .eq. 1) write(6,*) 'omega3 case'
      end if


c..various exponents, kamm equations 42-47
c..in terms of book's notation:
c..a0=beta6 a1=beta1  a2=-beta2 a3=beta3 a4=beta4 a5=-beta5

      a0  = 2.0q0/xg2
      a2  = -gamm1/denom2
      a1  =  xg2*gamma/(2.0q0 + xgeom*gamm1) * 
     1      (((2.0q0*(xgeom*(2.0q0-gamma) - omega))/(gamma*xg2*xg2))-a2)
      a3  = (xgeom - omega) / denom2
      a4  = xg2 * (xgeom - omega) * a1 /denom3
      a5  = (omega*gamp1 - 2.0q0*xgeom)/denom3


c..frequent combinations, kamm equations 33-37
      a_val = 0.25q0 * xg2 * gamp1
      b_val = gpogm
      c_val = 0.5q0 * xg2 * gamma
      d_val = (xg2 * gamp1)/(xg2*gamp1 - 2.0q0*(2.0q0 + xgeom*gamm1))
      e_val = 0.5q0 * (2.0q0 + xgeom * gamm1)



c..evaluate the energy integrals
c..the singular case can be done by hand; save some cpu cycles
c..kamm equations 80, 81, and 85

      if (lsingular) then

       eval2 = gamp1/(xgeom*(gamm1*xgeom + 2.0q0)**2)      
       eval1 = 2.0q0/gamm1 * eval2
       alpha = gpogm * 2**(xgeom)/(xgeom*(gamm1*xgeom + 2.0q0)**2) 
       if (int(xgeom) .ne. 1) alpha = pi * alpha



c..for the standard or vacuum cases
c..v0 = post-shock origin v0 and vv = vacuum boundary vv
c..set the radius corespondin to vv to zero for now
c..kamm equations 18, and 28.

      else
       v0  = 2.0q0 / (xg2 * gamma)
       vv  = 2.0q0 / xg2
       rvv = 0.0d0
       if (lstandard) vmin = v0
       if (lvacuum)   vmin  = vv



c..the first energy integral
c..in the standard case the term (c_val*v - 1) might be singular at v=vmin

       if (lstandard) then 
        gam_int = a3 - a2*xg2 - 1.0q0
        if (gam_int .ge. 0) then
         call qromo(efun01,vmin,v2,eps,eval1,midpnt)
        else
         gam_int = abs(gam_int)
         call qromo(efun01,vmin,v2,eps,eval1,midpowl)
        end if

c..in the vacuum case the term (1 - c_val/gamma*v) might be singular at v=vmin

       else if (lvacuum) then
        gam_int = a5
        if (gam_int .ge. 0) then
         call qromo(efun01,vmin,v2,eps,eval1,midpnt)
        else
         gam_int = abs(gam_int)
         call qromo(efun01,vmin,v2,eps,eval1,midpowl2)
        end if
       end if



c..the second energy integral
c..in the standard case the term (c_val*v - 1) might be singular at v=vmin

       if (lstandard) then 
        gam_int = a3 - a2*xg2 - 2.0q0
        if (gam_int .ge. 0) then
         call qromo(efun02,vmin,v2,eps,eval2,midpnt)
        else
         gam_int = abs(gam_int)
         call qromo(efun02,vmin,v2,eps,eval2,midpowl)
        end if

c..in the vacuum case the term (1 - c_val/gamma*v) might be singular at v=vmin

       else if (lvacuum) then
        gam_int = a5
        if (gam_int .ge. 0) then
         call qromo(efun02,vmin,v2,eps,eval2,midpnt)
        else
         gam_int = abs(gam_int)
         call qromo(efun02,vmin,v2,eps,eval2,midpowl2)
        end if
       end if



c..kamm equations 57 and 58 for alpha, in a slightly different form.
       if (int(xgeom) .eq. 1) then
        alpha = 0.5q0*eval1 + eval2/gamm1
       else
        alpha = (xgeom - 1.0q0) * pi * (eval1 + 2.0q0 * eval2/gamm1)
       end if
      end if


c..write what we have for the energy integrals
      if (iprint .eq. 1) 
     1    write(6,88) 'xgeom =',xgeom,'eblast=',eblast,
     2                'omega =',omega,'alpha =',alpha,
     3                'j1    =',eval1,'j2    =',eval2

c      write(6,87) omega,alpha





c..immediate post-shock values 
c..kamm page 14 or equations 14, 16, 5, 13
c..r2 = shock position, u2 = shock speed, rho1 = pre-shock density,
c..u2 = post-shock material speed, rho2 = post-shock density, 
c..p2 = post-shock pressure, e2 = post-shoock specific internal energy, 
c..and cs2 = post-shock sound speed

      r2   = (eblast/(alpha*rho0))**(1.0q0/xg2) * time**(2.0q0/xg2)
      us   = (2.0q0/xg2) * r2 / time
      rho1 = rho0 * r2**(-omega)
      u2   = 2.0q0 * us / gamp1
      rho2 = gpogm * rho1
      p2   = 2.0q0 * rho1 * us**2 / gamp1
      e2   = p2/(gamm1*rho2)
      cs2  = sqrt(gamma*p2/rho2)


c..find the radius corresponding to vv
       if (lvacuum)   then
        vwant = vv
        rvv = zeroin(0.0q0,r2,sed_r_find,eps2)
       end if


      if (lstandard .and. iprint .eq. 1) 
     1     write(6,88) 'r2    =',r2,'rho2  =',rho2,
     2                 'u2    =',u2,'e2    =',e2,
     3                 'p2    =',p2,'cs2   =',cs2

      if (lvacuum .and. iprint .eq. 1) 
     1     write(6,88) 
     2                 'rv    =',rvv,
     3                 'r2    =',r2,'rho2  =',rho2,
     4                 'u2    =',u2,'e2    =',e2,
     5                 'p2    =',p2,'cs2   =',cs2




c..now start the loop over spatial positions
      do i=1,nstep
       rwant  = xpos(i)


c..if we are upstream from the shock front
       if (rwant .gt. r2) then
        den(i)  = rho0 * rwant**(-omega)
        vel(i)  = vel0
        pres(i) = pres0
        ener(i) = ener0
        cs(i)   = cs0


c..if we are between the origin and the shock front
c..find the correct similarity value for this radius in the standard or vacuum cases
       else  
        if (lstandard) then 
         vat = zeroin(0.90q0*v0,v2,sed_v_find,eps2)
        else if (lvacuum) then
         vat = zeroin(v2,1.2q0*vv,sed_v_find,eps2)
        end if

c..the physical solution
        call sedov_funcs(vat,l_fun,dlamdv,f_fun,g_fun,h_fun)
        den(i)   = rho2 * g_fun
        vel(i)   = u2   * f_fun
        pres(i)  = p2   * h_fun
        ener(i)  = 0.0q0
        cs(i)    = 0.0q0
        if (den(i) .ne. 0.0) then
         ener(i)  = pres(i) / (gamm1 * den(i))
         cs(i)    = sqrt(gamma * pres(i)/den(i))  
        end if
       end if

c..end of loop over positions
      enddo

      return
      end






      real*16 function efun01(v)
      implicit none
      save

c..evaluates the first energy integrand, kamm equations 67 and 10.
c..the (c_val*v - 1) term might be singular at v=vmin in the standard case.
c..the (1 - c_val/gamma * v) term might be singular at v=vmin in the vacuum case.
c..due care should be taken for these removable singularities by the integrator.

c..declare the pass
      real*16  v


c..common block communication
      logical          lsingular,lstandard,lvacuum,lomega2,lomega3
      real*16          gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv
      common /slap/    gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv,
     3                 lsingular,lstandard,lvacuum,lomega2,lomega3

c..local variables
      real*16 l_fun,dlamdv,f_fun,g_fun,h_fun

c..go
      call sedov_funcs(v,l_fun,dlamdv,f_fun,g_fun,h_fun)
      efun01 = dlamdv * l_fun**(xgeom + 1.0q0) * gpogm * g_fun * v**2

      return
      end






      real*16 function efun02(v)
      implicit none
      save

c..evaluates the second energy integrand, kamm equations 68 and 11.
c..the (c_val*v - 1) term might be singular at v=vmin in the standard case.
c..the (1 - c_val/gamma * v) term might be singular at v=vmin in the vacuum case.
c..due care should be taken for these removable singularities by the integrator.

c..declare the pass
      real*16  v


c..common block communication
      logical          lsingular,lstandard,lvacuum,lomega2,lomega3
      real*16          gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv
      common /slap/    gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv,
     3                 lsingular,lstandard,lvacuum,lomega2,lomega3

c..local variables
      real*16 l_fun,dlamdv,f_fun,g_fun,h_fun,z

c..go
      call sedov_funcs(v,l_fun,dlamdv,f_fun,g_fun,h_fun)
      z = 8.0q0/( (xgeom + 2.0q0 - omega)**2 * gamp1)
      efun02 = dlamdv * l_fun**(xgeom - 1.0q0 ) * h_fun * z

      return
      end







      real*16 function sed_v_find(v)
      implicit none
      save

c..given corresponding physical distances, find the similarity variable v
c..kamm equation 38 as a root find
 
c..declare the pass
      real*16  v


c..common block communication
      logical          lsingular,lstandard,lvacuum,lomega2,lomega3
      real*16          gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv
      common /slap/    gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv,
     3                 lsingular,lstandard,lvacuum,lomega2,lomega3


c..local variables
      real*16 l_fun,dlamdv,f_fun,g_fun,h_fun


      call sedov_funcs(v,l_fun,dlamdv,f_fun,g_fun,h_fun)
      sed_v_find = r2*l_fun - rwant

      return
      end





      real*16 function sed_r_find(r)
      implicit none
      save

c..given the similarity variable v, find the corresponding physical distance
c..kamm equation 38 as a root find
 
c..declare the pass
      real*16  r


c..common block communication
      logical          lsingular,lstandard,lvacuum,lomega2,lomega3
      real*16          gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv
      common /slap/    gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv,
     3                 lsingular,lstandard,lvacuum,lomega2,lomega3


c..local variables
      real*16 l_fun,dlamdv,f_fun,g_fun,h_fun

      call sedov_funcs(vwant,l_fun,dlamdv,f_fun,g_fun,h_fun)
      sed_r_find = r2*l_fun - r

      return
      end








      real*16 function sed_lam_find(v)
      implicit none
      save

c..given the similarity variable v, find the corresponding physical distance
c..kamm equation 38 as a root find
 
c..declare the pass
      real*16  v


c..common block communication
      logical          lsingular,lstandard,lvacuum,lomega2,lomega3
      real*16          gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv
      common /slap/    gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv,
     3                 lsingular,lstandard,lvacuum,lomega2,lomega3


c..local variables
      real*16 l_fun,dlamdv,f_fun,g_fun,h_fun


      call sedov_funcs(v,l_fun,dlamdv,f_fun,g_fun,h_fun)
      sed_lam_find = l_fun - xlam_want

      return
      end







      subroutine sedov_funcs(v,l_fun,dlamdv,f_fun,g_fun,h_fun)
      implicit none
      save

c..given the similarity variable v, returns functions 
c..lambda, f, g, and h and the derivative of lambda with v dlamdv

c..although the ordinary differential equations are analytic, 
c..the sedov expressions appear to become singular for various 
c..combinations of parameters and at the lower limits of the integration
c..range. all these singularies are removable and done so by this routine.


c..declare the pass
      real*16          v,l_fun,dlamdv,f_fun,g_fun,h_fun
 

c..common block communication
      logical          lsingular,lstandard,lvacuum,lomega2,lomega3
      real*16          gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv
      common /slap/    gamma,gamm1,gamp1,gpogm,xgeom,xg2,rwant,r2,
     1                 a0,a1,a2,a3,a4,a5,a_val,b_val,c_val,d_val,e_val,
     2                 omega,vv,xlam_want,vwant,rvv,
     3                 lsingular,lstandard,lvacuum,lomega2,lomega3


c..local variables
      real*16          x1,x2,x3,x4,dx1dv,dx2dv,dx3dv,dx4dv,
     1                 cbag,ebag,beta0,pp1,pp2,pp3,pp4,c2,c6,y,z,
     2                 dpp2dv,eps
      parameter        (eps = 1.0q-30)


c..frequent combinations and their derivative with v
c..kamm equation 29-32, x4 a bit different to save a divide 
c..x1 is book's F

       x1 = a_val * v
       dx1dv = a_val

       cbag = max(eps, c_val * v - 1.0q0)
       x2 = b_val * cbag
       dx2dv = b_val * c_val
      
       ebag = 1.0q0 - e_val * v
       x3 = d_val * ebag
       dx3dv = -d_val * e_val

       x4 = b_val * (1.0q0 - 0.5q0 * xg2 *v)
       dx4dv = -b_val * 0.5q0 * xg2


c..transition region between standard and vacuum cases
c..kamm page 15 or equations 88-92
c..lambda = l_fun is book's zeta
c..f_fun is books V, g_fun is book's D, h_fun is book's P

      if (lsingular) then
       l_fun  = rwant/r2
       dlamdv = 0.0q0
       f_fun  = l_fun
       g_fun  = l_fun**(xgeom - 2.0q0)
       h_fun  = l_fun**xgeom 



c..for the vacuum case in the hole
      else if (lvacuum .and. rwant .lt. rvv) then

       l_fun  = 0.0q0
       dlamdv = 0.0q0
       f_fun  = 0.0q0
       g_fun  = 0.0q0
       h_fun  = 0.0q0



c..omega = omega2 = (2*(gamma -1) + xgeom)/gamma case, denom2 = 0
c..book expressions 20-22

      else if (lomega2) then

       beta0 = 1.0q0/(2.0q0 * e_val)
       pp1   = gamm1 * beta0
       c6    = 0.5q0 * gamp1
       c2    = c6/gamma
       y     = 1.0q0/(x1 - c2)
       z     = (1.0q0 - x1)*y
       pp2   = gamp1 * beta0 * z
       dpp2dv = -gamp1 * beta0 * dx1dv * y * (1.0q0 + z)
       pp3   = (4.0q0 - xgeom - 2.0q0*gamma) * beta0
       pp4   = -xgeom * gamma * beta0

       l_fun = x1**(-a0) * x2**(pp1) * exp(pp2)
       dlamdv = (-a0*dx1dv/x1 + pp1*dx2dv/x2 + dpp2dv) * l_fun
       f_fun = x1 * l_fun
       g_fun = x1**(a0*omega) * x2**pp3 * x4**a5 * exp(-2.0q0*pp2)
       h_fun = x1**(a0*xgeom) * x2**pp4 * x4**(1.0q0 + a5)



c..omega = omega3 = xgeom*(2 - gamma) case, denom3 = 0
c..book expressions 23-25

      else if (lomega3) then

       beta0 = 1.0q0/(2.0q0 * e_val)
       pp1   = a3 + omega * a2
       pp2   = 1.0q0 - 4.0q0 * beta0
       c6    = 0.5q0 * gamp1
       pp3   = -xgeom * gamma * gamp1 * beta0 * (1.0q0 - x1)/(c6 - x1)
       pp4   = 2.0q0 * (xgeom * gamm1 - gamma) * beta0

       l_fun = x1**(-a0) * x2**(-a2) * x4**(-a1)
       dlamdv = -(a0*dx1dv/x1 + a2*dx2dv/x2 + a1*dx4dv/x4) * l_fun
       f_fun = x1 * l_fun
       g_fun = x1**(a0*omega) * x2**pp1 * x4**pp2 * exp(pp3)
       h_fun = x1**(a0*xgeom) * x4**pp4 * exp(pp3)


c..for the standard or vacuum case not in the hole
c..kamm equations 38-41 

      else      
       l_fun = x1**(-a0) * x2**(-a2) * x3**(-a1)
       dlamdv = -(a0*dx1dv/x1 + a2*dx2dv/x2 + a1*dx3dv/x3) * l_fun
       f_fun = x1 * l_fun
       g_fun = x1**(a0*omega)*x2**(a3+a2*omega)*x3**(a4+a1*omega)*x4**a5
       h_fun = x1**(a0*xgeom)*x3**(a4+a1*(omega-2.0q0))*x4**(1.0q0 + a5)

      end if

      return
      end
 







      subroutine midpnt(func,a,b,s,n)    
      implicit none
      save

c..this routine computes the n'th stage of refinement of an extended midpoint 
c..rule. func is input as the name of the function to be integrated between  
c..limits a and b. when called with n=1, the routine returns as s the crudest 
c..estimate of the integralof func from a to b. subsequent calls with n=2,3... 
c..improve the accuracy of s by adding 2/3*3**(n-1) addtional interior points. 

c..declare 
      external          func 
      integer           n,it,j   
      real*16           func,a,b,s,tnm,del,ddel,x,sum 

      if (n.eq.1) then   
       s  = (b-a) * func(0.5q0*(a+b)) 
      else   
       it   = 3**(n-2) 
       tnm  = it  
       del  = (b-a)/(3.0q0*tnm) 
       ddel = del + del  
       x    = a + (0.5q0 * del)   
       sum  = 0.0q0 
       do j=1,it  
        sum = sum + func(x)  
        x   = x + ddel   
        sum = sum + func(x)  
        x   = x + del    
       enddo
       s  = (s + ((b-a) * sum/tnm)) / 3.0q0 
      end if 
      return 
      end    







      subroutine midpowl(funk,aa,bb,s,n)  
      implicit none
      save

c..this routine is an exact replacement for midpnt, except that it allows for 
c..an integrable power-law singularity of the form (x - a)**(-gam_int)
c..at the lower limit aa for 0 < gam_int < 1.

c..declare 
      external          funk  
      integer           n,it,j   
      real*16           func,funk,a,aa,b,bb,s,tnm,del,ddel,x,sum 


c..common block communication
      real*16          gam_int
      common /cmidp/   gam_int


c..a little conversion, recipe equation 4.4.3
      func(x) = 1.0q0/(1.0q0 - gam_int) * x**(gam_int/(1.0q0 - gam_int)) 
     &          * funk(x**(1.0q0/(1.0q0 - gam_int)) + aa)
      b = (bb - aa)**(1.0q0 - gam_int)    
      a = 0.0q0 

c..now exactly as midpnt 
      if (n .eq. 1) then   
       s = (b-a) * func(0.5q0*(a+b)) 
      else   
       it   = 3**(n-2) 
       tnm = it  
       del = (b-a)/(3.0q0*tnm) 
       ddel = del + del  
       x = a + (0.5q0 * del)   
       sum = 0.0q0 
       do j=1,it  
        sum = sum + func(x)  
        x   = x + ddel   
        sum = sum + func(x)  
        x   = x + del    
       enddo
       s = (s + ((b-a) * sum/tnm)) / 3.0q0 
      end if 
      return 
      end    




      subroutine midpowl2(funk,aa,bb,s,n)  
      implicit none
      save

c..this routine is an exact replacement for midpnt, except that it allows for 
c..an integrable power-law singularity of the form (a - x)**(-gam_int) 
c..at the lower limit aa for 0 < gam_int < 1.

c..declare 
      external          funk  
      integer           n,it,j   
      real*16           func,funk,a,aa,b,bb,s,tnm,del,ddel,x,sum 


c..common block communication
      real*16          gam_int
      common /cmidp/   gam_int


c..a little conversion, modulo recipe equation 4.4.3
      func(x) = 1.0q0/(gam_int - 1.0q0) * x**(gam_int/(1.0q0 - gam_int)) 
     &          * funk(aa - x**(1.0q0/(1.0q0 - gam_int)))
      b = (aa - bb)**(1.0q0 - gam_int)    
      a = 0.0q0 

c..now exactly as midpnt 
      if (n .eq. 1) then   
       s = (b-a) * func(0.5q0*(a+b)) 
      else   
       it   = 3**(n-2) 
       tnm = it  
       del = (b-a)/(3.0q0*tnm) 
       ddel = del + del  
       x = a + (0.5q0 * del)   
       sum = 0.0q0 
       do j=1,it  
        sum = sum + func(x)  
        x   = x + ddel   
        sum = sum + func(x)  
        x   = x + del    
       enddo
       s = (s + ((b-a) * sum/tnm)) / 3.0q0 
      end if 
      return 
      end    






      subroutine qromo(func,a,b,eps,ss,choose)  
      implicit none
      save

c..this routine returns as s the integral of the function func from a to b 
c..with fractional accuracy eps. 
c..jmax limits the number of steps; nsteps = 3**(jmax-1)  
c..integration is done via romberg algorithm.    

c..it is assumed the call to choose triples the number of steps on each call  
c..and that its error series contains only even powers of the number of steps. 
c..the external choose may be any of the above drivers, i.e midpnt,midinf... 


c..declare 
      external          choose,func 
      integer           j,jmax,jmaxp,k,km,i
      parameter         (jmax=14, jmaxp=jmax+1, k=5, km=k-1) 
      real*16           a,b,ss,s(jmaxp),h(jmaxp),eps,dss,func


      h(1) = 1.0q0 
      do j=1,jmax 
       call choose(func,a,b,s(j),j)  
       if (j .ge. k) then    
        call polint(h(j-km),s(j-km),k,0.0q0,ss,dss)    
        if (abs(dss) .le. eps*abs(ss)) return    
       end if    
       s(j+1) = s(j) 
       h(j+1) = h(j)/9.0q0 
      enddo
      write(6,*)  'too many steps in qromo' 
      return 
      end    








      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      save

c..given arrays xa and ya of length n and a value x, this routine returns a 
c..value y and an error estimate dy. if p(x) is the polynomial of degree n-1
c..such that ya = p(xa) ya then the returned value is y = p(x) 

c..declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=20)
      real*16          xa(n),ya(n),x,y,dy,c(nmax),d(nmax),dif,dift,
     1                 ho,hp,w,den

c..find the index ns of the closest table entry; initialize the c and d tables
      ns  = 1
      dif = abs(x - xa(1))
      do i=1,n
       dift = abs(x - xa(i))
       if (dift .lt. dif) then
        ns  = i
        dif = dift
       end if
       c(i)  = ya(i)
       d(i)  = ya(i)
      enddo

c..first guess for y
      y = ya(ns)

c..for each column of the table, loop over the c's and d's and update them
      ns = ns - 1
      do m=1,n-1
       do i=1,n-m
        ho   = xa(i) - x
        hp   = xa(i+m) - x
        w    = c(i+1) - d(i)
        den  = ho - hp
        if (den .eq. 0.0) stop ' 2 xa entries are the same in polint'
        den  = w/den
        d(i) = hp * den
        c(i) = ho * den
       enddo

c..after each column is completed, decide which correction c or d, to add
c..to the accumulating value of y, that is, which path to take in the table
c..by forking up or down. ns is updated as we go to keep track of where we
c..are. the last dy added is the error indicator.
       if (2*ns .lt. n-m) then
        dy = c(ns+1)
       else
        dy = d(ns)
        ns = ns - 1
       end if
       y = y + dy
      enddo
      return
      end







      real*16 function zeroin( ax, bx, f, tol)
      implicit real*16 (a-h,o-z)

c-----------------------------------------------------------------------
c
c This subroutine solves for a zero of the function  f(x)  in the
c interval ax,bx.
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0q0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c-----------------------------------------------------------------------

c.... call list variables

      real*16  ax
      real*16  bx
      real*16  f
      real*16  tol
c
      real*16  a
      real*16  b
      real*16  c
      real*16  d
      real*16  e
      real*16  eps
      real*16  fa
      real*16  fb
      real*16  fc
      real*16  tol1
      real*16  xm
      real*16  p
      real*16  q
      real*16  r
      real*16  s

      external f

c----------------------------------------------------------------------

c
c  compute eps, the relative machine precision
c
      eps = 1.0q0
   10 eps = eps/2.0q0
      tol1 = 1.0q0 + eps
      if (tol1 .gt. 1.0q0) go to 10
c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if ( abs(fc) .ge.  abs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0q0*eps*abs(b) + 0.5q0*tol
      xm = 0.5q0*(c - b)
      if (abs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0q0) go to 90
c
c is bisection necessary?
c
      if (abs(e) .lt. tol1) go to 70
      if (abs(fa) .le. abs(fb)) go to 70
c
c is quadratic interpolation possible?
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0q0*xm*s
      q = 1.0q0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0q0*xm*q*(q - r) - (b - a)*(r - 1.0q0))
      q = (q - 1.0q0)*(r - 1.0q0)*(s - 1.0q0)
c
c adjust signs
c
   60 if (p .gt. 0.0q0) q = -q
      p = abs(p)
c
c is interpolation acceptable?
c
      if ((2.0q0*p) .ge. (3.0q0*xm*q - abs(tol1*q))) go to 70
      if (p .ge. abs(0.5q0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (abs(d) .gt. tol1) b = b + d
      if (abs(d) .le. tol1) b = b + Sign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/abs(fc))) .gt. 0.0q0) go to 20
      go to 30
c
c done
c
   90 zeroin = b

      return
      end









      real*16 function value(string)
      implicit none
      save


c..this routine takes a character string and converts it to a real number. 
c..on error during the conversion, a fortran stop is issued

c..declare
      logical          pflag
      character*(*)    string
      character*1      plus,minus,decmal,blank,se,sd,sq,se1,sd1,sq1
      integer          noblnk,long,ipoint,power,psign,iten,j,z,i
      real*16           x,sign,factor,rten,temp
      parameter        (plus = '+'  , minus = '-' , decmal = '.'   ,
     1                  blank = ' ' , 
     2                  se  = 'e'   , sd  = 'd'   , sq = 'q'       ,
     3                  se1 = 'E'   , sd1 = 'D'   , sq1 ='Q'       ,
     4                  rten =  10.0, iten = 10 )

c..initialize
      x      =  0.0q0
      sign   =  1.0q0
      factor =  rten
      pflag  =  .false.
      noblnk =  0
      power  =  0
      psign  =  1
      long   =  len(string)


c..remove any leading blanks and get the sign of the number
      do z = 1,7
       noblnk = noblnk + 1
       if ( string(noblnk:noblnk) .eq. blank) then
        if (noblnk .gt. 6 ) goto  30
       else
        if (string(noblnk:noblnk) .eq. plus) then
         noblnk = noblnk + 1
        else if (string(noblnk:noblnk) .eq. minus) then
         noblnk = noblnk + 1
         sign =  -1.0q0
        end if
        goto 10
       end if
      enddo


c..main number conversion loop
 10   continue
      do i = noblnk,long
       ipoint = i + 1


c..if a blank character then we are done
       if ( string(i:i) .eq. blank ) then
        x     = x * sign
        value = x 
        return


c..if an exponent character, process the whole exponent, and return
       else if (string(i:i) .eq. se .or. 
     1          string(i:i) .eq. sd .or. 
     2          string(i:i) .eq. sq .or.
     3          string(i:i) .eq. se1 .or. 
     4          string(i:i) .eq. sd1 .or. 
     5          string(i:i) .eq. sq1   ) then
        if (x .eq. 0.0 .and. ipoint.eq.2)     x = 1.0q0
        if (sign .eq. -1.0 .and. ipoint.eq.3) x = 1.0q0
        if (string(ipoint:ipoint) .eq. plus) ipoint = ipoint + 1
        if (string(ipoint:ipoint) .eq. minus) then
         ipoint = ipoint + 1
         psign = -1
        end if
        do z = ipoint,long
         if (string(z:z) .eq. blank)  then
          x = sign * x * rten**(power*psign)
          value = x
          return
         else
          j = ichar(string(z:z)) - 48
          if ( (j.lt.0) .or. (j.gt.9) ) goto 30
          power= (power * iten)  + j
         end if
        enddo


c..if an ascii number character, process ie
       else if (string(i:i) .ne. decmal) then
        j = ichar(string(i:i)) - 48
        if ( (j.lt.0) .or. (j.gt.9) ) goto 30
        if (.not.(pflag) ) then
         x = (x*rten) + j
        else
         temp   = j
         x      = x + (temp/factor)
         factor = factor * rten
         goto 20
        end if

c..must be a decimal point if none of the above
c..check that there are not two decimal points!
       else
        if (pflag) goto 30
        pflag = .true.
       end if
 20   continue
      end do

c..if we got through the do loop ok, then we must be done
      x     = x * sign
      value = x 
      return
      

c..error processing the number
 30   write(6,40) long,string(1:long)
 40   format(' error converting the ',i4,' characters ',/,
     1       ' >',a,'< ',/,
     2       ' into a real number in function value')
      stop ' error in routine value'
      end

