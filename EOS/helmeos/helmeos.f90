module helmeos_module

!..for the tables, in general
      integer          imax,jmax,itmax,jtmax
      parameter        (imax = 271, jmax = 101)
      double precision d(imax),t(jmax)
      
      double precision tlo, thi, tstp, tstpi
      double precision dlo, dhi, dstp, dstpi

!..for the helmholtz free energy tables
      double precision f(imax,jmax),fd(imax,jmax),                     &
                       ft(imax,jmax),fdd(imax,jmax),ftt(imax,jmax),    &
                       fdt(imax,jmax),fddt(imax,jmax),fdtt(imax,jmax), &
                       fddtt(imax,jmax)

!..for the pressure derivative with density ables
      double precision dpdf(imax,jmax),dpdfd(imax,jmax),                &
                       dpdft(imax,jmax),dpdfdt(imax,jmax)

!..for chemical potential tables
      double precision ef(imax,jmax),efd(imax,jmax),                    &
                       eft(imax,jmax),efdt(imax,jmax)

!..for the number density tables
      double precision xf(imax,jmax),xfd(imax,jmax),                    &
                       xft(imax,jmax),xfdt(imax,jmax)

!..for storing the differences
      double precision dt_sav(jmax),dt2_sav(jmax),                      &
                       dti_sav(jmax),dt2i_sav(jmax),                    &
                       dd_sav(imax),dd2_sav(imax),                      &
                       ddi_sav(imax),dd2i_sav(imax)


      integer,          save, private :: max_newton = 100

      double precision, save, private :: ttol = 1.0d-8
      double precision, save, private :: dtol = 1.0d-8

      ! 2006 CODATA physical constants                                                                             

      ! Math constants                                                                                             
      double precision :: pi       = 3.1415926535897932384d0
      double precision :: eulercon = 0.577215664901532861d0
      double precision :: a2rad    
      double precision :: rad2a

      ! Physical constants                                                                                         
      double precision :: g       = 6.6742867d-8
      double precision :: h       = 6.6260689633d-27
      double precision :: hbar    
      double precision :: qe      = 4.8032042712d-10
      double precision :: avo     = 6.0221417930d23
      double precision :: clight  = 2.99792458d10
      double precision :: kerg    = 1.380650424d-16
      double precision :: ev2erg  = 1.60217648740d-12
      double precision :: kev     
      double precision :: amu     = 1.66053878283d-24
      double precision :: mn      = 1.67492721184d-24
      double precision :: mp      = 1.67262163783d-24
      double precision :: me      = 9.1093821545d-28
      double precision :: rbohr   
      double precision :: fine    
      double precision :: hion    = 13.605698140d0

      double precision :: ssol    = 5.67051d-5
      double precision :: asol
      double precision :: weinlam
      double precision :: weinfre 
      double precision :: rhonuc  = 2.342d14

      ! Astronomical constants                                                                                     
      double precision :: msol    = 1.9892d33
      double precision :: rsol    = 6.95997d10
      double precision :: lsol    = 3.8268d33
      double precision :: mearth  = 5.9764d27
      double precision :: rearth  = 6.37d8
      double precision :: ly      = 9.460528d17
      double precision :: pc
      double precision :: au      = 1.495978921d13
      double precision :: secyer  = 3.1558149984d7

      ! Some other useful combinations of the constants
      double precision :: sioncon
      double precision :: forth  
      double precision :: forpi  
      double precision :: kergavo
      double precision :: ikavo  
      double precision :: asoli3 
      double precision :: light2

      ! Constants used for the Coulomb corrections
      double precision :: a1    = -0.898004d0
      double precision :: b1    =  0.96786d0
      double precision :: c1    =  0.220703d0
      double precision :: d1    = -0.86097d0
      double precision :: e1    =  2.5269d0
      double precision :: a2    =  0.29561d0
      double precision :: b2    =  1.9885d0
      double precision :: c2    =  0.288675d0
      double precision :: onethird = 1.0d0/3.0d0
      double precision :: esqu

contains

      subroutine helmeos(do_coulomb, eosfail, state, N, input, input_is_constant, acc_cutoff)

      use meth_params_module, only: do_acc
      use bl_error_module
      use bl_types
      use bl_constants_module
      use eos_type_module
      use eos_data_module

      implicit none

!  Frank Timmes Helmholtz based Equation of State
!  http://cococubed.asu.edu/


!..given a temperature temp [K], density den [g/cm**3], and a composition 
!..characterized by abar and zbar, this routine returns most of the other 
!..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
!..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
!..their derivatives with respect to temperature, density, abar, and zbar.
!..other quantites such the normalized chemical potential eta (plus its
!..derivatives), number density of electrons and positron pair (along 
!..with their derivatives), adiabatic indices, specific heats, and 
!..relativistically correct sound speed are also returned.
!..
!..this routine assumes planckian photons, an ideal gas of ions,
!..and an electron-positron gas with an arbitrary degree of relativity
!..and degeneracy. interpolation in a table of the helmholtz free energy
!..is used to return the electron-positron thermodynamic quantities.
!..all other derivatives are analytic.
!..
!..references: cox & giuli chapter 24 ; timmes & swesty apj 1999

!..input arguments
      logical :: do_coulomb, eosfail
      integer :: N
      type (eos_t), dimension(N) :: state
      integer :: input
      logical :: input_is_constant
      integer :: acc_cutoff

!..outputs
      double precision temp_row(N), den_row(N), abar_row(N), &
                       zbar_row(N), etot_row(N), ptot_row(N), &
                       cv_row(N), cp_row(N),  &
                       xne_row(N), xnp_row(N), etaele_row(N), &
                       pele_row(N), ppos_row(N), dpd_row(N),  &
                       dpt_row(N), dpa_row(N), dpz_row(N),  &
                       ded_row(N), det_row(N), dea_row(N),  &
                       dez_row(N),  &
                       stot_row(N), dsd_row(N), dst_row(N), &
                       htot_row(N), dhd_row(N), dht_row(N), &
                       dpe_row(N), dpdr_e_row(N), &
                       gam1_row(N), cs_row(N)

! these directives are used by f2py to generate a python wrapper around this
! fortran code
!
!f2py intent(in)    :: do_coulomb
!f2py intent(out)   :: eosfail
!f2py intent(inout) :: state

!..declare local variables

      integer j

      logical :: single_iter, double_iter, converged
      logical :: use_acc
      integer :: var, dvar, var1, var2, dvar1, dvar2, iter
      double precision :: v_want(N), v1_want(N), v2_want(N)
      double precision :: xnew, xtol, dvdx, smallx, error, v
      double precision :: v1, v2, dv1dt, dv1dr, dv2dt,dv2dr, delr, error1, error2, told, rold, tnew, rnew, v1i, v2i

      double precision x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
                       dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                       dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                       deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                       dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
                       sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
                       dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
                       gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
                       detadt,detadd,xnefer,dxnedt,dxnedd,s, &
                       temp,den,abar,zbar,ytot1,ye

!..for the abar derivatives
      double precision dpradda,deradda,dsradda, &
                       dpionda,deionda,dsionda, &
                       dpepda,deepda,dsepda,    &
                       dpresda,denerda,dentrda, &
                       detada,dxneda


!..for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz, &
                       dpiondz,deiondz,dsiondz, &
                       dpepdz,deepdz,dsepdz,    &
                       dpresdz,denerdz,dentrdz ,&
                       detadz,dxnedz



!..for the interpolations
      integer          iat,jat
      double precision free,df_d,df_t,df_tt,df_dt
      double precision xt,xd,mxt,mxd, &
                       si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                       si0d,si1d,si2d,si0md,si1md,si2md, &
                       dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                       dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                       ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                       z,din,fi(36)

!..for the coulomb corrections
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd,     &
                       plasg,plasgdd,plasgdt,plasgda,plasgdz,     &
                       ecoul,decouldd,decouldt,decoulda,decouldz, &
                       pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                       scoul,dscouldd,dscouldt,dscoulda,dscouldz

      temp_row = state(:) % T
      den_row  = state(:) % rho
      abar_row = state(:) % abar
      zbar_row = state(:) % zbar

      ! Initial setup for iterations

      if (input .eq. eos_input_rt) then

        ! Nothing to do here.

      elseif (input .eq. eos_input_rh) then

        single_iter = .true.
        v_want(:) = state(:) % h
        var  = ienth
        dvar = itemp

      elseif (input .eq. eos_input_tp) then

        single_iter = .true.
        v_want(:) = state(:) % p
        var  = ipres
        dvar = idens

      elseif (input .eq. eos_input_rp) then

        single_iter = .true.
        v_want(:) = state(:) % p
        var  = ipres
        dvar = itemp

      elseif (input .eq. eos_input_re) then

        single_iter = .true.
        v_want(:) = state(:) % e
        var  = iener
        dvar = itemp

      elseif (input .eq. eos_input_ps) then

        double_iter = .true.
        v1_want(:) = state(:) % p
        v2_want(:) = state(:) % s
        var1 = ipres
        var2 = ientr

      elseif (input .eq. eos_input_ph) then

        double_iter = .true.
        v1_want(:) = state(:) % p
        v2_want(:) = state(:) % h
        var1 = ipres
        var2 = ienth

      elseif (input .eq. eos_input_th) then

        single_iter = .true.
        v_want(:) = state(:) % h
        var  = ienth
        dvar = idens

      endif

      ptot_row = 0.0d0
      dpt_row = 0.0d0
      dpd_row = 0.0d0
      dpa_row = 0.0d0
      dpz_row = 0.0d0
      dpe_row = 0.0d0
      dpdr_e_row = 0.0d0

      etot_row = 0.0d0
      det_row = 0.0d0
      ded_row = 0.0d0
      dea_row = 0.0d0
      dez_row = 0.0d0

      stot_row = 0.0d0
      dst_row = 0.0d0
      dsd_row = 0.0d0

      htot_row = 0.0d0
      dhd_row = 0.0d0
      dht_row = 0.0d0

      pele_row = 0.0d0
      ppos_row = 0.0d0

      xne_row = 0.0d0
      xnp_row = 0.0d0

      etaele_row = 0.0d0

      cv_row = 0.0d0
      cp_row = 0.0d0
      cs_row = 0.0d0
      gam1_row = 0.0d0

      if (N .gt. acc_cutoff .and. do_acc .eq. 1) then
         use_acc = .true.
      else
         use_acc = .false.
      endif

      converged = .false.

      ! Note that for OpenACC, we do not need to have a present clause
      ! for the various constants and arrays -- the enter data constructs in helm_init
      ! is enough for the compiler to recognize they are present at runtime.

      !$acc parallel loop if(use_acc) &
      !$acc copy(den_row,temp_row) &
      !$acc copyin(abar_row,zbar_row) &
      !$acc copyin(input,var,dvar,v_want) &
      !$acc copyout(ptot_row,dpt_row,dpd_row,dpe_row,dpdr_e_row) &
      !$acc copyout(etot_row,det_row,ded_row,dea_row,dez_row) &
      !$acc copyout(stot_row,dst_row,dsd_row) &
      !$acc copyout(htot_row,dht_row,dhd_row) &
      !$acc copyout(cs_row,cv_row,cp_row,gam1_row) &
      !$acc copyout(etaele_row,pele_row,ppos_row,xne_row,xnp_row)
      include 'helm_loop.dek'
      !$acc end parallel loop

      state(:) % T    = temp_row
      state(:) % rho  = den_row

      state(:) % p    = ptot_row(:)
      state(:) % dpdT = dpt_row(:)
      state(:) % dpdr = dpd_row(:)

      state(:) % dpdA = dpa_row(:)
      state(:) % dpdZ = dpz_row(:)
      state(:) % dpde = dpe_row(:)
      state(:) % dpdr_e = dpdr_e_row(:)

      state(:) % e    = etot_row(:)
      state(:) % dedT = det_row(:)
      state(:) % dedr = ded_row(:)
      state(:) % dedA = dea_row(:)   
      state(:) % dedZ = dez_row(:)

      state(:) % s    = stot_row(:)
      state(:) % dsdT = dst_row(:)
      state(:) % dsdr = dsd_row(:)

      state(:) % h    = htot_row(:)
      state(:) % dhdR = dhd_row(:)
      state(:) % dhdT = dht_row(:)

      state(:) % pele = pele_row(:)
      state(:) % ppos = ppos_row(:)

      state(:) % xne = xne_row(:)
      state(:) % xnp = xnp_row(:)

      state(:) % eta = etaele_row(:)

      state(:) % cv   = cv_row(:)
      state(:) % cp   = cp_row(:)
      state(:) % gam1 = gam1_row(:)
!      state(:) % cs   = cs_row

      ! Take care of final housekeeping.

      ! Count the positron contribution in the electron quantities.
      state(:) % xne  = state(:) % xne  + state(:) % xnp
      state(:) % pele = state(:) % pele + state(:) % ppos
 
      ! Use the non-relativistic version of the sound speed, cs = sqrt(gam_1 * P / rho).
      ! This replaces the relativistic version that comes out of helmeos.

      state(:) % cs = sqrt(state(:) % gam1 * state(:) % p / state(:) % rho)

      if (input_is_constant) then

        if (input .eq. eos_input_rh) then

          state(:) % h = v_want(:)

        elseif (input .eq. eos_input_tp) then

          state(:) % p = v_want(:)

        elseif (input .eq. eos_input_rp) then

          state(:) % p = v_want(:)

        elseif (input .eq. eos_input_re) then

          state(:) % e = v_want(:)

        elseif (input .eq. eos_input_ps) then

          state(:) % p = v1_want(:)
          state(:) % s = v2_want(:)

        elseif (input .eq. eos_input_ph) then

          state(:) % p = v1_want(:)
          state(:) % h = v2_want(:)

        elseif (input .eq. eos_input_th) then

          state(:) % h = v_want(:)

        endif

      endif

      end subroutine helmeos

!     -----------------------------------------------------------------
      subroutine helmeos_init

      use bl_error_module
      use eos_data_module

      implicit none
      
      double precision dth, dt2, dti, dt2i
      double precision dd, dd2, ddi, dd2i
      double precision tsav, dsav
      integer i, j

      if ( initialized ) return

!..   open the table
      open(unit=2,file='helm_table.dat',status='old',err=1010)
      GO TO 1011
 1010 CONTINUE
      call bl_error('EOSINIT: Failed to open helm_table.dat')
 1011 CONTINUE

      itmax = imax
      jtmax = jmax

      !$acc enter data copyin(itmax,jtmax)

!..   read the helmholtz free energy table
      tlo   = 3.0d0
      thi   = 13.0d0
      tstp  = (thi - tlo)/float(jmax-1)
      tstpi = 1.0d0/tstp
      dlo   = -12.0d0
      dhi   = 15.0d0
      dstp  = (dhi - dlo)/float(imax-1)
      dstpi = 1.0d0/dstp
      do j=1,jmax
         tsav = tlo + (j-1)*tstp
         t(j) = 10.0d0**(tsav)
         do i=1,imax
            dsav = dlo + (i-1)*dstp
            d(i) = 10.0d0**(dsav)
            read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                 fddt(i,j),fdtt(i,j),fddtt(i,j)
         end do
      end do

      !$acc enter data copyin(tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi)
      !$acc enter data copyin(d,t)
      !$acc enter data copyin(f,fd,ft,fdd,ftt,fdt,fdtt,fddt,fddtt)

!..   read the pressure derivative with density table
      do j = 1, jmax
         do i = 1, imax
            read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
         end do
      end do

      !$acc enter data copyin(dpdf,dpdfd,dpdft,dpdfdt)

!..   read the electron chemical potential table
      do j = 1, jmax
         do i = 1, imax
            read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
         end do
      end do

      !$acc enter data copyin(ef,efd,eft,efdt)

!..   read the number density table
      do j = 1, jmax
         do i = 1, imax
            read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
         end do
      end do

      !$acc enter data copyin(xf,xfd,xft,xfdt)

!..   construct the temperature and density deltas and their inverses 
      do j = 1, jmax-1
         dth         = t(j+1) - t(j)
         dt2         = dth * dth
         dti         = 1.0d0/dth
         dt2i        = 1.0d0/dt2
         dt_sav(j)   = dth
         dt2_sav(j)  = dt2
         dti_sav(j)  = dti
         dt2i_sav(j) = dt2i
      end do
      do i = 1, imax-1
         dd          = d(i+1) - d(i)
         dd2         = dd * dd
         ddi         = 1.0d0/dd
         dd2i        = 1.0d0/dd2
         dd_sav(i)   = dd
         dd2_sav(i)  = dd2
         ddi_sav(i)  = ddi
         dd2i_sav(i) = dd2i
      end do

      !$acc enter data copyin(dt_sav,dt2_sav,dti_sav,dt2i_sav)
      !$acc enter data copyin(dd_sav,dd2_sav,ddi_sav,dd2i_sav)

      close(unit=2)

      ! Some initialization of constants

      esqu = qe * qe

      a2rad   = pi/180.0d0
      rad2a   = 180.0d0/pi

      hbar    = 0.5d0 * h/pi
      kev     = kerg/ev2erg
      rbohr   = hbar*hbar/(me * qe * qe)
      fine    = qe*qe/(hbar*clight)

      asol    = 4.0d0 * ssol / clight
      weinlam = h*clight/(kerg * 4.965114232d0)
      weinfre = 2.821439372d0*kerg/h
      pc      = 3.261633d0 * ly

      sioncon = (2.0d0 * pi * amu * kerg)/(h*h)
      forth   = 4.0d0/3.0d0
      forpi   = 4.0d0 * pi
      kergavo = kerg * avo
      ikavo   = 1.0d0/kergavo
      asoli3  = asol/3.0d0
      light2  = clight * clight


      !$acc enter data &                                                                           
      !$acc copyin(msol,rsol,lsol,mearth,rearth,ly,pc,au,secyer) &
      !$acc copyin(ssol,asol,weinlam,weinfre,rhonuc) &
      !$acc copyin(pi,eulercon,a2rad,rad2a) &
      !$acc copyin(g,h,hbar,qe,avo,clight,kerg,ev2erg,kev,amu,mn,mp,me,rbohr,fine,hion) &
      !$acc copyin(sioncon,forth,forpi,kergavo,ikavo,asoli3,light2) &
      !$acc copyin(a1,b1,c1,d1,e1,a2,b2,c2,onethird,esqu)

      !$acc enter data &
      !$acc copyin(eos_input_rt,eos_input_rh,eos_input_tp,eos_input_rp) &
      !$acc copyin(eos_input_re,eos_input_ps,eos_input_ph,eos_input_th) &
      !$acc copyin(iener,ienth,itemp,idens,ientr,ipres) &
      !$acc copyin(smallt,smalld,ttol,dtol)

      initialized = .true.
    
      end subroutine helmeos_init



      ! quintic hermite polynomial functions                                                             
      ! psi0 and its derivatives                                                                                   
      function psi0(z)
      !$acc routine vector                                                                                         
            double precision :: z, psi0
            psi0 = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      end function

      function dpsi0(z)
      !$acc routine vector                                                                                         
            double precision :: z, dpsi0
            dpsi0 = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      end function

      function ddpsi0(z)
      !$acc routine vector                                                                                         
            double precision :: z, ddpsi0
            ddpsi0 = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)
      end function

      ! psi1 and its derivatives                                                                                   
      function psi1(z)
      !$acc routine vector                                                                                         
            double precision :: z, psi1
            psi1 = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      end function

      function dpsi1(z)
      !$acc routine vector                                                                                         
            double precision :: z, dpsi1
            dpsi1 = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      end function

      function ddpsi1(z)
      !$acc routine vector                                                                                         
            double precision :: z, ddpsi1
            ddpsi1 = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)
      end function

      ! psi2  and its derivatives                                                                                  
      function psi2(z)
      !$acc routine vector                                                                                         
            double precision :: z, psi2
            psi2 = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      end function

      function dpsi2(z)
      !$acc routine vector                                                                                         
            double precision :: z, dpsi2
            dpsi2 = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      end function

      function ddpsi2(z)
      !$acc routine vector                                                                                         
            double precision :: z, ddpsi2
            ddpsi2 = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)
      end function


      ! biquintic hermite polynomial function                                                            
      function h5(fi,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)
      !$acc routine vector                                                                                         
            double precision :: fi(36)
            double precision :: w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md,h5

            h5 =   fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
                 + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
                 + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
                 + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
                 + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
                 + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
                 + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
                 + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
                 + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
                 + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
                 + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
                 + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
                 + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
                 + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
                 + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
                 + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
                 + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
                 + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt
      end function


      ! cubic hermite polynomial functions                                                               
      ! psi0 & derivatives                                                                                         
      function xpsi0(z)
      !$acc routine vector                                                                                         
            double precision :: z, xpsi0
            xpsi0 = z * z * (2.0d0*z - 3.0d0) + 1.0
      end function

      function xdpsi0(z)
      !$acc routine vector                                                                                         
            double precision :: z, xdpsi0
            xdpsi0 = z * (6.0d0*z - 6.0d0)
      end function


      ! psi1 & derivatives                                                                                         
      function xpsi1(z)
      !$acc routine vector                                                                                         
            double precision :: z, xpsi1
            xpsi1 = z * ( z * (z - 2.0d0) + 1.0d0)
      end function

      function xdpsi1(z)
      !$acc routine vector                                                                                         
            double precision :: z, xdpsi1
            xdpsi1 = z * (3.0d0*z - 4.0d0) + 1.0d0
      end function

      ! bicubic hermite polynomial function                                                              
      function h3(fi,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md)
      !$acc routine vector                                                                                         
            double precision :: fi(36)
            double precision :: w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md,h3
            h3 =   fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
                 + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
                 + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
                 + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
                 + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
                 + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
                 + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
                 + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt
      end function

end module helmeos_module
