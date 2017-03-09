
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

      use probdata_module
      use network, only : network_init
      use amrex_fort_module, only : rt => c_real
      implicit none

      integer :: init, namlen
      integer :: name(namlen)
      real(rt)         :: problo(1), probhi(1)

      integer :: untin,i

      namelist /fortin/ rho_0, T_0, kappa_0, x_jump, R, wref_l1, wref_l2

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
      integer maxlen
      parameter (maxlen=256)
      character probin*(maxlen)

      call network_init()

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
      end if

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do
         
      ! set namelist defaults

      rho_0 = 1.8212111e-5_rt
      T_0 = 0.1e0_rt           ! keV
      kappa_0 = 0.1e0_rt
      x_jump = 0.5e0_rt
      R = 1.e0_rt

      wref_l1 = 0.e0_rt
      wref_l2 = 0.e0_rt

!     Read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

!     domain extrema
      xmin = problo(1)
      xmax = probhi(1)

      end subroutine amrex_probinit

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
      subroutine ca_initdata(level,time,lo,hi,nscal, &
                             state,state_l1,state_h1,delta,xlo,xhi)

      use fundamental_constants_module, only : hplanck, k_B, ev2erg, c_light
      use probdata_module
      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UFS, UFX, UTEMP
      use network, only : nspec, naux
      
      use amrex_fort_module, only : rt => c_real
      implicit none
      
      integer level, nscal
      integer lo(1), hi(1)
      integer state_l1,state_h1
      real(rt)         state(state_l1:state_h1,NVAR)
      real(rt)         time, delta(1)
      real(rt)         xlo(1), xhi(1)
      
      real(rt)         xcen
      integer i
      real(rt)         Tcgs, B0, nu0, l0, x0, u0
      real(rt)         pi, rhoe_0

      Tcgs = T_0 * 1.e3_rt * ev2erg / k_B
      
      pi = 4.0e0_rt*atan(1.0e0_rt)
      B0 = 8.*pi*hplanck/c_light**3
      nu0 = k_B*Tcgs/hplanck
      l0 = nu0**3/kappa_0
      x0 = l0/sqrt(3.e0_rt)
      u0 = B0*nu0**3
      
!      rhoe_0 = R * k_B * Tcgs * u0 / hplanck
      rhoe_0 = 99968636.6828e0_rt * Tcgs * rho_0
      !          cv            * T    * rho

      do i = lo(1), hi(1)
         xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5e0_rt)
         
         state(i,URHO ) = rho_0
         state(i,UMX  ) = 0.e0_rt
         
         ! set the composition to be all in the first species
         state(i,UFS:UFS-1+nspec) = 0.e0_rt
         state(i,UFS) = state(i,URHO)
         if (naux > 0) then
            state(i,UFX) = state(i,URHO)
         end if
         
         if (abs(xcen)/x0 .lt. x_jump) then
            state(i,UEDEN) = rhoe_0 
            state(i,UEINT) = rhoe_0
            state(i,UTEMP) = Tcgs
         else
            state(i,UEDEN) = 0.e0_rt
            state(i,UEINT) = 0.e0_rt
            state(i,UTEMP) = 0.e0_rt
         end if
         
      enddo
     
      end subroutine ca_initdata

! ::: 
! ::: -----------------------------------------------------------
! :::
      subroutine ca_initrad(level,time,lo,hi,nrad, &
           rad_state,rad_state_l1,rad_state_h1, &
           delta,xlo,xhi)

        use probdata_module

        use amrex_fort_module, only : rt => c_real
        implicit none
        integer level, nrad
        integer lo(1), hi(1)
        integer rad_state_l1,rad_state_h1

        real(kind=8) xlo(1), xhi(1), time, delta(1)
        real(kind=8) rad_state(rad_state_l1:rad_state_h1,nrad)

        ! local variables
        real(rt)         xcen
        integer i

        do i = lo(1), hi(1)  
           rad_state(i, :) = 0.e0_rt
        enddo

      end subroutine ca_initrad

