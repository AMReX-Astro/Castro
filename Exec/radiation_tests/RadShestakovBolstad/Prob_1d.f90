
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

      use probdata_module
      use network, only : network_init
      implicit none

      integer :: init, namlen
      integer :: name(namlen)
      double precision :: problo(1), probhi(1)

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

      rho_0 = 1.8212111d-5
      T_0 = 0.1d0           ! keV
      kappa_0 = 0.1d0
      x_jump = 0.5d0
      R = 1.d0

      wref_l1 = 0.d0
      wref_l2 = 0.d0

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
      
      implicit none
      
      integer level, nscal
      integer lo(1), hi(1)
      integer state_l1,state_h1
      double precision state(state_l1:state_h1,NVAR)
      double precision time, delta(1)
      double precision xlo(1), xhi(1)
      
      double precision xcen
      integer i
      double precision Tcgs, B0, nu0, l0, x0, u0
      double precision pi, rhoe_0

      Tcgs = T_0 * 1.d3 * ev2erg / k_B
      
      pi = 4.0d0*atan(1.0d0)
      B0 = 8.*pi*hplanck/c_light**3
      nu0 = k_B*Tcgs/hplanck
      l0 = nu0**3/kappa_0
      x0 = l0/sqrt(3.d0)
      u0 = B0*nu0**3
      
!      rhoe_0 = R * k_B * Tcgs * u0 / hplanck
      rhoe_0 = 99968636.6828d0 * Tcgs * rho_0
      !          cv            * T    * rho

      do i = lo(1), hi(1)
         xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)
         
         state(i,URHO ) = rho_0
         state(i,UMX  ) = 0.d0
         
         ! set the composition to be all in the first species
         state(i,UFS:UFS-1+nspec) = 0.d0
         state(i,UFS) = state(i,URHO)
         if (naux > 0) then
            state(i,UFX) = state(i,URHO)
         end if
         
         if (abs(xcen)/x0 .lt. x_jump) then
            state(i,UEDEN) = rhoe_0 
            state(i,UEINT) = rhoe_0
            state(i,UTEMP) = Tcgs
         else
            state(i,UEDEN) = 0.d0
            state(i,UEINT) = 0.d0
            state(i,UTEMP) = 0.d0
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

        implicit none
        integer level, nrad
        integer lo(1), hi(1)
        integer rad_state_l1,rad_state_h1

        real(kind=8) xlo(1), xhi(1), time, delta(1)
        real(kind=8) rad_state(rad_state_l1:rad_state_h1,nrad)

        ! local variables
        double precision xcen
        integer i

        do i = lo(1), hi(1)  
           rad_state(i, :) = 0.d0
        enddo

      end subroutine ca_initrad

