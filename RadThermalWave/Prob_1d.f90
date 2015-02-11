
      subroutine PROBINIT (init,name,namlen,problo,probhi)
      use probdata_module
      use network, only : network_init
      implicit none
      integer init, namlen
      integer name(namlen)
      double precision problo(1), probhi(1)

      integer untin,i
      double precision rn

      namelist /fortin/ rhocv, T0, Eexp, rexp

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

!     Set defaults
      rhocv = -1.0d50
      T0 = -1.0d50
      Eexp = -1.d-50
      rexp = -1.d50

      ! domain extrema and center
      xmin = problo(1)
      xmax = probhi(1)

!     Read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      end

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
       state,state_l1,state_h1, &
       delta,xlo,xhi)

      use probdata_module
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS, UFX

      use network, only : nspec

      implicit none
      integer level, nscal
      integer lo(1), hi(1)
      integer state_l1,state_h1
      double precision xlo(1), xhi(1), time, delta(1)
      double precision state(state_l1:state_h1, nscal)

      integer i
      type(eos_t) :: eos_state
      double precision :: rho, cv, T, p, eint, xcell
      double precision :: rhoeexp, Vexp, rhoe0

      Vexp = 4.d0/3.d0*Pi * rexp**3

      rhoeexp = Eexp / Vexp

      rho = 1.0d0
      T = 1.0d0

      eos_state % rho = rho
      eos_state % T   = T
      eos_state % xn  = 0.d0
      eos_state % xn(1)  = 1.d0

      call eos(eos_input_rt, eos_state)

      cv = eos_state % cv

      rho = rhocv / cv

      eint = rhoeexp / rho
      T = eint / cv

      rhoe0 = rho * cv * T0

      do i = lo(1), hi(1)   
         xcell = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)

         state(i,URHO) = rho
         state(i,UMX)  = 0.d0

         if (abs(xcell) .lt. rexp) then
            state(i,UTEMP) = T
            state(i,UEDEN) = rhoeexp
            state(i,UEINT) = rhoeexp
         else
            state(i,UTEMP) = T0
            state(i,UEDEN) = rhoe0
            state(i,UEINT) = rhoe0
         endif
         state(i,UFS) = state(i,URHO)
      enddo

      end subroutine ca_initdata


! ::: 
! ::: -----------------------------------------------------------
! :::
      subroutine ca_initrad(level,time,lo,hi,nrad, &
           rad_state,rad_state_l1, rad_state_h1, &
           delta,xlo,xhi)

        use probdata_module

        implicit none
        integer level, nrad
        integer lo(1), hi(1)
        integer rad_state_l1
        integer rad_state_h1
        real(kind=8) xlo(1), xhi(1), time, delta(1)
        real(kind=8) rad_state(rad_state_l1:rad_state_h1, nrad)

        integer i

        do i = lo(1), hi(1)
           rad_state(i,:) = 0.d0
        end do
      end subroutine ca_initrad


! ::: -----------------------------------------------------------

      subroutine ca_hypfill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)

      use meth_params_module, only : NVAR

      implicit none

      include 'bc_types.fi'
      integer adv_l1,adv_h1
      integer bc(1,2,*)
      integer domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision adv(adv_l1:adv_h1,NVAR)

      integer n

      do n = 1,NVAR
         call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
             domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      end

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_denfill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_h1
      integer bc(1,2,*)
      integer domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision adv(adv_l1:adv_h1)

      call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

      end

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_radfill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_h1
      integer bc(1,2,*)
      integer domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision adv(adv_l1:adv_h1)

      call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

      end

!-----------------------------------------------------------------------
