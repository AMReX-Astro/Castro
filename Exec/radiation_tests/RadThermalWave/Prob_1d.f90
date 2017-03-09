
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
      use probdata_module
      use network, only : network_init
      use bl_fort_module, only : rt => c_real
      implicit none
      integer init, namlen
      integer name(namlen)
      real(rt)         problo(1), probhi(1)

      integer untin,i
      real(rt)         rn

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
      rhocv = -1.0e50_rt
      T0 = -1.0e50_rt
      Eexp = -1.e-50_rt
      rexp = -1.e50_rt

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

      use bl_fort_module, only : rt => c_real
      implicit none
      integer level, nscal
      integer lo(1), hi(1)
      integer state_l1,state_h1
      real(rt)         xlo(1), xhi(1), time, delta(1)
      real(rt)         state(state_l1:state_h1, nscal)

      integer i
      type(eos_t) :: eos_state
      real(rt)         :: rho, cv, T, p, eint, xcell
      real(rt)         :: rhoeexp, Vexp, rhoe0

      Vexp = 4.e0_rt/3.e0_rt*Pi * rexp**3

      rhoeexp = Eexp / Vexp

      rho = 1.0e0_rt
      T = 1.0e0_rt

      eos_state % rho = rho
      eos_state % T   = T
      eos_state % xn  = 0.e0_rt
      eos_state % xn(1)  = 1.e0_rt

      call eos(eos_input_rt, eos_state)

      cv = eos_state % cv

      rho = rhocv / cv

      eint = rhoeexp / rho
      T = eint / cv

      rhoe0 = rho * cv * T0

      do i = lo(1), hi(1)   
         xcell = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5e0_rt)

         state(i,URHO) = rho
         state(i,UMX)  = 0.e0_rt

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

        use bl_fort_module, only : rt => c_real
        implicit none
        integer level, nrad
        integer lo(1), hi(1)
        integer rad_state_l1
        integer rad_state_h1
        real(kind=8) xlo(1), xhi(1), time, delta(1)
        real(kind=8) rad_state(rad_state_l1:rad_state_h1, nrad)

        integer i

        do i = lo(1), hi(1)
           rad_state(i,:) = 0.e0_rt
        end do
      end subroutine ca_initrad


