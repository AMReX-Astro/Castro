
      subroutine PROBINIT (init,name,namlen,problo,probhi)
      use probdata_module
      use network, only : network_init
      implicit none
      integer init, namlen
      integer name(namlen)
      double precision problo(2), probhi(2)

      integer untin,i,j,k,dir
      double precision rn

      namelist /fortin/ rhocv, T0, Eexp, rexp, &
           denerr,   dengrad,  max_denerr_lev,  max_dengrad_lev, &
           velerr,   velgrad,  max_velerr_lev,  max_pressgrad_lev, &
           presserr, pressgrad,max_presserr_lev,max_pressgrad_lev, &
           temperr,  tempgrad, max_temperr_lev, max_tempgrad_lev, &
           raderr,   radgrad,  max_raderr_lev,  max_radgrad_lev

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

      denerr = 1.d20
      dengrad = 1.d20
      max_denerr_lev = -1
      max_dengrad_lev = -1
      
      velerr = 1.d20
      velgrad = 1.d20
      max_velerr_lev = -1
      max_velgrad_lev = -1
      
      presserr = 1.d20
      pressgrad = 1.d20
      max_presserr_lev = -1
      max_pressgrad_lev = -1
      
      temperr = 1.d20
      tempgrad = 1.d20
      max_temperr_lev = -1
      max_tempgrad_lev = -1
      
      raderr = 1.d20
      radgrad = 1.d20
      max_raderr_lev = -1
      max_radgrad_lev = -1
      
      ! domain extrema and center
      xmin = problo(1)
      xmax = probhi(1)
      center(1) = 0.5d0*(xmin + xmax)

      ymin = problo(2)
      ymax = probhi(2)
      center(2) = 0.5d0*(ymin + ymax)

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
       state,state_l1,state_l2,state_h1,state_h2, &
       delta,xlo,xhi)

      use probdata_module
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS, UFX

      use network, only : nspec

      implicit none
      integer level, nscal
      integer lo(2), hi(2)
      integer state_l1,state_l2,state_h1,state_h2
      double precision xlo(2), xhi(2), time, delta(2)
      double precision state(state_l1:state_h1,state_l2:state_h2, &
          nscal)

      integer i,j
      double precision :: X(nspec)
      double precision :: rho, cv, T, p, eint, xcell, ycell
      double precision :: rhoeexp, Vexp, rhoe0

      Vexp = 2.0d0 * Pi * rexp**3

      rhoeexp = Eexp / Vexp

      rho = 1.0d0
      T = 1.0d0
      X = 1.0d0
      call eos_get_cv(cv, rho, T, X)
      rho = rhocv / cv

      eint = rhoeexp / rho
      T = eint / cv

      rhoe0 = rho * cv * T0

      do j = lo(2), hi(2)     
         ycell = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)
         do i = lo(1), hi(1)   
            xcell = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)

            state(i,j,URHO) = rho
            state(i,j,UMX)  = 0.d0
            state(i,j,UMY)  = 0.d0

            if (abs(xcell) .lt. rexp .and. abs(ycell) .lt. rexp) then
               state(i,j,UTEMP) = T
               state(i,j,UEDEN) = rhoeexp
               state(i,j,UEINT) = rhoeexp
            else
               state(i,j,UTEMP) = T0
               state(i,j,UEDEN) = rhoe0
               state(i,j,UEINT) = rhoe0
            endif
            state(i,j,UFS) = state(i,j,URHO)
         enddo
      enddo

      end subroutine ca_initdata


! ::: 
! ::: -----------------------------------------------------------
! :::
      subroutine ca_initrad(level,time,lo,hi,nrad, &
           rad_state,rad_state_l1,rad_state_l2, &
           rad_state_h1,rad_state_h2, &
           delta,xlo,xhi)

        use probdata_module

        implicit none
        integer level, nrad
        integer lo(2), hi(2)
        integer rad_state_l1,rad_state_l2
        integer rad_state_h1,rad_state_h2
        real(kind=8) xlo(2), xhi(2), time, delta(2)
        real(kind=8) rad_state(rad_state_l1:rad_state_h1, &
             rad_state_l2:rad_state_h2, nrad)

        integer i,j

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              rad_state(i,j,:) = 0.d0
           end do
        end do
      end subroutine ca_initrad

! ::: 
! ::: -----------------------------------------------------------
! :::
     subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
          domlo,domhi,delta,xlo,time,bc)

      use meth_params_module, only : NVAR

      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_h1,adv_h2
      integer bc(2,2,*)
      integer domlo(2), domhi(2)
      double precision delta(2), xlo(2), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

      integer n

      do n = 1,NVAR
         call filcc(adv(adv_l1,adv_l2,n), &
             adv_l1,adv_l2,adv_h1,adv_h2, &
             domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      end

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
          domlo,domhi,delta,xlo,time,bc)
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_h1,adv_h2
      integer bc(2,2,*)
      integer domlo(2), domhi(2)
      double precision delta(2), xlo(2), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

      call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
          domlo,domhi,delta,xlo,bc)

      end

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_radfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
          domlo,domhi,delta,xlo,time,bc)
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_h1,adv_h2
      integer bc(2,2,*)
      integer domlo(2), domhi(2)
      double precision delta(2), xlo(2), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

      call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
          domlo,domhi,delta,xlo,bc)

      end

!-----------------------------------------------------------------------
