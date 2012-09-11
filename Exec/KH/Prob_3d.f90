
      subroutine PROBINIT (init,name,namlen,problo,probhi)

      use probdata_module
      use eos_module, only : gamma_const
      implicit none

      integer init, namlen
      integer name(namlen)
      double precision problo(3), probhi(3)

      integer untin,i

      namelist /fortin/ u_amb, rho_amb, T_amb, u_pert, rho_pert, T_pert, dengrad,  velgrad, &
                        probtype, idir

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
      integer maxlen
      parameter (maxlen=256)
      character probin*(maxlen)

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
      end if

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do
         
! set namelist defaults

      rho_amb  = 1.d0
        u_amb  = 10.0d0
        T_amb  = 2.e-4

      rho_pert = 2.0d0
        u_pert = -10.0d0
        T_pert = 1.d-4

      probtype = 1
      idir     = 1

      dengrad = 1.d20
      velgrad = 1.d20

      idir = 1     

      center(1) = 0.5d0 * (problo(1) + probhi(1))
      center(2) = 0.5d0 * (problo(2) + probhi(2))
      center(3) = 0.5d0 * (problo(3) + probhi(3))

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
                            state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                            delta,xlo,xhi)
     use probdata_module
     use eos_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP

     implicit none

     integer level, nscal
     integer lo(3), hi(3)
     integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     double precision xlo(3), xhi(3), time, delta(3)
     double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

     double precision xcen,ycen,zcen
     double precision top, bot
     double precision e,p
     integer i,j,k

      do k = lo(3), hi(3)
         zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)
      do j = lo(2), hi(2)
         ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)
         
         do i = lo(1), hi(1)
            xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)

            top = 2.5d0 + 0.05 * sin(xcen / 0.05) 
            bot = 1.5d0 - 0.05 * sin(xcen / 0.05) 
            if (ycen > top .or. ycen < bot) then
               state(i,j,k,URHO)  = rho_amb
               state(i,j,k,UMX)   = rho_amb * u_amb
               state(i,j,k,UTEMP) = T_amb
            else
               state(i,j,k,URHO)  = rho_pert
               state(i,j,k,UMX)   = rho_pert * u_pert
               state(i,j,k,UTEMP) = T_pert
            endif

            state(i,j,k,UMY)   = 0.d0
            state(i,j,k,UMZ)   = 0.d0

            ! This just has a single species ("X")
            state(i,j,k,UFS) = 1.d0

            ! rho, T, X are the inputs
            ! e, p    are the outputs
            call eos_given_RTX(e,p,state(i,j,k,URHO), &
                               state(i,j,k,UTEMP),state(i,j,k,UFS:))

            !  rho X = rho * X
            state(i,j,k,UFS  ) = state(i,j,k,URHO) * state(i,j,k,URHO)

            !  rho e = rho * e
            state(i,j,k,UEINT) = state(i,j,k,URHO) * e

            !  rho E = rho e + 1/2 rho u^2 
            state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
               0.5d0 * ( state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + &
                         state(i,j,k,UMZ)**2 ) / state(i,j,k,URHO)

         enddo
      enddo
      enddo

      end subroutine ca_initdata

! ::: -----------------------------------------------------------

     subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                           domlo,domhi,delta,xlo,time,bc)
 
     use meth_params_module, only : NVAR

     implicit none
     integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
     integer bc(3,2,*)
     integer domlo(3), domhi(3)
     double precision delta(3), xlo(3), time
     double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

     integer n

      do n = 1,NVAR
         call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
              adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
              domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      end subroutine ca_hypfill

! ::: -----------------------------------------------------------

      subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
           domlo,domhi,delta,xlo,time,bc)

      implicit none
      integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

      call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,delta,xlo,bc)

      end subroutine ca_denfill
