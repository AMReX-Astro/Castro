subroutine PROBINIT (init,name,namlen,problo,probhi)

  use bl_types
  use prob_params_module, only: center
  use probdata_module
  use bl_error_module

  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ thermal_conductivity

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
  ! Set namelist defaults
  thermal_conductivity = 1.0_dp_t

  ! set center, domain extrema
  center(1) = (problo(1)+probhi(1))/2.d0
  center(2) = (problo(2)+probhi(2))/2.d0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

end subroutine PROBINIT


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
  use network, only: nspec
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only : problo
  
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)

  double precision xc, yc
  double precision dens, eint, xvel, X(nspec), temp
  
  integer i,j

  type (eos_t) :: eos_state

  do j = lo(2), hi(2)
     yc = problo(2) + delta(2)*(dble(j) + HALF)

     do i = lo(1), hi(1)
        xc = problo(1) + delta(1)*(dble(i) + HALF)

        state(i,j,URHO) = ONE

        state(i,j,UMX) = ZERO
        state(i,j,UMY) = ZERO

        ! set the composition
        X(:) = 0.d0
        X(1) = 1.d0
        
        ! compute the internal energy and temperature
        eos_state%T = 1.d0 ! initial guess
        eos_state%rho = state(i,j,URHO)
        eos_state%xn(:) = X

        call eos(eos_input_rt, eos_state)

        temp = eos_state%T
        eint = eos_state%e

        state(i,j,UEDEN) = dens*eint +  &
             0.5d0*(state(i,j,UMX)**2/state(i,j,URHO) + &
                    state(i,j,UMY)**2/state(i,j,URHO))

        state(i,j,UEINT) = dens*eint
        state(i,j,UTEMP) = temp

        state(i,j,UFS:UFS-1+nspec) = dens*X(:)

     enddo
  enddo

end subroutine ca_initdata

