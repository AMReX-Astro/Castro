subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use prob_params_module, only: center
  use probdata_module
  use amrex_constants_module
  use amrex_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ p_ambient, dens_ambient, dens_pert_factor, vel_pert

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) call amrex_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
  ! set namelist defaults

  ! set center, domain extrema
  center(1) = (problo(1)+probhi(1))/2.e0_rt
  center(2) = (problo(2)+probhi(2))/2.e0_rt
  
  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(unit=untin)

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
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use amrex_constants_module
  use probdata_module
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rp
  use network, only: nspec
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only: center
  
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), time, delta(2)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)

  real(rt)         xcen, ycen
  real(rt)         dens, eint, xvel, X(nspec), temp
  
  integer i,j, icen, jcen

  type (eos_t) :: eos_state

  ! compute the integer location of the center of the domain
  icen = center(1)/delta(1)
  jcen = center(2)/delta(2)

    
  do j = lo(2), hi(2)
     ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

     do i = lo(1), hi(1)
        xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

        if (i == icen .and. j == jcen) then
           dens = dens_ambient*dens_pert_factor
        else
           dens = dens_ambient
        endif

        state(i,j,URHO) = dens

        ! velocity perturbation
        if (xcen < center(1)) then 
           xvel = vel_pert
        else if (xcen > center(1)) then
           xvel = -vel_pert
        else
           xvel = ZERO
        endif

        state(i,j,UMX) = dens*xvel
        state(i,j,UMY) = ZERO
        state(i,j,UMZ) = ZERO

        ! set the composition
        X(:) = ZERO
        X(1) = ONE
        
        
        ! compute the internal energy and temperature
        eos_state%T = ONE ! initial guess
        eos_state%rho = dens
        eos_state%p = p_ambient
        eos_state%xn(:) = X

        call eos(eos_input_rp, eos_state)

        temp = eos_state%T
        eint = eos_state%e

        state(i,j,UEDEN) = dens*eint +  &
             HALF*(state(i,j,UMX)**2/state(i,j,URHO) + &
                    state(i,j,UMY)**2/state(i,j,URHO) + &
                    state(i,j,UMZ)**2/state(i,j,URHO))

        state(i,j,UEINT) = dens*eint
        state(i,j,UTEMP) = temp

        state(i,j,UFS:UFS-1+nspec) = dens*X(:)

     enddo
  enddo

end subroutine ca_initdata

