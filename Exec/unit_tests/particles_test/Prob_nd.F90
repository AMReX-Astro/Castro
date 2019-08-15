subroutine amrex_probinit(init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module, only: ZERO, HALF
  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real

  use prob_params_module, only: center, coord_type
  use probdata_module

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer untin, i

  namelist /fortin/ vel_amp

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen > maxlen) call castro_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Set namelist defaults
  vel_amp = 1.0d0

  ! set center, domain extrema
  if (coord_type == 0) then
     center(1) = HALF*(problo(1)+probhi(1))
  elseif (coord_type >= 1) then
     center(1) = ZERO
  endif

#if AMREX_SPACEDIM == 2
  center(2) = HALF*(problo(2)+probhi(2))
#endif
#if AMREX_SPACEDIM == 3
  center(3) = HALF*(problo(3)+probhi(3))
#endif

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
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
subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, state_lo, state_hi, &
                       delta, xlo, xhi)

  use probdata_module, only : vel_amp
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rp
  use network, only: nspec
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only : problo, probhi
  use amrex_constants_module, only : ZERO, ONE, HALF

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt) :: x, y, xc, yc
  real(rt) :: Xn(nspec)

  integer :: i, j, k

  type (eos_t) :: eos_state

  ! set the composition
  Xn(:) = 0.e0_rt
  Xn(1) = 1.e0_rt

  xc = HALF * (probhi(1) - problo(1))
  yc = HALF * (probhi(2) - problo(2))

  state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:NVAR) = ZERO

  do k = lo(3), hi(3)

     do j = lo(2), hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF)

        do i = lo(1), hi(1)
           x = problo(1) + delta(1)*(dble(i) + HALF)

           state(i,j,k,URHO) = ONE

           state(i,j,k,UMX) = -vel_amp * (y - yc)
           state(i,j,k,UMY) =  vel_amp * (x - xc)

           ! compute the internal energy and temperature
           eos_state%p = ONE
           eos_state%rho = state(i,j,k,URHO)
           eos_state%xn(:) = Xn

           call eos(eos_input_rp, eos_state)

           state(i,j,k,UEDEN) = state(i,j,k,URHO)*eos_state%e +  &
                HALF*sum(state(i,j,k,UMX:UMY)**2)/state(i,j,k,URHO)

           state(i,j,k,UEINT) = state(i,j,k,URHO)*eos_state%e

           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO)*Xn(:)

        end do
     end do
  end do

end subroutine ca_initdata
