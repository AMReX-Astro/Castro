subroutine amrex_probinit(init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module, only: ZERO, HALF
  use amrex_error_module, only: amrex_error
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only: center, coord_type
  use probdata_module
  use eos_type_module, only: eos_t, eos_input_rt
  use eos_module, only : eos
  use network, only : nspec
  use extern_probin_module, only: const_conductivity

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer untin, i

  namelist /fortin/ diff_coeff, T1, T2, t_0

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen = 256
  character probin*(maxlen)
  real(rt) :: X(nspec)

  type (eos_t) :: eos_state

  if (namlen > maxlen) call amrex_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Set namelist defaults
  T1 = 1.0_rt
  T2 = 2.0_rt
  t_0 = 0.001_rt
  diff_coeff = 1.0_rt

  ! set center, domain extrema
  if (coord_type == 0) then
     center(1) = HALF*(problo(1)+probhi(1))
  elseif (coord_type >= 1) then
     center(1) = ZERO
  endif

#if AMREX_SPACEDIM >= 2
  center(2) = HALF*(problo(2)+probhi(2))
#endif
#if AMREX_SPACEDIM == 3
  center(3) = HALF*(problo(3)+probhi(3))
#endif

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(unit=untin)

  ! the conductivity is the physical quantity that appears in the
  ! diffusion term of the energy equation.  It is set via
  ! diffusion.conductivity in the inputs file.  For this test problem,
  ! we want to set the diffusion coefficient, D = k/(rho c_v), so the
  ! free parameter we have to play with is rho.  Note that for an
  ! ideal gas, c_v does not depend on rho, so we can call it the EOS
  ! with any density.
  X(:) = 0.e0_rt
  X(1) = 1.e0_rt

  eos_state%T = T1
  eos_state%rho = 1.0
  eos_state%xn(:) = X(:)

  call eos(eos_input_rt, eos_state)

  ! diffusion coefficient is D = k/(rho c_v). we are doing an ideal
  ! gas, so c_v is constant, so find the rho that combines with
  ! the conductivity
  rho0 = const_conductivity/(diff_coeff*eos_state%cv)

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

  use probdata_module, only : T1, T2, diff_coeff, t_0, rho0
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use network, only: nspec
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only : problo
  use prob_util_module, only : analytic
  use amrex_constants_module, only : ZERO, HALF

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt) :: r(3)
  real(rt) :: X(nspec), temp

  integer :: i, j, k

  type (eos_t) :: eos_state

  ! set the composition
  X(:) = 0.e0_rt
  X(1) = 1.e0_rt

  do k = lo(3), hi(3)
     r(3) = problo(3) + delta(3) * (dble(k) + HALF)

     do j = lo(2), hi(2)
        r(3) = problo(2) + delta(2) * (dble(j) + HALF)

        do i = lo(1), hi(1)
           r(3) = problo(1) + delta(1) * (dble(i) + HALF)

           state(i,j,k,URHO) = rho0

           call analytic(r, ZERO, temp)

           state(i,j,k,UTEMP) = temp

           ! compute the internal energy and temperature
           eos_state%T = temp
           eos_state%rho = state(i,j,k,URHO)
           eos_state%xn(:) = X

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UMX) = ZERO
           state(i,j,k,UMY) = ZERO
           state(i,j,k,UMZ) = ZERO

           state(i,j,k,UEDEN) = rho0*eos_state%e +  &
                HALF*sum(state(i,j,k,UMX:UMZ)**2)/state(i,j,k,URHO)

           state(i,j,k,UEINT) = rho0*eos_state%e

           state(i,j,k,UFS:UFS-1+nspec) = rho0*X(:)

        end do
     end do
  end do

end subroutine ca_initdata
