subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(C, name="amrex_probinit")

  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rt
  use castro_error_module
  use network, only: nspec
  use probdata_module
  use extern_probin_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module, only: ZERO, ONE

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  real(rt) :: xn(nspec)

  type (eos_t) :: eos_state

  real(rt) :: lambda_f, v_f

  ! output flame speed and width estimates
  eos_state%rho = rho_burn_ref
  eos_state%T = T_burn_ref
  eos_state%xn(:) = ZERO
  eos_state%xn(1) = ONE

  call eos(eos_input_rt, eos_state)

  lambda_f = sqrt(const_conductivity*T_burn_ref/ &
       (rho_burn_ref*specific_q_burn*nu*rtilde))

  v_f = sqrt(const_conductivity*specific_q_burn*nu*rtilde/ &
       (rho_burn_ref*eos_state%cp**2*T_burn_ref))

  print *, 'flame width = ', lambda_f
  print *, 'flame speed = ', v_f

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
                       state, s_lo, s_hi, &
                       delta, xlo, xhi)

  use network, only: nspec, network_species_index
  use probdata_module
  use prob_params_module, only: problo, probhi
  use extern_probin_module, only : specific_q_burn
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UTEMP, UFS
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t
  use eos_module, only: eos
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module, only: ZERO, ONE

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: xx, x_int, L, f, pert_width
  integer :: i, j, k

  real(rt) :: e_fuel, p_fuel
  real(rt) :: rho_ash, T_ash, e_ash
  real(rt) :: xn_fuel(nspec), xn_ash(nspec)

  type (eos_t) :: eos_state

  integer :: ifuel, iash

  ifuel = network_species_index("fuel")
  iash = network_species_index("ash")

  L = probhi(1) - problo(1)
  x_int = problo(1) + pert_frac*L

  pert_width = pert_delta*L

  ! fuel state
  xn_fuel(:) = ZERO
  xn_fuel(ifuel) = ONE

  eos_state%rho = rho_fuel
  eos_state%T = T_fuel
  eos_state%xn(:) = xn_fuel(:)

  call eos(eos_input_rt, eos_state)

  e_fuel = eos_state%e
  p_fuel = eos_state%p

  ! compute the ash state
  rho_ash = rho_fuel / (ONE + 0.6_rt*specific_q_burn/e_fuel)
  e_ash = e_fuel - p_fuel*(ONE/rho_ash - ONE/rho_fuel) + specific_q_burn
  xn_ash(:) = ZERO
  xn_ash(iash) = ONE

  eos_state%rho = rho_ash
  eos_state%e = e_ash
  eos_state%xn(:) = xn_ash(:)

  call eos(eos_input_re, eos_state)

  T_ash = eos_state%T

  print *, 'fuel: ', rho_fuel, T_fuel, xn_fuel
  print *, 'ash: ', rho_ash, T_ash, xn_ash

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           xx = problo(1) + delta(1)*(dble(i) + 0.5e0_rt)

           if (xx <= x_int) then

              ! ash
              state(i,j,k,URHO ) = rho_ash
              state(i,j,k,UMX:UMZ) = ZERO
              state(i,j,k,UEDEN) = rho_ash*e_ash
              state(i,j,k,UEINT) = rho_ash*e_ash
              state(i,j,k,UTEMP) = T_ash
              state(i,j,k,UFS:UFS-1+nspec) = rho_ash*xn_ash(:)

           elseif (xx > x_int .and. xx < x_int + pert_width) then

              ! linearly interpolate
              f = (xx - x_int)/pert_width
              eos_state%e = (ONE-f)*e_ash + f*e_fuel
              eos_state%rho = (ONE-f)*rho_ash + f*rho_fuel
              eos_state%xn(:) = (ONE-f)*xn_ash(:) + f*xn_fuel(:)

              call eos(eos_input_re, eos_state)

              state(i,j,k,URHO ) = eos_state%rho
              state(i,j,k,UMX:UMZ) = ZERO
              state(i,j,k,UEDEN) = eos_state%rho*eos_state%e
              state(i,j,k,UEINT) = eos_state%rho*eos_state%e
              state(i,j,k,UTEMP) = eos_state%T
              state(i,j,k,UFS:UFS-1+nspec) = eos_state%rho*eos_state%xn(:)

           else

              ! fuel
              state(i,j,k,URHO ) = rho_fuel
              state(i,j,k,UMX:UMZ) = ZERO
              state(i,j,k,UEDEN) = rho_fuel*e_fuel
              state(i,j,k,UEINT) = rho_fuel*e_fuel
              state(i,j,k,UTEMP) = T_fuel
              state(i,j,k,UFS:UFS-1+nspec) = rho_fuel*xn_fuel(:)
           end if

        end do
     end do
  end do

end subroutine ca_initdata
