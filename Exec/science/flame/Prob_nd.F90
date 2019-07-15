subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(C, name="amrex_probinit")

  use eos_module
  use eos_type_module
  use castro_error_module
  use network
  use probdata_module
  use extern_probin_module
  use amrex_constants_module
  use conservative_map_module, only : read_conserved_model_file

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt) :: problo(3), probhi(3)
  real(rt) :: xn(nspec)

  integer untin, i

  type (eos_t) :: eos_state

  real(rt) :: lambda_f, v_f

  integer :: ifuel1, iash1, ifuel2, iash2, ifuel3, iash3, ifuel4, iash4

  namelist /fortin/ pert_frac, pert_delta, rho_fuel, T_fuel, T_ash, v_inflow, &
                    fuel1_name, fuel2_name, fuel3_name, fuel4_name, &
                    ash1_name, ash2_name, ash3_name, ash4_name, &
                    X_fuel1, X_fuel2, X_fuel3, X_fuel4, X_ash1, X_ash2, X_ash3, X_ash4, &
                    interp_model, model_file, smallx_init

  ! Build "probin" filename -- the name of file containing
  ! fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call castro_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults
  pert_frac = 0.2e0_rt
  pert_delta = 0.02e0_rt
  rho_fuel = ONE
  T_fuel = ONE
  T_ash = ONE
  v_inflow = ZERO

  interp_model = 0
  model_file = ""

  ! defaults
  fuel1_name = "helium-4"
  X_fuel1 = 1.0

  ash1_name = "oxygen-16"
  X_ash1 = 1.0

  fuel2_name = ""
  X_fuel2 = 0.0

  ash2_name = ""
  X_ash2 = 0.0

  fuel3_name = ""
  X_fuel3 = 0.0

  ash3_name = ""
  X_ash3 = 0.0

  fuel4_name = ""
  X_fuel4 = 0.0

  ash4_name = ""
  X_ash4 = 0.0

  smallx_init = 0.0

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(untin)


  ifuel1 = network_species_index(fuel1_name)
  iash1 = network_species_index(ash1_name)

  ifuel2 = network_species_index(fuel2_name)
  iash2 = network_species_index(ash2_name)

  ifuel3 = network_species_index(fuel3_name)
  iash3 = network_species_index(ash3_name)

  ifuel4 = network_species_index(fuel4_name)
  iash4 = network_species_index(ash4_name)

  if (iash1 < 0 .and. iash2 < 0 .and. iash3 < 0 .and. iash4 < 0) then
     call castro_error("no valid ash state defined")
  endif

  ! fuel state
  xn_fuel(:) = smallx_init
  xn_ash(:) = smallx_init

  if (ifuel1 > 0) then
     xn_fuel(ifuel1) = X_fuel1
  endif

  if (ifuel2 > 0) then
     xn_fuel(ifuel2) = X_fuel2
  endif

  if (ifuel3 > 0) then
     xn_fuel(ifuel3) = X_fuel3
  endif

  if (ifuel4 > 0) then
     xn_fuel(ifuel4) = X_fuel4
  endif

  if (iash1 > 0) then
     xn_ash(iash1) = X_ash1
  endif

  if (iash2 > 0) then
     xn_ash(iash2) = X_ash2
  endif

  if (iash3 > 0) then
     xn_ash(iash3) = X_ash3
  endif

  if (iash4 > 0) then
     xn_ash(iash4) = X_ash4
  endif

  ! normalize
  xn_fuel(:) = xn_fuel(:)/sum(xn_fuel(:))
  xn_ash(:) = xn_ash(:)/sum(xn_ash(:))

  eos_state % rho = rho_fuel
  eos_state % T = T_fuel
  eos_state % xn(:) = xn_fuel(:)

  call eos(eos_input_rt, eos_state)

  e_fuel = eos_state % e
  p_fuel = eos_state % p

  ! compute the ash state -- this should be hot but in pressure equilibrium
  eos_state % rho = rho_fuel  ! initial guess
  eos_state % p = p_fuel      ! pressure equilibrum
  eos_state % T = T_ash
  eos_state % xn(:) = xn_ash(:)

  call eos(eos_input_tp, eos_state)

  rho_ash = eos_state % rho
  e_ash = eos_state % e

  ! mass flux will be constant across the flame
  mass_flux = rho_fuel * v_inflow

  ! if we are going to conservatively interpolate the model from a
  ! model file instead of initializing from start, let's set that up
  ! now
  if (interp_model > 0) then
     call read_conserved_model_file(model_file)
  end if

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

  use probdata_module
  use prob_params_module, only: problo, probhi
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UTEMP, UFS
  use eos_type_module
  use eos_module
  use amrex_constants_module
  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  use conservative_map_module, only : interpolate_conservative, interpolate_avg_to_center

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: xx, xl, xr, x_int, L, f, pert_width
  real(rt) :: val
  integer :: i, j, k, n

  type (eos_t) :: eos_state

  if (interp_model == 0) then

     L = probhi(1) - problo(1)
     x_int = problo(1) + pert_frac*L

     pert_width = pert_delta*L

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              xx = problo(1) + delta(1)*(dble(i) + HALF)

              ! blend the fuel and ash state, keeping the pressure constant
              eos_state % rho = (rho_ash - rho_fuel) * HALF * (ONE - tanh((xx - x_int)/pert_width)) + rho_fuel
              eos_state % T = (T_ash - T_fuel) * HALF * (ONE - tanh((xx - x_int)/pert_width)) + T_fuel
              eos_state % xn(:) = (xn_ash(:) - xn_fuel(:)) * HALF * (ONE - tanh((xx - x_int)/pert_width)) + xn_fuel(:)
              eos_state % p = p_fuel

              call eos(eos_input_tp, eos_state)

              state(i,j,k,URHO ) = eos_state % rho
              state(i,j,k,UMX:UMZ) = ZERO
              state(i,j,k,UMX) = mass_flux
              state(i,j,k,UEDEN) = eos_state % rho * eos_state % e + &
                   HALF * sum(state(i,j,k,UMX:UMZ)**2)/state(i,j,k,URHO)
              state(i,j,k,UEINT) = eos_state % rho * eos_state % e
              state(i,j,k,UTEMP) = eos_state % T
              state(i,j,k,UFS:UFS-1+nspec) = eos_state % rho * eos_state % xn(:)
           end do
        end do
     end do

  else if (interp_model == 1) then

     ! we are going to do a conservative interpolation of a
     ! (presumably higher-resolution) model onto our grid.

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              xl = problo(1) + delta(1)*(dble(i))
              xr = problo(1) + delta(1)*(dble(i) + ONE)

              do n = 1, NVAR
                 call interpolate_conservative(val, xl, xr, n)
                 state(i,j,k,n) = val
              end do

           end do
        end do
     end do

  else

     ! we are going to use a conservative interpolant to convert from
     ! the cell-averages in the model to the cell-center on our grid.
     ! We will then make everything thermodynamically consistent and
     ! leave it to the post init to convert to cell-averages.

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              xl = problo(1) + delta(1)*(dble(i))
              xr = problo(1) + delta(1)*(dble(i) + ONE)

              do n = 1, NVAR
                 call interpolate_avg_to_center(val, xl, xr, n)
                 state(i,j,k,n) = val
              end do

              ! we will respect the thermodynamics, so lets use rho,
              ! X, and T to define the energy.
              eos_state % rho = state(i,j,k,URHO)
              eos_state % T = state(i,j,k,UTEMP)
              eos_state % xn(:) = state(i,j,k,UFS:UFS-1+nspec) / eos_state % rho

              call eos(eos_input_rt, eos_state)

              state(i,j,k,UEINT) = eos_state % rho * eos_state % e
              state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * state(i,j,k,UMX)**2 / eos_state % rho

           end do
        end do
     end do

  end if

end subroutine ca_initdata
