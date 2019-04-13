subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use eos_module
  use eos_type_module
  use amrex_error_module
  use network
  use probdata_module
  use extern_probin_module
  use amrex_constants_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(1), probhi(1)
  real(rt)         xn(nspec)

  integer untin,i

  type (eos_t) :: eos_state

  real(rt)         :: lambda_f, v_f

  namelist /fortin/ pert_frac, pert_delta, rho_fuel, T_fuel, T_ash, &
       fuel1_name, fuel2_name, fuel3_name, fuel4_name, &
       ash1_name, ash2_name, ash3_name, ash4_name, &
       X_fuel1, X_fuel2, X_fuel3, X_fuel4, X_ash1, X_ash2, X_ash3, X_ash4

  ! Build "probin" filename -- the name of file containing 
  ! fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call amrex_error("probin file name too long")
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
  
  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(untin)

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
                      state,state_l1,state_h1,delta,xlo,xhi)

  use network, only: nspec, network_species_index
  use probdata_module
  use prob_params_module, only: problo, probhi
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UTEMP, UFS
  use eos_type_module
  use eos_module
  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(1), hi(1)
  integer state_l1,state_h1
  real(rt)         state(state_l1:state_h1,NVAR)
  real(rt)         time, delta(1)
  real(rt)         xlo(1), xhi(1)

  real(rt)         xx, x_int, L, f, pert_width
  integer i

  real(rt)         :: e_fuel, p_fuel
  real(rt)         :: rho_ash, e_ash
  real(rt)         :: xn_fuel(nspec), xn_ash(nspec)

  type (eos_t) :: eos_state

  integer :: ifuel1, iash1, ifuel2, iash2, ifuel3, iash3, ifuel4, iash4

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

  ifuel1 = network_species_index(fuel1_name)
  iash1 = network_species_index(ash1_name)

  ifuel2 = network_species_index(fuel2_name)
  iash2 = network_species_index(ash2_name)

  ifuel3 = network_species_index(fuel4_name)
  iash3 = network_species_index(ash4_name)

  ifuel4 = network_species_index(fuel4_name)
  iash4 = network_species_index(ash4_name)

  if (iash1 < 0 .and. iash2 < 0 .and. iash3 < 0 .and. iash4 < 0) then
     call amrex_error("no valid ash state defined")
  endif
  

  L = probhi(1) - problo(1)
  x_int = problo(1) + pert_frac*L

  pert_width = pert_delta*L

  ! fuel state
  xn_fuel(:) = ZERO
  xn_ash(:) = ZERO

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

  !print *, 'fuel: ', rho_fuel, T_fuel, xn_fuel
  !print *, 'ash: ', rho_ash, T_ash, xn_ash

  do i = lo(1), hi(1)
     xx = problo(1) + delta(1)*(dble(i) + 0.5e0_rt)

     if (xx <= x_int) then

        ! ash
        state(i,URHO ) = rho_ash
        state(i,UMX:UMZ) = ZERO
        state(i,UEDEN) = rho_ash*e_ash
        state(i,UEINT) = rho_ash*e_ash
        state(i,UTEMP) = T_ash
        state(i,UFS:UFS-1+nspec) = rho_ash*xn_ash(:)

     else if (xx > x_int .and. xx < x_int + pert_width) then

        ! linearly interpolate T, X, and find the rho that keeps us isobaric

        f = (xx - x_int)/pert_width
        eos_state % T = (ONE-f)*T_ash + f*T_fuel
        eos_state % rho = (ONE-f)*rho_ash + f*rho_fuel
        eos_state % xn(:) = (ONE-f)*xn_ash(:) + f*xn_fuel(:)
        eos_state % p = p_fuel

        call eos(eos_input_tp, eos_state)

        state(i,URHO ) = eos_state % rho
        state(i,UMX:UMZ) = ZERO
        state(i,UEDEN) = eos_state % rho * eos_state % e
        state(i,UEINT) = eos_state % rho * eos_state % e
        state(i,UTEMP) = eos_state % T
        state(i,UFS:UFS-1+nspec) = eos_state % rho * eos_state % xn(:)

     else

        ! fuel
        state(i,URHO ) = rho_fuel
        state(i,UMX:UMZ) = ZERO
        state(i,UEDEN) = rho_fuel*e_fuel
        state(i,UEINT) = rho_fuel*e_fuel
        state(i,UTEMP) = T_fuel
        state(i,UFS:UFS-1+nspec) = rho_fuel*xn_fuel(:)
     end if

  end do

end subroutine ca_initdata
