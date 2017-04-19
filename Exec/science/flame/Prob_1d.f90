subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use eos_module
  use eos_type_module
  use bl_error_module
  use network
  use probdata_module
  use extern_probin_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(1), probhi(1)
  real(rt)         xn(nspec)

  integer untin,i

  type (eos_t) :: eos_state

  real(rt)         :: lambda_f, v_f

  namelist /fortin/ pert_frac, pert_delta, rho_fuel, T_fuel, T_ash

  ! Build "probin" filename -- the name of file containing 
  ! fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error("probin file name too long")
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

  integer :: ifuel, iash

  ifuel = network_species_index("helium-4")
  iash = network_species_index("oxygen-16")

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

  ! compute the ash state -- this should be hot but in pressure equilibrium
  xn_ash(:) = ZERO
  xn_ash(iash) = ONE

  eos_state%rho = rho_fuel  ! initial guess
  eos_state%p = p_fuel      ! pressure equilibrum
  eos_state%T = T_ash
  eos_state%xn(:) = xn_ash(:)

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

     elseif (xx > x_int .and. xx < x_int + pert_width) then

        ! linearly interpolate
        f = (xx - x_int)/pert_width
        eos_state%T = (ONE-f)*T_ash + f*T_fuel
        eos_state%rho = (ONE-f)*rho_ash + f*rho_fuel
        eos_state%xn(:) = (ONE-f)*xn_ash(:) + f*xn_fuel(:)
        
        call eos(eos_input_rt, eos_state)

        state(i,URHO ) = eos_state%rho
        state(i,UMX:UMZ) = ZERO
        state(i,UEDEN) = eos_state%rho*eos_state%e
        state(i,UEINT) = eos_state%rho*eos_state%e
        state(i,UTEMP) = eos_state%T
        state(i,UFS:UFS-1+nspec) = eos_state%rho*eos_state%xn(:)

     else

        ! fuel
        state(i,URHO ) = rho_fuel
        state(i,UMX:UMZ) = ZERO
        state(i,UEDEN) = rho_fuel*e_fuel
        state(i,UEINT) = rho_fuel*e_fuel
        state(i,UTEMP) = T_fuel
        state(i,UFS:UFS-1+nspec) = rho_fuel*xn_fuel(:)
     endif

  enddo

end subroutine ca_initdata
