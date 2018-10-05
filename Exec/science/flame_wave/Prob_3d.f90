subroutine amrex_probinit (init, name, namlen, problo, probhi) bind(c)

  use amrex_constants_module
  use amrex_error_module
  use initial_model_module
  use model_parser_module, only : model_r, model_state, npts_model, model_initialized

  use probdata_module

  use amrex_fort_module, only : rt => amrex_real

  use eos_type_module, only : eos_t, eos_input_rt, eos_input_tp
  use eos_module, only : eos
  use network

  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt) :: problo(3), probhi(3)

  type (eos_t) :: eos_state

  integer :: untin, i

  namelist /fortin/ nx_model, interp_BC, zero_vels, &
                    dtemp, x_half_max, x_half_width, &
                    X_min, cutoff_density, hot_ash, &
                    dens_base, T_star, T_hi, T_lo, H_star, atm_delta, &
                    fuel1_name, fuel2_name, fuel3_name, &
                    ash1_name, ash2_name, ash3_name, &
                    fuel1_frac, fuel2_frac, fuel3_frac, &
                    ash1_frac, ash2_frac, ash3_frac, &
                    low_density_cutoff, index_base_from_temp, smallx, &
                    max_hse_tagging_level, max_base_tagging_level

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  type(model_t) :: model_params

  integer :: iash1, iash2, iash3, ifuel1, ifuel2, ifuel3
  logical :: species_defined

  real(rt) :: dx_model
  integer :: ng

  if (namlen > maxlen) call amrex_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults here
  X_min = 1.e-4_rt
  cutoff_density = 500.e0_rt

  dtemp = 3.81e8_rt
  x_half_max = 1.2e5_rt
  x_half_width = 3.6e4_rt

  interp_BC = .false.
  zero_vels = .false.

  dens_base = 2.d6

  T_star = 1.d8
  T_hi = 5.d8
  T_lo   = 5.e7

  H_star = 500.d0
  atm_delta  = 25.d0

  fuel1_name = "helium-4"
  fuel2_name = ""
  fuel3_name = ""

  ash1_name  = "iron-56"
  ash2_name  = ""
  ash3_name  = ""

  fuel1_frac = ONE
  fuel2_frac = ZERO
  fuel3_frac = ZERO

  ash1_frac = ONE
  ash2_frac = ZERO
  ash3_frac = ZERO

  index_base_from_temp = .false.

  low_density_cutoff = 1.d-4

  smallx = 1.d-10

  max_hse_tagging_level = 2
  max_base_tagging_level = 2

  open(newunit=untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! get the species indices
  species_defined = .true.
  ifuel1 = network_species_index(trim(fuel1_name))
  if (ifuel1 < 0) species_defined = .false.

  if (fuel2_name /= "") then
     ifuel2 = network_species_index(trim(fuel2_name))
     if (ifuel2 < 0) species_defined = .false.
  endif

  if (fuel3_name /= "") then
     ifuel3 = network_species_index(trim(fuel3_name))
     if (ifuel3 < 0) species_defined = .false.
  endif

  iash1 = network_species_index(trim(ash1_name))
  if (iash1 < 0) species_defined = .false.

  if (ash2_name /= "") then
     iash2 = network_species_index(trim(ash2_name))
     if (iash2 < 0) species_defined = .false.
  endif

  if (ash3_name /= "") then
     iash3 = network_species_index(trim(ash3_name))
     if (iash3 < 0) species_defined = .false.
  endif

  if (.not. species_defined) then
     print *, ifuel1, ifuel2, ifuel3
     print *, iash1, iash2, iash3
     call amrex_error("ERROR: species not defined")
  endif


  ! set the composition of the underlying star
  model_params % xn_star(:) = smallx
  model_params % xn_star(iash1) = ash1_frac
  if (ash2_name /= "") model_params % xn_star(iash2) = ash2_frac
  if (ash3_name /= "") model_params % xn_star(iash3) = ash3_frac

  ! and the composition of the accreted layer
  model_params % xn_base(:) = smallx
  model_params % xn_base(ifuel1) = fuel1_frac
  if (fuel2_name /= "") model_params % xn_base(ifuel2) = fuel2_frac
  if (fuel3_name /= "") model_params % xn_base(ifuel3) = fuel3_frac

  ! check if they sum to 1
  if (abs(sum(model_params % xn_star) - ONE) > nspec*smallx) then
     call amrex_error("ERROR: ash mass fractions don't sum to 1")
  endif

  if (abs(sum(model_params % xn_base) - ONE) > nspec*smallx) then
     call amrex_error("ERROR: fuel mass fractions don't sum to 1")
  endif

  ! we are going to generate an initial model from problo(2) to
  ! probhi(2) with nx_model zones.  But to allow for a interpolated
  ! lower boundary, we'll add 4 ghostcells to this, so we need to
  ! compute dx
  dx_model = (probhi(3) - problo(3))/nx_model
  ng = 4

  ! now generate the initial models
  call init_model_data(nx_model+ng, 2)

  model_params % dens_base = dens_base
  model_params % T_star = T_star
  model_params % T_hi = T_hi
  model_params % T_lo = T_lo

  model_params % H_star = H_star
  model_params % atm_delta = atm_delta

  model_params % low_density_cutoff = low_density_cutoff

  model_params % index_base_from_temp = index_base_from_temp

  call init_1d_tanh(nx_model+ng, problo(3)-ng*dx_model, probhi(3), model_params, 1)

  ! store the model in the model_parser_module since that is used in
  ! the boundary conditions
  allocate(model_r(nx_model+ng))
  model_r(:) = gen_model_r(:, 1)

  allocate(model_state(nx_model+ng, nvars_model))
  model_state(:, :) = gen_model_state(:, :, 1)

  allocate(npts_model)

  npts_model = nx_model+ng
  model_initialized = .true.

  ! now create a perturbed model -- we want the same base conditions
  ! a hotter temperature
  model_params % T_hi = model_params % T_hi + dtemp

  call init_1d_tanh(nx_model+ng, problo(3)-ng*dx_model, probhi(3), model_params, 2)

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
! ::: -----------------------------------------------------------
subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, state_l1, state_l2, state_l3, state_h1, state_h2, state_h3, &
                       delta, xlo, xhi)

  use amrex_constants_module
  use probdata_module
  use interpolate_module
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt, eos_input_tp, eos_input_rp
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only: problo
  use network, only: nspec
  use initial_model_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(3), hi(3)
  integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  real(rt)         xlo(3), xhi(3), time, delta(3)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

  real(rt)         dist,x,y,z, r
  integer i,j,k,n

  real(rt)         t0,x1,y1,r1,temp

  real(rt)         temppres(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3)

  type (eos_t) :: eos_state
  real(rt) :: sum_excess, sum_excess2, current_fuel, f


  do k = lo(3), hi(3)
     z = problo(3) + (dble(k)+HALF)*delta(3)

     do j = lo(2), hi(2)
        y = problo(2) + (dble(j)+HALF)*delta(2)

        do i = lo(1), hi(1)
           x = problo(1) + (dble(i)+HALF)*delta(1)

           r = sqrt(x**2 + y**2)

           if (r < x_half_max) then
              f = 1.0_rt
           else if (r > x_half_max + x_half_width) then
              f = 0.0_rt
           else
              f = -(r - x_half_max)/x_half_width + ONE
           endif

           state(i,j,k,URHO)  = f * interpolate(z,gen_npts_model,gen_model_r(:,2), &
                                                gen_model_state(:,idens_model,2)) + &
                     (1.0_rt - f) * interpolate(z,gen_npts_model,gen_model_r(:,1), &
                                                gen_model_state(:,idens_model,1))

           state(i,j,k,UTEMP) = f * interpolate(z,gen_npts_model,gen_model_r(:,2), &
                                                gen_model_state(:,itemp_model,2)) + &
                     (1.0_rt - f) * interpolate(z,gen_npts_model,gen_model_r(:,1), &
                                                gen_model_state(:,itemp_model,1))

           temppres(i,j,k) = f * interpolate(z,gen_npts_model,gen_model_r(:,2), &
                                             gen_model_state(:,ipres_model,2)) + &
                  (1.0_rt - f) * interpolate(z,gen_npts_model,gen_model_r(:,1), &
                                             gen_model_state(:,ipres_model,1))

           state(i,j,k,UFS:UFS-1+nspec) = ZERO

           do n = 1, nspec
              state(i,j,k,UFS-1+n) = f * interpolate(z,gen_npts_model,gen_model_r(:,2), &
                                                     gen_model_state(:,ispec_model-1+n,2)) + &
                          (1.0_rt - f) * interpolate(z,gen_npts_model,gen_model_r(:,1), &
                                                     gen_model_state(:,ispec_model-1+n,1))
           enddo

           eos_state%rho = state(i,j,k,URHO)
           eos_state%T = state(i,j,k,UTEMP)
           eos_state%p = temppres(i,j,k)
           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_rp, eos_state)

           state(i,j,k,UTEMP) = eos_state % T
           state(i,j,k,UEINT) = eos_state % rho * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEDEN)

           ! Initial velocities = 0
           state(i,j,k,UMX:UMZ) = 0.e0_rt

           ! convert to partial densities
           do n = 1, nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do
        end do
     enddo
  enddo

end subroutine ca_initdata
