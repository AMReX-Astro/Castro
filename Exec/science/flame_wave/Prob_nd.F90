subroutine amrex_probinit (init, name, namlen, problo, probhi) bind(c)

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: ZERO, ONE, HALF
  use castro_error_module, only: castro_error
  use model_parser_module, only: model_parser_init
  use initial_model_module, only: model_t, init_model_data, gen_model_r, gen_model_state, init_1d_tanh
  use probdata_module, only: dx_model, dtemp, x_half_max, x_half_width, &
                             X_min, cutoff_density, dens_base, T_star, &
                             T_hi, T_lo, H_star, atm_delta, &
                             fuel1_name, fuel2_name, fuel3_name, fuel4_name, &
                             ash1_name, ash2_name, ash3_name, &
                             fuel1_frac, fuel2_frac, fuel3_frac, fuel4_frac, &
                             ash1_frac, ash2_frac, ash3_frac, &
                             low_density_cutoff, smallx, &
                             max_hse_tagging_level, max_base_tagging_level, x_refine_distance
  use network, only: nspec, network_species_index
  use prob_params_module, only : center
  use meth_params_module, only : small_dens

  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt) :: problo(3), probhi(3)

  ! Build "probin" filename -- the name of file containing fortin namelist.
  type(model_t) :: model_params

  integer :: iash1, iash2, iash3, ifuel1, ifuel2, ifuel3, ifuel4
  logical :: species_defined

  integer :: nx_model
  integer :: ng

  ! get the problm parameters
  call probdata_init(name, namlen)


  ! check to make sure that small_dens is less than low_density_cutoff
  ! if not, funny things can happen above the atmosphere
  if (small_dens >= 0.99_rt * low_density_cutoff) then
     call castro_error("ERROR: small_dens should be set lower than low_density_cutoff")
  end if

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

  if (fuel4_name /= "") then
     ifuel4 = network_species_index(trim(fuel4_name))
     if (ifuel4 < 0) species_defined = .false.
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
     print *, ifuel1, ifuel2, ifuel3, ifuel4
     print *, iash1, iash2, iash3
     call castro_error("ERROR: species not defined")
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
  if (fuel4_name /= "") model_params % xn_base(ifuel4) = fuel4_frac

  ! check if they sum to 1
  if (abs(sum(model_params % xn_star) - ONE) > nspec*smallx) then
     call castro_error("ERROR: ash mass fractions don't sum to 1")
  endif

  if (abs(sum(model_params % xn_base) - ONE) > nspec*smallx) then
     call castro_error("ERROR: fuel mass fractions don't sum to 1")
  endif

  ! we are going to generate an initial model from problo(2) to
  ! probhi(2) with nx_model zones.  But to allow for a interpolated
  ! lower boundary, we'll add 4 ghostcells to this, so we need to
  ! compute dx
  nx_model = int((probhi(AMREX_SPACEDIM) - problo(AMREX_SPACEDIM))/dx_model)

  !dx_model = (probhi(AMREX_SPACEDIM) - problo(AMREX_SPACEDIM))/nx_model
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

  call init_1d_tanh(nx_model+ng, &
                    problo(AMREX_SPACEDIM)-ng*dx_model, probhi(AMREX_SPACEDIM), &
                    model_params, 1)

  ! store the model in the model_parser_module since that is used in
  ! the boundary conditions
  call model_parser_init(nx_model+ng, gen_model_r(:,1), gen_model_state(:,:,1))

  ! now create a perturbed model -- we want the same base conditions
  ! a hotter temperature
  model_params % T_hi = model_params % T_hi + dtemp

  call init_1d_tanh(nx_model+ng, &
                    problo(AMREX_SPACEDIM)-ng*dx_model, probhi(AMREX_SPACEDIM), &
                    model_params, 2)

  ! set center
  center(:) = HALF * (problo(:) + probhi(:))

#if AMREX_SPACEDIM == 2
  ! for axisymmetry, put the x-center on the x-axis
  center(1) = ZERO
#endif

end subroutine amrex_probinit



subroutine ca_initdata(lo, hi, &
                       state, s_lo, s_hi, &
                       dx, problo) bind(C, name='ca_initdata')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: ZERO, HALF, ONE
#ifndef AMREX_USE_CUDA
  use castro_error_module, only: castro_error
#endif
  use probdata_module, only: x_half_width, x_half_max
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rt, eos_input_tp, eos_input_rp
  use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec
  use initial_model_module, only: gen_npts_model, gen_model_r, gen_model_state, &
                                  idens_model, itemp_model, ipres_model, ispec_model
  use interpolate_module, only: interpolate ! function

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
  real(rt), intent(in   ) :: dx(3), problo(3)

  real(rt) :: x, y, z, r, height
  integer :: i, j, k, n

  real(rt) :: temppres

  type (eos_t) :: eos_state
  real(rt) :: f

  !$gpu

  do k = lo(3), hi(3)
     z = problo(3) + (dble(k) + HALF) * dx(3)

     do j = lo(2), hi(2)
        y = problo(2) + (dble(j) + HALF) * dx(2)

        do i = lo(1), hi(1)
           x = problo(1) + (dble(i) + HALF) * dx(1)

           ! lateral distance
           if (AMREX_SPACEDIM == 2) then
              r = x
              height = y
           else if (AMREX_SPACEDIM == 3) then
              r = sqrt(x**2 + y**2)
              height = z
#ifndef AMREX_USE_CUDA
           else
              call castro_error("ERROR: problem not setup for 1D")
#endif
           end if

           if (r < x_half_max) then
              f = 1.0_rt
           else if (r > x_half_max + x_half_width) then
              f = 0.0_rt
           else
              f = -(r - x_half_max)/x_half_width + ONE
           endif

           state(i,j,k,URHO)  = f * interpolate(height,gen_npts_model,gen_model_r(:,2), &
                                                gen_model_state(:,idens_model,2)) + &
                     (1.0_rt - f) * interpolate(height,gen_npts_model,gen_model_r(:,1), &
                                                gen_model_state(:,idens_model,1))

           state(i,j,k,UTEMP) = f * interpolate(height,gen_npts_model,gen_model_r(:,2), &
                                                gen_model_state(:,itemp_model,2)) + &
                     (1.0_rt - f) * interpolate(height,gen_npts_model,gen_model_r(:,1), &
                                                gen_model_state(:,itemp_model,1))

           temppres = f * interpolate(height,gen_npts_model,gen_model_r(:,2), &
                                             gen_model_state(:,ipres_model,2)) + &
                  (1.0_rt - f) * interpolate(height,gen_npts_model,gen_model_r(:,1), &
                                             gen_model_state(:,ipres_model,1))

           state(i,j,k,UFS:UFS-1+nspec) = ZERO

           do n = 1, nspec
              state(i,j,k,UFS-1+n) = f * interpolate(height,gen_npts_model,gen_model_r(:,2), &
                                                     gen_model_state(:,ispec_model-1+n,2)) + &
                          (1.0_rt - f) * interpolate(height,gen_npts_model,gen_model_r(:,1), &
                                                     gen_model_state(:,ispec_model-1+n,1))
           enddo

           eos_state%rho = state(i,j,k,URHO)
           eos_state%T = state(i,j,k,UTEMP)
           eos_state%p = temppres
           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_rp, eos_state)

           state(i,j,k,UTEMP) = eos_state % T
           state(i,j,k,UEINT) = eos_state % rho * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT)

           ! Initial velocities = 0
           state(i,j,k,UMX:UMZ) = 0.e0_rt

           ! convert to partial densities
           do n = 1, nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

        end do
     end do
  end do

end subroutine ca_initdata
