!!  Create a 1-d hydrostatic, atmosphere with an isothermal region
!!  (T_star) representing the NS, a hyperbolic tangent rise to a
!!  peak temperature (T_base) representing the base of an accreted
!!  layer, an isoentropic profile down to a lower temperature (T_lo),
!!  and then isothermal. This can serve as an initial model for a
!!  nova or XRB.
!!
!!  The temperature profile is:
!!
!!         ^
!!         |
!!  T_base +        /\
!!         |       /  \
!!         |      /  . \
!!  T_star +-----+      \
!!         |     .   .   \
!!         |              \
!!         |     .   .     \
!!  T_lo   +                +-----------
!!         |     .   .
!!         +-----+---+---------------> r
!!         |      \  /
!!         |      atm_delta
!!         |< H_star>|
!!
!!  We take dens_base, the density at the base of the isentropic layer
!!  as input.  The composition is "ash" in the lower isothermal region
!!  and "fuel" in the isentropic and upper isothermal regions.  In the
!!  linear transition region, we linearly interpolate the composition.
!!
!!  The fuel and ash compositions are specified by the fuel?_name,
!!  fuel?_frac and ash?_name, ash?_frac parameters (name of the species
!!  and mass fraction).  Where ? = 1,2,3.
!!
!!  The model is placed into HSE by the following differencing:
!!
!!   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
!!
!!  This will be iterated over in tandem with the EOS call,
!!  P(i-1) = P_eos(rho(i-1), T(i-1), X(i-1)
!!

!!
!! this uses the model_parser_module to store the result

subroutine init_1d_tanh(nx, xmin, xmax)

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rt
  use network, only : nspec, network_species_index, spec_names
  use fundamental_constants_module, only: Gconst
  use model_parser_module
  use probdata_module
  use meth_params_module, only : const_grav

  implicit none

  integer, intent(in) :: nx
  real(kind=dp_t) :: xmin, xmax

  integer :: i, n

  real (kind=dp_t) :: slope_T, slope_xn(nspec)

  real (kind=dp_t) :: pres_base, entropy_base
  real (kind=dp_t), DIMENSION(nspec) :: xn_base, xn_star

  real :: A, B

  ! we'll get the composition from the network module
  ! we allow for 3 different species separately in the fuel and ash
  integer :: ifuel1, ifuel2, ifuel3
  integer :: iash1, iash2, iash3
  logical :: species_defined

  real (kind=dp_t) :: dCoord

  real (kind=dp_t) :: dens_zone, temp_zone, pres_zone, entropy
  real (kind=dp_t) :: dpd, dpt, dsd, dst

  real (kind=dp_t) :: p_want, drho, dtemp, delx

  real (kind=dp_t), parameter :: TOL = 1.e-10

  integer, parameter :: MAX_ITER = 250

  integer :: iter

  logical :: converged_hse, fluff

  real (kind=dp_t), dimension(nspec) :: xn

  integer :: index_base

  logical :: isentropic

  integer :: narg

  type (eos_t) :: eos_state

  real(kind=dp_t) :: sumX

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
     call bl_error("ERROR: species not defined")
  endif



  ! set the composition of the underlying star
  xn_star(:) = smallx
  xn_star(iash1) = ash1_frac
  if (ash2_name /= "") xn_star(iash2) = ash2_frac
  if (ash3_name /= "") xn_star(iash3) = ash3_frac

  ! and the composition of the accreted layer
  xn_base(:) = smallx
  xn_base(ifuel1) = fuel1_frac
  if (fuel2_name /= "") xn_base(ifuel2) = fuel2_frac
  if (fuel3_name /= "") xn_base(ifuel3) = fuel3_frac

  ! check if they sum to 1
  if (abs(sum(xn_star) - ONE) > nspec*smallx) then
     call bl_error("ERROR: ash mass fractions don't sum to 1")
  endif

  if (abs(sum(xn_base) - ONE) > nspec*smallx) then
     call bl_error("ERROR: fuel mass fractions don't sum to 1")
  endif



!-----------------------------------------------------------------------------
! Create a 1-d uniform grid that is identical to the mesh that we are
! mapping onto, and then we want to force it into HSE on that mesh.
!-----------------------------------------------------------------------------

  ! allocate the storage in the model_parser module
  allocate(model_r(nx))
  allocate(model_state(nx,nvars_model))
  npts_model = nx

  ! compute the coordinates of the new gridded function
  dCoord = (xmax - xmin) / dble(nx)

  do i = 1, nx
     model_r(i)  = xmin + (dble(i) - HALF)*dCoord
  enddo


  ! find the index of the base height
  index_base = -1
  do i = 1, nx
     if (model_r(i) >= xmin + H_star + atm_delta) then
        index_base = i+1
        exit
     endif
  enddo

  if (index_base == -1) then
     print *, 'ERROR: base_height not found on grid'
     call bl_error('ERROR: invalid base_height')
  endif


!-----------------------------------------------------------------------------
! put the model onto our new uniform grid
!-----------------------------------------------------------------------------

  fluff = .false.

  ! determine the conditions at the base
  eos_state%T     = T_base
  eos_state%rho   = dens_base
  eos_state%xn(:) = xn_base(:)

  call eos(eos_input_rt, eos_state)

  ! store the conditions at the base -- we'll use the entropy later
  ! to constrain the isentropic layer
  pres_base = eos_state%p
  entropy_base = eos_state%s

  print *, 'entropy_base = ', entropy_base
  print *, 'pres_base = ', pres_base

  ! set an initial temperature profile and composition
  do i = 1, nx

     ! hyperbolic tangent transition:
     model_state(i,ispec_model:ispec_model-1+nspec) = xn_star(1:nspec) + &
          HALF*(xn_base(1:nspec) - xn_star(1:nspec))* &
          (ONE + tanh((model_r(i) - (xmin + H_star - atm_delta) + atm_delta)/atm_delta))

     ! force them to sum to 1
     sumX = sum(model_state(i,ispec_model:ispec_model-1+nspec))
     model_state(i,ispec_model:ispec_model-1+nspec) = model_state(i,ispec_model:ispec_model-1+nspec) / sumX

     model_state(i,itemp_model) = T_star + HALF*(T_base - T_star)* &
          (ONE + tanh((model_r(i) - (xmin + H_star - atm_delta) + atm_delta)/atm_delta))


     ! the density and pressure will be determined via HSE,
     ! for now, set them to the base conditions
     model_state(i,idens_model) = dens_base
     model_state(i,ipres_model) = pres_base

  enddo


  if (index_base_from_temp) then
     ! find the index of the base height -- look at the temperature for this
     index_base = -1
     do i = 1, nx
        !if (model_r(i) >= xmin + H_star + atm_delta) then
        if (model_state(i,itemp_model) > 0.9995*T_base) then
           index_base = i+1
           exit
        endif
     enddo

     if (index_base == -1) then
        print *, 'ERROR: base_height not found on grid'
        call bl_error('ERROR: invalid base_height')
     endif
  endif

  print *, 'index_base = ', index_base

  ! make the base thermodynamics consistent for this base point -- that is
  ! what we will integrate from!
  eos_state%rho = model_state(index_base,idens_model)
  eos_state%T = model_state(index_base,itemp_model)
  eos_state%xn(:) = model_state(index_base,ispec_model:ispec_model-1+nspec)

  call eos(eos_input_rt, eos_state)

  model_state(index_base,ipres_model) = eos_state%p


!-----------------------------------------------------------------------------
! HSE + entropy solve
!-----------------------------------------------------------------------------

! the HSE state will be done putting creating an isentropic state until
! the temperature goes below T_lo -- then we will do isothermal.
! also, once the density goes below low_density_cutoff, we stop HSE

  isentropic = .true.

  !---------------------------------------------------------------------------
  ! integrate up
  !---------------------------------------------------------------------------
  do i = index_base+1, nx

     delx = model_r(i) - model_r(i-1)

     ! we've already set initial guesses for density, temperature, and
     ! composition
     dens_zone = model_state(i,idens_model)
     temp_zone = model_state(i,itemp_model)
     xn(:) = model_state(i,ispec_model:nvars_model)


     !-----------------------------------------------------------------------
     ! iteration loop
     !-----------------------------------------------------------------------

     ! start off the Newton loop by saying that the zone has not converged
     converged_hse = .FALSE.

     if (.not. fluff) then

        do iter = 1, MAX_ITER

           if (isentropic) then

              ! get the pressure we want from the HSE equation, just the
              ! zone below the current.  Note, we are using an average of
              ! the density of the two zones as an approximation of the
              ! interface value -- this means that we need to iterate for
              ! find the density and pressure that are consistent

              ! furthermore, we need to get the entropy that we need,
              ! which will come from adjusting the temperature in
              ! addition to the density.

              ! HSE differencing
              p_want = model_state(i-1,ipres_model) + &
                   delx*0.5*(dens_zone + model_state(i-1,idens_model))*const_grav



              ! now we have two functions to zero:
              !   A = p_want - p(rho,T)
              !   B = entropy_base - s(rho,T)
              ! We use a two dimensional Taylor expansion and find the deltas
              ! for both density and temperature


              ! now we know the pressure and the entropy that we want, so we
              ! need to find the temperature and density through a two
              ! dimensional root find

              ! (t, rho) -> (p, s)
              eos_state%T     = temp_zone
              eos_state%rho   = dens_zone
              eos_state%xn(:) = xn(:)

              call eos(eos_input_rt, eos_state)

              entropy = eos_state%s
              pres_zone = eos_state%p

              dpt = eos_state%dpdt
              dpd = eos_state%dpdr
              dst = eos_state%dsdt
              dsd = eos_state%dsdr

              A = p_want - pres_zone
              B = entropy_base - entropy

              dtemp = ((dsd/(dpd-0.5*delx*const_grav))*A - B)/ &
                   (dsd*dpt/(dpd -0.5*delx*const_grav) - dst)

              drho = (A - dpt*dtemp)/(dpd - 0.5*delx*const_grav)

              dens_zone = max(0.9_dp_t*dens_zone, &
                   min(dens_zone + drho, 1.1_dp_t*dens_zone))

              temp_zone = max(0.9_dp_t*temp_zone, &
                   min(temp_zone + dtemp, 1.1_dp_t*temp_zone))


              ! check if the density falls below our minimum cut-off --
              ! if so, floor it
              if (dens_zone < low_density_cutoff) then

                 dens_zone = low_density_cutoff
                 temp_zone = T_lo
                 converged_hse = .TRUE.
                 fluff = .TRUE.
                 exit

              endif

              ! if (A < TOL .and. B < ETOL) then
              if (abs(drho) < TOL*dens_zone .and. &
                  abs(dtemp) < TOL*temp_zone) then
                 converged_hse = .TRUE.
                 exit
              endif

           else

              ! do isothermal
              p_want = model_state(i-1,ipres_model) + &
                   delx*0.5*(dens_zone + model_state(i-1,idens_model))*const_grav

              temp_zone = T_lo

              ! (t, rho) -> (p)
              eos_state%T   = temp_zone
              eos_state%rho = dens_zone
              eos_state%xn(:) = xn(:)

              call eos(eos_input_rt, eos_state)

              entropy = eos_state%s
              pres_zone = eos_state%p

              dpd = eos_state%dpdr

              drho = (p_want - pres_zone)/(dpd - 0.5*delx*const_grav)

              dens_zone = max(0.9*dens_zone, &
                   min(dens_zone + drho, 1.1*dens_zone))

              if (abs(drho) < TOL*dens_zone) then
                 converged_hse = .TRUE.
                 exit
              endif

              if (dens_zone < low_density_cutoff) then

                 dens_zone = low_density_cutoff
                 temp_zone = T_lo
                 converged_hse = .TRUE.
                 fluff = .TRUE.
                 exit

              endif

           endif

           if (temp_zone < T_lo) then
              temp_zone = T_lo
              isentropic = .false.
           endif

        enddo


        if (.NOT. converged_hse) then

           print *, 'Error zone', i, ' did not converge in init_1d'
           print *, 'integrate up'
           print *, dens_zone, temp_zone
           print *, p_want, entropy_base, entropy
           print *, drho, dtemp
           call bl_error('Error: HSE non-convergence')

        endif

     else
        dens_zone = low_density_cutoff
        temp_zone = T_lo
     endif


     ! call the EOS one more time for this zone and then go on to the next
     ! (t, rho) -> (p)
     eos_state%T     = temp_zone
     eos_state%rho   = dens_zone
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)

     pres_zone = eos_state%p

     ! update the thermodynamics in this zone
     model_state(i,idens_model) = dens_zone
     model_state(i,itemp_model) = temp_zone
     model_state(i,ipres_model) = pres_zone


     ! to make this process converge faster, set the density in the
     ! next zone to the density in this zone
     ! model_state(i+1,idens_model) = dens_zone

  enddo


  !---------------------------------------------------------------------------
  ! integrate down -- using the temperature profile defined above
  !---------------------------------------------------------------------------
  do i = index_base-1, 1, -1

     delx = model_r(i+1) - model_r(i)

     ! we already set the temperature and composition profiles
     temp_zone = model_state(i,itemp_model)
     xn(:) = model_state(i,ispec_model:nvars_model)

     ! use our previous initial guess for density
     dens_zone = model_state(i+1,idens_model)


     !-----------------------------------------------------------------------
     ! iteration loop
     !-----------------------------------------------------------------------

     ! start off the Newton loop by saying that the zone has not converged
     converged_hse = .FALSE.

     do iter = 1, MAX_ITER

        ! get the pressure we want from the HSE equation, just the
        ! zone below the current.  Note, we are using an average of
        ! the density of the two zones as an approximation of the
        ! interface value -- this means that we need to iterate for
        ! find the density and pressure that are consistent

        ! HSE differencing
        p_want = model_state(i+1,ipres_model) - &
             delx*0.5*(dens_zone + model_state(i+1,idens_model))*const_grav


        ! we will take the temperature already defined in model_state
        ! so we only need to zero:
        !   A = p_want - p(rho)

        ! (t, rho) -> (p)
        eos_state%T     = temp_zone
        eos_state%rho   = dens_zone
        eos_state%xn(:) = xn(:)

        call eos(eos_input_rt, eos_state)

        pres_zone = eos_state%p

        dpd = eos_state%dpdr

        A = p_want - pres_zone

        drho = A/(dpd + 0.5*delx*const_grav)

        dens_zone = max(0.9_dp_t*dens_zone, &
             min(dens_zone + drho, 1.1_dp_t*dens_zone))


        if (abs(drho) < TOL*dens_zone) then
           converged_hse = .TRUE.
           exit
        endif

     enddo

     if (.NOT. converged_hse) then

        print *, 'Error zone', i, ' did not converge in init_1d'
        print *, 'integrate down'
        print *, dens_zone, temp_zone
        print *, p_want
        print *, drho
        call bl_error('Error: HSE non-convergence')

     endif


     ! call the EOS one more time for this zone and then go on to the next
     ! (t, rho) -> (p)
     eos_state%T     = temp_zone
     eos_state%rho   = dens_zone
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)

     pres_zone = eos_state%p

     ! update the thermodynamics in this zone
     model_state(i,idens_model) = dens_zone
     model_state(i,itemp_model) = temp_zone
     model_state(i,ipres_model) = pres_zone

  enddo

  do i = 1, nx
     print *, model_r(i), model_state(i,idens_model), model_state(i,ispec_model:ispec_model-1+nspec)
  enddo

end subroutine init_1d_tanh
