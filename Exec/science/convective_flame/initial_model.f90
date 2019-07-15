!!  Create a 1-d hydrostatic, atmosphere with an isothermal region
!!  (T_star) representing the NS, a hyperbolic tangent rise to a
!!  peak temperature (T_hi) representing the base of an accreted
!!  layer, an isoentropic profile down to a lower temperature (T_lo),
!!  and then isothermal. This can serve as an initial model for a
!!  nova or XRB.
!!
!!  The temperature profile is:
!!
!!         ^
!!         |
!!  T_hi   +            /\
!!         |           /  \
!!         |          /  . \
!!  T_star +---------+      \
!!         |         .   .   \
!!         |                  \
!!         |         .   .     \
!!  T_lo   +                    +-----------
!!         |         .   .
!!         +---------+---+---------------> r
!!         |         \  /
!!         |       atm_delta
!!         |< H_star>|
!!
!!
!!                   ^
!!                   |
!!                   +-- dens_base
!!
!!  dens_base is the density at a height H_star -- just below the rise
!!  in T up to the peak T_hi.  The composition is "ash" in the lower
!!  isothermal region and "fuel" in the isentropic and upper
!!  isothermal regions.  In the transition region, we apply the same
!!  hyperbolic tangent profile to interpolate the composition.
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

! this version allows for mutiple initial models


module initial_model_module

  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec

  implicit none

  ! integer keys for indexing the model_state array
  integer, parameter :: nvars_model = 3 + nspec
  integer, parameter :: idens_model = 1
  integer, parameter :: itemp_model = 2
  integer, parameter :: ipres_model = 3
  integer, parameter :: ispec_model = 4

  ! number of points in the model file
  integer, save :: gen_npts_model, num_models

  ! arrays for storing the model data -- we have an extra index here
  ! which is the model number
  real(rt), allocatable, save :: gen_model_state(:,:,:)
  real(rt), allocatable, save :: gen_model_r(:,:)

  type :: model_t

     real(rt) :: xn_base(nspec)
     real(rt) :: xn_star(nspec)
     real(rt) :: xn_perturb(nspec)

     real(rt) :: dens_base
     real(rt) :: T_star
     real(rt) :: T_hi
     real(rt) :: T_lo

     real(rt) :: H_star
     real(rt) :: atm_delta

     real(rt) :: low_density_cutoff

     logical :: index_base_from_temp

  end type model_t

contains

  subroutine init_model_data(nx, num_models_in)

    implicit none

    integer, intent(in) :: nx, num_models_in

    ! allocate storage for the model data
    allocate (gen_model_state(nx, nvars_model, num_models_in))
    allocate (gen_model_r(nx, num_models_in))

    gen_npts_model = nx
    num_models = num_models_in

  end subroutine init_model_data


  subroutine init_1d_tanh(nx, xmin, xmax, model_params, model_num)

    use amrex_constants_module
    use castro_error_module
    use amrex_fort_module, only : rt => amrex_real

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only : nspec, network_species_index, spec_names
    use fundamental_constants_module, only: Gconst
    use meth_params_module, only : const_grav

    use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor

    implicit none

    integer, intent(in) :: nx
    type(model_t), intent(in) :: model_params
    integer, intent(in) :: model_num
    real(rt) :: xmin, xmax

    integer :: i, n

    real (rt) :: slope_T, slope_xn(nspec)

    real (rt) :: pres_base, entropy_base

    real :: A, B

    real (rt) :: dCoord

    real (rt) :: dens_zone, temp_zone, pres_zone, entropy
    real (rt) :: dpd, dpt, dsd, dst

    real (rt) :: p_want, drho, dtemp, delx

    real (rt), parameter :: TOL = 1.e-10

    integer, parameter :: MAX_ITER = 250

    integer :: iter

    logical :: converged_hse, fluff

    real (rt), dimension(nspec) :: xn

    integer :: index_base

    logical :: isentropic, flipped

    integer :: narg

    type (eos_t) :: eos_state

    real(rt) :: sumX, xc


    !-----------------------------------------------------------------------------
    ! Create a 1-d uniform grid that is identical to the mesh that we are
    ! mapping onto, and then we want to force it into HSE on that mesh.
    !-----------------------------------------------------------------------------

    gen_npts_model = nx

    ! compute the coordinates of the new gridded function
    dCoord = (xmax - xmin) / dble(nx)

    do i = 1, nx
       gen_model_r(i, model_num)  = xmin + (dble(i) - HALF)*dCoord
    enddo


    ! find the index of the base height
    index_base = -1
    do i = 1, nx
       if (gen_model_r(i, model_num) >= xmin + model_params % H_star) then
          index_base = i+1
          exit
       endif
    enddo

    if (index_base == -1) then
       print *, 'ERROR: base_height not found on grid'
       call castro_error('ERROR: invalid base_height')
    endif


    !-----------------------------------------------------------------------------
    ! put the model onto our new uniform grid
    !-----------------------------------------------------------------------------
    fluff = .false.

    ! determine the conditions at the base -- this is below the atmosphere
    eos_state%T     = model_params % T_star
    eos_state%rho   = model_params % dens_base
    eos_state%xn(:) = model_params % xn_star(:)

    call eos(eos_input_rt, eos_state)

    ! store the conditions at the base -- we'll use the entropy later
    ! to constrain the isentropic layer
    pres_base = eos_state%p

    ! set an initial temperature profile and composition
    print *, "in model stuff", model_params % T_hi

    do i = 1, nx

       xc = gen_model_r(i,model_num) - (xmin + model_params % H_star) - 1.5_rt * model_params % atm_delta

       ! hyperbolic tangent transition:
       gen_model_state(i,ispec_model:ispec_model-1+nspec,model_num) = model_params % xn_star(1:nspec) + &
            HALF*(model_params % xn_base(1:nspec) - model_params % xn_star(1:nspec))* &
            (ONE + tanh(xc/(HALF*model_params % atm_delta)))

       ! force them to sum to 1
       sumX = sum(gen_model_state(i,ispec_model:ispec_model-1+nspec,model_num))
       gen_model_state(i,ispec_model:ispec_model-1+nspec,model_num) = gen_model_state(i,ispec_model:ispec_model-1+nspec,model_num) / sumX

       gen_model_state(i,itemp_model,model_num) = model_params % T_star + &
            HALF*(model_params % T_hi - model_params % T_star)* &
            (ONE + tanh(xc/(HALF*model_params % atm_delta)))

       gen_model_state(1:index_base,itemp_model,model_num) = model_params % T_star

       ! the density and pressure will be determined via HSE,
       ! for now, set them to the base conditions
       gen_model_state(i,idens_model,model_num) = model_params % dens_base
       gen_model_state(i,ipres_model,model_num) = pres_base

    enddo

    ! make the base thermodynamics consistent for this base point -- that is
    ! what we will integrate from!
    eos_state%rho = gen_model_state(index_base,idens_model,model_num)
    eos_state%T = gen_model_state(index_base,itemp_model,model_num)
    eos_state%xn(:) = gen_model_state(index_base,ispec_model:ispec_model-1+nspec,model_num)

    call eos(eos_input_rt, eos_state)

    gen_model_state(index_base,ipres_model,model_num) = eos_state%p


    !-----------------------------------------------------------------------------
    ! HSE + entropy solve
    !-----------------------------------------------------------------------------

    ! the HSE state will be done putting creating an isentropic state until
    ! the temperature goes below T_lo -- then we will do isothermal.
    ! also, once the density goes below low_density_cutoff, we stop HSE

    isentropic = .false.
    flipped = .false.   ! we start out isothermal and then 'flip' to isentropic

    !---------------------------------------------------------------------------
    ! integrate up
    !---------------------------------------------------------------------------
    do i = index_base+1, nx

       if ((gen_model_r(i,model_num) > xmin + model_params % H_star + 3.0_rt * model_params % atm_delta) .and. .not. flipped) then
          isentropic = .true.
          flipped = .true.

          ! now we need to know the entropy we are confining ourselves to
          eos_state % rho = gen_model_state(i-1,idens_model,model_num)
          eos_state % T = gen_model_state(i-1,itemp_model,model_num)
          eos_state % xn(:) = gen_model_state(i-1,ispec_model:ispec_model-1+nspec,model_num)

          call eos(eos_input_rt, eos_state)

          entropy_base = eos_state % s


          if (parallel_IOProcessor()) then
             print *, "base density = ", eos_state % rho, eos_state % T
          endif
       endif

       delx = gen_model_r(i,model_num) - gen_model_r(i-1,model_num)

       ! we've already set initial guesses for density, temperature, and
       ! composition
       dens_zone = gen_model_state(i,idens_model,model_num)
       temp_zone = gen_model_state(i,itemp_model,model_num)
       xn(:) = gen_model_state(i,ispec_model:nvars_model,model_num)


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
                p_want = gen_model_state(i-1,ipres_model,model_num) + &
                     delx*0.5*(dens_zone + gen_model_state(i-1,idens_model,model_num))*const_grav


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

                dens_zone = max(0.9_rt*dens_zone, &
                     min(dens_zone + drho, 1.1_rt*dens_zone))

                temp_zone = max(0.9_rt*temp_zone, &
                     min(temp_zone + dtemp, 1.1_rt*temp_zone))


                ! check if the density falls below our minimum cut-off --
                ! if so, floor it
                if (dens_zone < model_params % low_density_cutoff) then

                   dens_zone = model_params % low_density_cutoff
                   temp_zone = model_params % T_lo
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
                p_want = gen_model_state(i-1,ipres_model,model_num) + &
                     delx*0.5*(dens_zone + gen_model_state(i-1,idens_model,model_num))*const_grav

                if (gen_model_r(i,model_num) > xmin + model_params % H_star + 3.0_rt * model_params % atm_delta) then
                   temp_zone = model_params % T_lo
                endif

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

                if (dens_zone < model_params % low_density_cutoff) then

                   dens_zone = model_params % low_density_cutoff
                   temp_zone = model_params % T_lo
                   converged_hse = .TRUE.
                   fluff = .TRUE.
                   exit

                endif

             endif

             if (temp_zone < model_params % T_lo) then
                temp_zone = model_params % T_lo
                isentropic = .false.
             endif

          enddo


          if (.NOT. converged_hse) then

             print *, 'Error zone', i, ' did not converge in init_1d'
             print *, 'integrate up'
             print *, dens_zone, temp_zone
             print *, p_want, entropy_base, entropy
             print *, drho, dtemp
             call castro_error('Error: HSE non-convergence')

          endif

       else
          dens_zone = model_params % low_density_cutoff
          temp_zone = model_params % T_lo
       endif


       ! call the EOS one more time for this zone and then go on to the next
       ! (t, rho) -> (p)
       eos_state%T     = temp_zone
       eos_state%rho   = dens_zone
       eos_state%xn(:) = xn(:)

       call eos(eos_input_rt, eos_state)

       pres_zone = eos_state%p

       ! update the thermodynamics in this zone
       gen_model_state(i,idens_model,model_num) = dens_zone
       gen_model_state(i,itemp_model,model_num) = temp_zone
       gen_model_state(i,ipres_model,model_num) = pres_zone


       ! to make this process converge faster, set the density in the
       ! next zone to the density in this zone
       ! gen_model_state(i+1,idens_model) = dens_zone

    enddo


    !---------------------------------------------------------------------------
    ! integrate down -- using the temperature profile defined above
    !---------------------------------------------------------------------------
    do i = index_base-1, 1, -1

       delx = gen_model_r(i+1,model_num) - gen_model_r(i,model_num)

       ! we already set the temperature and composition profiles
       temp_zone = gen_model_state(i,itemp_model,model_num)
       xn(:) = gen_model_state(i,ispec_model:nvars_model,model_num)

       ! use our previous initial guess for density
       dens_zone = gen_model_state(i+1,idens_model,model_num)


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
          p_want = gen_model_state(i+1,ipres_model,model_num) - &
               delx*0.5*(dens_zone + gen_model_state(i+1,idens_model,model_num))*const_grav


          ! we will take the temperature already defined in gen_model_state
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

          dens_zone = max(0.9_rt*dens_zone, &
               min(dens_zone + drho, 1.1_rt*dens_zone))


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
          call castro_error('Error: HSE non-convergence')

       endif


       ! call the EOS one more time for this zone and then go on to the next
       ! (t, rho) -> (p)
       eos_state%T     = temp_zone
       eos_state%rho   = dens_zone
       eos_state%xn(:) = xn(:)

       call eos(eos_input_rt, eos_state)

       pres_zone = eos_state%p

       ! update the thermodynamics in this zone
       gen_model_state(i,idens_model,model_num) = dens_zone
       gen_model_state(i,itemp_model,model_num) = temp_zone
       gen_model_state(i,ipres_model,model_num) = pres_zone

    enddo

    !do i = 1, nx
    !   print *, gen_model_r(i), gen_model_state(i,idens_model), gen_model_state(i,ispec_model:ispec_model-1+nspec)
    !enddo

  end subroutine init_1d_tanh

end module initial_model_module
