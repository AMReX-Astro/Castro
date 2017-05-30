! This model contains routines to:
! Generate a 1D isothermal white dwarf in hydrostatic equilibrium
! Interpolate a 1D WD model onto a 3D Cartesian grid

module initial_model_module

  use bl_types
  use bl_constants_module
  use bl_error_module, only: bl_error
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rt
  use network, only: nspec
  use model_parser_module, only: itemp_model, idens_model, ipres_model, ispec_model
  use fundamental_constants_module, only: Gconst, M_solar
  use interpolate_module, only: interpolate
  use meth_params_module, only: small_temp

  type :: initial_model

     ! Physical characteristics

     double precision :: mass = ZERO
     double precision :: envelope_mass = ZERO
     double precision :: central_density = ZERO
     double precision :: central_temp = ZERO
     double precision :: min_density = ZERO
     double precision :: radius = ZERO

     ! Composition

     double precision :: core_comp(nspec) = ZERO
     double precision :: envelope_comp(nspec) = ZERO

     ! Model storage

     double precision :: dx
     integer          :: npts
     double precision :: mass_tol, hse_tol

     double precision, allocatable :: r(:), rl(:), rr(:)
     double precision, allocatable :: M_enclosed(:), g(:)
     type (eos_t), allocatable :: state(:)

  end type initial_model
  
contains

  subroutine initialize_model(model, dx, npts, mass_tol, hse_tol)

    implicit none

    type (initial_model) :: model
    integer :: npts
    double precision :: dx, mass_tol, hse_tol

    integer :: i

    model % dx = dx
    model % npts = npts
    model % mass_tol = mass_tol
    model % hse_tol = hse_tol

    allocate(model % r(npts))
    allocate(model % rl(npts))
    allocate(model % rr(npts))
    allocate(model % M_enclosed(npts))
    allocate(model % g(npts))
    allocate(model % state(npts))

    do i = 1, npts

       model % rl(i) = (dble(i) - ONE)*dx
       model % rr(i) = (dble(i)      )*dx
       model % r(i)  = HALF*(model % rl(i) + model % rr(i))

    enddo

  end subroutine initialize_model



  subroutine establish_hse(model)

    implicit none

    ! Arguments

    type (initial_model) :: model

    ! Local variables

    integer :: i, icutoff

    double precision :: rho_c, rho_c_old, drho_c, mass, mass_old, radius

    double precision :: p_want, rho_avg, rho, drho

    integer :: max_hse_iter = 250, max_mass_iter

    integer :: hse_iter, mass_iter

    logical :: converged_hse, fluff, mass_converged

    ! Note that if central_density > 0, then this initial model generator will use it in calculating
    ! the model. If mass is also provided in this case, we assume it is an estimate used for the purpose of 
    ! determining the envelope mass boundary. 

    ! Check to make sure we've specified at least one of them.

    if (model % mass < ZERO .and. model % central_density < ZERO) then
       call bl_error('Error: Must specify either mass or central density in the initial model generator.')
    endif

    ! If we are specifying the mass, then we don't know what WD central density 
    ! will give the desired total mass, so we need to do a secant iteration 
    ! over central density. rho_c_old is the 'old' guess for the central 
    ! density and rho_c is the current guess.  After two loops, we can 
    ! start estimating the density required to yield our desired mass.

    ! If instead we are specifying the central density, then we only need to do a 
    ! single HSE integration.

    if (model % central_density > ZERO) then

       max_mass_iter = 1

       rho_c_old = model % central_density
       rho_c     = model % central_density

    else

       max_mass_iter = max_hse_iter

       rho_c_old = -ONE
       rho_c     = 1.d7     ! A reasonable starting guess for moderate-mass WDs

    endif

    ! Check to make sure the initial temperature makes sense.

    if (model % central_temp < small_temp) then
       call bl_error("Error: WD central temperature is less than small_temp. Aborting.")
    endif

    mass_converged = .false.

    do mass_iter = 1, max_mass_iter

       fluff = .false.

       ! We start at the center of the WD and integrate outward.  Initialize
       ! the central conditions.

       model % state(1) % T    = model % central_temp
       model % state(1) % rho  = rho_c
       model % state(1) % xn   = model % core_comp

       call eos(eos_input_rt, model % state(1))

       ! Make the initial guess be completely uniform.

       model % state(:) = model % state(1)

       ! Keep track of the mass enclosed below the current zone.

       model % M_enclosed(1) = FOUR3RD * M_PI * (model % rr(1)**3 - model % rl(1)**3) * model % state(1) % rho

       !-------------------------------------------------------------------------
       ! HSE solve
       !-------------------------------------------------------------------------
       do i = 2, model % npts

          ! As the initial guess for the density, use the underlying zone.

          model % state(i) % rho = model % state(i-1) % rho

          if (model % mass > ZERO .and. model % M_enclosed(i-1) .ge. model % mass - model % envelope_mass) then
             model % state(i) % xn = model % envelope_comp
          else
             model % state(i) % xn = model % core_comp
          endif

          model % g(i) = -Gconst * model % M_enclosed(i-1) / model % rl(i)**2


          !----------------------------------------------------------------------
          ! Iteration loop
          !----------------------------------------------------------------------

          ! Start off the Newton loop by assuming that the zone has not converged.

          converged_hse = .FALSE.

          do hse_iter = 1, max_hse_iter

             if (fluff) then
                model % state(i) % rho = model % min_density
                exit
             endif

             ! The core is isothermal, so we just need to constrain
             ! the density and pressure to agree with the EOS and HSE.

             ! We difference HSE about the interface between the current
             ! zone and the one just inside.

             rho_avg = HALF * (model % state(i) % rho + model % state(i-1) % rho)
             p_want = model % state(i-1) % p + model % dx * rho_avg * model % g(i)

             call eos(eos_input_rt, model % state(i))

             drho = (p_want - model % state(i) % p) / (model % state(i) % dpdr - HALF * model % dx * model % g(i))
             rho = model % state(i) % rho

             rho = max(0.9 * rho, min(rho + drho, 1.1 * rho))

             model % state(i) % rho = rho

             if (rho < model % min_density) then
                icutoff = i
                fluff = .TRUE.
             endif

             if (abs(drho) < model % hse_tol * rho) then
                converged_hse = .TRUE.
                exit
             endif

          enddo

          if (.NOT. converged_hse .and. (.NOT. fluff)) then

             print *, 'Error: zone', i, ' did not converge in init_hse().'
             print *, model % state(i) % rho, model % state (i) % T
             print *, p_want, model % state(i) % p
             print *, drho, model % hse_tol * model % state(i) % rho
             call bl_error('Error: HSE non-convergence.')

          endif

          ! Call the EOS to establish the final properties of this zone.

          call eos(eos_input_rt, model % state(i))

          ! Discretize the mass enclose as (4 pi / 3) * rho * dr * (rl**2 + rl * rr + rr**2).

          model % M_enclosed(i) = model % M_enclosed(i-1) + &
                                  FOUR3RD * M_PI * model % state(i) % rho * model % dx * &
                                  (model % rr(i)**2 + model % rl(i) * model % rr(i) + model % rl(i)**2)

       enddo  ! End loop over zones

       mass = model % M_enclosed(icutoff)
       radius = model % r(icutoff)

       if (rho_c_old < ZERO) then

          ! Not enough iterations yet -- use an arbitrary guess for the next iteration.

          rho_c_old = rho_c
          rho_c = HALF * rho_c_old

       else

          ! Check if we have converged.

          if ( abs(mass - model % mass) / model % mass < model % mass_tol ) then
             mass_converged = .true.
             exit
          endif

          ! Do a secant iteration:
          ! M_tot = M(rho_c) + dM/drho |_rho_c x drho + ...        

          drho_c = (model % mass - mass) / ( (mass  - mass_old) / (rho_c - rho_c_old) )

          rho_c_old = rho_c
          rho_c = min(1.1d0 * rho_c_old, max((rho_c + drho_c), 0.9d0 * rho_c_old))

       endif     

       mass_old = mass

    enddo  ! End mass constraint loop

    if (.not. mass_converged .and. max_mass_iter .gt. 1) then
       call bl_error("ERROR: WD mass did not converge.")
    endif

    model % central_density = model % state(1) % rho
    model % radius = radius
    model % mass = mass

  end subroutine establish_hse



  ! Takes a one-dimensional stellar model and interpolates it to a point in
  ! 3D Cartesian grid space. Optionally does a sub-grid-scale interpolation
  ! if nsub > 1 (set in the probin file).

  subroutine interpolate_3d_from_1d(rho, T, xn, r, npts, loc, dx, state, nsub_in)

    implicit none

    double precision,  intent(in   ) :: rho(npts)
    double precision,  intent(in   ) :: T(npts)
    double precision,  intent(in   ) :: xn(npts, nspec)
    double precision,  intent(in   ) :: r(npts)
    integer,           intent(in   ) :: npts
    double precision,  intent(in   ) :: loc(3), dx(3)
    type (eos_t),      intent(inout) :: state
    integer, optional, intent(in   ) :: nsub_in
    
    integer :: i, j, k, n
    integer :: nsub
    double precision :: x, y, z, dist

    if (present(nsub_in)) then
       nsub = nsub_in
    else
       nsub = 1
    endif
    
    state % rho = ZERO 
    state % p   = ZERO
    state % T   = ZERO
    state % xn  = ZERO

    ! We perform a sub-grid-scale interpolation, where 
    ! nsub determines the number of intervals we split the zone into.
    ! If nsub = 1, we simply interpolate using the cell-center location.
    ! As an example, if nsub = 3, then the sampled locations will be
    ! k = 0 --> z = loc(3) - dx(3) / 3   (1/6 of way from left edge of zone)
    ! k = 1 --> z = loc(3)               (halfway between left and right edge)
    ! k = 2 --> z = loc(3) + dx(3) / 3   (1/6 of way from right edge of zone)

    do k = 0, nsub-1
       z = loc(3) + dble(k + HALF * (1 - nsub)) * dx(3) / nsub

       do j = 0, nsub-1
          y = loc(2) + dble(j + HALF * (1 - nsub)) * dx(2) / nsub

          do i = 0, nsub-1
             x = loc(1) + dble(i + HALF * (1 - nsub)) * dx(1) / nsub

             dist = (x**2 + y**2 + z**2)**HALF

             state % rho = state % rho + interpolate(dist, npts, r, rho)
             state % T   = state % T   + interpolate(dist, npts, r, T)

             do n = 1, nspec
                state % xn(n) = state % xn(n) + interpolate(dist, npts, r, xn(:,n))
             enddo

          enddo
       enddo
    enddo

    ! Now normalize by the number of intervals, and complete the thermodynamics.

    state % rho = state % rho / (nsub**3)
    state % T   = state % T   / (nsub**3)
    state % xn  = state % xn  / (nsub**3)

    call eos(eos_input_rt, state)
                    
  end subroutine interpolate_3d_from_1d
  
end module initial_model_module
