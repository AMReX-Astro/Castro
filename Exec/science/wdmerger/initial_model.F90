! This model contains routines to:
! Generate a 1D isothermal white dwarf in hydrostatic equilibrium
! Interpolate a 1D WD model onto a 3D Cartesian grid

module initial_model_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use castro_error_module, only: castro_error
  use eos_type_module, only: eos_t, eos_input_rt
  use network, only: nspec
  use model_parser_module, only: itemp_model, idens_model, ipres_model, ispec_model
  use fundamental_constants_module, only: Gconst, M_solar
  use meth_params_module, only: small_temp
  use probdata_module

contains

  subroutine initialize_model(primary, dx, npts, mass_tol, hse_tol)

    implicit none

    logical,  intent(in   ) :: primary
    integer,  intent(in   ) :: npts
    real(rt), intent(in   ) :: dx, mass_tol, hse_tol

    integer :: i

    if (primary) then

       model_P % mass = ZERO
       model_P % envelope_mass = ZERO
       model_P % central_density = ZERO
       model_P % central_temp = ZERO
       model_P % min_density = ZERO
       model_P % radius = ZERO

       allocate(model_P % core_comp(nspec))
       allocate(model_P % envelope_comp(nspec))

       model_P % core_comp(:) = ZERO
       model_P % envelope_comp(:) = ZERO

       model_P % dx = dx
       model_P % npts = npts
       model_P % mass_tol = mass_tol
       model_P % hse_tol = hse_tol

       allocate(model_P % r(npts))
       allocate(model_P % rl(npts))
       allocate(model_P % rr(npts))
       allocate(model_P % M_enclosed(npts))
       allocate(model_P % g(npts))
       allocate(model_P % state(npts))

       allocate(rho_P(npts))
       allocate(T_P(npts))
       allocate(xn_P(npts,nspec))
       allocate(r_P(npts))

       do i = 1, npts

          model_P % rl(i) = (dble(i) - ONE)*dx
          model_P % rr(i) = (dble(i)      )*dx
          model_P % r(i)  = HALF*(model_P % rl(i) + model_P % rr(i))
          r_P(i) = model_P % r(i)

       end do

    else

       model_S % mass = ZERO
       model_S % envelope_mass = ZERO
       model_S % central_density = ZERO
       model_S % central_temp = ZERO
       model_S % min_density = ZERO
       model_S % radius = ZERO

       allocate(model_S % core_comp(nspec))
       allocate(model_S % envelope_comp(nspec))

       model_S % core_comp(:) = ZERO
       model_S % envelope_comp(:) = ZERO

       model_S % dx = dx
       model_S % npts = npts
       model_S % mass_tol = mass_tol
       model_S % hse_tol = hse_tol

       allocate(model_S % r(npts))
       allocate(model_S % rl(npts))
       allocate(model_S % rr(npts))
       allocate(model_S % M_enclosed(npts))
       allocate(model_S % g(npts))
       allocate(model_S % state(npts))

       allocate(rho_S(npts))
       allocate(T_S(npts))
       allocate(xn_S(npts,nspec))
       allocate(r_S(npts))

       do i = 1, npts

          model_S % rl(i) = (dble(i) - ONE)*dx
          model_S % rr(i) = (dble(i)      )*dx
          model_S % r(i)  = HALF*(model_S % rl(i) + model_S % rr(i))
          r_S(i) = model_S % r(i)

       end do

    end if

  end subroutine initialize_model



  subroutine establish_hse(model, rho, T, xn, r)

    use eos_module, only: eos_on_host

    implicit none

    ! Arguments

    type (initial_model), intent(inout) :: model
    real(rt), intent(inout) :: rho(model % npts)
    real(rt), intent(inout) :: T(model % npts)
    real(rt), intent(inout) :: xn(model % npts, nspec)
    real(rt), intent(inout) :: r(model % npts)

    ! Local variables

    integer :: i, icutoff, n

    real(rt) :: rho_c, rho_c_old, drho_c, mass, mass_old, radius

    real(rt) :: p_want, rho_avg, drho

    integer :: max_hse_iter = 250, max_mass_iter

    integer :: hse_iter, mass_iter

    logical :: converged_hse, fluff, mass_converged

    ! Note that if central_density > 0, then this initial model generator will use it in calculating
    ! the model. If mass is also provided in this case, we assume it is an estimate used for the purpose of 
    ! determining the envelope mass boundary.

    ! Check to make sure we've specified at least one of them.

    if (model % mass < ZERO .and. model % central_density < ZERO) then
       call castro_error('Error: Must specify either mass or central density in the initial model generator.')
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
       rho_c     = 1.e7_rt     ! A reasonable starting guess for moderate-mass WDs

    endif

    ! Check to make sure the initial temperature makes sense.

    if (model % central_temp < small_temp) then
       call castro_error("Error: WD central temperature is less than small_temp. Aborting.")
    endif

    mass_converged = .false.

    do mass_iter = 1, max_mass_iter

       fluff = .false.

       ! We start at the center of the WD and integrate outward.  Initialize
       ! the central conditions.

       T(1)    = model % central_temp
       rho(1)  = rho_c
       xn(1,:) = model % core_comp

       model % state(1) % rho  = rho(1)
       model % state(1) % T    = T(1)
       model % state(1) % xn   = xn(1,:)

       call eos_on_host(eos_input_rt, model % state(1))

       ! Make the initial guess be completely uniform.

       model % state(:) = model % state(1)
       rho(:) = rho(1)
       T(:)   = T(1)
       do n = 1, nspec
          xn(:,n) = xn(1,n)
       end do

       ! Keep track of the mass enclosed below the current zone.

       model % M_enclosed(1) = FOUR3RD * M_PI * (model % rr(1)**3 - model % rl(1)**3) * rho(1)

       !-------------------------------------------------------------------------
       ! HSE solve
       !-------------------------------------------------------------------------
       do i = 2, model % npts

          ! As the initial guess for the density, use the underlying zone.

          rho(i) = rho(i-1)
          model % state(i) % rho = model % state(i-1) % rho

          if (model % mass > ZERO .and. model % M_enclosed(i-1) .ge. model % mass - model % envelope_mass) then
             xn(i,:) = model % envelope_comp
             model % state(i) % xn = xn(i,:)
          else
             xn(i,:) = model % core_comp
             model % state(i) % xn = xn(i,:)
          endif

          model % g(i) = -Gconst * model % M_enclosed(i-1) / model % rl(i)**2


          !----------------------------------------------------------------------
          ! Iteration loop
          !----------------------------------------------------------------------

          ! Start off the Newton loop by assuming that the zone has not converged.

          converged_hse = .FALSE.

          do hse_iter = 1, max_hse_iter

             if (fluff) then
                rho(i) = model % min_density
                model % state(i) % rho = model % min_density
                exit
             endif

             ! The core is isothermal, so we just need to constrain
             ! the density and pressure to agree with the EOS and HSE.

             ! We difference HSE about the interface between the current
             ! zone and the one just inside.

             rho_avg = HALF * (rho(i) + rho(i-1))
             p_want = model % state(i-1) % p + model % dx * rho_avg * model % g(i)

             call eos_on_host(eos_input_rt, model % state(i))

             drho = (p_want - model % state(i) % p) / (model % state(i) % dpdr - HALF * model % dx * model % g(i))

             rho(i) = max(0.9 * rho(i), min(rho(i) + drho, 1.1 * rho(i)))
             model % state(i) % rho = rho(i)

             if (rho(i) < model % min_density) then
                icutoff = i
                fluff = .TRUE.
             endif

             if (abs(drho) < model % hse_tol * rho(i)) then
                converged_hse = .TRUE.
                exit
             endif

          enddo

          if (.NOT. converged_hse .and. (.NOT. fluff)) then

             print *, 'Error: zone', i, ' did not converge in init_hse().'
             print *, rho(i), T(i)
             print *, p_want, model % state(i) % p
             print *, drho, model % hse_tol * rho(i)
             call castro_error('Error: HSE non-convergence.')

          endif

          ! Call the EOS to establish the final properties of this zone.

          call eos_on_host(eos_input_rt, model % state(i))

          ! Discretize the mass enclose as (4 pi / 3) * rho * dr * (rl**2 + rl * rr + rr**2).

          model % M_enclosed(i) = model % M_enclosed(i-1) + &
                                  FOUR3RD * M_PI * rho(i) * model % dx * &
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
          rho_c = min(1.1e0_rt * rho_c_old, max((rho_c + drho_c), 0.9e0_rt * rho_c_old))

       endif

       mass_old = mass

    enddo  ! End mass constraint loop

    if (.not. mass_converged .and. max_mass_iter .gt. 1) then
       call castro_error("ERROR: WD mass did not converge.")
    endif

    model % central_density = rho(1)
    model % radius = radius
    model % mass = mass

  end subroutine establish_hse



  ! Takes a one-dimensional stellar model and interpolates it to a point in
  ! 3D Cartesian grid space. Optionally does a sub-grid-scale interpolation
  ! if nsub > 1 (set in the probin file).

  subroutine interpolate_3d_from_1d(rho, T, xn, r, npts, loc, star_radius, dx, state, nsub_in)

    use interpolate_module, only: interpolate ! function
    use eos_module, only: eos

    implicit none

    real(rt),     intent(in   ) :: rho(npts)
    real(rt),     intent(in   ) :: T(npts)
    real(rt),     intent(in   ) :: xn(npts, nspec)
    real(rt),     intent(in   ) :: r(npts), star_radius
    integer,      intent(in   ) :: npts
    real(rt),     intent(in   ) :: loc(3), dx(3)
    type (eos_t), intent(inout) :: state
    integer,      intent(in   ), optional :: nsub_in

    integer  :: i, j, k, n
    integer  :: nsub
    real(rt) :: x, y, z, dist

    !$gpu

    if (present(nsub_in)) then
       nsub = nsub_in
    else
       nsub = 1
    endif

    state % rho = ZERO
    state % p   = ZERO
    state % T   = ZERO
    state % xn  = ZERO

    ! If the model radius is smaller than the zone size, just use the center of the model.

    dist = (loc(1)**2 + loc(2)**2 + loc(3)**2)**HALF

    if (star_radius <= maxval(dx) .and. dist < maxval(dx)) then

       state % rho = rho(1)
       state % T   = T(1)
       state % xn  = xn(1,:)

    else

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

       ! Now normalize by the number of intervals.

       state % rho = state % rho / (nsub**3)
       state % T   = state % T   / (nsub**3)
       state % xn  = state % xn  / (nsub**3)

    end if

    ! Complete the thermodynamics.

    call eos(eos_input_rt, state)

  end subroutine interpolate_3d_from_1d

end module initial_model_module
