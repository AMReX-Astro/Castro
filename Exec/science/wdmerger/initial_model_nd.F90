! This model contains routines to:
! Generate a 1D isothermal white dwarf in hydrostatic equilibrium
! Interpolate a 1D WD model onto a 3D Cartesian grid

module initial_model_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: ZERO, HALF, ONE, FOUR3RD, M_PI
  use eos_type_module, only: eos_t, eos_input_rt
  use network, only: nspec
  use model_parser_module, only: itemp_model, idens_model, ipres_model, ispec_model
  use fundamental_constants_module, only: Gconst, M_solar
  use meth_params_module, only: small_temp
  use iso_c_binding

  integer, parameter :: initial_model_max_npts = 10000

  type, bind(C) :: initial_model

     ! Physical characteristics

     real(rt) :: mass
     real(rt) :: envelope_mass
     real(rt) :: central_density
     real(rt) :: central_temp
     real(rt) :: min_density
     real(rt) :: radius

     ! Composition

     real(rt) :: core_comp(nspec)
     real(rt) :: envelope_comp(nspec)

     ! Model storage

     real(rt) :: dx
     integer  :: npts
     real(rt) :: mass_tol, hse_tol

     real(rt) :: r(initial_model_max_npts), rl(initial_model_max_npts), rr(initial_model_max_npts)
     real(rt) :: M_enclosed(initial_model_max_npts), g(initial_model_max_npts)

  end type initial_model

  ! 1D initial models

  type (initial_model) :: model_P, model_S
  real(rt) :: rho_P(initial_model_max_npts), rho_S(initial_model_max_npts)
  real(rt) :: T_P(initial_model_max_npts), T_S(initial_model_max_npts)
  real(rt) :: xn_P(initial_model_max_npts,nspec), xn_S(initial_model_max_npts,nspec)
  real(rt) :: r_P(initial_model_max_npts), r_S(initial_model_max_npts)

contains

  subroutine initialize_model(primary, dx, npts, mass_tol, hse_tol)

    use castro_error_module, only: castro_error

    implicit none

    logical,  intent(in   ) :: primary
    integer,  intent(in   ) :: npts
    real(rt), intent(in   ) :: dx, mass_tol, hse_tol

    integer :: i

    if (npts > initial_model_max_npts) then
       call castro_error("npts too large, please increase initial_model_max_npts")
    end if

    if (primary) then

       model_P % mass = ZERO
       model_P % envelope_mass = ZERO
       model_P % central_density = ZERO
       model_P % central_temp = ZERO
       model_P % min_density = ZERO
       model_P % radius = ZERO

       model_P % core_comp(:) = ZERO
       model_P % envelope_comp(:) = ZERO

       model_P % dx = dx
       model_P % npts = npts
       model_P % mass_tol = mass_tol
       model_P % hse_tol = hse_tol

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

       model_S % core_comp(:) = ZERO
       model_S % envelope_comp(:) = ZERO

       model_S % dx = dx
       model_S % npts = npts
       model_S % mass_tol = mass_tol
       model_S % hse_tol = hse_tol

       do i = 1, npts

          model_S % rl(i) = (dble(i) - ONE)*dx
          model_S % rr(i) = (dble(i)      )*dx
          model_S % r(i)  = HALF*(model_S % rl(i) + model_S % rr(i))
          r_S(i) = model_S % r(i)

       end do

    end if

  end subroutine initialize_model



  ! Takes a one-dimensional stellar model and interpolates it to a point in
  ! 3D Cartesian grid space. Optionally does a sub-grid-scale interpolation
  ! if nsub > 1 (set in the probin file).

  subroutine interpolate_3d_from_1d(rho, T, xn, r, npts, loc, star_radius, dx, state, nsub_in)

    use interpolate_module, only: interpolate ! function
    use eos_module, only: eos

    implicit none

    real(rt),     intent(in   ) :: rho(initial_model_max_npts)
    real(rt),     intent(in   ) :: T(initial_model_max_npts)
    real(rt),     intent(in   ) :: xn(initial_model_max_npts, nspec)
    real(rt),     intent(in   ) :: r(initial_model_max_npts), star_radius
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
