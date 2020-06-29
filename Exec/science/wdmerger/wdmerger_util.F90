module wdmerger_util_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  ! This routine calls all of the other subroutines at the beginning
  ! of a simulation to fill in the basic problem data.

  subroutine initialize_problem(init_in)

    use castro_error_module, only: castro_error
    use prob_params_module, only: dim
    use probdata_module, only: init

    implicit none

    integer :: init_in

    init = init_in

    ! We have already read in the namelist; do some parameter postprocessing.

    call finalize_probdata()

    ! Establish binary parameters and create initial models.

    call binary_setup()

    ! Set small_pres and small_ener.

    call set_small()

  end subroutine initialize_problem



  ! This routine reads in the namelist

  subroutine finalize_probdata ()

    use amrex_constants_module, only: ZERO, THIRD, HALF, ONE, TWO, THREE, M_PI, FOUR, SIX, EIGHT
    use meth_params_module, only: point_mass, do_rotation, rot_axis
    use prob_params_module, only: dim, coord_type
    use castro_error_module, only: castro_error
    use network, only: nspec
    use fundamental_constants_module, only: M_solar
    use initial_model_module, only: model_P, model_S
    use probdata_module, only: mass_P, mass_S, central_density_S, single_star, &
                               problem, axis_1, axis_2, &
                               max_he_wd_mass, max_hybrid_wd_mass, max_co_wd_mass, &
                               hybrid_wd_he_shell_mass, co_wd_he_shell_mass, &
                               collision_impact_parameter, orbital_angle, orbital_eccentricity, &
                               tde_beta, tde_separation, &
                               relaxation_damping_factor, relaxation_is_done

    implicit none

    ! Allocate parameters and set defaults.

    allocate(model_P)
    allocate(model_S)

    ! Convert masses from solar masses to grams.

    mass_P = mass_P * M_solar
    mass_S = mass_S * M_solar

    max_he_wd_mass = max_he_wd_mass * M_solar
    max_hybrid_wd_mass = max_hybrid_wd_mass * M_solar
    max_co_wd_mass = max_co_wd_mass * M_solar

    hybrid_wd_he_shell_mass = hybrid_wd_he_shell_mass * M_solar
    co_wd_he_shell_mass = co_wd_he_shell_mass * M_solar

    if (mass_S < ZERO .and. central_density_S < ZERO) single_star = .true.

    ! Make sure that the primary mass is really the larger mass

    call ensure_primary_mass_larger()

    ! We enforce that the orbital plane is (z, phi) for two-dimensional problems,
    ! with rotation about z = 0 along the radial axis.

    if (dim .eq. 2) then

       if (coord_type .ne. 1) then
          call castro_error("We only support cylindrical coordinates in two dimensions. Set coord_type == 1.")
       endif

       axis_1 = 2
       axis_2 = 3
       rot_axis = 1

    endif

    ! Make sure we have a sensible collision impact parameter.

    if (collision_impact_parameter > 1.0) then
       call castro_error("Impact parameter must be less than one in our specified units.")
    endif

    ! Safety check: we can't run most problems in one dimension.

    if (dim .eq. 1 .and. problem /= 0) then
       call castro_error("Can only run a collision or freefall in 1D. Exiting.")
    endif

    ! Don't do a collision, free-fall, or TDE in a rotating reference frame.

    if (problem .eq. 0 .and. do_rotation .eq. 1) then
       call castro_error("The free-fall/collision problem does not make sense in a rotating reference frame.")
    endif

    if (problem .eq. 2 .and. do_rotation .eq. 1) then
       call castro_error("The TDE problem does not make sense in a rotating reference frame.")
    end if

    ! Make sure we have a sensible eccentricity.

    if (orbital_eccentricity >= 1.0) then
       call castro_error("Orbital eccentricity cannot be larger than one.")
    endif

    ! Make sure we have a sensible angle. Then convert it to radians.

    if (orbital_angle < 0.0 .or. orbital_angle > 360.0) then
       call castro_error("Orbital angle must be between 0 and 360 degrees.")
    endif

    orbital_angle = orbital_angle * M_PI / 180.0

    ! If we're doing a relaxation, we need to reset the relaxation_is_done parameter.
    ! This will be reset as appropriate from the checkpoint if we're performing a restart.

    if (problem .eq. 1 .and. relaxation_damping_factor > ZERO) then
       relaxation_is_done = 0
    end if

    ! TDE sanity checks

    if (problem .eq. 2) then

       ! We must have a BH point mass defined.

       if (point_mass <= ZERO) then

          call castro_error("No point mass specified for the TDE problem.")

       end if

       ! This problem cannot have a secondary WD.

       if (mass_S >= ZERO) then

          call castro_error("TDE problem cannot have a secondary WD.")

       end if

       ! Beta parameter must be positive.

       if (tde_beta <= ZERO) then

          call castro_error("TDE beta must be positive.")

       end if

       ! Initial distance must be positive.

       if (tde_separation <= ZERO) then

          call castro_error("TDE separation must be positive.")

       end if

    end if

  end subroutine finalize_probdata



  ! Calculate small_pres and small_ener

  subroutine set_small

    use meth_params_module, only: small_temp, small_pres, small_dens, small_ener, &
                                  UFS, URHO
    use ambient_module, only: ambient_state
    use network, only: nspec
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_module, only: eos

    implicit none

    type (eos_t) :: eos_state

    ! Given the inputs of small_dens and small_temp, figure out small_pres.

    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = ambient_state(UFS:UFS+nspec-1) / ambient_state(URHO)

    call eos(eos_input_rt, eos_state)

    small_pres = eos_state % p
    small_ener = eos_state % e

  end subroutine set_small



  ! Given a WD mass, set its core and envelope composition.

  subroutine set_wd_composition(model)

    use extern_probin_module, only: small_x
    use network, only: network_species_index
    use castro_error_module, only: castro_error
    use amrex_constants_module, only: ZERO, ONE
    use initial_model_module, only: initial_model
    use probdata_module, only: max_he_wd_mass, max_hybrid_wd_mass, max_co_wd_mass, &
                               hybrid_wd_c_frac, hybrid_wd_o_frac, hybrid_wd_he_shell_mass, &
                               co_wd_c_frac, co_wd_o_frac, co_wd_he_shell_mass, &
                               onemg_wd_o_frac, onemg_wd_ne_frac, onemg_wd_mg_frac

    implicit none

    type (initial_model), intent(inout) :: model

    integer :: iHe4, iC12, iO16, iNe20, iMg24

    iHe4 = network_species_index("helium-4")
    iC12 = network_species_index("carbon-12")
    iO16 = network_species_index("oxygen-16")
    iNe20 = network_species_index("neon-20")
    iMg24 = network_species_index("magnesium-24")

    model % core_comp = small_x
    model % envelope_comp = small_x

    model % envelope_mass = ZERO

    ! Here we follow the prescription of Dan et al. 2012.

    if (model % mass > ZERO .and. model % mass < max_he_wd_mass) then

       if (iHe4 < 0) call castro_error("Must have He4 in the nuclear network.")

       model % core_comp(iHe4) = ONE

       model % envelope_comp = model % core_comp

    else if (model % mass >= max_he_wd_mass .and. model % mass < max_hybrid_wd_mass) then

       if (iC12 < 0) call castro_error("Must have C12 in the nuclear network.")
       if (iO16 < 0) call castro_error("Must have O16 in the nuclear network.")

       model % core_comp(iC12) = hybrid_wd_c_frac
       model % core_comp(iO16) = hybrid_wd_o_frac

       model % envelope_mass = hybrid_wd_he_shell_mass

       if (model % envelope_mass > ZERO) then
          if (iHe4 < 0) call castro_error("Must have He4 in the nuclear network.")
          model % envelope_comp(iHe4) = ONE
       else
          model % envelope_comp = model % core_comp
       endif

    else if (model % mass >= max_hybrid_wd_mass .and. model % mass < max_co_wd_mass) then

       if (iC12 < 0) call castro_error("Must have C12 in the nuclear network.")
       if (iO16 < 0) call castro_error("Must have O16 in the nuclear network.")

       model % core_comp(iC12) = co_wd_c_frac
       model % core_comp(iO16) = co_wd_o_frac

       model % envelope_mass = co_wd_he_shell_mass

       if (model % envelope_mass > ZERO) then
          if (iHe4 < 0) call castro_error("Must have He4 in the nuclear network.")
          model % envelope_comp(iHe4) = ONE
       else
          model % envelope_comp = model % core_comp
       endif

    else if (model % mass > max_co_wd_mass) then

       if (iO16 < 0) call castro_error("Must have O16 in the nuclear network.")
       if (iNe20 < 0) call castro_error("Must have Ne20 in the nuclear network.")
       if (iMg24 < 0) call castro_error("Must have Mg24 in the nuclear network.")

       model % core_comp(iO16)  = onemg_wd_o_frac
       model % core_comp(iNe20) = onemg_wd_ne_frac
       model % core_comp(iMg24) = onemg_wd_mg_frac

       model % envelope_comp = model % core_comp

    endif

    ! Normalize compositions so that they sum to one.

     model % core_comp = model % core_comp / sum(model % core_comp)
     model % envelope_comp = model % envelope_comp / sum(model % envelope_comp)

  end subroutine set_wd_composition



  ! This routine checks to see if the primary mass is actually larger
  ! than the secondary mass, and switches them if not.

  subroutine ensure_primary_mass_larger ()

    use problem_io_module, only: ioproc
    use probdata_module, only: mass_P, mass_S

    implicit none

    real(rt) :: temp_mass

    ! We want the primary WD to be more massive. If what we're calling
    ! the primary is less massive, switch the stars.

    if (mass_P < mass_S) then

      if (ioproc) then
        print *, "Primary mass is less than secondary mass; switching the stars so that the primary is more massive."
      endif

      temp_mass = mass_P
      mass_P = mass_S
      mass_S = temp_mass

    endif

  end subroutine ensure_primary_mass_larger



  ! Return the locations of the stellar centers of mass

  subroutine get_star_data(P_com, S_com, P_vel, S_vel, P_mass, S_mass, P_t_ff, S_t_ff) bind(C,name='get_star_data')

    use probdata_module, only: mass_P, mass_S, t_ff_P, t_ff_S, com_P, com_S, vel_P, vel_S

    implicit none

    real(rt), intent(inout) :: P_com(3), S_com(3)
    real(rt), intent(inout) :: P_vel(3), S_vel(3)
    real(rt), intent(inout) :: P_mass, S_mass
    real(rt), intent(inout) :: P_t_ff, S_t_ff

    P_com = com_P
    S_com = com_S

    P_vel = vel_P
    S_vel = vel_S

    P_mass = mass_P
    S_mass = mass_S

    P_t_ff = t_ff_P
    S_t_ff = t_ff_S

  end subroutine get_star_data



  ! Set the locations of the stellar centers of mass

  subroutine set_star_data(P_com, S_com, P_vel, S_vel, P_mass, S_mass, P_t_ff, S_t_ff) bind(C,name='set_star_data')

    use amrex_constants_module, only: TENTH, ZERO
    use prob_params_module, only: center
    use binary_module, only: get_roche_radii
    use probdata_module, only: mass_P, mass_S, t_ff_P, t_ff_S, com_P, com_S, vel_P, vel_S, roche_rad_P, roche_rad_S

    implicit none

    real(rt), intent(in) :: P_com(3), S_com(3)
    real(rt), intent(in) :: P_vel(3), S_vel(3)
    real(rt), intent(in) :: P_mass, S_mass
    real(rt), intent(in) :: P_t_ff, S_t_ff

    real(rt) :: r

    r = ZERO

    if (mass_P > ZERO) then

       com_P  = P_com
       vel_P  = P_vel
       mass_P = P_mass
       t_ff_P = P_t_ff

    endif

    if (mass_S > ZERO) then

       com_S  = S_com
       vel_S  = S_vel
       mass_S = S_mass
       t_ff_S = S_t_ff

    endif

    if (mass_P > ZERO .and. mass_S > ZERO) then

       r = sum((com_P-com_S)**2)**(0.5)

       call get_roche_radii(mass_S/mass_P, roche_rad_S, roche_rad_P, r)

       ! Beyond a certain point, it doesn't make sense to track the stars separately
       ! anymore. We'll set the secondary to a fixed constant and keep it there
       ! if its Roche radius becomes smaller than 10% of the primary's. Also, for exactly 
       ! equal mass systems sometimes it is the primary that disrupts, perhaps
       ! just due to numerical noise, so do the same check for the primary.

       if (roche_rad_S < TENTH * roche_rad_P) then
          com_S = center
          vel_S = ZERO
          mass_S = ZERO
          roche_rad_S = ZERO
          t_ff_S = ZERO
       else if (roche_rad_P < TENTH * roche_rad_S) then
          com_P = center
          vel_P = ZERO
          mass_S = ZERO
          roche_rad_P = ZERO
          t_ff_P = ZERO
       endif

    endif

  end subroutine set_star_data



  ! Set up a binary simulation

  subroutine binary_setup

    use meth_params_module, only: rot_period, point_mass, URHO, UTEMP, UEINT, UEDEN, UFS, UFX
    use network, only: nspec, naux
    use initial_model_module
    use prob_params_module, only: center, problo, probhi, dim, max_level, dx_level, physbc_lo, Symmetry
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product
    use binary_module, only: get_roche_radii
    use problem_io_module, only: ioproc
    use castro_error_module, only: castro_error
    use ambient_module, only: ambient_state
    use eos_type_module, only: eos_input_rt, eos_t
    use eos_module, only: eos
    use fundamental_constants_module, only: Gconst, c_light, AU, M_solar
    use amrex_constants_module, only: ZERO, THIRD, HALF, ONE, TWO, M_PI
    use probdata_module, only: center_P_initial, center_S_initial, r_P_initial, r_S_initial, a, orbital_eccentricity, orbital_angle, &
                               com_P, com_S, mass_P, mass_S, vel_P, vel_S, roche_rad_P, roche_rad_S, central_density_P, central_density_S, &
                               initial_model_dx, initial_model_npts, initial_model_mass_tol, initial_model_hse_tol, &
                               single_star, axis_1, axis_2, &
                               stellar_temp, &
                               collision_separation, collision_velocity, collision_impact_parameter, &
                               tde_separation, tde_beta, tde_initial_velocity, &
                               tde_tidal_radius, tde_schwarzschild_radius, tde_pericenter_radius, &
                               roche_radius_factor, &
                               init, problem

    implicit none

    real(rt) :: v_ff, collision_offset
    real(rt) :: omega(3)

    integer :: lev

    real(rt) :: v_P_r, v_S_r, v_P_phi, v_S_phi, v_P, v_S

    type(eos_t) :: eos_state

    call get_omega(omega)

    ! Safety check: ensure that if we have a symmetric lower boundary, that the
    ! domain center (and thus the stars) are on that boundary.

    if (physbc_lo(1) .eq. Symmetry .and. center(1) .ne. problo(1)) then
       call castro_error("Symmetric lower x-boundary but the center is not on this boundary.")
    end if

    if (physbc_lo(2) .eq. Symmetry .and. center(2) .ne. problo(2)) then
       call castro_error("Symmetric lower y-boundary but the center is not on this boundary.")
    end if

    if (physbc_lo(3) .eq. Symmetry .and. center(3) .ne. problo(3)) then
       call castro_error("Symmetric lower z-boundary but the center is not on this boundary.")
    end if

    ! Set some default values for these quantities;
    ! we'll update them soon.

    center_P_initial = center
    center_S_initial = center

    com_P = center
    com_S = center

    vel_P = ZERO
    vel_S = ZERO

    ! Allocate arrays to hold the stellar models.

    call initialize_model(.true.,  initial_model_dx, initial_model_npts, initial_model_mass_tol, initial_model_hse_tol)
    call initialize_model(.false., initial_model_dx, initial_model_npts, initial_model_mass_tol, initial_model_hse_tol)

    model_P % min_density = ambient_state(URHO)
    model_S % min_density = ambient_state(URHO)

    model_P % central_temp = stellar_temp
    model_S % central_temp = stellar_temp



    ! Fill in the model's physical details.
    ! If we're integrating to reach a desired mass, set the composition accordingly.
    ! If instead we're fixing the central density, then first we'll assume the composition is
    ! that of a solar mass WD as a initial guess, and get the corresponding mass. 
    ! Then we set the composition to match this preliminary mass, and we'll get a final mass later.

    if (mass_P > ZERO) then

       model_P % mass = mass_P

       call set_wd_composition(model_P)

    elseif (central_density_P > ZERO) then

       model_P % mass = M_solar

       call set_wd_composition(model_P)

       model_P % central_density = central_density_P

       call establish_hse(model_P, rho_P, T_P, xn_P, r_P)

       call set_wd_composition(model_P)

    else

       call castro_error("Must specify either a positive primary mass or a positive primary central density.")

    endif



    if (.not. single_star) then

       if (mass_S > ZERO) then

          model_S % mass = mass_S

          call set_wd_composition(model_S)

       elseif (central_density_S > ZERO) then

          model_S % mass = M_solar

          call set_wd_composition(model_S)

          model_S % central_density = central_density_S

          call establish_hse(model_S, rho_S, T_S, xn_S, r_S)

          call set_wd_composition(model_S)

       else

          call castro_error("If we are doing a binary calculation, we must specify either a " // &
                           "positive secondary mass or a positive secondary central density.")

       endif

       ambient_state(UFS:UFS+nspec-1) = ambient_state(URHO) * (model_P % envelope_comp + model_S % envelope_comp) / 2

    else

       ambient_state(UFS:UFS+nspec-1) = ambient_state(URHO) * model_P % envelope_comp

    endif

    ! We have completed (rho, T, xn) for the ambient state, so we can call the EOS.

    eos_state % rho = ambient_state(URHO)
    eos_state % T   = ambient_state(UTEMP)
    eos_state % xn  = ambient_state(UFS:UFS+nspec-1) / ambient_state(URHO)
    eos_state % aux = ambient_state(UFX:UFX+naux-1) / ambient_state(URHO)

    call eos(eos_input_rt, eos_state)

    ambient_state(UEINT) = ambient_state(URHO) * eos_state % e
    ambient_state(UEDEN) = ambient_state(UEINT)



    roche_rad_P = ZERO
    roche_rad_S = ZERO

    ! Generate primary and secondary WD models.

    call establish_hse(model_P, rho_P, T_P, xn_P, r_P)

    if (ioproc .and. init == 1) then

       ! Set the color to bold green for printing to terminal in this section. See:
       ! http://stackoverflow.com/questions/6402700/coloured-terminal-output-from-fortran

       print *, ''//achar(27)//'[1;32m'

       write (*,1001) model_P % mass / M_solar, model_P % central_density, model_P % radius
       1001 format ("Generated initial model for primary WD of mass ", f4.2, &
                    " solar masses, central density ", ES8.2, " g cm**-3, and radius ", ES8.2, " cm.")
       print *, ""
    endif

    mass_P = model_P % mass
    roche_rad_P = model_P % radius

    if (.not. single_star) then

       call establish_hse(model_S, rho_S, T_S, xn_S, r_S)

       if (ioproc .and. init == 1) then
          write (*,1002) model_S % mass / M_solar, model_S % central_density, model_S % radius
          1002 format ("Generated initial model for secondary WD of mass ", f4.2, &
                       " solar masses, central density ", ES8.2, " g cm**-3, and radius ", ES8.2, " cm.")
          print *, ""
       endif

       mass_S = model_S % mass
       roche_rad_S = model_S % radius

       ! Compute initial Roche radii

       call get_roche_radii(mass_S / mass_P, roche_rad_S, roche_rad_P)

       ! Set up the stellar distances and velocities according to the problem choice

       if (problem == 0) then

          collision_separation = collision_separation * model_S % radius

          if (collision_velocity < 0.0e0_rt) then

             call freefall_velocity(mass_P + mass_S, collision_separation, v_ff)

             vel_P(axis_1) =  (mass_P / (mass_S + mass_P)) * v_ff
             vel_S(axis_1) = -(mass_S / (mass_S + mass_P)) * v_ff

          else

             vel_P(axis_1) =  collision_velocity
             vel_S(axis_1) = -collision_velocity

          end if

          r_P_initial = -(mass_P / (mass_S + mass_P)) * collision_separation
          r_S_initial =  (mass_S / (mass_S + mass_P)) * collision_separation

          a = r_S_initial - r_P_initial

          center_P_initial(axis_1) = center_P_initial(axis_1) + r_P_initial
          center_S_initial(axis_1) = center_S_initial(axis_1) + r_S_initial

          ! We also permit a non-zero impact parameter b in the direction perpendicular
          ! to the motion of the stars. This is measured in units of the radius of the
          ! primary, so that b > 1 doesn't make any sense as the stars won't collide.
          ! Since the secondary's radius is greater than the primary's, measuring in the
          ! units of the primary's radius will guarantee contact.

          collision_offset = collision_impact_parameter * model_P % radius

          center_P_initial(axis_2) = center_P_initial(axis_2) - collision_offset
          center_S_initial(axis_2) = center_S_initial(axis_2) + collision_offset

       else if (problem == 1) then

          if (roche_radius_factor < ZERO) then

             ! Determine the orbital distance based on the rotational period.

             a = -ONE

          else

             ! Set the orbital distance, then calculate the rotational period.

             a = roche_radius_factor * (model_S % radius / roche_rad_S)

             rot_period = -ONE

          endif

          call kepler_third_law(model_P % radius, model_P % mass, model_S % radius, model_S % mass, &
                                rot_period, orbital_eccentricity, orbital_angle, &
                                a, r_P_initial, r_S_initial, v_P_r, v_S_r, v_P_phi, v_S_phi)

          if (ioproc .and. init == 1) then
             write (*,1003) a, a / AU
             write (*,1004) r_P_initial, r_P_initial / AU
             write (*,1005) r_S_initial, r_S_initial / AU
             write (*,1006) rot_period
             write (*,1007) TWO * M_PI * r_P_initial / rot_period
             write (*,1008) TWO * M_PI * r_S_initial / rot_period
1003         format ("Generated binary orbit of distance ", ES8.2, " cm = ", ES8.2, " AU.")
1004         format ("The primary orbits the center of mass at distance ", ES9.2, " cm = ", ES9.2, " AU.")
1005         format ("The secondary orbits the center of mass at distance ", ES9.2, " cm = ", ES9.2, " AU.")
1006         format ("The initial orbital period is ", F6.2 " s.")
1007         format ("The initial orbital speed of the primary is ", ES9.2 " cm/s.")
1008         format ("The initial orbital speed of the secondary is ", ES9.2 " cm/s.")
          endif

          ! Star center positions -- we'll put them in the midplane, with the center of mass at the center of the domain.

          center_P_initial(axis_1) = center_P_initial(axis_1) + r_P_initial * cos(orbital_angle)
          center_P_initial(axis_2) = center_P_initial(axis_2) + r_P_initial * sin(orbital_angle)

          center_S_initial(axis_1) = center_S_initial(axis_1) + r_S_initial * cos(orbital_angle)
          center_S_initial(axis_2) = center_S_initial(axis_2) + r_S_initial * sin(orbital_angle)           

          ! Star velocities, from Kepler's third law. Note that these are the velocities in the inertial frame.

          vel_P(axis_1) = v_P_r   * cos(orbital_angle) - v_P_phi * sin(orbital_angle)
          vel_P(axis_2) = v_P_phi * cos(orbital_angle) + v_P_r   * sin(orbital_angle)

          vel_S(axis_1) = v_S_r   * cos(orbital_angle) - v_S_phi * sin(orbital_angle)
          vel_S(axis_2) = v_S_phi * cos(orbital_angle) + v_S_r   * sin(orbital_angle)

       else

          call castro_error("Error: Unknown problem choice.")

       endif

       ! Scale the Roche radii by the initial distance.

       roche_rad_P = roche_rad_P * a
       roche_rad_S = roche_rad_S * a

    else

       if (problem == 2) then

          ! The tidal radius is given by (M_BH / M_WD)^(1/3) * R_WD.

          tde_tidal_radius = (point_mass / mass_P)**THIRD * model_P % radius

          ! The usual definition for the Schwarzschild radius.

          tde_schwarzschild_radius = TWO * Gconst * point_mass / c_light**TWO

          ! The pericenter radius is the distance of closest approach,
          ! for a point mass on a parabolic orbit.

          tde_pericenter_radius = tde_tidal_radius / tde_beta

          ! Given the pericenter distance, we can calculate the parameters of
          ! the parabolic orbit. A parabolic orbit has E = 0, so KE = PE,
          ! or (1/2) M_WD v**2 = G * M_WD * M_BH / r. This simplifies to
          ! v = sqrt(2 * G * M_BH / r). The initial distance is set at runtime.

          r_P_initial = tde_separation * tde_tidal_radius

          v_P = sqrt(TWO * Gconst * point_mass / r_P_initial)

          ! Now we need to convert this into angular and radial components. To
          ! do this, we need the orbital angle, which comes from the orbit equation,
          ! r = r_0 / (1 - eccentricity * cos(phi)), where for a parabolic orbit,
          ! r_0 is twice the pericenter distance.

          orbital_eccentricity = ONE
          orbital_angle = acos(ONE - (TWO * tde_pericenter_radius) / r_P_initial)

          ! Now set the x and y components of the position and velocity. The position is
          ! straightforward: the orbital angle is the usual angle phi such that x = r cos(phi)
          ! and y = r sin(phi). The velocity is a little more involved and depends on the orbit
          ! equation. The Cartesian form of the parabolic orbit equation is x = y**2 / (2 * r_0) + r_0 / 2,
          ! so dx/dt = (y / r_0) * dy/dt. Given that v**2 = v_x**2 + v_y**2, we have
          ! v_x = v / sqrt( 1 + (r_0 / (r * sin(phi))) )**2 , and
          ! v_y = v / sqrt( 1 + (r * sin(phi) / r_0)**2 ).

          center_P_initial(axis_1) = center_P_initial(axis_1) - r_P_initial * cos(orbital_angle)
          center_P_initial(axis_2) = center_P_initial(axis_2) - r_P_initial * sin(orbital_angle)

          if (tde_initial_velocity == 1) then

             vel_P(axis_1) = v_P / sqrt(ONE + (TWO * tde_pericenter_radius / (r_P_initial * sin(orbital_angle)))**2)
             vel_P(axis_2) = v_P / sqrt(ONE + (r_P_initial * sin(orbital_angle) / (TWO * tde_pericenter_radius))**2)

          end if

       end if

    endif

    ! Reset the terminal color to its previous state.

    if (ioproc .and. init == 1) then
       print *, ''//achar(27)//'[0m'
    endif

    com_P = center_P_initial
    com_S = center_S_initial

    ! Safety check: make sure the stars are actually inside the computational domain.

    if (.not. (dim .eq. 2 .and. physbc_lo(2) .eq. Symmetry)) then

       if ( (HALF * (probhi(1) - problo(1)) < model_P % radius) .or. &
            (HALF * (probhi(2) - problo(2)) < model_P % radius) .or. &
            (HALF * (probhi(3) - problo(3)) < model_P % radius .and. dim .eq. 3) ) then
          call castro_error("Primary does not fit inside the domain.")
       endif

       if ( (HALF * (probhi(1) - problo(1)) < model_S % radius) .or. &
            (HALF * (probhi(2) - problo(2)) < model_S % radius) .or. &
            (HALF * (probhi(3) - problo(3)) < model_S % radius .and. dim .eq. 3) ) then
          call castro_error("Secondary does not fit inside the domain.")
       endif

    else

       if ( (probhi(1) - problo(1) < model_S % radius) .or. &
            (probhi(2) - problo(2) < 2 * model_S % radius) ) then
          call castro_error("Secondary does not fit inside the domain.")
       end if

    end if

  end subroutine binary_setup



  ! Accepts the masses of two stars (in solar masses)
  ! and the orbital period of a system,
  ! and returns the semimajor axis of the orbit (in cm),
  ! as well as the distances a_1 and a_2 from the center of mass.

  subroutine kepler_third_law(radius_1, mass_1, radius_2, mass_2, period, eccentricity, phi, a, r_1, r_2, v_1r, v_2r, v_1p, v_2p)

    use amrex_constants_module, only: ZERO, THIRD, HALF, ONE, TWO, M_PI, FOUR
    use prob_params_module, only: problo, probhi, physbc_lo, Symmetry
    use sponge_module, only: sponge_lower_radius
    use meth_params_module, only: do_sponge
    use fundamental_constants_module, only: Gconst
    use castro_error_module, only: castro_error
    use probdata_module, only: axis_1

    implicit none

    real(rt), intent(in   ) :: mass_1, mass_2, eccentricity, phi, radius_1, radius_2
    real(rt), intent(inout) :: period, a, r_1, r_2, v_1r, v_2r, v_1p, v_2p

    real(rt) :: length

    real(rt) :: mu, M ! Reduced mass, total mass
    real(rt) :: r     ! Position
    real(rt) :: v_r, v_phi ! Radial and azimuthal velocity

    ! Definitions of total and reduced mass

    M  = mass_1 + mass_2
    mu = mass_1 * mass_2 / M

    ! First, solve for the orbit in the reduced one-body problem, where
    ! an object of mass mu orbits an object with mass M located at r = 0.
    ! For this we follow Carroll and Ostlie, Chapter 2, but many texts discuss this.
    ! Note that we use the convention that phi measures angle from aphelion,
    ! which is opposite to the convention they use.

    if (period > ZERO .and. a < ZERO) then

       a = (Gconst * M * period**2 / (FOUR * M_PI**2))**THIRD ! C + O, Equation 2.37

    else if (period < ZERO .and. a > ZERO) then

       period = (a**3 * FOUR * M_PI**2 / (Gconst * M))**HALF

    else

       call castro_error("Error: overspecified Kepler's third law calculation.")

    endif

    r = a * (ONE - eccentricity**2) / (ONE - eccentricity * cos(phi)) ! C + O, Equation 2.3

    ! To get the radial and azimuthal velocity, we take the appropriate derivatives of the above.
    ! v_r = dr / dt = dr / d(phi) * d(phi) / dt, with d(phi) / dt being derived from
    ! C + O, Equation 2.30 for the angular momentum, and the fact that L = mu * r**2 * d(phi) / dt.

    v_r   = -TWO * M_PI * a * eccentricity * sin(phi) / (period * (ONE - eccentricity**2)**HALF)
    v_phi =  TWO * M_PI * a * (ONE - eccentricity * cos(phi)) / (period * (ONE - eccentricity**2)**HALF)

    ! Now convert everything back to the binary frame, using C+O, Equation 2.23 and 2.24. This applies
    ! to the velocities as well as the positions because the factor in front of r_1 and r_2 is constant.

    r_1  = -(mu / mass_1) * r
    r_2  =  (mu / mass_2) * r

    v_1r = -(mu / mass_1) * v_r
    v_2r =  (mu / mass_2) * v_r

    v_1p = -(mu / mass_1) * v_phi
    v_2p =  (mu / mass_2) * v_phi

    ! Make sure the domain is big enough to hold stars in an orbit this size.

    if (physbc_lo(axis_1) .eq. Symmetry) then

       ! In this case we're only modelling the secondary.
       length = r_2 + radius_2

    else

       length = (r_2 - r_1) + radius_1 + radius_2

    end if

    if (length > (probhi(axis_1)-problo(axis_1))) then
       call castro_error("ERROR: The domain width is too small to include the binary orbit.")
    endif

    ! We want to do a similar check to make sure that no part of the stars
    ! land in the sponge region.

    if (do_sponge .eq. 1 .and. sponge_lower_radius > ZERO) then

       if (abs(r_1) + radius_1 .ge. sponge_lower_radius) then
          call castro_error("ERROR: Primary contains material inside the sponge region.")
       endif

       if (abs(r_2) + radius_2 .ge. sponge_lower_radius) then
          call castro_error("ERROR: Secondary contains material inside the sponge region.")
       endif

    endif

    ! Make sure the stars are not touching.
    if (radius_1 + radius_2 > a) then
       call castro_error("ERROR: Stars are touching!")
    endif

  end subroutine kepler_third_law



  ! Given total mass of a binary system and the initial separation of
  ! two point particles, obtain the velocity at this separation 
  ! assuming the point masses fell in from infinity. This will
  ! be the velocity in the frame where the center of mass is stationary.

  subroutine freefall_velocity(mass, distance, vel)

    use amrex_constants_module, only: HALF, TWO
    use fundamental_constants_module, only: Gconst

    implicit none

    real(rt), intent(in   ) :: mass, distance
    real(rt), intent(inout) :: vel

    vel = (TWO * Gconst * mass / distance)**HALF

  end subroutine freefall_velocity



  ! Given a rotating frame velocity, get the inertial frame velocity.
  ! Note that we simply return the original velocity if we're
  ! already in the inertial frame.

  function inertial_velocity(loc, vel, time) result(vel_i)

    use meth_params_module, only: do_rotation, state_in_rotating_frame
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product ! function

    implicit none

    real(rt), intent(in   ) :: loc(3), vel(3), time
    real(rt) :: omega(3)

    real(rt) :: vel_i(3)

    !$gpu

    call get_omega(omega)

    vel_i = vel

    if (do_rotation .eq. 1 .and. state_in_rotating_frame .eq. 1) then
       vel_i = vel_i + cross_product(omega, loc)
    endif

  end function inertial_velocity



  ! C++ interface for inertial_velocity.

  subroutine get_inertial_velocity(loc, vel, time, inertial_vel) bind(C,name='get_inertial_velocity')

    implicit none

    real(rt), intent(in   ) :: loc(3), vel(3), time
    real(rt), intent(inout) :: inertial_vel(3)

    inertial_vel = inertial_velocity(loc, vel, time)

  end subroutine get_inertial_velocity



  subroutine set_relaxation_damping_factor(factor) bind(C,name='set_relaxation_damping_factor')

    use probdata_module, only: relaxation_damping_factor

    implicit none

    real(rt), intent(in), value :: factor

    relaxation_damping_factor = factor

  end subroutine set_relaxation_damping_factor



  subroutine get_axes(axis_1_in, axis_2_in, axis_3_in) bind(C,name='get_axes')

    use probdata_module, only: axis_1, axis_2, axis_3

    implicit none

    integer, intent(inout) :: axis_1_in, axis_2_in, axis_3_in

    axis_1_in = axis_1
    axis_2_in = axis_2
    axis_3_in = axis_3

  end subroutine get_axes



  ! Return whether we're doing a single star simulation or not.

  subroutine get_single_star(flag) bind(C,name='get_single_star')

    use probdata_module, only: single_star

    implicit none

    integer, intent(inout) :: flag

    flag = 0

    if (single_star) flag = 1

  end subroutine get_single_star



  ! Return the problem type.

  subroutine get_problem_number(problem_out) bind(C,name='get_problem_number')

    use probdata_module, only: problem

    implicit none

    integer :: problem_out

    problem_out = problem

  end subroutine get_problem_number



  ! Return the mass-weighted center of mass and velocity
  ! for the primary and secondary, for a given FAB.
  ! Since this will rely on a sum over processors,
  ! we should only add to the relevant variables
  ! in anticipation of a MPI reduction, and not
  ! overwrite them. Note that ultimately what we
  ! are doing here is to use an old guess at the
  ! effective potential of the primary and secondary
  ! to generate a new estimate.

  subroutine wdcom(rho,  r_lo, r_hi, &
                   xmom, px_lo, px_hi, &
                   ymom, py_lo, py_hi, &
                   zmom, pz_lo, pz_hi, &
                   pmask, pm_lo, pm_hi, &
                   smask, sm_lo, sm_hi, &
                   vol,  vo_lo, vo_hi, &
                   lo, hi, dx, time, &
                   com_p_x, com_p_y, com_p_z, &
                   com_s_x, com_s_y, com_s_z, &
                   vel_p_x, vel_p_y, vel_p_z, &
                   vel_s_x, vel_s_y, vel_s_z, &
                   m_p, m_s) bind(C,name='wdcom')

    use amrex_constants_module, only: HALF, ZERO, ONE, TWO
    use prob_params_module, only: problo, probhi, physbc_lo, physbc_hi, Symmetry, coord_type
    use castro_util_module, only: position ! function
    use reduction_module, only: reduce_add

    implicit none

    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: px_lo(3), px_hi(3)
    integer,  intent(in   ) :: py_lo(3), py_hi(3)
    integer,  intent(in   ) :: pz_lo(3), pz_hi(3)
    integer,  intent(in   ) :: pm_lo(3), pm_hi(3)
    integer,  intent(in   ) :: sm_lo(3), sm_hi(3)
    integer,  intent(in   ) :: vo_lo(3), vo_hi(3)

    real(rt), intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: xmom(px_lo(1):px_hi(1),px_lo(2):px_hi(2),px_lo(3):px_hi(3))
    real(rt), intent(in   ) :: ymom(py_lo(1):py_hi(1),py_lo(2):py_hi(2),py_lo(3):py_hi(3))
    real(rt), intent(in   ) :: zmom(pz_lo(1):pz_hi(1),pz_lo(2):pz_hi(2),pz_lo(3):pz_hi(3))
    real(rt), intent(in   ) :: pmask(pm_lo(1):pm_hi(1),pm_lo(2):pm_hi(2),pm_lo(3):pm_hi(3))
    real(rt), intent(in   ) :: smask(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
    real(rt), intent(in   ) :: vol(vo_lo(1):vo_hi(1),vo_lo(2):vo_hi(2),vo_lo(3):vo_hi(3))

    integer,  intent(in   ) :: lo(3), hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: com_p_x, com_p_y, com_p_z
    real(rt), intent(inout) :: com_s_x, com_s_y, com_s_z
    real(rt), intent(inout) :: vel_p_x, vel_p_y, vel_p_z
    real(rt), intent(inout) :: vel_s_x, vel_s_y, vel_s_z
    real(rt), intent(inout) :: m_p, m_s
    real(rt), intent(in   ), value :: time

    integer  :: i, j, k
    real(rt) :: r(3), rSymmetric(3), dm, dmSymmetric, momSymmetric(3)
    real(rt) :: primary_factor, secondary_factor

    !$gpu

    ! Add to the COM locations and velocities of the primary and secondary
    ! depending on which potential dominates, ignoring unbound material.
    ! Note that in this routine we actually are
    ! summing mass-weighted quantities for the COM and the velocity; 
    ! we will account for this at the end of the calculation in 
    ! post_timestep() by dividing by the mass.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! Our convention is that the COM locations for the WDs are 
             ! absolute positions on the grid, not relative to the center.

             r = position(i,j,k)

             ! We account for symmetric boundaries in this sum as usual,
             ! by adding to the position the locations that would exist
             ! on the opposite side of the symmetric boundary. Note that
             ! in axisymmetric coordinates, some of this work is already
             ! done for us in the definition of the zone volume.

             rSymmetric = r
             rSymmetric = merge(rSymmetric + (problo - rSymmetric), rSymmetric, physbc_lo(:) .eq. Symmetry)
             rSymmetric = merge(rSymmetric + (rSymmetric - probhi), rSymmetric, physbc_hi(:) .eq. Symmetry)

             dm = rho(i,j,k) * vol(i,j,k)

             dmSymmetric = dm
             momSymmetric(1) = xmom(i,j,k)
             momSymmetric(2) = ymom(i,j,k)
             momSymmetric(3) = zmom(i,j,k)

             if (coord_type .eq. 0) then

                if (physbc_lo(1) .eq. Symmetry) then
                   dmSymmetric = TWO * dmSymmetric
                   momSymmetric = TWO * momSymmetric
                end if

                if (physbc_lo(2) .eq. Symmetry) then
                   dmSymmetric = TWO * dmSymmetric
                   momSymmetric = TWO * momSymmetric
                end if

                if (physbc_lo(3) .eq. Symmetry) then
                   dmSymmetric = TWO * dmSymmetric
                   momSymmetric = TWO * momSymmetric
                end if

             end if

             primary_factor = ZERO
             secondary_factor = ZERO

             if (pmask(i,j,k) > ZERO) then

                primary_factor = ONE

             else if (smask(i,j,k) > ZERO) then

                secondary_factor = ONE

             endif

             call reduce_add(m_p, dmSymmetric * primary_factor)

             call reduce_add(com_p_x, dmSymmetric * rSymmetric(1) * primary_factor)
             call reduce_add(com_p_y, dmSymmetric * rSymmetric(2) * primary_factor)
             call reduce_add(com_p_z, dmSymmetric * rSymmetric(3) * primary_factor)

             call reduce_add(vel_p_x, momSymmetric(1) * vol(i,j,k) * primary_factor)
             call reduce_add(vel_p_y, momSymmetric(2) * vol(i,j,k) * primary_factor)
             call reduce_add(vel_p_z, momSymmetric(3) * vol(i,j,k) * primary_factor)

             call reduce_add(m_s, dmSymmetric * secondary_factor)

             call reduce_add(com_s_x, dmSymmetric * rSymmetric(1) * secondary_factor)
             call reduce_add(com_s_y, dmSymmetric * rSymmetric(2) * secondary_factor)
             call reduce_add(com_s_z, dmSymmetric * rSymmetric(3) * secondary_factor)

             call reduce_add(vel_s_x, momSymmetric(1) * vol(i,j,k) * secondary_factor)
             call reduce_add(vel_s_y, momSymmetric(2) * vol(i,j,k) * secondary_factor)
             call reduce_add(vel_s_z, momSymmetric(3) * vol(i,j,k) * secondary_factor)

          enddo
       enddo
    enddo

  end subroutine wdcom



  ! This function uses the known center of mass of the two white dwarfs,
  ! and given a density cutoff, computes the total volume of all zones
  ! whose density is greater or equal to that density cutoff.
  ! We also impose a distance requirement so that we only look
  ! at zones within the Roche lobe of the white dwarf.

  subroutine ca_volumeindensityboundary(rho, r_lo, r_hi, &
                                        pmask, pm_lo, pm_hi, &
                                        smask, sm_lo, sm_hi, &
                                        vol, v_lo, v_hi, &
                                        lo, hi, dx, &
                                        volp, vols, rho_cutoff) &
                                        bind(C, name='ca_volumeindensityboundary')

    use amrex_constants_module, only: ZERO, ONE
    use reduction_module, only: reduce_add

    implicit none

    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: pm_lo(3), pm_hi(3)
    integer,  intent(in   ) :: sm_lo(3), sm_hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    integer,  intent(in   ) :: lo(3), hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: pmask(pm_lo(1):pm_hi(1),pm_lo(2):pm_hi(2),pm_lo(3):pm_hi(3))
    real(rt), intent(in   ) :: smask(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
    real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(inout) :: volp, vols
    real(rt), intent(in   ), value :: rho_cutoff

    integer :: i, j, k
    real(rt) :: primary_factor, secondary_factor

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             primary_factor = ZERO
             secondary_factor = ZERO

             if (rho(i,j,k) > rho_cutoff) then

                if (pmask(i,j,k) > ZERO) then

                   primary_factor = ONE

                else if (smask(i,j,k) > ZERO) then

                   secondary_factor = ONE

                endif

             endif

             call reduce_add(volp, vol(i,j,k) * primary_factor)
             call reduce_add(vols, vol(i,j,k) * secondary_factor)

          enddo
       enddo
    enddo

  end subroutine ca_volumeindensityboundary



  ! Given state data in the rotating frame, transform it to the inertial frame.

  subroutine transform_to_inertial_frame(state, s_lo, s_hi, lo, hi, time) &
                                         bind(C,name='transform_to_inertial_frame')

    use meth_params_module, only: NVAR, URHO, UMX, UMZ
    use castro_util_module, only: position

    implicit none

    integer  :: lo(3), hi(3)
    integer  :: s_lo(3), s_hi(3)
    real(rt) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt) :: time

    real(rt) :: loc(3), vel(3)
    integer  :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k)
             vel = state(i,j,k,UMX:UMZ) / state(i,j,k,URHO)

             state(i,j,k,UMX:UMZ) = state(i,j,k,URHO) * inertial_velocity(loc, vel, time)

          enddo
       enddo
    enddo

  end subroutine transform_to_inertial_frame



  ! Given the above quadrupole tensor, calculate the strain tensor.

  subroutine gw_strain_tensor(h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3, Qtt, time) &
                              bind(C,name='gw_strain_tensor')

    use amrex_constants_module, only: ZERO, HALF, ONE, TWO
    use fundamental_constants_module, only: Gconst, c_light, parsec
    use prob_params_module, only: dim
    use probdata_module, only: axis_1, axis_2, axis_3, gw_dist

    implicit none

    real(rt), intent(inout) :: h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3
    real(rt), intent(in   ) :: Qtt(3,3)
    real(rt), intent(in   ) :: time

    integer  :: i, j, k, l, dir
    real(rt) :: h(3,3), proj(3,3,3,3), delta(3,3), n(3), r
    real(rt) :: dist(3)

    ! Standard Kronecker delta.

    delta(:,:) = ZERO

    do i = 1, 3
       delta(i,i) = ONE
    enddo

    ! Unit vector for the wave is simply the distance
    ! vector to the observer normalized by the total distance.
    ! We are going to repeat this process by looking along
    ! all three coordinate axes.

    do dir = 1, 3

       dist(:) = ZERO
       dist(dir) = gw_dist

       r = sqrt(sum(dist**2))

       n(:) = dist(:) / r

       h = ZERO

       ! Projection operator onto the unit vector n.

       do l = 1, 3
          do k = 1, 3
             do j = 1, 3
                do i = 1, 3
                   proj(i,j,k,l) = (delta(i,k) - n(i) * n(k)) * (delta(j,l) - n(j) * n(l)) &
                                 - HALF * (delta(i,j) - n(i) * n(j)) * (delta(k,l) - n(k) * n(l))
                enddo
             enddo
          enddo
       enddo

       ! Now we can calculate the strain tensor.

       do l = 1, 3
          do k = 1, 3
             do j = 1, 3
                do i = 1, 3
                   h(i,j) = h(i,j) + proj(i,j,k,l) * Qtt(k,l)
                enddo
             enddo
          enddo
       enddo

       ! Finally multiply by the coefficients.

       r = r * parsec * 1d3 ! Convert from kpc to cm.

       h(:,:) = h(:,:) * TWO * Gconst / (c_light**4 * r)

       if (dim .eq. 3) then

          ! If rot_axis == 3, then h_+ = h_{11} = -h_{22} and h_x = h_{12} = h_{21}.
          ! Analogous statements hold along the other axes.

          ! We are adding here so that this calculation makes sense on multiple levels.

          if (dir .eq. axis_1) then

             h_plus_1  = h_plus_1  + h(axis_2,axis_2)
             h_cross_1 = h_cross_1 + h(axis_2,axis_3)

          else if (dir .eq. axis_2) then

             h_plus_2  = h_plus_2  + h(axis_3,axis_3)
             h_cross_2 = h_cross_2 + h(axis_3,axis_1)

          else if (dir .eq. axis_3) then

             h_plus_3  = h_plus_3  + h(axis_1,axis_1)
             h_cross_3 = h_cross_3 + h(axis_1,axis_2)

          endif

       else

          ! In 2D axisymmetric coordinates, enforce that axis_1 is the x-axis,
          ! axis_2 is the y-axis, and axis_3 is the z-axis.

          if (dir .eq. 1) then

             h_plus_1  = h_plus_1  + h(2,2)
             h_cross_1 = h_cross_1 + h(2,3)

          else if (dir .eq. 2) then

             h_plus_2  = h_plus_2  + h(3,3)
             h_cross_2 = h_cross_2 + h(3,1)

          else if (dir .eq. 3) then

             h_plus_3  = h_plus_3  + h(1,1)
             h_cross_3 = h_cross_3 + h(1,2)

          endif

       endif

    enddo

  end subroutine gw_strain_tensor



  subroutine update_center(time) bind(C,name='update_center')

    use amrex_constants_module, only: ZERO
    use castro_error_module, only: castro_error
    use prob_params_module, only: center, problo, probhi, dim
    use probdata_module, only: center_fracx, center_fracy, center_fracz, &
                               bulk_velx, bulk_vely, bulk_velz

    implicit none

    real(rt), intent(in) :: time

    ! Determine the original location of the center.

    if (dim .eq. 3) then

       center(1) = problo(1) + center_fracx * (probhi(1) - problo(1))
       center(2) = problo(2) + center_fracy * (probhi(2) - problo(2))
       center(3) = problo(3) + center_fracz * (probhi(3) - problo(3))

    else if (dim .eq. 2) then

       center(1) = problo(1)
       center(2) = problo(2) + center_fracz * (probhi(2) - problo(2))
       center(3) = ZERO

    else if (dim .eq. 1) then

       center(1) = problo(1) + center_fracx * (probhi(1) - problo(1))
       center(2) = ZERO
       center(3) = ZERO

    else

       call castro_error("Error: unknown dim in subroutine update_center.")

    endif

    ! Now update using the time passed since the beginning of the simulation.

    center(1) = center(1) + bulk_velx * time
    center(2) = center(2) + bulk_vely * time
    center(3) = center(3) + bulk_velz * time

  end subroutine update_center



  ! Updates the CASTRO rotational period.

  subroutine set_period(period) bind(C,name='set_period')

    use meth_params_module, only: rot_period

    implicit none

    real(rt) :: period

    rot_period = period

  end subroutine set_period



  ! Returns the CASTRO rotational period.

  subroutine get_period(period) bind(C,name='get_period')

    use meth_params_module, only: rot_period

    implicit none

    real(rt) :: period

    period = rot_period

  end subroutine get_period



  ! Returns the CASTRO rotation frequency vector.

  subroutine get_omega_vec(omega_in) bind(C,name='get_omega_vec')

    use rotation_frequency_module, only: get_omega

    implicit none

    real(rt), intent(inout) :: omega_in(3)

    call get_omega(omega_in)

  end subroutine get_omega_vec



  ! Updates the global extrema.

  subroutine set_extrema(T_max, rho_max, ts_te_max) bind(C,name='set_extrema')

    use probdata_module, only: T_global_max, rho_global_max, ts_te_global_max

    implicit none

    real(rt), intent(in) :: T_max, rho_max, ts_te_max

    T_global_max     = T_max
    rho_global_max   = rho_max
    ts_te_global_max = ts_te_max

  end subroutine set_extrema



  ! Retrieves the global extrema.

  subroutine get_extrema(T_max, rho_max, ts_te_max) bind(C,name='get_extrema')

    use probdata_module, only: T_global_max, rho_global_max, ts_te_global_max

    implicit none

    real(rt), intent(inout) :: T_max, rho_max, ts_te_max

    T_max     = T_global_max
    rho_max   = rho_global_max
    ts_te_max = ts_te_global_max

  end subroutine get_extrema



  ! Returns whether the simulation is done.

  subroutine get_job_status(jobDoneStatus) bind(C,name='get_job_status')

    use probdata_module, only: jobIsDone

    implicit none

    integer, intent(inout) :: jobDoneStatus

    if (jobIsDone) then
       jobDoneStatus = 1
    else
       jobDoneStatus = 0
    endif

  end subroutine get_job_status



  ! Sets whether the simulation is done.

  subroutine set_job_status(jobDoneStatus) bind(C,name='set_job_status')

    use probdata_module, only: jobIsDone

    implicit none

    integer, intent(in) :: jobDoneStatus

    if (jobDoneStatus == 1) then
       jobIsDone = .true.
    else
       jobIsDone = .false.
    endif

  end subroutine set_job_status



  ! Get the relaxation_cutoff_density parameter.

  subroutine get_relaxation_density_cutoff(relaxation_density_cutoff_in) bind(C, name='get_relaxation_density_cutoff')

    use probdata_module, only: relaxation_density_cutoff
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(inout) :: relaxation_density_cutoff_in

    relaxation_density_cutoff_in = relaxation_density_cutoff

  end subroutine get_relaxation_density_cutoff




  ! Get the relaxation_cutoff_time parameter.

  subroutine get_relaxation_cutoff_time(relaxation_cutoff_time_in) bind(C,name='get_relaxation_cutoff_time')

    use probdata_module, only: relaxation_cutoff_time

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(inout) :: relaxation_cutoff_time_in

    relaxation_cutoff_time_in = relaxation_cutoff_time

  end subroutine get_relaxation_cutoff_time



  ! Gets whether the relaxation is done.

  subroutine get_relaxation_status(relaxation_status) bind(C,name='get_relaxation_status')

    use probdata_module, only: relaxation_is_done

    implicit none

    integer, intent(inout) :: relaxation_status

    relaxation_status = relaxation_is_done

  end subroutine get_relaxation_status



  ! Sets whether the relaxation is done.

  subroutine set_relaxation_status(relaxation_status) bind(C,name='set_relaxation_status')

    use probdata_module, only: relaxation_is_done

    implicit none

    integer, intent(in) :: relaxation_status

    relaxation_is_done = relaxation_status

  end subroutine set_relaxation_status



  ! Retrieve the total energy array.

  subroutine get_total_ener_array(ener_array_in) bind(C,name='get_total_ener_array')

    use probdata_module, only: num_previous_ener_timesteps, total_ener_array

    implicit none

    real(rt), intent(inout) :: ener_array_in(num_previous_ener_timesteps)

    ener_array_in(:) = total_ener_array(:)

  end subroutine get_total_ener_array



  ! Set the total energy array.

  subroutine set_total_ener_array(ener_array_in) bind(C,name='set_total_ener_array')

    use probdata_module, only: num_previous_ener_timesteps, total_ener_array

    implicit none

    real(rt), intent(in) :: ener_array_in(num_previous_ener_timesteps)

    total_ener_array(:) = ener_array_in(:)

  end subroutine set_total_ener_array

end module wdmerger_util_module
