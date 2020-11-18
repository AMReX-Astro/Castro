module wdmerger_util_module

  use iso_c_binding
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  interface
     subroutine kepler_third_law (radius_1, mass_1, radius_2, mass_2, period, eccentricity, &
                                  phi, a, r_1, r_2, v_1r, v_2r, v_1p, v_2p) bind(C)
       use amrex_fort_module, only: rt => amrex_real
       implicit none
       real(rt), intent(in   ), value :: mass_1, mass_2, eccentricity, phi, radius_1, radius_2
       real(rt), intent(inout) :: period, a, r_1, r_2, v_1r, v_2r, v_1p, v_2p
     end subroutine kepler_third_law

     subroutine freefall_velocity(mass, distance, vel) bind(C)
       use amrex_fort_module, only: rt => amrex_real
       implicit none
       real(rt), intent(in   ), value :: mass, distance
       real(rt), intent(inout) :: vel
     end subroutine freefall_velocity
  end interface

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

    use meth_params_module, only: rotational_period, point_mass, URHO, UTEMP, UEINT, UEDEN, UFS, UFX, do_sponge
    use sponge_module, only: sponge_lower_radius
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
    real(rt) :: length

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

             rotational_period = -ONE

          endif

          call kepler_third_law(model_P % radius, model_P % mass, model_S % radius, model_S % mass, &
                                rotational_period, orbital_eccentricity, orbital_angle, &
                                a, r_P_initial, r_S_initial, v_P_r, v_S_r, v_P_phi, v_S_phi)

          ! Make sure the domain is big enough to hold stars in an orbit this size.

          if (physbc_lo(axis_1) .eq. Symmetry) then

             ! In this case we're only modelling the secondary.
             length = r_P_initial + model_P % radius

          else

             length = (r_S_initial - r_P_initial) + model_P % radius + model_S % radius

          end if

          if (length > (probhi(axis_1) - problo(axis_1))) then
             call castro_error("ERROR: The domain width is too small to include the binary orbit.")
          endif

          ! We want to do a similar check to make sure that no part of the stars
          ! land in the sponge region.

          if (do_sponge .eq. 1 .and. sponge_lower_radius > ZERO) then

             if (abs(r_P_initial) + model_P % radius .ge. sponge_lower_radius) then
                call castro_error("ERROR: Primary contains material inside the sponge region.")
             endif

             if (abs(r_S_initial) + model_S % radius .ge. sponge_lower_radius) then
                call castro_error("ERROR: Secondary contains material inside the sponge region.")
             endif

          end if

          ! Make sure the stars are not touching.

          if (model_P % radius + model_S % radius > a) then
             call castro_error("ERROR: Stars are touching!")
          endif

          if (ioproc .and. init == 1) then
             write (*,1003) a, a / AU
             write (*,1004) r_P_initial, r_P_initial / AU
             write (*,1005) r_S_initial, r_S_initial / AU
             write (*,1006) rotational_period
             write (*,1007) TWO * M_PI * r_P_initial / rotational_period
             write (*,1008) TWO * M_PI * r_S_initial / rotational_period
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

    use meth_params_module, only: rotational_period

    implicit none

    real(rt) :: period

    rotational_period = period

  end subroutine set_period



  ! Returns the CASTRO rotational period.

  subroutine get_period(period) bind(C,name='get_period')

    use meth_params_module, only: rotational_period

    implicit none

    real(rt) :: period

    period = rotational_period

  end subroutine get_period



  ! Returns the CASTRO rotational axis.

  subroutine get_rot_axis(rot_axis_out) bind(C,name='get_rot_axis')

    use meth_params_module, only: rot_axis

    implicit none

    integer, intent(inout) :: rot_axis_out

    rot_axis_out = rot_axis

  end subroutine get_rot_axis



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

end module wdmerger_util_module
