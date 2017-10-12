module wdmerger_util_module

  use probdata_module

  implicit none

contains

  ! This routine calls all of the other subroutines at the beginning
  ! of a simulation to fill in the basic problem data.

  subroutine initialize_problem(init_in)

    use bl_error_module, only: bl_error
    use prob_params_module, only: dim

    implicit none

    integer :: init_in

    init = init_in

    ! Read in the namelist to set problem parameters.

    call read_namelist

    ! Establish binary parameters and create initial models.

    call binary_setup

    ! Set small_pres and small_ener.

    call set_small

  end subroutine initialize_problem



  ! This routine reads in the namelist

  subroutine read_namelist

    use bl_constants_module, only: ZERO
    use meth_params_module
    use prob_params_module, only: dim, coord_type
    use problem_io_module, only: probin

    implicit none

    integer :: untin

    ! Read namelist to override the module defaults.

    untin = 9 
    open(untin,file=probin,form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

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

    call ensure_primary_mass_larger

    ! We enforce that the orbital plane is (z, phi) for two-dimensional problems,
    ! with rotation about z = 0 along the radial axis.

    if (dim .eq. 2) then

       if (coord_type .ne. 1) then
          call bl_error("We only support cylindrical coordinates in two dimensions. Set coord_type == 1.")
       endif

       axis_1 = 2
       axis_2 = 3
       rot_axis = 1

    endif

    ! Make sure we have a sensible collision impact parameter.

    if (collision_impact_parameter > 1.0) then
       call bl_error("Impact parameter must be less than one in our specified units.")
    endif

    ! Safety check: we can't run most problems in one dimension.

    if (dim .eq. 1 .and. (.not. (problem .eq. 0 .or. problem .eq. 4))) then
       call bl_error("Can only run a collision or freefall in 1D. Exiting.")
    endif

    ! Don't do a collision or a free-fall in a rotating reference frame.

    if (problem .eq. 0 .and. do_rotation .eq. 1) then
       call bl_error("The collision problem does not make sense in a rotating reference frame.")
    endif

    if (problem .eq. 4 .and. do_rotation .eq. 1) then
       call bl_error("The free-fall problem does not make sense in a rotating reference frame.")
    endif

    ! Make sure we have a sensible eccentricity.

    if (orbital_eccentricity >= 1.0) then
       call bl_error("Orbital eccentricity cannot be larger than one.")
    endif

    ! Make sure we have a sensible angle. Then convert it to radians.

    if (orbital_angle < 0.0 .or. orbital_angle > 360.0) then
       call bl_error("Orbital angle must be between 0 and 360 degrees.")
    endif

    orbital_angle = orbital_angle * M_PI / 180.0

    ! If we're doing problem 3, which has an initial relaxation, ensure
    ! that the damping timescale is positive. Be aware that this safety
    ! check implies that relaxation_damping_timescale must be positive even
    ! if we're restarting from a checkpoint that indicates that the relaxation
    ! has already completed. There's no way for us to test that here because
    ! the checkpoint is read after the namelist initialization step we are
    ! currently in. If the checkpoint does indicate that, the relaxation
    ! damping timescale will just be set negative, effectively disabling
    ! relaxation, at the time the checkpoint is read.

    if (problem .eq. 3 .and. relaxation_damping_timescale <= ZERO) then
       call bl_error("The relaxation step requires a positive relaxation_damping_timescale.")
    endif

    ! Disable the Coriolis term if we're doing a relaxation.

    if (problem .eq. 3) then
       rotation_include_coriolis = 0
    endif

  end subroutine read_namelist



  ! Calculate small_pres and small_ener

  subroutine set_small

    use meth_params_module, only: small_temp, small_pres, small_dens, small_ener

    implicit none

    type (eos_t) :: eos_state

    ! Given the inputs of small_dens and small_temp, figure out small_pres.

    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = ambient_comp

    call eos(eos_input_rt, eos_state)

    small_pres = eos_state % p
    small_ener = eos_state % e

  end subroutine set_small



  ! Returns the ambient state

  subroutine get_ambient(ambient_state)

    implicit none

    type (eos_t) :: ambient_state

    ! Define ambient state, using a composition that is an 
    ! even mixture of the primary and secondary composition, 
    ! and then call the EOS to get internal energy and pressure.

    ambient_state % rho = ambient_density
    ambient_state % T   = ambient_temp
    ambient_state % xn  = ambient_comp

    call eos(eos_input_rt, ambient_state)

  end subroutine get_ambient



  ! Given a WD mass, set its core and envelope composition.

  subroutine set_wd_composition(model)

    use extern_probin_module, only: small_x
    use network, only: network_species_index

    implicit none

    type (initial_model), intent(inout) :: model

    integer :: iHe4, iC12, iO16, iNe20, iMg24

    iHe4 = network_species_index("helium-4")
    iC12 = network_species_index("carbon-12")
    iO16 = network_species_index("oxygen-16")
    iNe20 = network_species_index("neon-20")
    iMg24 = network_species_index("magnesium-24")

    if (iHe4 < 0) call bl_error("Must have He4 in the nuclear network.")
    if (iC12 < 0) call bl_error("Must have C12 in the nuclear network.")
    if (iO16 < 0) call bl_error("Must have O16 in the nuclear network.")
    if (iNe20 < 0) call bl_error("Must have Ne20 in the nuclear network.")
    if (iMg24 < 0) call bl_error("Must have Mg24 in the nuclear network.")

    model % core_comp = small_x
    model % envelope_comp = small_x

    model % envelope_mass = ZERO

    ! Here we follow the prescription of Dan et al. 2012.

    if (model % mass > ZERO .and. model % mass < max_he_wd_mass) then

       model % core_comp(iHe4) = ONE

       model % envelope_comp = model % core_comp

    else if (model % mass >= max_he_wd_mass .and. model % mass < max_hybrid_wd_mass) then

       model % core_comp(iC12) = hybrid_wd_c_frac
       model % core_comp(iO16) = hybrid_wd_o_frac

       model % envelope_mass = hybrid_wd_he_shell_mass

       if (model % envelope_mass > ZERO) then
          model % envelope_comp(iHe4) = ONE
       else
          model % envelope_comp = model % core_comp
       endif

    else if (model % mass >= max_hybrid_wd_mass .and. model % mass < max_co_wd_mass) then

       model % core_comp(iC12) = co_wd_c_frac
       model % core_comp(iO16) = co_wd_o_frac

       model % envelope_mass = co_wd_he_shell_mass

       if (model % envelope_mass > ZERO) then
          model % envelope_comp(iHe4) = ONE
       else
          model % envelope_comp = model % core_comp
       endif

    else if (model % mass > max_co_wd_mass) then

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

  subroutine ensure_primary_mass_larger

    use problem_io_module, only: ioproc

    implicit none

    double precision :: temp_mass

    ! We want the primary WD to be more massive. If what we're calling
    ! the primary is less massive, switch the stars.

    if ( mass_P < mass_S ) then

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

    implicit none

    double precision, intent(inout) :: P_com(3), S_com(3)
    double precision, intent(inout) :: P_vel(3), S_vel(3)
    double precision, intent(inout) :: P_mass, S_mass
    double precision, intent(inout) :: P_t_ff, S_t_ff

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

    use bl_constants_module, only: TENTH, ZERO
    use prob_params_module, only: center
    use binary_module, only: get_roche_radii

    implicit none

    double precision, intent(in) :: P_com(3), S_com(3)
    double precision, intent(in) :: P_vel(3), S_vel(3)
    double precision, intent(in) :: P_mass, S_mass
    double precision, intent(in) :: P_t_ff, S_t_ff

    double precision :: r

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

    use meth_params_module, only: rot_period
    use initial_model_module, only: initialize_model, establish_hse
    use prob_params_module, only: center, problo, probhi, dim
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product
    use binary_module, only: get_roche_radii
    use problem_io_module, only: ioproc
    use bl_error_module, only: bl_error

    implicit none

    double precision :: v_ff, collision_offset
    double precision :: omega(3)

    omega = get_omega(ZERO)

    ! Set up the center variable. We want it to be at 
    ! problo + center_frac * domain_width in each direction.
    ! center_frac is 1/2 by default, so the problem
    ! would be set up exactly in the center of the domain.
    ! Note that we override this for 2D axisymmetric, as the
    ! radial coordinate must be centered at zero for the problem to make sense.

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

       call bl_error("Error: unknown value for dim in subroutine binary_setup.")

    endif

    ! Set some default values for these quantities;
    ! we'll update them soon.

    center_P_initial = center
    center_S_initial = center

    com_P = center
    com_S = center

    vel_P = ZERO
    vel_S = ZERO

    ! Allocate arrays to hold the stellar models.

    call initialize_model(model_P, initial_model_dx, initial_model_npts, initial_model_mass_tol, initial_model_hse_tol)
    call initialize_model(model_S, initial_model_dx, initial_model_npts, initial_model_mass_tol, initial_model_hse_tol)

    model_P % min_density = ambient_density
    model_S % min_density = ambient_density

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

       call establish_hse(model_P)

       call set_wd_composition(model_P)

    else

       call bl_error("Must specify either a positive primary mass or a positive primary central density.")

    endif



    if (.not. single_star) then

       if (mass_S > ZERO) then

          model_S % mass = mass_S

          call set_wd_composition(model_S)

       elseif (central_density_S > ZERO) then

          model_S % mass = M_solar

          call set_wd_composition(model_S)

          model_S % central_density = central_density_S

          call establish_hse(model_S)

          call set_wd_composition(model_S)

       else

          call bl_error("If we are doing a binary calculation, we must specify either a ", &
                        "positive secondary mass or a positive secondary central density.")

       endif

       ambient_comp = (model_P % envelope_comp + model_S % envelope_comp) / 2

    else

       ambient_comp = model_P % envelope_comp

    endif



    roche_rad_P = ZERO
    roche_rad_S = ZERO

    ! Generate primary and secondary WD models.

    call establish_hse(model_P)

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

       call establish_hse(model_S)

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

       if (problem == 0 .or. problem == 4) then

          collision_separation = collision_separation * model_S % radius

          if (problem == 0) then

             call freefall_velocity(mass_P + mass_S, collision_separation, v_ff)

             vel_P(axis_1) =  (mass_P / (mass_S + mass_P)) * v_ff
             vel_S(axis_1) = -(mass_S / (mass_S + mass_P)) * v_ff

          endif

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

       else if (problem == 1 .or. problem == 2 .or. problem == 3) then

          if (problem == 1) then

             ! Determine the orbital distance based on the rotational period.

             a = -ONE

          else if (problem == 2 .or. problem == 3) then

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
1003         format ("Generated binary orbit of distance ", ES8.2, " cm = ", ES8.2, " AU.")
1004         format ("The primary orbits the center of mass at distance ", ES9.2, " cm = ", ES9.2, " AU.")
1005         format ("The secondary orbits the center of mass at distance ", ES9.2, " cm = ", ES9.2, " AU.")
1006         format ("The initial orbital period is ", F6.2 " s.")
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

          call bl_error("Error: Unknown problem choice.")

       endif

       ! Scale the Roche radii by the initial distance.

       roche_rad_P = roche_rad_P * a
       roche_rad_S = roche_rad_S * a

    endif

    ! Reset the terminal color to its previous state.

    if (ioproc .and. init == 1) then
       print *, ''//achar(27)//'[0m'
    endif

    com_P = center_P_initial
    com_S = center_S_initial

    ! Safety check: make sure the stars are actually inside the computational domain.

    if ( ( center_P_initial(1) - model_P % radius .lt. problo(1) .and. dim .ne. 2 ) .or. &
         center_P_initial(1) + model_P % radius .gt. probhi(1) .or. &
         ( center_P_initial(2) - model_P % radius .lt. problo(2) .and. dim .ge. 2 ) .or. &
         ( center_P_initial(2) + model_P % radius .gt. probhi(2) .and. dim .ge. 2 ) .or. &
         ( center_P_initial(3) - model_P % radius .lt. problo(3) .and. dim .eq. 3 ) .or. &
         ( center_P_initial(3) + model_P % radius .gt. probhi(3) .and. dim .eq. 3 ) ) then
       call bl_error("Primary does not fit inside the domain.")
    endif

    if ( ( center_S_initial(1) - model_S % radius .lt. problo(1) .and. dim .ne. 2 ) .or. &
         center_S_initial(1) + model_S % radius .gt. probhi(1) .or. &
         ( center_S_initial(2) - model_S % radius .lt. problo(2) .and. dim .ge. 2 ) .or. &
         ( center_S_initial(2) + model_S % radius .gt. probhi(2) .and. dim .ge. 2 ) .or. &
         ( center_S_initial(3) - model_S % radius .lt. problo(3) .and. dim .eq. 3 ) .or. &
         ( center_S_initial(3) + model_S % radius .gt. probhi(3) .and. dim .eq. 3 ) ) then
       call bl_error("Secondary does not fit inside the domain.")
    endif


  end subroutine binary_setup



  ! Accepts the masses of two stars (in solar masses)
  ! and the orbital period of a system,
  ! and returns the semimajor axis of the orbit (in cm),
  ! as well as the distances a_1 and a_2 from the center of mass.

  subroutine kepler_third_law(radius_1, mass_1, radius_2, mass_2, period, eccentricity, phi, a, r_1, r_2, v_1r, v_2r, v_1p, v_2p)

    use bl_constants_module
    use prob_params_module, only: problo, probhi
    use sponge_module, only: sponge_lower_radius
    use meth_params_module, only: do_sponge
    use fundamental_constants_module, only: Gconst

    implicit none

    double precision, intent(in   ) :: mass_1, mass_2, eccentricity, phi, radius_1, radius_2
    double precision, intent(inout) :: period, a, r_1, r_2, v_1r, v_2r, v_1p, v_2p

    double precision :: length

    double precision :: mu, M ! Reduced mass, total mass
    double precision :: r     ! Position
    double precision :: v_r, v_phi ! Radial and azimuthal velocity

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

       call bl_error("Error: overspecified Kepler's third law calculation.")

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

    length = (r_2 - r_1) + radius_1 + radius_2

    if (length > (probhi(axis_1)-problo(axis_1))) then
       call bl_error("ERROR: The domain width is too small to include the binary orbit.")
    endif

    ! We want to do a similar check to make sure that no part of the stars
    ! land in the sponge region.

    if (do_sponge .eq. 1 .and. sponge_lower_radius > ZERO) then

       if (abs(r_1) + radius_1 .ge. sponge_lower_radius) then
          call bl_error("ERROR: Primary contains material inside the sponge region.")
       endif

       if (abs(r_2) + radius_2 .ge. sponge_lower_radius) then
          call bl_error("ERROR: Secondary contains material inside the sponge region.")
       endif

    endif

    ! Make sure the stars are not touching.
    if (radius_1 + radius_2 > a) then
       call bl_error("ERROR: Stars are touching!")
    endif

  end subroutine kepler_third_law



  ! Given total mass of a binary system and the initial separation of
  ! two point particles, obtain the velocity at this separation 
  ! assuming the point masses fell in from infinity. This will
  ! be the velocity in the frame where the center of mass is stationary.

  subroutine freefall_velocity(mass, distance, vel)

    use bl_constants_module, only: HALF, TWO
    use fundamental_constants_module, only: Gconst

    implicit none

    double precision, intent(in   ) :: mass, distance
    double precision, intent(inout) :: vel

    vel = (TWO * Gconst * mass / distance)**HALF

  end subroutine freefall_velocity



  ! Given a zone state, fill it with ambient material.

  subroutine fill_ambient(state, loc, time)

    use bl_constants_module, only: ZERO
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UTEMP, UEINT, UEDEN, UFS, do_rotation
    use network, only: nspec
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product

    implicit none

    double precision :: state(NVAR)
    double precision :: loc(3), time

    type (eos_t) :: ambient_state
    double precision :: omega(3)

    omega = get_omega(time)

    call get_ambient(ambient_state)

    state(URHO) = ambient_state % rho
    state(UTEMP) = ambient_state % T
    state(UFS:UFS-1+nspec) = ambient_state % rho * ambient_state % xn(:)                 

    ! If we're in the inertial frame, give the material the rigid-body rotation speed.
    ! Otherwise set it to zero.

    if ( (do_rotation .ne. 1) .and. (problem .eq. 1 .or. problem .eq. 2 .or. problem .eq. 3) ) then

       state(UMX:UMZ) = state(URHO) * cross_product(omega, loc)

    else

       state(UMX:UMZ) = ZERO

    endif

    state(UEINT) = ambient_state % rho * ambient_state % e
    state(UEDEN) = state(UEINT) + HALF * sum(state(UMX:UMZ)**2) / state(URHO)

  end subroutine fill_ambient



  ! If we are in a rotating reference frame, then rotate a vector
  ! by an amount corresponding to the time that has passed
  ! since the beginning of the simulation.

  function inertial_rotation(vec, time) result(vec_i)

    use bl_constants_module, only: ZERO
    use rotation_frequency_module, only: get_omega
    use meth_params_module, only: do_rotation, rot_period, rot_period_dot

    implicit none

    double precision :: vec(3), time

    double precision :: vec_i(3)

    double precision :: omega(3), theta(3), rot_matrix(3,3)


    ! To get the angle, we integrate omega over the time of the
    ! simulation. Since the time rate of change is linear in the
    ! period, let's work that variable. At some time t the current
    ! period P is given by P = P_0 + Pdot * t. Then:
    !
    ! theta(t) = int( omega(t) * dt )
    !      theta = int( omega(t) dt )
    !            = (2 * pi / P_0) * int( dt / (1 + (dPdt / P_0) * t) )
    !
    ! if dPdt = 0, then theta = 2 * pi * t / P_0 = omega_0 * t, as expected.
    ! if dPdt > 0, then theta = (2 * pi / P_0) * (P_0 / dPdt) * ln| (dPdt / P_0) * t + 1 |
    ! Note that if dPdt << P_0, then we have ln(1 + x) = x, and we again
    ! recover the original expression as expected.

    if (do_rotation .eq. 1) then

       if (abs(rot_period_dot) > ZERO .and. time > ZERO) then
          theta = get_omega(ZERO) * (rot_period / rot_period_dot) * &
                  log( abs( (rot_period_dot / rot_period) * time + 1 ) )
       else
          theta = get_omega(ZERO) * time
       endif

       omega = get_omega(time)

    else

       omega = ZERO
       theta = ZERO

    endif

    ! This is the 3D rotation matrix for converting between reference frames.
    ! It is the composition of rotations along the x, y, and z axes. Therefore 
    ! it allows for the case where we are rotating about multiple axes. Normally 
    ! we use the right-hand convention for constructing the usual rotation matrix, 
    ! but this is the transpose of that rotation matrix to account for the fact 
    ! that we are rotating *back* to the inertial frame, rather than from the 
    ! inertial frame to the rotating frame.

    rot_matrix(1,1) =  cos(theta(2)) * cos(theta(3))
    rot_matrix(1,2) = -cos(theta(2)) * sin(theta(3))
    rot_matrix(1,3) =  sin(theta(2))
    rot_matrix(2,1) =  cos(theta(1)) * sin(theta(3)) + sin(theta(1)) * sin(theta(2)) * cos(theta(3))
    rot_matrix(2,2) =  cos(theta(1)) * cos(theta(3)) - sin(theta(1)) * sin(theta(2)) * sin(theta(3))
    rot_matrix(2,3) = -sin(theta(1)) * cos(theta(2))
    rot_matrix(3,1) =  sin(theta(1)) * sin(theta(3)) - cos(theta(1)) * sin(theta(2)) * cos(theta(3))
    rot_matrix(3,2) =  sin(theta(1)) * cos(theta(3)) + cos(theta(1)) * sin(theta(2)) * sin(theta(3))
    rot_matrix(3,3) =  cos(theta(1)) * cos(theta(2))

    vec_i = matmul(rot_matrix, vec)

  end function inertial_rotation



  ! Given a rotating frame velocity, get the inertial frame velocity.
  ! Note that we simply return the original velocity if we're
  ! already in the inertial frame.

  function inertial_velocity(loc, vel, time) result (vel_i)

    use meth_params_module, only: do_rotation, state_in_rotating_frame
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product

    implicit none

    double precision :: loc(3), vel(3), time
    double precision :: omega(3)

    double precision :: vel_i(3)

    omega = get_omega(time)

    vel_i = vel

    if (do_rotation .eq. 1 .and. state_in_rotating_frame .eq. 1) then
       vel_i = vel_i + cross_product(omega, loc)
    endif

  end function inertial_velocity



  ! Check whether we should stop the initial relaxation.
  ! The criterion is that we're outside the critical Roche surface
  ! and the density is greater than a specified threshold.
  ! If so, set do_initial_relaxation to false, which will effectively
  ! turn off the external source terms.

  subroutine check_relaxation(state, s_lo, s_hi, &
                              phiEff, p_lo, p_hi, &
                              lo, hi, potential, is_done) bind(C,name='check_relaxation')

    use meth_params_module, only: URHO, NVAR
    use castro_util_module, only: position_to_index

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    integer          :: p_lo(3), p_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: phiEff(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision :: potential
    integer          :: is_done

    integer          :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (phiEff(i,j,k) > potential .and. state(i,j,k,URHO) > relaxation_density_cutoff) then

                is_done = 1

             endif

          enddo
       enddo
    enddo

  end subroutine check_relaxation



  ! This routine is called when we've satisfied our criterion
  ! for disabling the initial relaxation phase. We set the
  ! relaxation timescale to a negative number, which disables
  ! the damping, and we set the sponge timescale to a negative
  ! number, which disables the sponging.

  subroutine turn_off_relaxation(time) bind(C,name='turn_off_relaxation')

    use problem_io_module, only: ioproc
    use sponge_module, only: sponge_timescale
    use meth_params_module, only: rotation_include_coriolis

    implicit none

    double precision :: time

    relaxation_damping_timescale = -ONE
    sponge_timescale = -ONE
    rotation_include_coriolis = 1

    ! If we got a valid simulation time, print to the log when we stopped.

    if (ioproc .and. time >= 0.0d0) then
       print *, ""
       print *, "Initial relaxation phase terminated at t = ", time
       print *, ""
    endif

  end subroutine turn_off_relaxation



  subroutine get_axes(axis_1_in, axis_2_in, axis_3_in) bind(C,name='get_axes')

    implicit none

    integer :: axis_1_in, axis_2_in, axis_3_in

    axis_1_in = axis_1
    axis_2_in = axis_2
    axis_3_in = axis_3

  end subroutine get_axes



  ! Return whether we're doing a single star simulation or not.

  subroutine get_single_star(flag) bind(C,name='get_single_star')

    implicit none

    integer :: flag

    flag = 0

    if (single_star) flag = 1

  end subroutine get_single_star



  ! Return the problem type.

  subroutine get_problem_number(problem_out) bind(C,name='get_problem_number')

    implicit none

    integer :: problem_out

    problem_out = problem

  end subroutine get_problem_number



  ! Computes the sum of the hydrodynamic and gravitational forces acting on the WDs.

  subroutine sum_force_on_stars(lo, hi, &
                                force, f_lo, f_hi, &
                                state, s_lo, s_hi, &
                                vol, v_lo, v_hi, &
                                pmask, pm_lo, pm_hi, &
                                smask, sm_lo, sm_hi, &
                                fpx, fpy, fpz, fsx, fsy, fsz) &
                                bind(C,name='sum_force_on_stars')

    use bl_constants_module, only: ZERO
    use prob_params_module, only: center
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ
    use castro_util_module, only: position

    implicit none

    integer :: lo(3), hi(3)
    integer :: f_lo(3), f_hi(3)
    integer :: s_lo(3), s_hi(3)
    integer :: v_lo(3), v_hi(3)
    integer :: pm_lo(3), pm_hi(3)
    integer :: sm_lo(3), sm_hi(3)

    double precision :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3), NVAR)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3), NVAR)
    double precision :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    double precision :: pmask(pm_lo(1):pm_hi(1),pm_lo(2):pm_hi(2),pm_lo(3):pm_hi(3))
    double precision :: smask(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))

    double precision :: fpx, fpy, fpz, fsx, fsy, fsz
    double precision :: dt

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (pmask(i,j,k) > ZERO) then

                fpx = fpx + vol(i,j,k) * force(i,j,k,UMX)
                fpy = fpy + vol(i,j,k) * force(i,j,k,UMY)
                fpz = fpz + vol(i,j,k) * force(i,j,k,UMZ)

             else if (smask(i,j,k) > ZERO) then

                fsx = fsx + vol(i,j,k) * force(i,j,k,UMX)
                fsy = fsy + vol(i,j,k) * force(i,j,k,UMY)
                fsz = fsz + vol(i,j,k) * force(i,j,k,UMZ)

             endif

          enddo
       enddo
    enddo

  end subroutine sum_force_on_stars



  subroutine find_ignited_zones(lo, hi, &
                                ignition_radius, ir_lo, ir_hi, &
                                state, s_lo, s_hi, &
                                dx, num_zones, level) &
                                bind(C,name='find_ignited_zones')

    use meth_params_module, only: NVAR, URHO, UTEMP, UFS
    use prob_params_module, only: dim
    use network, only: network_species_index

    implicit none

    integer :: lo(3), hi(3)
    integer :: ir_lo(3), ir_hi(3)
    integer :: s_lo(3), s_hi(3)

    double precision :: ignition_radius(ir_lo(1):ir_hi(1),ir_lo(2):ir_hi(2),ir_lo(3):ir_hi(3))
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    double precision :: dx(3)

    integer :: num_zones, level

    integer :: i, j, k

    double precision, parameter :: safety_factor = 0.1d0

    integer :: iC12

    iC12 = network_species_index("carbon-12")

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (maxval(dx(1:dim)) < safety_factor * ignition_radius(i,j,k)) then

                num_zones = num_zones + 1

                print *, "  Zone ignited on level ", level, " with ignition radius ", ignition_radius(i,j,k), " at zone "
                print *, "  Zone ignited at (i,j,k) = ", i, j, k
                print *, "  Zone ignited with T = ", state(i,j,k,UTEMP)
                print *, "  Zone ignited with rho = ", state(i,j,k,URHO)
                print *, "  Zone ignited with X_C = ", state(i,j,k,UFS+iC12-1) / state(i,j,k,URHO)

             end if

          enddo
       enddo
    enddo

  end subroutine find_ignited_zones

end module wdmerger_util_module
