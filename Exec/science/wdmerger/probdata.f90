module probdata_module

  use network, only: nspec, network_species_index
  use eos_type_module, only: eos_t, eos_input_rt
  use eos_module, only: eos
  use bl_constants_module, only: ZERO, THIRD, HALF, ONE, TWO, THREE, M_PI, FOUR
  use fundamental_constants_module, only: Gconst, M_solar, AU
  use initial_model_module, only: initial_model

  ! Initial stellar properties
  ! Note that the envelope mass is included within the total mass of the star

  double precision, save :: mass_P = ONE
  double precision, save :: mass_S = ONE
  double precision, save :: central_density_P = -ONE
  double precision, save :: central_density_S = -ONE
  double precision, save :: stellar_temp = 1.0d7
  double precision, save :: primary_envelope_mass, secondary_envelope_mass
  double precision, save :: primary_envelope_comp(nspec), secondary_envelope_comp(nspec)



  ! Ambient medium

  double precision, save :: ambient_density = 1.0d-4
  double precision, save :: ambient_temp = 1.0d7
  double precision, save :: ambient_comp(nspec)



  ! Smallest allowed velocity on the grid

  double precision, save :: smallu = ZERO



  ! Parameters for interpolation from 1D model to 3D model:

  ! Number of sub-grid-scale zones to use

  integer, save :: nsub = 1

  ! Default to interpolation that preserves temperature; otherwise, use pressure

  logical, save :: interp_temp = .true.



  ! Method for determining the initial problem setup.
  !
  ! 0 = Collision; distance determined by a multiple of the secondary WD radius
  ! 1 = Keplerian orbit; distance determined by the rotation period
  ! 2 = Keplerian orbit; distance set so that the secondary exactly fills its Roche lobe radius
  ! 3 = Problem 2 with an initial relaxation step
  ! 4 = Free-fall; distance determined by a multiple of the secondary WD radius

  integer, save :: problem = 2



  ! If we're automatically determining the initial distance based on the Roche lobe
  ! radii for problem 3, this is the sizing factor we use. We set the default value
  ! to be a large enough distance so that the system is close to stable.

  double precision, save :: roche_radius_factor = 2.0d0



  ! Collision parameters

  ! For a collision, number of (secondary) WD radii to 
  ! separate the WDs by.

  double precision, save :: collision_separation = 4.0d0

  ! For a collision, the impact parameter measured in
  ! units of the primary's initial radius.

  double precision, save :: collision_impact_parameter = 0.0d0



  ! Binary orbit properties

  double precision, save :: r_P_initial, r_S_initial, a_P_initial, a_S_initial, a  
  double precision, save :: v_P_r, v_S_r, v_P_phi, v_S_phi
  double precision, save :: center_P_initial(3), center_S_initial(3)
  double precision, save :: orbital_eccentricity = 0.0d0
  double precision, save :: orbital_angle = 0.0d0



  ! Axis is in orbital plane; we measure angle with respect to this axis. Normally the x axis.

  integer, save :: axis_1 = 1

  ! Perpendicular axis in the orbital plane. Normally the y axis.

  integer, save :: axis_2 = 2

  ! Perpendicular to both other axes. Normally the z axis and also the rotation axis.

  integer, save :: axis_3 = 3

  ! Location of the physical center of the problem, as a fraction of domain size

  double precision, save :: center_fracx = HALF
  double precision, save :: center_fracy = HALF
  double precision, save :: center_fracz = HALF

  ! Bulk system motion

  double precision, save :: bulk_velx = ZERO
  double precision, save :: bulk_vely = ZERO
  double precision, save :: bulk_velz = ZERO

  ! Whether we're doing an initialization or a restart

  integer, save :: init

  ! Are we doing a single star simulation?

  logical, save :: single_star

  ! Should we override the domain boundary conditions with
  ! ambient material?

  logical, save :: fill_ambient_bc = .false.



  ! 1D initial models

  type (initial_model) :: model_P, model_S

  ! For the grid spacing for our model, we'll use 
  ! 6.25 km. No simulation we do is likely to have a resolution
  ! higher than that inside the stars (it represents
  ! three jumps by a factor of four compared to our 
  ! normal coarse grid resolution). By using 4096
  ! grid points, the size of the 1D domain will be 2.56e9 cm,
  ! which is larger than any reasonable mass white dwarf.

  double precision, save :: initial_model_dx = 6.25d5
  integer, save          :: initial_model_npts = 4096

  ! initial_model_mass_tol is tolerance used for getting the total WD mass 
  ! equal to the desired mass. It can be reasonably small, since there
  ! will always be a central density value that can give the desired
  ! WD mass on the grid we use.

  double precision, save :: initial_model_mass_tol = 1.d-6

  ! hse_tol is the tolerance used when iterating over a zone to force
  ! it into HSE by adjusting the current density (and possibly
  ! temperature).  hse_tol should be very small (~ 1.e-10).

  double precision, save :: initial_model_hse_tol = 1.d-10



  ! Composition properties of initial models.
  ! We follow the prescription of Dan et al. 2012 for determining
  ! the composition of the WDs. In this approach, below a certain 
  ! mass there are pure He WDs; just above that are hybrid WDs
  ! with pure CO cores and a He mantle; above that are pure CO WDs
  ! with slightly more oxygen than carbon; and above that are 
  ! ONeMg WDs. All masses are in solar masses.

  double precision, save :: max_he_wd_mass = 0.45d0
  double precision, save :: max_hybrid_wd_mass = 0.6d0
  double precision, save :: hybrid_wd_he_shell_mass = 0.1d0
  double precision, save :: max_co_wd_mass = 1.05d0
  double precision, save :: co_wd_he_shell_mass = 0.0d0

  double precision, save :: hybrid_wd_c_frac = 0.50d0
  double precision, save :: hybrid_wd_o_frac = 0.50d0

  double precision, save :: co_wd_c_frac = 0.40d0
  double precision, save :: co_wd_o_frac = 0.60d0

  double precision, save :: onemg_wd_o_frac  = 0.60d0
  double precision, save :: onemg_wd_ne_frac = 0.35d0
  double precision, save :: onemg_wd_mg_frac = 0.05d0


  ! Tagging criteria

  integer,          save :: max_stellar_tagging_level = 20
  integer,          save :: max_temperature_tagging_level = 20
  integer,          save :: max_center_tagging_level = 20
  double precision, save :: stellar_density_threshold = 1.0d0
  double precision, save :: temperature_tagging_threshold = 5.0d8
  double precision, save :: center_tagging_radius = 0.0d0
  double precision, save :: max_tagging_radius = 0.75d0
  double precision, save :: roche_tagging_factor = 2.0d0



  ! Stores the center of mass location of the stars throughout the run

  double precision, save :: com_P(3), com_S(3)
  double precision, save :: vel_P(3), vel_S(3)

  ! Stores the effective Roche radii

  double precision, save :: roche_rad_P, roche_rad_S



  ! Number of zones that have passed the critical radius for thermonuclear ignition
  integer,          save :: num_zones_ignited = 0

  ! Level on which the initial ignition happened; negative means it hasn't happened yet
  integer,          save :: ignition_level = -1

  ! Relaxation parameters for problem 3

  double precision, save :: relaxation_damping_timescale = -1.0d0
  double precision, save :: relaxation_density_cutoff = 1.0d3
  logical,          save :: relaxation_implicit = .false.
  integer,          save :: relaxation_is_done = 0

  ! Radial damping parameters for problem 3

  double precision, save :: radial_damping_factor = 1.0d3
  double precision, save :: initial_radial_velocity_factor = -1.0d-3

  ! Distance (in kpc) used for calculation of the gravitational wave amplitude
  ! (this wil be calculated along all three coordinate axes).

  double precision, save :: gw_dist = 10.0d0

  ! Current value of the dynamical timescale for each star

  double precision, save :: t_ff_P, t_ff_S

  ! Global extrema

  double precision, save :: T_global_max, rho_global_max, ts_te_global_max



  namelist /fortin/ &
       mass_P, mass_S, &
       central_density_P, central_density_S, &
       nsub, &
       roche_radius_factor, &
       problem, &
       collision_separation, &
       collision_impact_parameter, &
       interp_temp, &
       relaxation_damping_timescale, &
       relaxation_density_cutoff, &
       relaxation_implicit, &
       initial_radial_velocity_factor, &
       radial_damping_factor, &
       ambient_density, &
       stellar_temp, &
       ambient_temp, &
       max_he_wd_mass, &
       max_hybrid_wd_mass, hybrid_wd_he_shell_mass, &
       max_co_wd_mass, &
       co_wd_he_shell_mass, &
       hybrid_wd_c_frac, hybrid_wd_o_frac, &
       co_wd_c_frac, co_wd_o_frac, &
       onemg_wd_o_frac, onemg_wd_ne_frac, onemg_wd_mg_frac, &
       orbital_eccentricity, orbital_angle, &
       axis_1, axis_2, axis_3, &
       max_stellar_tagging_level, &
       max_temperature_tagging_level, &
       max_center_tagging_level, &
       stellar_density_threshold, &
       temperature_tagging_threshold, &
       center_tagging_radius, &
       max_tagging_radius, &
       roche_tagging_factor, &
       bulk_velx, bulk_vely, bulk_velz, &
       smallu, &
       center_fracx, center_fracy, center_fracz, &
       initial_model_dx, &
       initial_model_npts, &
       initial_model_mass_tol, &
       initial_model_hse_tol, &
       gw_dist, &
       fill_ambient_bc


  ! Stores whether we assert that the simulation has completed.

  logical, save :: jobIsDone = .false.
  logical, save :: signalJobIsNotDone = .false.

  ! Auxiliary data for determining whether the job is done.

  integer, parameter :: num_previous_ener_timesteps = 5
  double precision :: total_ener_array(num_previous_ener_timesteps)

end module probdata_module
