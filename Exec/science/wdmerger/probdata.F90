module probdata_module

  use amrex_fort_module, only: rt => amrex_real
  use network, only: nspec, network_species_index
  use eos_type_module, only: eos_t, eos_input_rt
  use eos_module, only: eos
  use amrex_constants_module, only: ZERO, THIRD, HALF, ONE, TWO, THREE, M_PI, FOUR, SIX, EIGHT
  use fundamental_constants_module, only: Gconst, M_solar, AU
  use initial_model_module, only: initial_model

  ! Initial stellar properties
  ! Note that the envelope mass is included within the total mass of the star

  real(rt), allocatable :: mass_P
  real(rt), allocatable :: mass_S
  real(rt), save :: central_density_P = -ONE
  real(rt), save :: central_density_S = -ONE
  real(rt), save :: stellar_temp = 1.0e7_rt
  real(rt), save :: primary_envelope_mass, secondary_envelope_mass
  real(rt), save :: primary_envelope_comp(nspec), secondary_envelope_comp(nspec)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: mass_P, mass_S
#endif


  ! Ambient medium

  real(rt), save :: ambient_density = 1.0e-4_rt
  real(rt), save :: ambient_temp = 1.0e7_rt
  real(rt), save :: ambient_comp(nspec)



  ! Smallest allowed velocity on the grid

  real(rt), save :: smallu = ZERO



  ! Parameters for interpolation from 1D model to 3D model:

  ! Number of sub-grid-scale zones to use

  integer, save :: nsub = 1

  ! Default to interpolation that preserves temperature; otherwise, use pressure

  logical, save :: interp_temp = .true.



  ! Method for determining the initial problem setup.
  !
  ! 0 = Collision; distance determined by a multiple of the secondary WD radius
  ! 1 = Keplerian orbit; distance determined by Roche radius (or rotation period)
  ! 5 = Tidal disruption event; distance determined by a multiple of the WD tidal radius

  integer, save :: problem = 1



  ! If we're automatically determining the initial distance based on the Roche lobe
  ! radii for the merger problem, this is the sizing factor we use. Negative means
  ! that we set the initial distance using the user-selected rotation period.

  real(rt), save :: roche_radius_factor = ONE



  ! Collision parameters

  ! For a collision, number of (secondary) WD radii to 
  ! separate the WDs by.

  real(rt), save :: collision_separation = FOUR

  ! For a collision, the impact parameter measured in
  ! units of the primary's initial radius.

  real(rt), save :: collision_impact_parameter = ZERO

  ! For a collision, the initial velocity of the WDs toward
  ! each other. If this is negative, the velocity will
  ! be set according to free-fall from an infinite distance.

  real(rt), save :: collision_velocity = -ONE



  ! TDE parameters

  ! For a TDE, number of WD tidal radii to separate the WD and BH.

  real(rt), save :: tde_separation = EIGHT

  ! For a TDE, the parameter beta: the ratio of the tidal radius to
  ! the Schwarzschild radius of the BH.

  real(rt), save :: tde_beta = SIX

  ! For a TDE, should we give the star an initial kick of velocity
  ! corresponding to its parabolic orbit? By default we will, but
  ! this option exists so we can test for HSE.

  integer, save :: tde_initial_velocity = 1

  real(rt), save :: tde_tidal_radius
  real(rt), save :: tde_schwarzschild_radius
  real(rt), save :: tde_pericenter_radius



  ! Binary orbit properties

  real(rt), save :: r_P_initial, r_S_initial, a
  real(rt), save :: center_P_initial(3), center_S_initial(3)
  real(rt), save :: orbital_eccentricity = ZERO
  real(rt), save :: orbital_angle = ZERO



  ! Axis is in orbital plane; we measure angle with respect to this axis. Normally the x axis.

  integer, save :: axis_1 = 1

  ! Perpendicular axis in the orbital plane. Normally the y axis.

  integer, save :: axis_2 = 2

  ! Perpendicular to both other axes. Normally the z axis and also the rotation axis.

  integer, save :: axis_3 = 3

  ! Location of the physical center of the problem, as a fraction of domain size

  real(rt), save :: center_fracx = HALF
  real(rt), save :: center_fracy = HALF
  real(rt), save :: center_fracz = HALF

  ! Bulk system motion

  real(rt), save :: bulk_velx = ZERO
  real(rt), save :: bulk_vely = ZERO
  real(rt), save :: bulk_velz = ZERO

  ! Whether we're doing an initialization or a restart

  integer, save :: init

  ! Are we doing a single star simulation?

  logical, save :: single_star



  ! 1D initial models

  type (initial_model) :: model_P, model_S

  ! For the grid spacing for our model, we'll use 
  ! 6.25 km. No simulation we do is likely to have a resolution
  ! higher than that inside the stars (it represents
  ! three jumps by a factor of four compared to our 
  ! normal coarse grid resolution). By using 4096
  ! grid points, the size of the 1D domain will be 2.56e9 cm,
  ! which is larger than any reasonable mass white dwarf.

  real(rt), save :: initial_model_dx = 6.25e5_rt
  integer,  save :: initial_model_npts = 4096

  ! initial_model_mass_tol is tolerance used for getting the total WD mass 
  ! equal to the desired mass. It can be reasonably small, since there
  ! will always be a central density value that can give the desired
  ! WD mass on the grid we use.

  real(rt), save :: initial_model_mass_tol = 1.e-6_rt

  ! hse_tol is the tolerance used when iterating over a zone to force
  ! it into HSE by adjusting the current density (and possibly
  ! temperature).  hse_tol should be very small (~ 1.e-10).

  real(rt), save :: initial_model_hse_tol = 1.e-10_rt



  ! Composition properties of initial models.
  ! We follow the prescription of Dan et al. 2012 for determining
  ! the composition of the WDs. In this approach, below a certain 
  ! mass there are pure He WDs; just above that are hybrid WDs
  ! with pure CO cores and a He mantle; above that are pure CO WDs
  ! with slightly more oxygen than carbon; and above that are 
  ! ONeMg WDs. All masses are in solar masses.

  real(rt), save :: max_he_wd_mass = 0.45e0_rt
  real(rt), save :: max_hybrid_wd_mass = 0.6e0_rt
  real(rt), save :: hybrid_wd_he_shell_mass = 0.1e0_rt
  real(rt), save :: max_co_wd_mass = 1.05e0_rt
  real(rt), save :: co_wd_he_shell_mass = 0.0e0_rt

  real(rt), save :: hybrid_wd_c_frac = 0.50e0_rt
  real(rt), save :: hybrid_wd_o_frac = 0.50e0_rt

  real(rt), save :: co_wd_c_frac = 0.40e0_rt
  real(rt), save :: co_wd_o_frac = 0.60e0_rt

  real(rt), save :: onemg_wd_o_frac  = 0.60e0_rt
  real(rt), save :: onemg_wd_ne_frac = 0.35e0_rt
  real(rt), save :: onemg_wd_mg_frac = 0.05e0_rt


  ! Tagging criteria

  integer,  save :: max_stellar_tagging_level = 20
  integer,  save :: max_temperature_tagging_level = 20
  integer,  save :: max_center_tagging_level = 20
  real(rt), allocatable :: stellar_density_threshold
  real(rt), save :: temperature_tagging_threshold = 5.0e8_rt
  real(rt), save :: center_tagging_radius = 0.0e0_rt
  real(rt), save :: max_tagging_radius = 0.75e0_rt
  real(rt), save :: roche_tagging_factor = 2.0e0_rt

#ifdef AMREX_USE_CUDA
  attributes(managed) :: stellar_density_threshold
#endif



  ! Stores the center of mass location of the stars throughout the run

  real(rt), allocatable :: com_P(:), com_S(:)
  real(rt), save :: vel_P(3), vel_S(3)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: com_P, com_S
#endif

  ! Stores the effective Roche radii

  real(rt), save :: roche_rad_P, roche_rad_S



  ! Relaxation parameters

  real(rt), save :: relaxation_damping_factor = -1.0e-1_rt
  real(rt), save :: relaxation_density_cutoff = 1.0e3_rt
  real(rt), save :: relaxation_cutoff_time = -1.e0_rt
  integer,  save :: relaxation_is_done = 1

  ! Radial damping parameters

  real(rt), save :: radial_damping_factor = -1.0e3_rt
  real(rt), save :: initial_radial_velocity_factor = -1.0e-3_rt

  ! Distance (in kpc) used for calculation of the gravitational wave amplitude
  ! (this wil be calculated along all three coordinate axes).

  real(rt), save :: gw_dist = 10.0e0_rt

  ! Current value of the dynamical timescale for each star

  real(rt), save :: t_ff_P, t_ff_S

  ! Global extrema

  real(rt), save :: T_global_max, rho_global_max, ts_te_global_max

  ! Stores whether we assert that the simulation has completed.

  logical, save :: jobIsDone = .false.
  logical, save :: signalJobIsNotDone = .false.

  ! Auxiliary data for determining whether the job is done.

  integer, parameter :: num_previous_ener_timesteps = 5
  real(rt) :: total_ener_array(num_previous_ener_timesteps)
  
end module probdata_module
