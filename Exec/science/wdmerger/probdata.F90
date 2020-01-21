module probdata_module

  use amrex_fort_module, only: rt => amrex_real
  use eos_type_module, only: eos_t

  type :: initial_model

     ! Physical characteristics

     real(rt) :: mass
     real(rt) :: envelope_mass
     real(rt) :: central_density
     real(rt) :: central_temp
     real(rt) :: min_density
     real(rt) :: radius

     ! Composition

     real(rt), allocatable :: core_comp(:)
     real(rt), allocatable :: envelope_comp(:)

     ! Model storage

     real(rt) :: dx
     integer  :: npts
     real(rt) :: mass_tol, hse_tol

     real(rt), allocatable :: r(:), rl(:), rr(:)
     real(rt), allocatable :: M_enclosed(:), g(:)
     type (eos_t), allocatable :: state(:)

  end type initial_model
  
  ! Initial stellar properties
  ! Note that the envelope mass is included within the total mass of the star

  real(rt), allocatable :: mass_P, mass_S
  real(rt), allocatable :: central_density_P, central_density_S
  real(rt), allocatable :: stellar_temp
  real(rt), allocatable :: primary_envelope_mass, secondary_envelope_mass
  real(rt), allocatable :: primary_envelope_comp(:), secondary_envelope_comp(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: mass_P, mass_S
  attributes(managed) :: central_density_P, central_density_S
  attributes(managed) :: stellar_temp
  attributes(managed) :: primary_envelope_mass, secondary_envelope_mass
  attributes(managed) :: primary_envelope_comp, secondary_envelope_comp
#endif



  ! Ambient medium

  real(rt), allocatable :: ambient_density, ambient_temp, ambient_comp(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: ambient_density, ambient_temp, ambient_comp
#endif



  ! Smallest allowed velocity on the grid

  real(rt), allocatable :: smallu

#ifdef AMREX_USE_CUDA
  attributes(managed) :: smallu
#endif



  ! Parameters for interpolation from 1D model to 3D model:

  ! Number of sub-grid-scale zones to use

  integer, allocatable :: nsub

  ! Default to interpolation that preserves temperature; otherwise, use pressure

  logical, allocatable :: interp_temp

#ifdef AMREX_USE_CUDA
  attributes(managed) :: nsub, interp_temp
#endif



  ! Method for determining the initial problem setup.
  !
  ! 0 = Collision; distance determined by a multiple of the secondary WD radius
  ! 1 = Keplerian orbit; distance determined by Roche radius (or rotation period)
  ! 5 = Tidal disruption event; distance determined by a multiple of the WD tidal radius

  integer, allocatable :: problem

#ifdef AMREX_USE_CUDA
  attributes(managed) :: problem
#endif



  ! If we're automatically determining the initial distance based on the Roche lobe
  ! radii for the merger problem, this is the sizing factor we use. Negative means
  ! that we set the initial distance using the user-selected rotation period.

  real(rt), allocatable :: roche_radius_factor

#ifdef AMREX_USE_CUDA
  attributes(managed) :: roche_radius_factor
#endif



  ! Collision parameters

  ! For a collision, number of (secondary) WD radii to 
  ! separate the WDs by.

  real(rt), allocatable :: collision_separation

  ! For a collision, the impact parameter measured in
  ! units of the primary's initial radius.

  real(rt), allocatable :: collision_impact_parameter

  ! For a collision, the initial velocity of the WDs toward
  ! each other. If this is negative, the velocity will
  ! be set according to free-fall from an infinite distance.

  real(rt), allocatable :: collision_velocity

#ifdef AMREX_USE_CUDA
  attributes(managed) :: collision_separation, collision_impact_parameter, collision_velocity
#endif



  ! TDE parameters

  ! For a TDE, number of WD tidal radii to separate the WD and BH.

  real(rt), allocatable :: tde_separation

  ! For a TDE, the parameter beta: the ratio of the tidal radius to
  ! the Schwarzschild radius of the BH.

  real(rt), allocatable :: tde_beta

  ! For a TDE, should we give the star an initial kick of velocity
  ! corresponding to its parabolic orbit? By default we will, but
  ! this option exists so we can test for HSE.

  integer, allocatable :: tde_initial_velocity

  real(rt), allocatable :: tde_tidal_radius
  real(rt), allocatable :: tde_schwarzschild_radius
  real(rt), allocatable :: tde_pericenter_radius

#ifdef AMREX_USE_CUDA
  attributes(managed) :: tde_separation, tde_beta, tde_initial_velocity
  attributes(managed) :: tde_tidal_radius, tde_schwarzschild_radius, tde_pericenter_radius
#endif



  ! Binary orbit properties

  real(rt), allocatable :: r_P_initial, r_S_initial, a
  real(rt), allocatable :: center_P_initial(:), center_S_initial(:)
  real(rt), allocatable :: orbital_eccentricity, orbital_angle

#ifdef AMREX_USE_CUDA
  attributes(managed) :: r_P_initial, r_S_initial, a
  attributes(managed) :: center_P_initial, center_S_initial
  attributes(managed) :: orbital_eccentricity, orbital_angle
#endif



  ! Axis is in orbital plane; we measure angle with respect to this axis. Normally the x axis.

  integer, allocatable :: axis_1

  ! Perpendicular axis in the orbital plane. Normally the y axis.

  integer, allocatable :: axis_2

  ! Perpendicular to both other axes. Normally the z axis and also the rotation axis.

  integer, allocatable :: axis_3

  ! Location of the physical center of the problem, as a fraction of domain size

  real(rt), allocatable :: center_fracx, center_fracy, center_fracz

  ! Lagrange points

  real(rt), allocatable :: L1(:), L2(:), L3(:)

  ! Bulk system motion

  real(rt), allocatable :: bulk_velx, bulk_vely, bulk_velz

  ! Whether we're doing an initialization or a restart

  integer, allocatable :: init

  ! Are we doing a single star simulation?

  logical, allocatable :: single_star

#ifdef AMREX_USE_CUDA
  attributes(managed) :: axis_1, axis_2, axis_3
  attributes(managed) :: center_fracx, center_fracy, center_fracz
  attributes(managed) :: L1, L2, L3
  attributes(managed) :: bulk_velx, bulk_vely, bulk_velz
  attributes(managed) :: init
  attributes(managed) :: single_star
#endif



  ! 1D initial models

  type (initial_model), allocatable :: model_P, model_S
  real(rt), allocatable :: rho_P(:), rho_S(:)
  real(rt), allocatable :: T_P(:), T_S(:)
  real(rt), allocatable :: xn_P(:,:), xn_S(:,:)
  real(rt), allocatable :: r_P(:), r_S(:)

  ! For the grid spacing for our model, we'll use 
  ! 6.25 km. No simulation we do is likely to have a resolution
  ! higher than that inside the stars (it represents
  ! three jumps by a factor of four compared to our 
  ! normal coarse grid resolution). By using 4096
  ! grid points, the size of the 1D domain will be 2.56e9 cm,
  ! which is larger than any reasonable mass white dwarf.

  real(rt), allocatable :: initial_model_dx
  integer,  allocatable :: initial_model_npts

  ! initial_model_mass_tol is tolerance used for getting the total WD mass 
  ! equal to the desired mass. It can be reasonably small, since there
  ! will always be a central density value that can give the desired
  ! WD mass on the grid we use.

  real(rt), allocatable :: initial_model_mass_tol

  ! hse_tol is the tolerance used when iterating over a zone to force
  ! it into HSE by adjusting the current density (and possibly
  ! temperature).  hse_tol should be very small (~ 1.e-10).

  real(rt), allocatable :: initial_model_hse_tol

#ifdef AMREX_USE_CUDA
  attributes(managed) :: model_P, model_S
  attributes(managed) :: rho_P, rho_S
  attributes(managed) :: T_P, T_S
  attributes(managed) :: xn_P, xn_S
  attributes(managed) :: r_P, r_S
  attributes(managed) :: initial_model_dx, initial_model_npts
  attributes(managed) :: initial_model_mass_tol, initial_model_hse_tol
#endif



  ! Composition properties of initial models.
  ! We follow the prescription of Dan et al. 2012 for determining
  ! the composition of the WDs. In this approach, below a certain 
  ! mass there are pure He WDs; just above that are hybrid WDs
  ! with pure CO cores and a He mantle; above that are pure CO WDs
  ! with slightly more oxygen than carbon; and above that are 
  ! ONeMg WDs. All masses are in solar masses.

  real(rt), allocatable :: max_he_wd_mass
  real(rt), allocatable :: max_hybrid_wd_mass
  real(rt), allocatable :: hybrid_wd_he_shell_mass
  real(rt), allocatable :: max_co_wd_mass
  real(rt), allocatable :: co_wd_he_shell_mass

  real(rt), allocatable :: hybrid_wd_c_frac
  real(rt), allocatable :: hybrid_wd_o_frac

  real(rt), allocatable :: co_wd_c_frac
  real(rt), allocatable :: co_wd_o_frac

  real(rt), allocatable :: onemg_wd_o_frac
  real(rt), allocatable :: onemg_wd_ne_frac
  real(rt), allocatable :: onemg_wd_mg_frac

#ifdef AMREX_USE_CUDA
  attributes(managed) :: max_he_wd_mass, max_hybrid_wd_mass, hybrid_wd_he_shell_mass
  attributes(managed) :: max_co_wd_mass, co_wd_he_shell_mass
  attributes(managed) :: hybrid_wd_c_frac, hybrid_wd_o_frac
  attributes(managed) :: co_wd_c_frac, co_wd_o_frac
  attributes(managed) :: onemg_wd_o_frac, onemg_wd_ne_frac, onemg_wd_mg_frac
#endif



  ! Tagging criteria

  integer,  allocatable :: max_stellar_tagging_level
  integer,  allocatable :: max_temperature_tagging_level
  integer,  allocatable :: max_center_tagging_level
  real(rt), allocatable :: stellar_density_threshold
  real(rt), allocatable :: temperature_tagging_threshold
  real(rt), allocatable :: center_tagging_radius
  real(rt), allocatable :: max_tagging_radius
  real(rt), allocatable :: roche_tagging_factor

#ifdef AMREX_USE_CUDA
  attributes(managed) :: max_stellar_tagging_level
  attributes(managed) :: max_temperature_tagging_level
  attributes(managed) :: max_center_tagging_level
  attributes(managed) :: stellar_density_threshold
  attributes(managed) :: temperature_tagging_threshold
  attributes(managed) :: center_tagging_radius
  attributes(managed) :: max_tagging_radius
  attributes(managed) :: roche_tagging_factor
#endif



  ! Stores the center of mass location of the stars throughout the run

  real(rt), allocatable :: com_P(:), com_S(:)
  real(rt), allocatable :: vel_P(:), vel_S(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: com_P, com_S
  attributes(managed) :: vel_P, vel_S
#endif

  ! Stores the effective Roche radii

  real(rt), allocatable :: roche_rad_P, roche_rad_S

#ifdef AMREX_USE_CUDA
  attributes(managed) :: roche_rad_P, roche_rad_S
#endif



  ! Relaxation parameters

  real(rt), allocatable :: relaxation_damping_factor
  real(rt), allocatable :: relaxation_density_cutoff
  real(rt), allocatable :: relaxation_cutoff_time
  integer,  allocatable :: relaxation_is_done

  ! Radial damping parameters

  real(rt), allocatable :: radial_damping_factor
  real(rt), allocatable :: initial_radial_velocity_factor

#ifdef AMREX_USE_CUDA
  attributes(managed) :: relaxation_damping_factor, relaxation_density_cutoff, relaxation_cutoff_time, relaxation_is_done
  attributes(managed) :: radial_damping_factor, initial_radial_velocity_factor
#endif



  ! Distance (in kpc) used for calculation of the gravitational wave amplitude
  ! (this wil be calculated along all three coordinate axes).

  real(rt), allocatable :: gw_dist

#ifdef AMREX_USE_CUDA
  attributes(managed) :: gw_dist
#endif



  ! Current value of the dynamical timescale for each star

  real(rt), allocatable :: t_ff_P, t_ff_S

#ifdef AMREX_USE_CUDA
  attributes(managed) :: t_ff_P, t_ff_S
#endif



  ! Global extrema

  real(rt), allocatable :: T_global_max, rho_global_max, ts_te_global_max

#ifdef AMREX_USE_CUDA
  attributes(managed) :: T_global_max, rho_global_max, ts_te_global_max
#endif



  ! Stores whether we assert that the simulation has completed.

  logical, allocatable :: jobIsDone
  logical, allocatable :: signalJobIsNotDone

#ifdef AMREX_USE_CUDA
  attributes(managed) :: jobIsDone
  attributes(managed) :: signalJobIsNotDone
#endif

  ! Auxiliary data for determining whether the job is done.

  integer, parameter :: num_previous_ener_timesteps = 5
  real(rt) :: total_ener_array(num_previous_ener_timesteps)

end module probdata_module
