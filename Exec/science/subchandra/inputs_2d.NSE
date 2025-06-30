# ------------------  INPUTS TO MAIN PROGRAM  -------------------
amr.plot_files_output = 1
amr.checkpoint_files_output = 1

max_step = 1000000
stop_time = 10.0

geometry.is_periodic = 0 0
geometry.coord_sys = 1         # r-z coordinates

geometry.prob_lo   =  0.    0.
geometry.prob_hi   =  5.12e9 1.024e10

amr.n_cell         = 1280 2560

amr.max_level      = 1       # maximum level number allowed

castro.lo_bc       =  3 2
castro.hi_bc       =  2 2

# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
# 0 = Interior           3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall
# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<

castro.do_hydro = 1
castro.do_grav  = 1
castro.do_react = 1
castro.do_sponge = 1

castro.react_rho_min = 1.0
castro.react_T_min = 5.e7

castro.disable_shock_burning = 1
castro.ppm_type = 1
castro.ppm_temp_fix = 0

castro.use_flattening = 1

castro.riemann_solver = 1

# Full self-gravity with the Poisson equation
gravity.gravity_type = PoissonGrav

# Multipole expansion includes terms up to r**(-max_multipole_order)
gravity.max_multipole_order = 6

# Tolerance for multigrid solver for phi solves
gravity.abs_tol = 1.e-10

# Use sync solve for gravity after refluxing
#gravity.no_sync = 0

# Disable the use of the lagged composite correction for the potential
gravity.do_composite_phi_correction = 0

castro.sponge_upper_density = 1.e4
castro.sponge_lower_density = 1.e2
castro.sponge_timescale     = 1.e-3

castro.cfl            = 0.2     # cfl number for hyperbolic system
castro.init_shrink    = 0.05     # scale back initial timestep by this factor
castro.change_max     = 1.025    # factor by which dt is allowed to change each timestep
castro.sum_interval   = 0       # timesteps between computing and printing volume averages
castro.update_sources_after_reflux = 0
castro.time_integration_method = 3

castro.use_retry = 1
castro.retry_subcycle_factor = 0.5
castro.max_subcycles = 32

castro.abundance_failure_rho_cutoff = 1.0

#castro.dtnuc_e = 0.25
#castro.dtnuc_X = 0.25

amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2   # how often to regrid
amr.n_error_buf     = 2 2 2 2 # number of buffer cells in error est
amr.grid_eff        = 0.7     # what constitutes an efficient grid

amr.check_file      = subch_chk     # root name of checkpoint file
amr.check_int       = 50     # number of timesteps between checkpoints
amr.plot_file       = subch_plt     # root name of plot file
amr.plot_int        = -1     # number of timesteps between plotfiles
amr.plot_per        = 2.e-3

amr.max_grid_size   = 256      # maximum grid size allowed -- used to control parallelism
amr.blocking_factor = 32       # block factor in grid generation

amr.v               = 1       # control verbosity in Amr.cpp
castro.v            = 1       # control verbosity in Castro.cpp

amr.derive_plot_vars = ALL
castro.store_burn_weights = 1

castro.small_dens   = 1.e-5
castro.small_temp   = 1.e5

# problem parameters

problem.pert_temp_factor = 20.0
problem.pert_rad_factor = 0.5
problem.R_pert = 1.e7

problem.M_WD = 1.1
problem.M_He  = 0.05

problem.delta = 5.e6

problem.temp_core = 1.e7
problem.temp_base = 1.75e8

problem.mixed_co_wd = 1

problem.X_N14 = 0.01

problem.low_density_cutoff = 1.e-4
problem.temp_fluff = 7.5e7

problem.model_min_res = 5.e5

problem.mu_p = -5.0
problem.mu_n = -12.0

# tagging

amr.refinement_indicators = tempgrad denerr temperr

amr.refine.tempgrad.relative_gradient = 2.0
amr.refine.tempgrad.field_name = Temp
amr.refine.tempgrad.max_level = 1

amr.refine.denerr.value_greater = 1.0
amr.refine.denerr.field_name = density
amr.refine.denerr.max_level = 1

amr.refine.temperr.value_greater = 1.e8
amr.refine.temperr.field_name = Temp
amr.refine.temperr.max_level = 3

# Microphysics

network.small_x = 1.e-10
integrator.SMALL_X_SAFE = 1.e-10

integrator.rtol_spec = 1.e-5
integrator.atol_spec = 1.e-5
integrator.rtol_enuc = 1.e-5
integrator.atol_enuc = 1.e-8
integrator.jacobian = 1

integrator.X_reject_buffer = 4.0

# disable jacobian caching in VODE
integrator.use_jacobian_caching = 0

integrator.ode_max_steps = 500000

nse.nse_dx_independent = 1
nse.nse_molar_independent = 0
nse.ase_tol = 0.1
nse.nse_skip_molar = 0