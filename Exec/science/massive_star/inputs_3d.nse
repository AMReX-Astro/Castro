# ------------------  INPUTS TO MAIN PROGRAM  -------------------

amr.plot_files_output = 1
amr.checkpoint_files_output = 1

max_step = 5000000
stop_time = 360000

geometry.is_periodic = 0 0 0
geometry.coord_sys = 0         # r-z coordinates

geometry.prob_lo   =  0.    0.    0.
geometry.prob_hi   =  1.6384e10 1.6384e10 1.6384e10

amr.n_cell         = 512 512 512

amr.max_level      = 2      # maximum level number allowed

castro.lo_bc       =  2 2 2
castro.hi_bc       =  2 2 2

# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
# 0 = Interior           3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall
# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<

# GPU options
castro.hydro_memory_footprint_ratio = 3


castro.do_hydro = 1
castro.do_grav  = 1
castro.do_react = 1
castro.do_sponge = 1

castro.ppm_type = 1
castro.ppm_temp_fix = 0
castro.use_pslope = 1

castro.use_flattening = 1

castro.riemann_solver = 1

gravity.gravity_type = MonopoleGrav
gravity.drdxfac = 2
castro.grav_source_type = 4

castro.sponge_upper_density = 1.e3
castro.sponge_lower_density = 1.e2
castro.sponge_timescale     = 1.e-3

castro.cfl            = 0.5     # cfl number for hyperbolic system
castro.init_shrink    = 0.1     # scale back initial timestep by this factor
castro.change_max     = 1.2    # factor by which dt is allowed to change each timestep
castro.sum_interval   = 10       # timesteps between computing and printing volume averages

#castro.dtnuc_e = 0.25
#castro.dtnuc_X = 0.25

amr.ref_ratio       = 4 4 2 2 2 # refinement ratio
amr.regrid_int      = 10000   # how often to regrid
amr.n_error_buf     = 4 4 2 2 2 # number of buffer cells in error est
amr.grid_eff        = 0.7     # what constitutes an efficient grid

amr.check_file      = massive_star_chk     # root name of checkpoint file
amr.check_int       = 10     # number of timesteps between checkpoints

amr.plot_file       = massive_star_plt     # root name of plot file
amr.plot_per = 5.0
amr.derive_plot_vars = ALL
castro.store_burn_weights = 0

amr.small_plot_file       = massive_star_smallplt     # root name of plot file
amr.small_plot_per = 0.5
amr.small_plot_vars = density Temp in_nse
amr.derive_small_plot_vars = abar Ye enuc MachNumber magvel magvort

fab.format = NATIVE_32

castro.plot_per_is_exact = 0


amr.max_grid_size   = 64       # maximum grid size allowed -- used to control parallelism
amr.blocking_factor = 32       # block factor in grid generation

amr.subcycling_mode = None

amr.v               = 1       # control verbosity in Amr.cpp
castro.v            = 1       # control verbosity in Castro.cpp


castro.small_dens   = 1.0
castro.small_temp   = 1.e6

castro.time_integration_method = 3
castro.use_retry = 1
castro.max_subcycles = 16

# problem initialization

problem.model_name =  "15m_500_sec.aprox19.hse.5.00km"

problem.perturb_model = 1
problem.velpert_amplitude = 5.e6

problem.interpolate_pres = 1

# convection

castro.drive_initial_convection = 1
castro.drive_initial_convection_reinit_period = 2
castro.drive_initial_convection_tmax = 50

# refinement

amr.refinement_indicators = denerr denerr3

amr.refine.denerr.max_level = 1
amr.refine.denerr.value_greater = 2.e3
amr.refine.denerr.field_name = density

amr.refine.denerr3.max_level = 2
amr.refine.denerr3.value_greater = 3.e5
amr.refine.denerr3.field_name = density


# Microphysics

integrator.rtol_spec = 1.e-5
integrator.atol_spec = 1.e-5
integrator.rtol_enuc = 1.e-5
integrator.atol_enuc = 1.e-5

integrator.jacobian = 1

network.rho_nse = 1.e7
network.T_nse = 3.e9
network.Si_nse = 0.02
network.C_nse = 1.0
network.O_nse = 1.0

integrator.ode_max_steps = 5000

integrator.use_burn_retry = 1
integrator.retry_swap_jacobian = 1

integrator.retry_rtol_spec = 1.e-5
integrator.retry_atol_spec = 1.e-5
integrator.retry_rtol_enuc = 1.e-5
integrator.retry_atol_enuc = 1.e-5

network.small_x = 1.e-10

network.nse_relax_factor = 0.9

# use cubic interpolation from the NSE table
network.nse_table_interp_linear = 0

# disable jacobian caching in VODE
integrator.use_jacobian_caching = 0

# do we include weak rate neutrino losses in the energy?
integrator.nse_include_enu_weak = 0
