# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 100000
stop_time =  6.e-5

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 0 0 0
geometry.coord_sys   = 0  # 0 => cart, 1 => RZ  2=>spherical
geometry.prob_lo     = 0        0    0
geometry.prob_hi     = 1.6384e5
amr.n_cell           = @@NZONES@@

# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
# 0 = Interior           3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall
# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
castro.lo_bc       =  3   4   4
castro.hi_bc       =  2   4   4

# WHICH PHYSICS
castro.do_hydro = 1
castro.do_react = 1

castro.ppm_type = 1
castro.ppm_temp_fix = 0

castro.transverse_reset_density = 1

castro.use_flattening = 1

castro.small_temp     = 1.e7

castro.riemann_solver = 0

castro.time_integration_method = @@method@@
castro.sdc_iters = @@SDC_ITERS@@

# TIME STEP CONTROL
castro.cfl            = @@CFL@@     # cfl number for hyperbolic system
castro.init_shrink    = 0.1     # scale back initial timestep
castro.change_max     = 1.05    # scale back initial timestep

castro.dtnuc_e = @@DTNUC_E@@
castro.use_retry = 1
castro.max_subcycles = 32

castro.plot_per_is_exact = 1

# DIAGNOSTICS & VERBOSITY
castro.sum_interval   = 1       # timesteps between computing mass
castro.v              = 1       # verbosity in Castro.cpp
amr.v                 = 1       # verbosity in Amr.cpp
#amr.grid_log        = grdlog  # name of grid logging file

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed
amr.ref_ratio       = 2 2 2 2 # refinement ratio
amr.regrid_int      = 2 2 2 2 # how often to regrid
amr.blocking_factor = 4       # block factor in grid generation
amr.max_grid_size   = 256
amr.n_error_buf     = 2 2 2 2 # number of buffer cells in error est

# CHECKPOINT FILES
amr.check_file      = det_x_chk  # root name of checkpoint file
amr.check_int       = 1000         # number of timesteps between checkpoints

# PLOTFILES
amr.plot_file       = det_x_plt  # root name of plotfile
amr.plot_per = 5.e-6
amr.derive_plot_vars = ALL

#PROBIN FILENAME
amr.probin_file = probin-det-x.nse_disabled
