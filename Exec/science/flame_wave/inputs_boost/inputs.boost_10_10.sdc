# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 9900000
stop_time =  3.0

# PROBLEM SIZE & GEOMETRY
geometry.is_periodic = 0       0
geometry.coord_sys   = 1                  # 0 => cart, 1 => RZ  2=>spherical
geometry.prob_lo     = 0       0
geometry.prob_hi     = 1.2288e5     3.072e4
amr.n_cell           = 6144         1536

# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
# 0 = Interior           3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall
# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
castro.lo_bc       =  3   1
castro.hi_bc       =  2   1

castro.yl_ext_bc_type = "hse"
castro.hse_interp_temp = 1
castro.hse_reflect_vels = 0

castro.yr_ext_bc_type = "interp"

# WHICH PHYSICS
castro.do_hydro = 1
castro.do_react = 1
castro.do_rotation = 1
castro.do_grav = 1
castro.do_sponge = 1

castro.small_temp = 1.e6
castro.small_dens = 1.e-5

castro.time_integration_method = 2
castro.sdc_order = 2
castro.ppm_type = 0
castro.plm_iorder = 2
castro.plm_well_balanced = 1

castro.sdc_solve_for_rhoe = 1
castro.sdc_solver_tol_dens = 1.e-10
castro.sdc_solver_tol_spec = 1.e-10
castro.sdc_solver_tol_ener = 1.e-6
castro.sdc_solver_atol = 1.e-10
castro.sdc_solver = 1

castro.grav_source_type = 2

gravity.gravity_type = ConstantGrav
gravity.const_grav   = -1.5e14

castro.rotational_period = 0.0005
castro.rotation_include_centrifugal = 0

castro.diffuse_temp = 1
castro.diffuse_cutoff_density_hi = 5.e4
#castro.diffuse_cutoff_density = 2.e4

castro.diffuse_cond_scale_fac = 10.0

castro.react_rho_min = 1.e3
castro.react_rho_max = 1.5e7

castro.react_T_min = 6.e7


# TIME STEP CONTROL
castro.cfl            = 0.5     # cfl number for hyperbolic system
castro.init_shrink    = 0.1     # scale back initial timestep
castro.change_max     = 1.1     # max time step growth

castro.dtnuc_e = 0.1

# DIAGNOSTICS & VERBOSITY
castro.sum_interval   = 1       # timesteps between computing mass
castro.v              = 1       # verbosity in Castro.cpp
amr.v                 = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed
amr.ref_ratio       = 4 2 2 2 # refinement ratio
amr.regrid_int      = 2 2 2 2 # how often to regrid
amr.blocking_factor = 16       # block factor in grid generation
amr.max_grid_size   = 128
amr.n_error_buf     = 2 2 2 2 # number of buffer cells in error est

# CHECKPOINT FILES
amr.check_file      = flame_wave_chk        # root name of checkpoint file
amr.check_int       = 100        # number of timesteps between checkpoints

# PLOTFILES
amr.plot_file        = flame_wave_plt        # root name of plotfile
amr.plot_per         = 2.e-3      # number of seconds between plotfiles
amr.derive_plot_vars = ALL

amr.small_plot_file        = flame_wave_smallplt        # root name of plotfile
#amr.small_plot_per         = 5.e-5      # number of seconds between plotfiles
amr.plot_int = 10
amr.small_plot_vars = density Temp enuc
amr.derive_small_plot_vars = abar x_velocity y_velocity z_velocity

#PROBIN FILENAME
amr.probin_file = probin.boost_10_10
