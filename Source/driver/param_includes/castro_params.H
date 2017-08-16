
// This file is automatically created by parse_castro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

static int state_interp_order;
static int lin_limit_state_interp;
static int state_nghost;
static int do_reflux;
static int update_sources_after_reflux;
static int use_custom_knapsack_weights;
static amrex::Real difmag;
static amrex::Real small_dens;
static amrex::Real small_temp;
static amrex::Real small_pres;
static amrex::Real small_ener;
static int do_hydro;
static int do_ctu;
static int add_ext_src;
static int hybrid_hydro;
static int ppm_type;
static int ppm_trace_sources;
static int ppm_temp_fix;
static int ppm_predict_gammae;
static int ppm_reference_eigenvectors;
static int plm_iorder;
static int hybrid_riemann;
static int use_colglaz;
static int riemann_solver;
static int cg_maxiter;
static amrex::Real cg_tol;
static int cg_blend;
static int use_eos_in_riemann;
static int use_flattening;
static int transverse_use_eos;
static int transverse_reset_density;
static int transverse_reset_rhoe;
static int dual_energy_update_E_from_e;
static amrex::Real dual_energy_eta1;
static amrex::Real dual_energy_eta2;
static amrex::Real dual_energy_eta3;
static int use_pslope;
static int fix_mass_flux;
static int limit_fluxes_on_small_dens;
static int density_reset_method;
static int allow_negative_energy;
static int allow_small_energy;
static int do_sponge;
static int sponge_implicit;
static int time_center_sponge;
static int update_state_between_sources;
static int keep_sources_until_end;
static int source_term_predictor;
static int first_order_hydro;
static std::string xl_ext_bc_type;
static std::string xr_ext_bc_type;
static std::string yl_ext_bc_type;
static std::string yr_ext_bc_type;
static std::string zl_ext_bc_type;
static std::string zr_ext_bc_type;
static int hse_zero_vels;
static int hse_interp_temp;
static int hse_reflect_vels;
static amrex::Real fixed_dt;
static amrex::Real initial_dt;
static amrex::Real dt_cutoff;
static amrex::Real max_dt;
static amrex::Real cfl;
static amrex::Real init_shrink;
static amrex::Real change_max;
static int use_retry;
static int use_post_step_regrid;
static int retry_max_subcycles;
static amrex::Real retry_neg_dens_factor;
static int sdc_iters;
static amrex::Real dtnuc_e;
static amrex::Real dtnuc_X;
static int dtnuc_mode;
static amrex::Real dxnuc;
static int do_react;
static amrex::Real react_T_min;
static amrex::Real react_T_max;
static amrex::Real react_rho_min;
static amrex::Real react_rho_max;
static int disable_shock_burning;
#ifdef DIFFUSION
static int diffuse_temp;
#endif
#ifdef DIFFUSION
static int diffuse_enth;
#endif
#ifdef DIFFUSION
static int diffuse_spec;
#endif
#ifdef DIFFUSION
static int diffuse_vel;
#endif
#ifdef DIFFUSION
static amrex::Real diffuse_cutoff_density;
#endif
#ifdef DIFFUSION
static amrex::Real diffuse_cond_scale_fac;
#endif
static int do_grav;
static int moving_center;
static int grav_source_type;
static int do_rotation;
#ifdef ROTATION
static amrex::Real rotational_period;
#endif
#ifdef ROTATION
static amrex::Real rotational_dPdt;
#endif
#ifdef ROTATION
static int rotation_include_centrifugal;
#endif
#ifdef ROTATION
static int rotation_include_coriolis;
#endif
#ifdef ROTATION
static int rotation_include_domegadt;
#endif
#ifdef ROTATION
static int state_in_rotating_frame;
#endif
#ifdef ROTATION
static int rot_source_type;
#endif
#ifdef ROTATION
static int implicit_rotation_update;
#endif
#ifdef ROTATION
static int rot_axis;
#endif
#ifdef POINTMASS
static amrex::Real point_mass;
#endif
#ifdef POINTMASS
static int point_mass_fix_solution;
#endif
static int do_acc;
static int bndry_func_thread_safe;
static int grown_factor;
static int star_at_center;
static int do_special_tagging;
static int spherical_star;
static int print_fortran_warnings;
static int print_update_diagnostics;
static int coalesce_update_diagnostics;
static int track_grid_losses;
static int sum_interval;
static amrex::Real sum_per;
static int show_center_of_mass;
static int hard_cfl_limit;
static std::string job_name;
static int output_at_completion;
