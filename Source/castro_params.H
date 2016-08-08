
// This file is automatically created by parse_castro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

static int state_interp_order;
static int lin_limit_state_interp;
static int state_nghost;
static Real difmag;
static Real small_dens;
static Real small_temp;
static Real small_pres;
static Real small_ener;
static int do_hydro;
static int add_ext_src;
static int hybrid_hydro;
static int ppm_type;
static int ppm_reference;
static int ppm_trace_sources;
static int ppm_temp_fix;
static int ppm_tau_in_tracing;
static int ppm_predict_gammae;
static int ppm_reference_edge_limit;
static int ppm_reference_eigenvectors;
static int plm_iorder;
static int hybrid_riemann;
static int use_colglaz;
static int riemann_solver;
static int cg_maxiter;
static Real cg_tol;
static int cg_blend;
static int use_flattening;
static int ppm_flatten_before_integrals;
static int transverse_use_eos;
static int transverse_reset_density;
static int transverse_reset_rhoe;
static int dual_energy_update_E_from_e;
static Real dual_energy_eta1;
static Real dual_energy_eta2;
static Real dual_energy_eta3;
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
static Real fixed_dt;
static Real initial_dt;
static Real dt_cutoff;
static Real max_dt;
static Real cfl;
static Real init_shrink;
static Real change_max;
static int use_retry;
static Real retry_neg_dens_factor;
static int sdc_iters;
static Real dtnuc_e;
static Real dtnuc_X;
static int dtnuc_mode;
static Real dxnuc;
static int do_react;
static Real react_T_min;
static Real react_T_max;
static Real react_rho_min;
static Real react_rho_max;
static int disable_shock_burning;
static int do_grav;
static int moving_center;
static int grav_source_type;
static int do_rotation;
#ifdef ROTATION
static Real rotational_period;
#endif
#ifdef ROTATION
static Real rotational_dPdt;
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
static Real point_mass;
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
static int print_energy_diagnostics;
static int track_grid_losses;
static int sum_interval;
static Real sum_per;
static int show_center_of_mass;
static int hard_cfl_limit;
static std::string job_name;
