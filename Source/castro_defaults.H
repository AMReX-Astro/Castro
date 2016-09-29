
// This file is automatically created by parse_castro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

int         Castro::state_interp_order = 1;
int         Castro::lin_limit_state_interp = 0;
int         Castro::state_nghost = 0;
int         Castro::do_reflux = 1;
int         Castro::use_custom_knapsack_weights = 0;
Real        Castro::difmag = 0.1;
Real        Castro::small_dens = -1.e200;
Real        Castro::small_temp = -1.e200;
Real        Castro::small_pres = -1.e200;
Real        Castro::small_ener = -1.e200;
int         Castro::do_hydro = -1;
int         Castro::add_ext_src = 0;
int         Castro::hybrid_hydro = 0;
int         Castro::ppm_type = 1;
int         Castro::ppm_trace_sources = 1;
int         Castro::ppm_temp_fix = 0;
int         Castro::ppm_predict_gammae = 0;
int         Castro::ppm_reference_eigenvectors = 0;
int         Castro::plm_iorder = 2;
int         Castro::hybrid_riemann = 0;
int         Castro::use_colglaz = -1;
int         Castro::riemann_solver = 0;
int         Castro::cg_maxiter = 12;
Real        Castro::cg_tol = 1.0e-5;
int         Castro::cg_blend = 2;
int         Castro::use_flattening = 1;
int         Castro::transverse_use_eos = 0;
int         Castro::transverse_reset_density = 1;
int         Castro::transverse_reset_rhoe = 0;
int         Castro::dual_energy_update_E_from_e = 1;
Real        Castro::dual_energy_eta1 = 1.0e0;
Real        Castro::dual_energy_eta2 = 1.0e-4;
Real        Castro::dual_energy_eta3 = 1.0e0;
int         Castro::use_pslope = 1;
int         Castro::fix_mass_flux = 0;
int         Castro::limit_fluxes_on_small_dens = 0;
int         Castro::density_reset_method = 1;
int         Castro::allow_negative_energy = 0;
int         Castro::allow_small_energy = 1;
int         Castro::do_sponge = 0;
int         Castro::sponge_implicit = 1;
int         Castro::time_center_sponge = 0;
int         Castro::update_state_between_sources = 1;
int         Castro::keep_sources_until_end = 0;
int         Castro::source_term_predictor = 0;
int         Castro::first_order_hydro = 0;
Real        Castro::fixed_dt = -1.0;
Real        Castro::initial_dt = -1.0;
Real        Castro::dt_cutoff = 0.0;
Real        Castro::max_dt = 1.e200;
Real        Castro::cfl = 0.8;
Real        Castro::init_shrink = 1.0;
Real        Castro::change_max = 1.1;
int         Castro::use_retry = 0;
Real        Castro::retry_neg_dens_factor = 1.e-1;
int         Castro::sdc_iters = 2;
Real        Castro::dtnuc_e = 1.e200;
Real        Castro::dtnuc_X = 1.e200;
int         Castro::dtnuc_mode = 1;
Real        Castro::dxnuc = 1.e200;
int         Castro::do_react = -1;
Real        Castro::react_T_min = 0.0;
Real        Castro::react_T_max = 1.e200;
Real        Castro::react_rho_min = 0.0;
Real        Castro::react_rho_max = 1.e200;
int         Castro::disable_shock_burning = 0;
int         Castro::do_grav = -1;
int         Castro::moving_center = 0;
int         Castro::grav_source_type = 4;
int         Castro::do_rotation = -1;
#ifdef ROTATION
Real        Castro::rotational_period = -1.e200;
#endif
#ifdef ROTATION
Real        Castro::rotational_dPdt = 0.0;
#endif
#ifdef ROTATION
int         Castro::rotation_include_centrifugal = 1;
#endif
#ifdef ROTATION
int         Castro::rotation_include_coriolis = 1;
#endif
#ifdef ROTATION
int         Castro::rotation_include_domegadt = 1;
#endif
#ifdef ROTATION
int         Castro::state_in_rotating_frame = 1;
#endif
#ifdef ROTATION
int         Castro::rot_source_type = 4;
#endif
#ifdef ROTATION
int         Castro::implicit_rotation_update = 1;
#endif
#ifdef ROTATION
int         Castro::rot_axis = 3;
#endif
#ifdef POINTMASS
Real        Castro::point_mass = 0.0;
#endif
#ifdef POINTMASS
int         Castro::point_mass_fix_solution = 0;
#endif
int         Castro::do_acc = -1;
int         Castro::bndry_func_thread_safe = 1;
int         Castro::grown_factor = 1;
int         Castro::star_at_center = -1;
int         Castro::do_special_tagging = 0;
int         Castro::spherical_star = 0;
#ifdef DEBUG
int         Castro::print_fortran_warnings = 1;
#else
int         Castro::print_fortran_warnings = 0;
#endif
#ifdef DEBUG
int         Castro::print_energy_diagnostics = 1;
#else
int         Castro::print_energy_diagnostics = 0;
#endif
int         Castro::track_grid_losses = 0;
int         Castro::sum_interval = -1;
Real        Castro::sum_per = -1.0e0;
int         Castro::show_center_of_mass = 0;
int         Castro::hard_cfl_limit = 1;
std::string Castro::job_name = "";
