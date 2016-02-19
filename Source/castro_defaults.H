Real        Castro::difmag = 0.1;
Real        Castro::small_dens = -1.e200;
Real        Castro::small_temp = -1.e200;
Real        Castro::small_pres = -1.e200;
Real        Castro::small_ener = -1.e200;
Real        Castro::small_x = 1.e-30;
int         Castro::do_hydro = -1;
int         Castro::add_ext_src = 0;
int         Castro::hybrid_hydro = 0;
int         Castro::ppm_type = 1;
int         Castro::ppm_reference = 1;
int         Castro::ppm_trace_sources = 0;
int         Castro::ppm_temp_fix = 0;
int         Castro::ppm_tau_in_tracing = 0;
int         Castro::ppm_predict_gammae = 0;
int         Castro::ppm_reference_edge_limit = 1;
int         Castro::ppm_reference_eigenvectors = 0;
int         Castro::hybrid_riemann = 0;
int         Castro::use_colglaz = 0;
int         Castro::riemann_solver = 0;
int         Castro::cg_maxiter = 12;
Real        Castro::cg_tol = 1.0e-5;
int         Castro::cg_blend = 0;
int         Castro::use_flattening = 1;
int         Castro::ppm_flatten_before_integrals = 0;
int         Castro::transverse_use_eos = 0;
int         Castro::transverse_reset_density = 0;
int         Castro::transverse_reset_rhoe = 0;
int         Castro::dual_energy_update_E_from_e = 1;
Real        Castro::dual_energy_eta1 = 1.0e0;
Real        Castro::dual_energy_eta2 = 1.0e-4;
Real        Castro::dual_energy_eta3 = 0.0e0;
int         Castro::use_pslope = 1;
int         Castro::fix_mass_flux = 0;
int         Castro::allow_negative_energy = 1;
int         Castro::do_sponge = 0;
#ifdef SGS
Real        Castro::sum_turb_src = 0.0;
#endif
Real        Castro::fixed_dt = -1.0;
Real        Castro::initial_dt = -1.0;
Real        Castro::dt_cutoff = 0.0;
Real        Castro::max_dt = 1.e200;
Real        Castro::cfl = 0.8;
Real        Castro::init_shrink = 1.0;
Real        Castro::change_max = 1.1;
int         Castro::use_retry = 0;
Real        Castro::dtnuc_e = 1.e200;
Real        Castro::dtnuc_X = 1.e200;
int         Castro::dtnuc_mode = 1;
Real        Castro::dxnuc = 1.e200;
int         Castro::do_react = -1;
Real        Castro::react_T_min = 0.0;
Real        Castro::react_T_max = 1.e200;
int         Castro::do_grav = -1;
int         Castro::moving_center = 0;
int         Castro::grav_source_type = 2;
int         Castro::source_term_predictor = 0;
int         Castro::do_rotation = -1;
#ifdef ROTATION
Real        Castro::rotational_period = -1.e200;
#endif
#ifdef ROTATION
Real        Castro::rotational_dPdt = 0.0;
#endif
#ifdef ROTATION
int         Castro::rot_source_type = 1;
#endif
#ifdef ROTATION
int         Castro::rot_axis = 3;
#endif
#ifdef POINTMASS
Real        Castro::point_mass = 0.0;
#endif
#ifdef POINTMASS
int         Castro::point_mass_fix_solution = 1;
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
int         Castro::sum_interval = -1;
Real        Castro::sum_per = -1.0e0;
int         Castro::show_center_of_mass = 0;
int         Castro::hard_cfl_limit = 1;
std::string Castro::job_name = "";
