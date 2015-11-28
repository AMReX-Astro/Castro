static Real burning_timestep_factor;
static Real difmag;
static Real small_dens;
static Real small_temp;
static Real small_pres;
static Real small_ener;
static int do_hydro;
static int do_react;
static int do_grav;
static int do_rotation;
static int add_ext_src;
static int do_acc;
#ifdef ROTATION
static Real rotational_period;
#endif
#ifdef ROTATION
static Real rotational_dPdt;
#endif
#ifdef ROTATION
static int rot_source_type;
#endif
#ifdef ROTATION
static int rot_axis;
#endif
static int deterministic;
static int bndry_func_thread_safe;
static int grown_factor;
static int star_at_center;
static int moving_center;
static int normalize_species;
static int fix_mass_flux;
static int allow_negative_energy;
static int do_special_tagging;
static int ppm_type;
static int ppm_reference;
static int ppm_trace_sources;
static int ppm_temp_fix;
static int ppm_tau_in_tracing;
static int ppm_predict_gammae;
static int ppm_reference_edge_limit;
static int ppm_reference_eigenvectors;
static int source_term_predictor;
static int hybrid_riemann;
static int use_colglaz;
static int riemann_solver;
static int cg_maxiter;
static Real cg_tol;
static int hard_cfl_limit;
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
static int grav_source_type;
static int spherical_star;
static int do_sponge;
static int print_fortran_warnings;
static int print_energy_diagnostics;
static int sum_interval;
static int show_center_of_mass;
#ifdef SGS
static Real sum_turb_src;
#endif
static std::string job_name;
