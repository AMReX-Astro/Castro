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
static int hybrid_riemann;
static int use_colglaz;
static int riemann_solver;
static int cg_maxiter;
static Real cg_tol;
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
static int allow_negative_energy;
static int do_sponge;
#ifdef SGS
static Real sum_turb_src;
#endif
static Real fixed_dt;
static Real initial_dt;
static Real dt_cutoff;
static Real max_dt;
static Real cfl;
static Real init_shrink;
static Real change_max;
static int use_retry;
static Real dtnuc_e;
static Real dtnuc_X;
static int dtnuc_mode;
static Real dxnuc;
static int do_react;
static Real react_T_min;
static Real react_T_max;
static int do_grav;
static int moving_center;
static int grav_source_type;
static int source_term_predictor;
static int do_rotation;
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
#ifdef POINTMASS
static Real point_mass;
#endif
static int do_acc;
static int bndry_func_thread_safe;
static int grown_factor;
static int star_at_center;
static int do_special_tagging;
static int spherical_star;
static int print_fortran_warnings;
static int print_energy_diagnostics;
static int sum_interval;
static Real sum_per;
static int show_center_of_mass;
static int hard_cfl_limit;
static std::string job_name;
