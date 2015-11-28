#ifndef ROTATION
static Real rotational_period = -1.e200;
#endif
#ifndef ROTATION
static Real rotational_dPdt = 0.0;
#endif
#ifndef ROTATION
static int rot_source_type = 1;
#endif
#ifndef ROTATION
static int rot_axis = 3;
#endif
BL_FORT_PROC_CALL(SET_CASTRO_METHOD_PARAMS, set_castro_method_params)
   (difmag, small_dens, 
    small_temp, small_pres, 
    small_ener, do_grav, 
    do_rotation, do_acc, 
    rotational_period, rotational_dPdt, 
    rot_source_type, rot_axis, 
    deterministic, normalize_species, 
    fix_mass_flux, allow_negative_energy, 
    ppm_type, ppm_reference, 
    ppm_trace_sources, ppm_temp_fix, 
    ppm_tau_in_tracing, ppm_predict_gammae, 
    ppm_reference_edge_limit, ppm_reference_eigenvectors, 
    hybrid_riemann, use_colglaz, 
    riemann_solver, cg_maxiter, 
    cg_tol, use_flattening, 
    ppm_flatten_before_integrals, transverse_use_eos, 
    transverse_reset_density, transverse_reset_rhoe, 
    dual_energy_update_E_from_e, dual_energy_eta1, 
    dual_energy_eta2, dual_energy_eta3, 
    use_pslope, grav_source_type, 
    do_sponge);
