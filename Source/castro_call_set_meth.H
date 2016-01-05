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
#ifndef POINTMASS
static Real point_mass = 0.0;
#endif
set_castro_method_params
   (difmag, small_dens, 
    small_temp, small_pres, 
    small_ener, hybrid_hydro, 
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
    use_pslope, normalize_species, 
    fix_mass_flux, allow_negative_energy, 
    do_sponge, burning_timestep_factor, 
    react_T_min, react_T_max, 
    do_grav, grav_source_type, 
    do_rotation, rotational_period, 
    rotational_dPdt, rot_source_type, 
    rot_axis, point_mass, 
    do_acc);
