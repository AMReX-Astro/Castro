pp.query("difmag", difmag);
pp.query("small_dens", small_dens);
pp.query("small_temp", small_temp);
pp.query("small_pres", small_pres);
pp.query("small_ener", small_ener);
pp.query("small_x", small_x);
pp.query("do_hydro", do_hydro);
pp.query("add_ext_src", add_ext_src);
pp.query("hybrid_hydro", hybrid_hydro);
pp.query("ppm_type", ppm_type);
pp.query("ppm_reference", ppm_reference);
pp.query("ppm_trace_sources", ppm_trace_sources);
pp.query("ppm_temp_fix", ppm_temp_fix);
pp.query("ppm_tau_in_tracing", ppm_tau_in_tracing);
pp.query("ppm_predict_gammae", ppm_predict_gammae);
pp.query("ppm_reference_edge_limit", ppm_reference_edge_limit);
pp.query("ppm_reference_eigenvectors", ppm_reference_eigenvectors);
pp.query("hybrid_riemann", hybrid_riemann);
pp.query("use_colglaz", use_colglaz);
pp.query("riemann_solver", riemann_solver);
pp.query("cg_maxiter", cg_maxiter);
pp.query("cg_tol", cg_tol);
pp.query("cg_blend", cg_blend);
pp.query("use_flattening", use_flattening);
pp.query("ppm_flatten_before_integrals", ppm_flatten_before_integrals);
pp.query("transverse_use_eos", transverse_use_eos);
pp.query("transverse_reset_density", transverse_reset_density);
pp.query("transverse_reset_rhoe", transverse_reset_rhoe);
pp.query("dual_energy_update_E_from_e", dual_energy_update_E_from_e);
pp.query("dual_energy_eta1", dual_energy_eta1);
pp.query("dual_energy_eta2", dual_energy_eta2);
pp.query("dual_energy_eta3", dual_energy_eta3);
pp.query("use_pslope", use_pslope);
pp.query("fix_mass_flux", fix_mass_flux);
pp.query("allow_negative_energy", allow_negative_energy);
pp.query("allow_small_energy", allow_small_energy);
pp.query("do_sponge", do_sponge);
#ifdef SGS
pp.query("sum_turb_src", sum_turb_src);
#endif
pp.query("fixed_dt", fixed_dt);
pp.query("initial_dt", initial_dt);
pp.query("dt_cutoff", dt_cutoff);
pp.query("max_dt", max_dt);
pp.query("cfl", cfl);
pp.query("init_shrink", init_shrink);
pp.query("change_max", change_max);
pp.query("use_retry", use_retry);
pp.query("retry_neg_dens_factor", retry_neg_dens_factor);
pp.query("dtnuc_e", dtnuc_e);
pp.query("dtnuc_X", dtnuc_X);
pp.query("dtnuc_mode", dtnuc_mode);
pp.query("dxnuc", dxnuc);
pp.query("do_react", do_react);
pp.query("react_T_min", react_T_min);
pp.query("react_T_max", react_T_max);
pp.query("react_rho_min", react_rho_min);
pp.query("react_rho_max", react_rho_max);
pp.query("disable_shock_burning", disable_shock_burning);
pp.query("do_grav", do_grav);
pp.query("moving_center", moving_center);
pp.query("grav_source_type", grav_source_type);
pp.query("source_term_predictor", source_term_predictor);
pp.query("do_rotation", do_rotation);
#ifdef ROTATION
pp.query("rotational_period", rotational_period);
#endif
#ifdef ROTATION
pp.query("rotational_dPdt", rotational_dPdt);
#endif
#ifdef ROTATION
pp.query("rotation_include_centrifugal", rotation_include_centrifugal);
#endif
#ifdef ROTATION
pp.query("rotation_include_coriolis", rotation_include_coriolis);
#endif
#ifdef ROTATION
pp.query("rotation_include_domegadt", rotation_include_domegadt);
#endif
#ifdef ROTATION
pp.query("rot_source_type", rot_source_type);
#endif
#ifdef ROTATION
pp.query("implicit_rotation_update", implicit_rotation_update);
#endif
#ifdef ROTATION
pp.query("rot_axis", rot_axis);
#endif
#ifdef POINTMASS
pp.query("point_mass", point_mass);
#endif
#ifdef POINTMASS
pp.query("point_mass_fix_solution", point_mass_fix_solution);
#endif
pp.query("do_acc", do_acc);
pp.query("bndry_func_thread_safe", bndry_func_thread_safe);
pp.query("grown_factor", grown_factor);
pp.query("star_at_center", star_at_center);
pp.query("do_special_tagging", do_special_tagging);
pp.query("spherical_star", spherical_star);
pp.query("print_fortran_warnings", print_fortran_warnings);
pp.query("print_energy_diagnostics", print_energy_diagnostics);
pp.query("track_grid_losses", track_grid_losses);
pp.query("sum_interval", sum_interval);
pp.query("sum_per", sum_per);
pp.query("show_center_of_mass", show_center_of_mass);
pp.query("hard_cfl_limit", hard_cfl_limit);
pp.query("job_name", job_name);
