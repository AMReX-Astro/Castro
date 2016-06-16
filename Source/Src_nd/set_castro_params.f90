subroutine set_castro_method_params( &
  difmag_in, small_dens_in, small_temp_in,  &
  small_pres_in, small_ener_in, do_hydro_in,  &
  hybrid_hydro_in, ppm_type_in, ppm_reference_in,  &
  ppm_trace_sources_in, ppm_temp_fix_in, ppm_tau_in_tracing_in,  &
  ppm_predict_gammae_in, ppm_reference_edge_limit_in, ppm_reference_eigenvectors_in,  &
  hybrid_riemann_in, use_colglaz_in, riemann_solver_in,  &
  cg_maxiter_in, cg_tol_in, cg_blend_in,  &
  use_flattening_in, ppm_flatten_before_integrals_in, transverse_use_eos_in,  &
  transverse_reset_density_in, transverse_reset_rhoe_in, dual_energy_update_E_from_e_in,  &
  dual_energy_eta1_in, dual_energy_eta2_in, dual_energy_eta3_in,  &
  use_pslope_in, fix_mass_flux_in, allow_negative_energy_in,  &
  limit_fluxes_on_small_dens_in, allow_small_energy_in, do_sponge_in,  &
  cfl_in, dtnuc_e_in, dtnuc_X_in,  &
  dtnuc_mode_in, dxnuc_in, do_react_in,  &
  react_T_min_in, react_T_max_in, react_rho_min_in,  &
  react_rho_max_in, disable_shock_burning_in, do_grav_in,  &
  grav_source_type_in, do_rotation_in, rot_period_in,  &
  rot_period_dot_in, rotation_include_centrifugal_in, rotation_include_coriolis_in,  &
  rotation_include_domegadt_in, state_in_rotating_frame_in, rot_source_type_in,  &
  implicit_rotation_update_in, rot_axis_in, point_mass_in,  &
  point_mass_fix_solution_in, do_acc_in, track_grid_losses_in) bind(C)

  use meth_params_module
  use network, only : nspec, naux
  use eos_module
  use parallel
  use bl_error_module

  implicit none

  integer :: ioproc

  double precision, intent(in) :: difmag_in
  double precision, intent(in) :: small_dens_in
  double precision, intent(in) :: small_temp_in
  double precision, intent(in) :: small_pres_in
  double precision, intent(in) :: small_ener_in
  integer,          intent(in) :: do_hydro_in
  integer,          intent(in) :: hybrid_hydro_in
  integer,          intent(in) :: ppm_type_in
  integer,          intent(in) :: ppm_reference_in
  integer,          intent(in) :: ppm_trace_sources_in
  integer,          intent(in) :: ppm_temp_fix_in
  integer,          intent(in) :: ppm_tau_in_tracing_in
  integer,          intent(in) :: ppm_predict_gammae_in
  integer,          intent(in) :: ppm_reference_edge_limit_in
  integer,          intent(in) :: ppm_reference_eigenvectors_in
  integer,          intent(in) :: hybrid_riemann_in
  integer,          intent(in) :: use_colglaz_in
  integer,          intent(in) :: riemann_solver_in
  integer,          intent(in) :: cg_maxiter_in
  double precision, intent(in) :: cg_tol_in
  integer,          intent(in) :: cg_blend_in
  integer,          intent(in) :: use_flattening_in
  integer,          intent(in) :: ppm_flatten_before_integrals_in
  integer,          intent(in) :: transverse_use_eos_in
  integer,          intent(in) :: transverse_reset_density_in
  integer,          intent(in) :: transverse_reset_rhoe_in
  integer,          intent(in) :: dual_energy_update_E_from_e_in
  double precision, intent(in) :: dual_energy_eta1_in
  double precision, intent(in) :: dual_energy_eta2_in
  double precision, intent(in) :: dual_energy_eta3_in
  integer,          intent(in) :: use_pslope_in
  integer,          intent(in) :: fix_mass_flux_in
  integer,          intent(in) :: allow_negative_energy_in
  integer,          intent(in) :: limit_fluxes_on_small_dens_in
  integer,          intent(in) :: allow_small_energy_in
  integer,          intent(in) :: do_sponge_in
  double precision, intent(in) :: cfl_in
  double precision, intent(in) :: dtnuc_e_in
  double precision, intent(in) :: dtnuc_X_in
  integer,          intent(in) :: dtnuc_mode_in
  double precision, intent(in) :: dxnuc_in
  integer,          intent(in) :: do_react_in
  double precision, intent(in) :: react_T_min_in
  double precision, intent(in) :: react_T_max_in
  double precision, intent(in) :: react_rho_min_in
  double precision, intent(in) :: react_rho_max_in
  integer,          intent(in) :: disable_shock_burning_in
  integer,          intent(in) :: do_grav_in
  integer,          intent(in) :: grav_source_type_in
  integer,          intent(in) :: do_rotation_in
  double precision, intent(in) :: rot_period_in
  double precision, intent(in) :: rot_period_dot_in
  integer,          intent(in) :: rotation_include_centrifugal_in
  integer,          intent(in) :: rotation_include_coriolis_in
  integer,          intent(in) :: rotation_include_domegadt_in
  integer,          intent(in) :: state_in_rotating_frame_in
  integer,          intent(in) :: rot_source_type_in
  integer,          intent(in) :: implicit_rotation_update_in
  integer,          intent(in) :: rot_axis_in
  double precision, intent(in) :: point_mass_in
  integer,          intent(in) :: point_mass_fix_solution_in
  integer,          intent(in) :: do_acc_in
  integer,          intent(in) :: track_grid_losses_in

  difmag = difmag_in
  small_dens = small_dens_in
  small_temp = small_temp_in
  small_pres = small_pres_in
  small_ener = small_ener_in
  do_hydro = do_hydro_in
  hybrid_hydro = hybrid_hydro_in
  ppm_type = ppm_type_in
  ppm_reference = ppm_reference_in
  ppm_trace_sources = ppm_trace_sources_in
  ppm_temp_fix = ppm_temp_fix_in
  ppm_tau_in_tracing = ppm_tau_in_tracing_in
  ppm_predict_gammae = ppm_predict_gammae_in
  ppm_reference_edge_limit = ppm_reference_edge_limit_in
  ppm_reference_eigenvectors = ppm_reference_eigenvectors_in
  hybrid_riemann = hybrid_riemann_in
  use_colglaz = use_colglaz_in
  riemann_solver = riemann_solver_in
  cg_maxiter = cg_maxiter_in
  cg_tol = cg_tol_in
  cg_blend = cg_blend_in
  use_flattening = use_flattening_in
  ppm_flatten_before_integrals = ppm_flatten_before_integrals_in
  transverse_use_eos = transverse_use_eos_in
  transverse_reset_density = transverse_reset_density_in
  transverse_reset_rhoe = transverse_reset_rhoe_in
  dual_energy_update_E_from_e = dual_energy_update_E_from_e_in .ne. 0
  dual_energy_eta1 = dual_energy_eta1_in
  dual_energy_eta2 = dual_energy_eta2_in
  dual_energy_eta3 = dual_energy_eta3_in
  use_pslope = use_pslope_in
  fix_mass_flux = fix_mass_flux_in
  allow_negative_energy = allow_negative_energy_in
  limit_fluxes_on_small_dens = limit_fluxes_on_small_dens_in
  allow_small_energy = allow_small_energy_in
  do_sponge = do_sponge_in
  cfl = cfl_in
  dtnuc_e = dtnuc_e_in
  dtnuc_X = dtnuc_X_in
  dtnuc_mode = dtnuc_mode_in
  dxnuc = dxnuc_in
  do_react = do_react_in
  react_T_min = react_T_min_in
  react_T_max = react_T_max_in
  react_rho_min = react_rho_min_in
  react_rho_max = react_rho_max_in
  disable_shock_burning = disable_shock_burning_in
  do_grav = do_grav_in
  grav_source_type = grav_source_type_in
  do_rotation = do_rotation_in
  rot_period = rot_period_in
  rot_period_dot = rot_period_dot_in
  rotation_include_centrifugal = rotation_include_centrifugal_in
  rotation_include_coriolis = rotation_include_coriolis_in
  rotation_include_domegadt = rotation_include_domegadt_in
  state_in_rotating_frame = state_in_rotating_frame_in
  rot_source_type = rot_source_type_in
  implicit_rotation_update = implicit_rotation_update_in
  rot_axis = rot_axis_in
  point_mass = point_mass_in
  point_mass_fix_solution = point_mass_fix_solution_in
  do_acc = do_acc_in
  track_grid_losses = track_grid_losses_in

  ! some checks
  call bl_pd_is_ioproc(ioproc)

  if (small_dens_in > 0.d0) then
     small_dens = small_dens_in
  else
     if (ioproc == 1) then
        call bl_warning("Warning:: small_dens has not been set, defaulting to 1.d-200.")
     endif
     small_dens = 1.d-200
  endif

  if (small_temp_in > 0.d0) then
     small_temp = small_temp_in
  else
     if (ioproc == 1) then
        call bl_warning("Warning:: small_temp has not been set, defaulting to 1.d-200.")
     endif
     small_temp = 1.d-200
  endif

  if (small_pres_in > 0.d0) then
     small_pres = small_pres_in
  else
     small_pres = 1.d-200
  endif
  
  if (small_ener_in > 0.d0) then
     small_ener = small_ener_in
  else
     small_ener = 1.d-200
  endif
  
  call eos_init(small_dens=small_dens, small_temp=small_temp)

  ! The EOS might have modified our choices because of its
  ! internal limitations, so let's get small_dens and small_temp
  ! again just to make sure we're consistent with the EOS.
  
  call eos_get_small_dens(small_dens)
  call eos_get_small_temp(small_temp)

end subroutine set_castro_method_params
