! This module stores the runtime parameters and integer names for 
! indexing arrays.
!
! The Fortran-specific parameters are initialized in set_method_params(),
! and the ones that we are mirroring from C++ and obtaining through the
! ParmParse module are initialized in set_castro_method_params().

module meth_params_module

  implicit none

  integer         , save :: iorder        ! used only in uslope 

  ! number of ghost cells for the hyperbolic solver
  integer, parameter     :: NHYP    = 4

  ! NTHERM: number of thermodynamic variables
  integer         , save :: NTHERM, NVAR
  integer         , save :: URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX
  integer         , save :: USHK

  ! QTHERM: number of primitive variables
  integer         , save :: QTHERM, QVAR
  integer         , save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP
  integer         , save :: QGAMC, QGAME
  integer         , save :: QFA, QFS, QFX

  ! These are only used when we use the SGS model.
  integer         , save :: UESGS,QESGS

  integer         , save :: nadv

  integer, save :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)

  ! These are used for the Godunov state
  ! Note that the velocity indices here are picked to be the same value
  ! as in the primitive variable array
  integer, save :: ngdnv, GDRHO, GDU, GDV, GDW, GDPRES, GDGAME, GDLAMS, GDERADS

  integer         , save :: numpts_1d

  double precision, save, allocatable :: outflow_data_old(:,:)
  double precision, save, allocatable :: outflow_data_new(:,:)
  double precision, save :: outflow_data_old_time
  double precision, save :: outflow_data_new_time
  logical         , save :: outflow_data_allocated
  double precision, save :: max_dist

  double precision, save :: diffuse_cutoff_density

  double precision, save :: const_grav

  logical, save :: get_g_from_phi
  
  character(len=:), allocatable :: gravity_type
  
  ! Begin the declarations of the ParmParse parameters

  double precision, save :: difmag
  double precision, save :: small_dens
  double precision, save :: small_temp
  double precision, save :: small_pres
  double precision, save :: small_ener
  integer         , save :: do_hydro
  integer         , save :: hybrid_hydro
  integer         , save :: ppm_type
  integer         , save :: ppm_reference
  integer         , save :: ppm_trace_sources
  integer         , save :: ppm_temp_fix
  integer         , save :: ppm_tau_in_tracing
  integer         , save :: ppm_predict_gammae
  integer         , save :: ppm_reference_edge_limit
  integer         , save :: ppm_reference_eigenvectors
  integer         , save :: hybrid_riemann
  integer         , save :: use_colglaz
  integer         , save :: riemann_solver
  integer         , save :: cg_maxiter
  double precision, save :: cg_tol
  integer         , save :: cg_blend
  integer         , save :: use_flattening
  integer         , save :: ppm_flatten_before_integrals
  integer         , save :: transverse_use_eos
  integer         , save :: transverse_reset_density
  integer         , save :: transverse_reset_rhoe
  logical         , save :: dual_energy_update_E_from_e
  double precision, save :: dual_energy_eta1
  double precision, save :: dual_energy_eta2
  double precision, save :: dual_energy_eta3
  integer         , save :: use_pslope
  integer         , save :: fix_mass_flux
  integer         , save :: limit_fluxes_on_small_dens
  integer         , save :: density_reset_method
  integer         , save :: allow_negative_energy
  integer         , save :: allow_small_energy
  integer         , save :: do_sponge
  double precision, save :: cfl
  double precision, save :: dtnuc_e
  double precision, save :: dtnuc_X
  integer         , save :: dtnuc_mode
  double precision, save :: dxnuc
  integer         , save :: do_react
  double precision, save :: react_T_min
  double precision, save :: react_T_max
  double precision, save :: react_rho_min
  double precision, save :: react_rho_max
  integer         , save :: disable_shock_burning
  integer         , save :: do_grav
  integer         , save :: grav_source_type
  integer         , save :: do_rotation
  double precision, save :: rot_period
  double precision, save :: rot_period_dot
  integer         , save :: rotation_include_centrifugal
  integer         , save :: rotation_include_coriolis
  integer         , save :: rotation_include_domegadt
  integer         , save :: state_in_rotating_frame
  integer         , save :: rot_source_type
  integer         , save :: implicit_rotation_update
  integer         , save :: rot_axis
  double precision, save :: point_mass
  integer         , save :: point_mass_fix_solution
  integer         , save :: do_acc
  integer         , save :: track_grid_losses

  ! End the declarations of the ParmParse parameters

  double precision, save :: rot_vec(3)

contains

  subroutine set_castro_method_params() bind(C,name="set_castro_method_params")

    use parmparse_module, only: parmparse_build, parmparse_destroy, ParmParse

    implicit none

    type (ParmParse) :: pp

    call parmparse_build(pp, "castro")

    difmag = 0.1;
    small_dens = -1.e200;
    small_temp = -1.e200;
    small_pres = -1.e200;
    small_ener = -1.e200;
    do_hydro = -1;
    hybrid_hydro = 0;
    ppm_type = 1;
    ppm_reference = 1;
    ppm_trace_sources = 0;
    ppm_temp_fix = 0;
    ppm_tau_in_tracing = 0;
    ppm_predict_gammae = 0;
    ppm_reference_edge_limit = 1;
    ppm_reference_eigenvectors = 0;
    hybrid_riemann = 0;
    use_colglaz = 0;
    riemann_solver = 0;
    cg_maxiter = 12;
    cg_tol = 1.0e-5;
    cg_blend = 0;
    use_flattening = 1;
    ppm_flatten_before_integrals = 1;
    transverse_use_eos = 0;
    transverse_reset_density = 1;
    transverse_reset_rhoe = 0;
    dual_energy_update_E_from_e = 1;
    dual_energy_eta1 = 1.0e0;
    dual_energy_eta2 = 1.0e-4;
    dual_energy_eta3 = 1.0e0;
    use_pslope = 1;
    fix_mass_flux = 0;
    limit_fluxes_on_small_dens = 0;
    density_reset_method = 1;
    allow_negative_energy = 1;
    allow_small_energy = 1;
    do_sponge = 0;
    cfl = 0.8;
    dtnuc_e = 1.e200;
    dtnuc_X = 1.e200;
    dtnuc_mode = 1;
    dxnuc = 1.e200;
    do_react = -1;
    react_T_min = 0.0;
    react_T_max = 1.e200;
    react_rho_min = 0.0;
    react_rho_max = 1.e200;
    disable_shock_burning = 0;
    do_grav = -1;
    grav_source_type = 2;
    do_rotation = -1;
#ifdef ROTATION
    rot_period = -1.e200;
#endif
#ifdef ROTATION
    rot_period_dot = 0.0;
#endif
#ifdef ROTATION
    rotation_include_centrifugal = 1;
#endif
#ifdef ROTATION
    rotation_include_coriolis = 1;
#endif
#ifdef ROTATION
    rotation_include_domegadt = 1;
#endif
#ifdef ROTATION
    state_in_rotating_frame = 1;
#endif
#ifdef ROTATION
    rot_source_type = 1;
#endif
#ifdef ROTATION
    implicit_rotation_update = 0;
#endif
#ifdef ROTATION
    rot_axis = 3;
#endif
#ifdef POINTMASS
    point_mass = 0.0;
#endif
#ifdef POINTMASS
    point_mass_fix_solution = 1;
#endif
    do_acc = -1;
    track_grid_losses = 0;

    call pp%query("difmag", difmag)
    call pp%query("small_dens", small_dens)
    call pp%query("small_temp", small_temp)
    call pp%query("small_pres", small_pres)
    call pp%query("small_ener", small_ener)
    call pp%query("do_hydro", do_hydro)
    call pp%query("hybrid_hydro", hybrid_hydro)
    call pp%query("ppm_type", ppm_type)
    call pp%query("ppm_reference", ppm_reference)
    call pp%query("ppm_trace_sources", ppm_trace_sources)
    call pp%query("ppm_temp_fix", ppm_temp_fix)
    call pp%query("ppm_tau_in_tracing", ppm_tau_in_tracing)
    call pp%query("ppm_predict_gammae", ppm_predict_gammae)
    call pp%query("ppm_reference_edge_limit", ppm_reference_edge_limit)
    call pp%query("ppm_reference_eigenvectors", ppm_reference_eigenvectors)
    call pp%query("hybrid_riemann", hybrid_riemann)
    call pp%query("use_colglaz", use_colglaz)
    call pp%query("riemann_solver", riemann_solver)
    call pp%query("cg_maxiter", cg_maxiter)
    call pp%query("cg_tol", cg_tol)
    call pp%query("cg_blend", cg_blend)
    call pp%query("use_flattening", use_flattening)
    call pp%query("ppm_flatten_before_integrals", ppm_flatten_before_integrals)
    call pp%query("transverse_use_eos", transverse_use_eos)
    call pp%query("transverse_reset_density", transverse_reset_density)
    call pp%query("transverse_reset_rhoe", transverse_reset_rhoe)
    call pp%query("dual_energy_update_E_from_e", dual_energy_update_E_from_e)
    call pp%query("dual_energy_eta1", dual_energy_eta1)
    call pp%query("dual_energy_eta2", dual_energy_eta2)
    call pp%query("dual_energy_eta3", dual_energy_eta3)
    call pp%query("use_pslope", use_pslope)
    call pp%query("fix_mass_flux", fix_mass_flux)
    call pp%query("limit_fluxes_on_small_dens", limit_fluxes_on_small_dens)
    call pp%query("density_reset_method", density_reset_method)
    call pp%query("allow_negative_energy", allow_negative_energy)
    call pp%query("allow_small_energy", allow_small_energy)
    call pp%query("do_sponge", do_sponge)
    call pp%query("cfl", cfl)
    call pp%query("dtnuc_e", dtnuc_e)
    call pp%query("dtnuc_X", dtnuc_X)
    call pp%query("dtnuc_mode", dtnuc_mode)
    call pp%query("dxnuc", dxnuc)
    call pp%query("do_react", do_react)
    call pp%query("react_T_min", react_T_min)
    call pp%query("react_T_max", react_T_max)
    call pp%query("react_rho_min", react_rho_min)
    call pp%query("react_rho_max", react_rho_max)
    call pp%query("disable_shock_burning", disable_shock_burning)
    call pp%query("do_grav", do_grav)
    call pp%query("grav_source_type", grav_source_type)
    call pp%query("do_rotation", do_rotation)
#ifdef ROTATION
    call pp%query("rotational_period", rot_period)
#endif
#ifdef ROTATION
    call pp%query("rotational_dPdt", rot_period_dot)
#endif
#ifdef ROTATION
    call pp%query("rotation_include_centrifugal", rotation_include_centrifugal)
#endif
#ifdef ROTATION
    call pp%query("rotation_include_coriolis", rotation_include_coriolis)
#endif
#ifdef ROTATION
    call pp%query("rotation_include_domegadt", rotation_include_domegadt)
#endif
#ifdef ROTATION
    call pp%query("state_in_rotating_frame", state_in_rotating_frame)
#endif
#ifdef ROTATION
    call pp%query("rot_source_type", rot_source_type)
#endif
#ifdef ROTATION
    call pp%query("implicit_rotation_update", implicit_rotation_update)
#endif
#ifdef ROTATION
    call pp%query("rot_axis", rot_axis)
#endif
#ifdef POINTMASS
    call pp%query("point_mass", point_mass)
#endif
#ifdef POINTMASS
    call pp%query("point_mass_fix_solution", point_mass_fix_solution)
#endif
    call pp%query("do_acc", do_acc)
    call pp%query("track_grid_losses", track_grid_losses)

    call parmparse_destroy(pp)

  end subroutine set_castro_method_params

end module meth_params_module
