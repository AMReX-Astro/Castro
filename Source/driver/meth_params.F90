
! This file is automatically created by parse_castro_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh

! This module stores the runtime parameters and integer names for 
! indexing arrays.
!
! The Fortran-specific parameters are initialized in set_method_params(),
! and the ones that we are mirroring from C++ and obtaining through the
! ParmParse module are initialized in ca_set_castro_method_params().

module meth_params_module

  use bl_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! number of ghost cells for the hyperbolic solver
  integer, parameter     :: NHYP    = 4

  ! NTHERM: number of thermodynamic variables
  integer, save :: NTHERM, NVAR
  integer, save :: URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX
  integer, save :: USHK

  ! QTHERM: number of primitive variables
  integer, save :: QTHERM, QVAR
  integer, save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME
  integer, save :: NQAUX, QGAMC, QC, QDPDR, QDPDE
#ifdef RADIATION
  integer, save :: QGAMCG, QCG, QLAMS
#endif
  integer, save :: QFA, QFS, QFX

  integer, save :: nadv

  ! NQ will be the total number of primitive variables, hydro + radiation
  integer, save :: NQ         

#ifdef RADIATION
  integer, save :: QRADVAR, QRAD, QRADHI, QPTOT, QREITOT
  integer, save :: fspace_type
  logical, save :: do_inelastic_scattering
  logical, save :: comoving

  real(rt)        , save :: flatten_pp_threshold = -1.e0_rt
#endif

  integer, save :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)

  ! These are used for the Godunov state
  ! Note that the velocity indices here are picked to be the same value
  ! as in the primitive variable array
  integer, save :: NGDNV, GDRHO, GDU, GDV, GDW, GDPRES, GDGAME
#ifdef RADIATION
  integer, save :: GDLAMS, GDERADS
#endif

  integer         , save :: numpts_1d

  real(rt)        , save, allocatable :: outflow_data_old(:,:)
  real(rt)        , save, allocatable :: outflow_data_new(:,:)
  real(rt)        , save :: outflow_data_old_time
  real(rt)        , save :: outflow_data_new_time
  logical         , save :: outflow_data_allocated
  real(rt)        , save :: max_dist

  character(len=:), allocatable :: gravity_type

  ! these flags are for interpreting the EXT_DIR BCs
  integer, parameter :: EXT_UNDEFINED = -1
  integer, parameter :: EXT_HSE = 1
  integer, parameter :: EXT_INTERP = 2 
  
  integer, save :: xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext

  ! Create versions of these variables on the GPU
  ! the device update is then done in Castro_nd.f90

  !$acc declare &
  !$acc create(NTHERM, NVAR) &
  !$acc create(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS,UFX) &
  !$acc create(USHK) &
  !$acc create(QTHERM, QVAR) &
  !$acc create(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP) &
  !$acc create(QGAMC, QGAME) &
  !$acc create(NQ) &
#ifdef RADIATION
  !$acc create(QGAMCG, QCG, QLAMS) &
  !$acc create(QRADVAR, QRAD, QRADHI, QPTOT, QREITOT) &
  !$acc create(fspace_type, do_inelastic_scattering, comoving) &
#endif
  !$acc create(QFA, QFS, QFX) &
  !$acc create(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)

  ! Begin the declarations of the ParmParse parameters

  real(rt), save :: difmag
  real(rt), save :: small_dens
  real(rt), save :: small_temp
  real(rt), save :: small_pres
  real(rt), save :: small_ener
  integer         , save :: do_hydro
  integer         , save :: do_ctu
  integer         , save :: hybrid_hydro
  integer         , save :: ppm_type
  integer         , save :: ppm_temp_fix
  integer         , save :: ppm_predict_gammae
  integer         , save :: ppm_reference_eigenvectors
  integer         , save :: plm_iorder
  integer         , save :: hybrid_riemann
  integer         , save :: riemann_solver
  integer         , save :: cg_maxiter
  real(rt), save :: cg_tol
  integer         , save :: cg_blend
  integer         , save :: use_eos_in_riemann
  integer         , save :: use_flattening
  integer         , save :: transverse_use_eos
  integer         , save :: transverse_reset_density
  integer         , save :: transverse_reset_rhoe
  integer         , save :: dual_energy_update_E_from_e
  real(rt), save :: dual_energy_eta1
  real(rt), save :: dual_energy_eta2
  real(rt), save :: dual_energy_eta3
  integer         , save :: use_pslope
  integer         , save :: fix_mass_flux
  integer         , save :: limit_fluxes_on_small_dens
  integer         , save :: density_reset_method
  integer         , save :: allow_negative_energy
  integer         , save :: allow_small_energy
  integer         , save :: do_sponge
  integer         , save :: sponge_implicit
  integer         , save :: first_order_hydro
  character (len=:), allocatable, save :: xl_ext_bc_type
  character (len=:), allocatable, save :: xr_ext_bc_type
  character (len=:), allocatable, save :: yl_ext_bc_type
  character (len=:), allocatable, save :: yr_ext_bc_type
  character (len=:), allocatable, save :: zl_ext_bc_type
  character (len=:), allocatable, save :: zr_ext_bc_type
  integer         , save :: hse_zero_vels
  integer         , save :: hse_interp_temp
  integer         , save :: hse_reflect_vels
  real(rt), save :: cfl
  real(rt), save :: dtnuc_e
  real(rt), save :: dtnuc_X
  real(rt), save :: dtnuc_X_threshold
  integer         , save :: dtnuc_mode
  real(rt), save :: dxnuc
  real(rt), save :: dxnuc_max
  integer         , save :: do_react
  real(rt), save :: react_T_min
  real(rt), save :: react_T_max
  real(rt), save :: react_rho_min
  real(rt), save :: react_rho_max
  integer         , save :: disable_shock_burning
  real(rt), save :: diffuse_cutoff_density
  real(rt), save :: diffuse_cond_scale_fac
  integer         , save :: do_grav
  integer         , save :: grav_source_type
  integer         , save :: do_rotation
  real(rt), save :: rot_period
  real(rt), save :: rot_period_dot
  integer         , save :: rotation_include_centrifugal
  integer         , save :: rotation_include_coriolis
  integer         , save :: rotation_include_domegadt
  integer         , save :: state_in_rotating_frame
  integer         , save :: rot_source_type
  integer         , save :: implicit_rotation_update
  integer         , save :: rot_axis
  real(rt), save :: point_mass
  integer         , save :: point_mass_fix_solution
  integer         , save :: do_acc
  integer         , save :: grown_factor
  integer         , save :: track_grid_losses
  real(rt), save :: const_grav
  integer         , save :: get_g_from_phi

  !$acc declare &
  !$acc create(difmag, small_dens, small_temp) &
  !$acc create(small_pres, small_ener, do_hydro) &
  !$acc create(do_ctu, hybrid_hydro, ppm_type) &
  !$acc create(ppm_temp_fix, ppm_predict_gammae, ppm_reference_eigenvectors) &
  !$acc create(plm_iorder, hybrid_riemann, riemann_solver) &
  !$acc create(cg_maxiter, cg_tol, cg_blend) &
  !$acc create(use_eos_in_riemann, use_flattening, transverse_use_eos) &
  !$acc create(transverse_reset_density, transverse_reset_rhoe, dual_energy_update_E_from_e) &
  !$acc create(dual_energy_eta1, dual_energy_eta2, dual_energy_eta3) &
  !$acc create(use_pslope, fix_mass_flux, limit_fluxes_on_small_dens) &
  !$acc create(density_reset_method, allow_negative_energy, allow_small_energy) &
  !$acc create(do_sponge, sponge_implicit, first_order_hydro) &
  !$acc create(hse_zero_vels, hse_interp_temp, hse_reflect_vels) &
  !$acc create(cfl, dtnuc_e, dtnuc_X) &
  !$acc create(dtnuc_X_threshold, dtnuc_mode, dxnuc) &
  !$acc create(dxnuc_max, do_react, react_T_min) &
  !$acc create(react_T_max, react_rho_min, react_rho_max) &
  !$acc create(disable_shock_burning, diffuse_cutoff_density, diffuse_cond_scale_fac) &
  !$acc create(do_grav, grav_source_type, do_rotation) &
  !$acc create(rot_period, rot_period_dot, rotation_include_centrifugal) &
  !$acc create(rotation_include_coriolis, rotation_include_domegadt, state_in_rotating_frame) &
  !$acc create(rot_source_type, implicit_rotation_update, rot_axis) &
  !$acc create(point_mass, point_mass_fix_solution, do_acc) &
  !$acc create(grown_factor, track_grid_losses, const_grav) &
  !$acc create(get_g_from_phi)

  ! End the declarations of the ParmParse parameters

  real(rt)        , save :: rot_vec(3)

contains

  subroutine ca_set_castro_method_params() bind(C, name="ca_set_castro_method_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, amrex_parmparse_destroy, amrex_parmparse

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (amrex_parmparse) :: pp


    const_grav = 0.0d0;
    get_g_from_phi = 0;

    call amrex_parmparse_build(pp, "gravity")
    call pp%query("const_grav", const_grav)
    call pp%query("get_g_from_phi", get_g_from_phi)
    call amrex_parmparse_destroy(pp)


    difmag = 0.1d0;
    small_dens = -1.d200;
    small_temp = -1.d200;
    small_pres = -1.d200;
    small_ener = -1.d200;
    do_hydro = -1;
    do_ctu = 1;
    hybrid_hydro = 0;
    ppm_type = 1;
    ppm_temp_fix = 0;
    ppm_predict_gammae = 0;
    ppm_reference_eigenvectors = 0;
    plm_iorder = 2;
    hybrid_riemann = 0;
    riemann_solver = 0;
    cg_maxiter = 12;
    cg_tol = 1.0d-5;
    cg_blend = 2;
    use_eos_in_riemann = 0;
    use_flattening = 1;
    transverse_use_eos = 0;
    transverse_reset_density = 1;
    transverse_reset_rhoe = 0;
    dual_energy_update_E_from_e = 1;
    dual_energy_eta1 = 1.0d0;
    dual_energy_eta2 = 1.0d-4;
    dual_energy_eta3 = 1.0d0;
    use_pslope = 1;
    fix_mass_flux = 0;
    limit_fluxes_on_small_dens = 0;
    density_reset_method = 1;
    allow_negative_energy = 0;
    allow_small_energy = 1;
    do_sponge = 0;
    sponge_implicit = 1;
    first_order_hydro = 0;
    allocate(character(len=1)::xl_ext_bc_type)
    xl_ext_bc_type = "";
    allocate(character(len=1)::xr_ext_bc_type)
    xr_ext_bc_type = "";
    allocate(character(len=1)::yl_ext_bc_type)
    yl_ext_bc_type = "";
    allocate(character(len=1)::yr_ext_bc_type)
    yr_ext_bc_type = "";
    allocate(character(len=1)::zl_ext_bc_type)
    zl_ext_bc_type = "";
    allocate(character(len=1)::zr_ext_bc_type)
    zr_ext_bc_type = "";
    hse_zero_vels = 0;
    hse_interp_temp = 0;
    hse_reflect_vels = 0;
    cfl = 0.8d0;
    dtnuc_e = 1.d200;
    dtnuc_X = 1.d200;
    dtnuc_X_threshold = 1.d-3;
    dtnuc_mode = 1;
    dxnuc = 1.d200;
    dxnuc_max = 1.d200;
    do_react = -1;
    react_T_min = 0.0d0;
    react_T_max = 1.d200;
    react_rho_min = 0.0d0;
    react_rho_max = 1.d200;
    disable_shock_burning = 0;
    diffuse_cutoff_density = -1.d200;
    diffuse_cond_scale_fac = 1.0d0;
    do_grav = -1;
    grav_source_type = 4;
    do_rotation = -1;
    rot_period = -1.d200;
    rot_period_dot = 0.0d0;
    rotation_include_centrifugal = 1;
    rotation_include_coriolis = 1;
    rotation_include_domegadt = 1;
    state_in_rotating_frame = 1;
    rot_source_type = 4;
    implicit_rotation_update = 1;
    rot_axis = 3;
    point_mass = 0.0d0;
    point_mass_fix_solution = 0;
    do_acc = -1;
    grown_factor = 1;
    track_grid_losses = 0;

    call amrex_parmparse_build(pp, "castro")
    call pp%query("difmag", difmag)
    call pp%query("small_dens", small_dens)
    call pp%query("small_temp", small_temp)
    call pp%query("small_pres", small_pres)
    call pp%query("small_ener", small_ener)
    call pp%query("do_hydro", do_hydro)
    call pp%query("do_ctu", do_ctu)
    call pp%query("hybrid_hydro", hybrid_hydro)
    call pp%query("ppm_type", ppm_type)
    call pp%query("ppm_temp_fix", ppm_temp_fix)
    call pp%query("ppm_predict_gammae", ppm_predict_gammae)
    call pp%query("ppm_reference_eigenvectors", ppm_reference_eigenvectors)
    call pp%query("plm_iorder", plm_iorder)
    call pp%query("hybrid_riemann", hybrid_riemann)
    call pp%query("riemann_solver", riemann_solver)
    call pp%query("cg_maxiter", cg_maxiter)
    call pp%query("cg_tol", cg_tol)
    call pp%query("cg_blend", cg_blend)
    call pp%query("use_eos_in_riemann", use_eos_in_riemann)
    call pp%query("use_flattening", use_flattening)
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
    call pp%query("sponge_implicit", sponge_implicit)
    call pp%query("first_order_hydro", first_order_hydro)
    call pp%query("xl_ext_bc_type", xl_ext_bc_type)
    call pp%query("xr_ext_bc_type", xr_ext_bc_type)
    call pp%query("yl_ext_bc_type", yl_ext_bc_type)
    call pp%query("yr_ext_bc_type", yr_ext_bc_type)
    call pp%query("zl_ext_bc_type", zl_ext_bc_type)
    call pp%query("zr_ext_bc_type", zr_ext_bc_type)
    call pp%query("hse_zero_vels", hse_zero_vels)
    call pp%query("hse_interp_temp", hse_interp_temp)
    call pp%query("hse_reflect_vels", hse_reflect_vels)
    call pp%query("cfl", cfl)
    call pp%query("dtnuc_e", dtnuc_e)
    call pp%query("dtnuc_X", dtnuc_X)
    call pp%query("dtnuc_X_threshold", dtnuc_X_threshold)
    call pp%query("dtnuc_mode", dtnuc_mode)
    call pp%query("dxnuc", dxnuc)
    call pp%query("dxnuc_max", dxnuc_max)
    call pp%query("do_react", do_react)
    call pp%query("react_T_min", react_T_min)
    call pp%query("react_T_max", react_T_max)
    call pp%query("react_rho_min", react_rho_min)
    call pp%query("react_rho_max", react_rho_max)
    call pp%query("disable_shock_burning", disable_shock_burning)
#ifdef DIFFUSION
    call pp%query("diffuse_cutoff_density", diffuse_cutoff_density)
#endif
#ifdef DIFFUSION
    call pp%query("diffuse_cond_scale_fac", diffuse_cond_scale_fac)
#endif
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
    call pp%query("grown_factor", grown_factor)
    call pp%query("track_grid_losses", track_grid_losses)
    call amrex_parmparse_destroy(pp)



    !$acc update &
    !$acc device(difmag, small_dens, small_temp) &
    !$acc device(small_pres, small_ener, do_hydro) &
    !$acc device(do_ctu, hybrid_hydro, ppm_type) &
    !$acc device(ppm_temp_fix, ppm_predict_gammae, ppm_reference_eigenvectors) &
    !$acc device(plm_iorder, hybrid_riemann, riemann_solver) &
    !$acc device(cg_maxiter, cg_tol, cg_blend) &
    !$acc device(use_eos_in_riemann, use_flattening, transverse_use_eos) &
    !$acc device(transverse_reset_density, transverse_reset_rhoe, dual_energy_update_E_from_e) &
    !$acc device(dual_energy_eta1, dual_energy_eta2, dual_energy_eta3) &
    !$acc device(use_pslope, fix_mass_flux, limit_fluxes_on_small_dens) &
    !$acc device(density_reset_method, allow_negative_energy, allow_small_energy) &
    !$acc device(do_sponge, sponge_implicit, first_order_hydro) &
    !$acc device(hse_zero_vels, hse_interp_temp, hse_reflect_vels) &
    !$acc device(cfl, dtnuc_e, dtnuc_X) &
    !$acc device(dtnuc_X_threshold, dtnuc_mode, dxnuc) &
    !$acc device(dxnuc_max, do_react, react_T_min) &
    !$acc device(react_T_max, react_rho_min, react_rho_max) &
    !$acc device(disable_shock_burning, diffuse_cutoff_density, diffuse_cond_scale_fac) &
    !$acc device(do_grav, grav_source_type, do_rotation) &
    !$acc device(rot_period, rot_period_dot, rotation_include_centrifugal) &
    !$acc device(rotation_include_coriolis, rotation_include_domegadt, state_in_rotating_frame) &
    !$acc device(rot_source_type, implicit_rotation_update, rot_axis) &
    !$acc device(point_mass, point_mass_fix_solution, do_acc) &
    !$acc device(grown_factor, track_grid_losses, const_grav) &
    !$acc device(get_g_from_phi)


    ! now set the external BC flags
    select case (xl_ext_bc_type)
    case ("hse", "HSE")
       xl_ext = EXT_HSE
    case ("interp", "INTERP")       
       xl_ext = EXT_INTERP
    case default
       xl_ext = EXT_UNDEFINED
    end select

    select case (yl_ext_bc_type)
    case ("hse", "HSE")
       yl_ext = EXT_HSE
    case ("interp", "INTERP")       
       yl_ext = EXT_INTERP
    case default
       yl_ext = EXT_UNDEFINED
    end select

    select case (zl_ext_bc_type)
    case ("hse", "HSE")
       zl_ext = EXT_HSE
    case ("interp", "INTERP")       
       zl_ext = EXT_INTERP
    case default
       zl_ext = EXT_UNDEFINED
    end select

    select case (xr_ext_bc_type)
    case ("hse", "HSE")
       xr_ext = EXT_HSE
    case ("interp", "INTERP")       
       xr_ext = EXT_INTERP
    case default
       xr_ext = EXT_UNDEFINED
    end select

    select case (yr_ext_bc_type)
    case ("hse", "HSE")
       yr_ext = EXT_HSE
    case ("interp", "INTERP")       
       yr_ext = EXT_INTERP
    case default
       yr_ext = EXT_UNDEFINED
    end select

    select case (zr_ext_bc_type)
    case ("hse", "HSE")
       zr_ext = EXT_HSE
    case ("interp", "INTERP")       
       zr_ext = EXT_INTERP
    case default
       zr_ext = EXT_UNDEFINED
    end select

    !$acc update device(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)


  end subroutine ca_set_castro_method_params


  subroutine ca_finalize_meth_params() bind(C, name="ca_finalize_meth_params")
    implicit none

    if (allocated(xl_ext_bc_type)) then
        deallocate(xl_ext_bc_type)
    end if
    if (allocated(xr_ext_bc_type)) then
        deallocate(xr_ext_bc_type)
    end if
    if (allocated(yl_ext_bc_type)) then
        deallocate(yl_ext_bc_type)
    end if
    if (allocated(yr_ext_bc_type)) then
        deallocate(yr_ext_bc_type)
    end if
    if (allocated(zl_ext_bc_type)) then
        deallocate(zl_ext_bc_type)
    end if
    if (allocated(zr_ext_bc_type)) then
        deallocate(zr_ext_bc_type)
    end if


    
  end subroutine ca_finalize_meth_params


#ifdef RADIATION
  subroutine ca_init_radhydro_pars(fsp_type_in, do_is_in, com_in,fppt) &
       bind(C, name="ca_init_radhydro_pars")

    use rad_params_module, only : ngroups

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: fsp_type_in, do_is_in, com_in
    real(rt)        , intent(in) :: fppt

    QPTOT  = QVAR+1
    QREITOT = QVAR+2
    QRAD = QVAR+3
    QRADHI = qrad+ngroups-1
  
    QRADVAR = QVAR + 2 + ngroups
  
    ! update NQ -- it was already initialized in the hydro
    NQ = QRADVAR

    ! NQAUX already knows about the hydro and the non-group-dependent
    ! rad variables, update it here
    NQAUX = NQAUX + ngroups

    if (ngroups .eq. 1) then
       fspace_type = 1
    else
       fspace_type = fsp_type_in
    end if
    
    if (fsp_type_in .ne. 1 .and. fsp_type_in .ne. 2) then
       call bl_error("Unknown fspace_type", fspace_type)
    end if
    
    do_inelastic_scattering = (do_is_in .ne. 0)
    
    if (com_in .eq. 1) then
       comoving = .true.
    else if (com_in .eq. 0) then
       comoving = .false.
    else
       call bl_error("Wrong value for comoving", fspace_type)
    end if
    
    flatten_pp_threshold = fppt
    
    !$acc update &
    !$acc device(NQ,NQAUX) &
    !$acc device(QRADVAR, QRAD, QRADHI, QPTOT, QREITOT) &
    !$acc device(fspace_type) &
    !$acc device(do_inelastic_scattering) &
    !$acc device(comoving)
    !$acc device(flatten_pp_threshold = -1.e0_rt)

  end subroutine ca_init_radhydro_pars
#endif

end module meth_params_module
