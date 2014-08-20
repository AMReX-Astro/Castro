
! This module stores the runtime parameters and integer names for 
! indexing arrays.
!
! These parameter are initialized in set_method_params().

module meth_params_module

  implicit none

  double precision, save :: difmag        ! used only in consup to weight the divu contribution
  integer         , save :: iorder        ! used only in uslope and uflaten

  integer, parameter     :: NHYP    = 4
  integer, parameter     :: MAXADV  = 2

  ! NTHERM: number of thermodynamic variables
  integer         , save :: NTHERM, NVAR
  integer         , save :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFA, UFS, UFX

  ! QTHERM: number of primitive variables
  integer         , save :: QTHERM, QVAR
  integer         , save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP
  integer         , save :: QGAMC, QGAME
  integer         , save :: QFA, QFS, QFX

  ! These are only used when we use the SGS model.
  integer         , save :: UESGS,QESGS

  integer         , save :: nadv

  double precision, save :: small_dens, small_temp, small_pres  

  integer         , save :: allow_negative_energy

  integer         , save :: ppm_type
  integer         , save :: ppm_reference
  integer         , save :: ppm_trace_grav
  integer         , save :: ppm_temp_fix
  integer         , save :: ppm_tau_in_tracing
  integer         , save :: ppm_predict_gammae
  integer         , save :: ppm_reference_edge_limit
  integer         , save :: ppm_flatten_before_integrals
  integer         , save :: ppm_reference_eigenvectors
  integer         , save :: use_colglaz
  integer         , save :: use_flattening
  integer         , save :: transverse_use_eos
  integer         , save :: transverse_reset_density
  integer         , save :: transverse_reset_rhoe

  integer         , save :: cg_maxiter
  double precision, save :: cg_tol
  integer         , save :: use_pslope
  integer         , save :: grav_source_type
  integer         , save :: do_sponge
  integer         , save :: normalize_species
  integer         , save :: fix_mass_flux

  integer         , save :: numpts_1d

  double precision, save, allocatable :: outflow_data_old(:,:)
  double precision, save, allocatable :: outflow_data_new(:,:)
  double precision, save :: outflow_data_old_time
  double precision, save :: outflow_data_new_time
  logical         , save :: outflow_data_allocated
  double precision, save :: max_dist

  double precision, save :: rot_period
  double precision, save :: const_grav

  integer, save :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)

end module meth_params_module
