module probdata_module

  integer         , save :: npts_model

  ! make the model name enter through the probin file
  character (len=80), save  :: model_name

  ! arrange storage for read_in model-- not worrying about efficiency, 
  ! since this will only be called once
  double precision, allocatable, save ::  hse_r(:), hse_rho(:)
  double precision, allocatable, save ::  hse_t(:), hse_p(:)
  double precision, allocatable, save ::  hse_s(:,:)
  
  ! hold the state at the top of the initial model for the boundary
  ! conditions
  double precision, save              :: hse_rho_top, hse_T_top
  double precision, save              :: hse_p_top, hse_eint_top
  double precision, allocatable, save :: hse_X_top(:)
  
  double precision, save :: xmin, xmax, ymin, ymax, zmin, zmax
  
  double precision, save :: heating_time, heating_rad, &
                            heating_peak, heating_sigma

  double precision, save :: sponge_weighting, sponge_start_density, &
                            sponge_width_factor


  ! the prob_type matches MAESTRO test_basestate
  integer         , save :: prob_type

end module probdata_module
