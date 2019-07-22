module probdata_module

  ! make the model name enter through the probin file
  use amrex_fort_module, only : rt => amrex_real
  character (len=80), save  :: model_name

  ! hold the state at the top of the initial model for the boundary
  ! conditions
  real(rt)        , save              :: hse_rho_top, hse_T_top
  real(rt)        , save              :: hse_p_top, hse_eint_top
  real(rt)        , allocatable, save :: hse_X_top(:)

  real(rt)        , save :: xmin, xmax, ymin, ymax, zmin, zmax

  real(rt)        , save :: heating_time, heating_rad, &
                            heating_peak, heating_sigma

  ! the prob_type matches MAESTRO test_basestate
  integer         , save :: prob_type

end module probdata_module
