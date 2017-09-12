module probdata_module

  use amrex_fort_module, only : rt => amrex_real
  character(len=80), save :: model_name

  real(rt)        , save :: dtemp, x_half_max, x_half_width

  real(rt)        , save :: X_min, cutoff_density
  ! lower boundary
  logical         , save :: interp_BC, zero_vels

end module probdata_module
