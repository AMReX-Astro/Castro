module probdata_module

  use amrex_fort_module, only : rt => c_real
  character(len=80), save :: model_name

  real(rt)        , save :: pert_temp_factor, pert_rad_factor

end module probdata_module
