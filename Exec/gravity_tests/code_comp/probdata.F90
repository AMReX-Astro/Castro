module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: heating_factor, g0, rho0, p0, gamma1
  logical, save :: do_pert
  integer, save :: ny

#ifdef AMREX_USE_CUDA
  attributes(managed) :: heating_factor, g0, rho0, p0, gamma1
  attributes(managed) :: do_pert
#endif

end module probdata_module
