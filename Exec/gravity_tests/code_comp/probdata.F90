module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: heating_factor, g0, rho0, p0, gamma1
  logical, save :: do_pert
  integer, save :: ny

end module probdata_module
