module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: T_min, T_max, rho_ambient, width
  real(rt), save :: cfrac, ofrac
  real(rt), allocatable, save :: xn(:)

  namelist /fortin/ T_min, T_max, width, rho_ambient, cfrac, ofrac

end module probdata_module
