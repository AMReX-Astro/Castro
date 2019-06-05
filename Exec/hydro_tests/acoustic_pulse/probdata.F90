module probdata_module

  use network, only: nspec
  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: xn_zone(nspec)
  real(rt), save :: rho0, drho0

end module probdata_module
