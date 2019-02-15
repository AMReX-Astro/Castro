module probdata_module

  use network, only: nspec
  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: xn_zone(nspec)
  real(rt), save :: rho0, T0, dp_fact, L_pert

  real(rt), save :: s0, p0

end module probdata_module
