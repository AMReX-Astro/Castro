module probdata_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , save :: pert_frac, pert_delta

  real(rt)        , save :: rho_fuel, T_fuel

end module probdata_module
