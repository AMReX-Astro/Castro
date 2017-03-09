module probdata_module

  use amrex_fort_module, only : rt => c_real
  implicit none

  real(rt)        , save :: pert_frac, pert_delta

  real(rt)        , save :: rho_fuel, T_fuel, T_ash

end module probdata_module
