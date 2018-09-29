module probdata_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , save :: pert_frac, pert_delta

  real(rt)        , save :: rho_fuel, T_fuel, T_ash

  character(len=32), save :: fuel1_name, fuel2_name, fuel3_name, fuel4_name
  real(rt), save :: X_fuel1, X_fuel2, X_fuel3, X_fuel4

  character(len=32), save :: ash1_name, ash2_name, ash3_name, ash4_name
  real(rt), save :: X_ash1, X_ash2, X_ash3, X_ash4

end module probdata_module
