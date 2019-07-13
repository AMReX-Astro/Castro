module probdata_module

  use amrex_fort_module, only : rt => amrex_real
  use network, only: nspec

  implicit none


  ! we take 3 options here.
  !   interp_model = 0: no interpolation
  !   interp_model = 1: conservative average from the model file
  !   interp_model = 2: interp cell-average in model to cell centers on grid
  integer, save :: interp_model

  character(len=64), save :: model_file

  real(rt), save :: pert_frac, pert_delta

  real(rt), save :: rho_fuel, T_fuel, e_fuel, p_fuel
  real(rt), save :: xn_fuel(nspec)

  real(rt), save :: rho_ash, T_ash, e_ash
  real(rt), save :: xn_ash(nspec)

  character(len=32), save :: fuel1_name, fuel2_name, fuel3_name, fuel4_name
  real(rt), save :: X_fuel1, X_fuel2, X_fuel3, X_fuel4

  character(len=32), save :: ash1_name, ash2_name, ash3_name, ash4_name
  real(rt), save :: X_ash1, X_ash2, X_ash3, X_ash4

  real(rt), save :: v_inflow, mass_flux

  real(rt), save :: smallx_init

end module probdata_module
