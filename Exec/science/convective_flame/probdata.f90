module probdata_module

  use network, only : nspec

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  character(len=80), save :: model_name

  real(rt)        , save :: pert_factor, x_pert_loc, pert_width
  real(rt)        , save :: cutoff_density

  real(rt)        , save :: thermal_conductivity

  logical         , save :: zero_vels

  real(rt)        , save :: rho_ambient, T_ambient, e_ambient, xn_ambient(nspec)

  real(rt)        , save :: refine_cutoff_height

  integer, save :: nx_model

  real(rt), save :: T_base, T_star, T_lo
  real(rt), save :: dens_base
  real(rt), save :: H_star, atm_delta

  character (len=32), save :: fuel1_name, fuel2_name, fuel3_name
  character (len=32), save :: ash1_name,  ash2_name,  ash3_name
  real (rt), save :: fuel1_frac, fuel2_frac, fuel3_frac
  real (rt), save :: ash1_frac, ash2_frac, ash3_frac

  real (rt), save :: low_density_cutoff, smallx

  logical, save :: index_base_from_temp


end module probdata_module
