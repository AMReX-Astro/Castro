module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), save :: dtemp, x_half_max, x_half_width

  real(rt), save :: X_min, cutoff_density

  real(rt), save :: dx_model

  real(rt), save :: T_hi, T_star, T_lo
  real(rt), save :: dens_base
  real(rt), save :: H_star, atm_delta

  character (len=32), save :: fuel1_name, fuel2_name, fuel3_name
  character (len=32), save :: ash1_name,  ash2_name,  ash3_name
  real (rt), save :: fuel1_frac, fuel2_frac, fuel3_frac
  real (rt), save :: ash1_frac, ash2_frac, ash3_frac

  real (rt), save :: low_density_cutoff, smallx
  real (rt), save :: x_refine_distance

  integer, save :: max_hse_tagging_level
  integer, save :: max_base_tagging_level


end module probdata_module
