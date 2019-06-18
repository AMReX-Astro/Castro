module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt) :: dtemp
  real(rt), allocatable :: x_half_max, x_half_width

  real(rt), allocatable :: X_min, cutoff_density

  real(rt), save :: dx_model

  real(rt), save :: T_hi, T_star, T_lo
  real(rt), save :: dens_base
  real(rt), save :: H_star, atm_delta

  character (len=32), save :: fuel1_name, fuel2_name, fuel3_name
  character (len=32), save :: ash1_name,  ash2_name,  ash3_name
  real (rt), save :: fuel1_frac, fuel2_frac, fuel3_frac
  real (rt), save :: ash1_frac, ash2_frac, ash3_frac

  real (rt), save :: low_density_cutoff, smallx
  real (rt), allocatable :: x_refine_distance

  integer, allocatable :: max_hse_tagging_level
  integer, allocatable :: max_base_tagging_level

#ifdef AMREX_USE_CUDA
  attributes(managed) :: x_half_max, x_half_width
  attributes(managed) :: X_min, cutoff_density
  attributes(managed) :: x_refine_distance
  attributes(managed) :: max_hse_tagging_level
  attributes(managed) :: max_base_tagging_level
#endif

end module probdata_module
