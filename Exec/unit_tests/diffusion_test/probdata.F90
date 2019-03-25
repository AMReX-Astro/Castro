module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), allocatable :: thermal_conductivity, diff_coeff
  real(rt), allocatable :: T1, T2, rho0
  real(rt), allocatable :: t_0

#ifdef AMREX_USE_CUDA
  attributes(managed) :: thermal_conductivity, diff_coeff
  attributes(managed) :: T1, T2, rho0
  attributes(managed) :: t_0
#endif

end module probdata_module
