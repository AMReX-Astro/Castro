module probdata_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , save :: thermal_conductivity, diff_coeff
  real(rt)        , save :: T1, T2, rho0
  real(rt)        , save :: t_0

end module probdata_module
