module probdata_module

  use amrex_fort_module, only : rt => amrex_real
  real(rt)         , save :: frac

  real(rt)         , save :: split(3)

  ! RT parameters
  real(rt)        , save :: rho_1, rho_2
  real(rt)        , save :: p0_base
  real(rt)        , save :: L_x

end module probdata_module
