module probdata_module

  ! Ensman test variables -- we set the defaults here
  use amrex_fort_module, only : rt => amrex_real
  real(rt), save :: rho0 = 5.4588672836e-13_rt
  real(rt), save :: T0 = 100.e0_rt
  real(rt), save :: v0 = 235435.371882e0_rt
  real(rt), save :: rho1 = 1.24793794736e-12_rt
  real(rt), save :: T1 = 207.756999533e0_rt
  real(rt), save :: v1 = 102986.727159e0_rt
  integer, save :: idir = 1   ! Direction integer x, y ,z = 1,2,3

  namelist /fortin/ rho0, T0, v0, rho1, T1, v1, idir

  ! for convenience
  real(rt), save :: xmin, xmax, ymin, ymax, zmin, zmax

  real(rt), save :: eint0, etot0, eint1, etot1

end module probdata_module
