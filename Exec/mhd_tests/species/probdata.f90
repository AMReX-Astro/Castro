module probdata_module

  use network, only: nspec
  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: xn_zone(nspec)

  real(rt), save :: p, u_x, u_y, u_z, rho, rhoe, frac, T, &
                    B_x, B_y, B_z
  !real(rt), save :: center(3)
  integer, save :: idir

end module probdata_module
