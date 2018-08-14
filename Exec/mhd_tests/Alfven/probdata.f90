module probdata_module

  use network, only: nspec
  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: xn_zone(nspec)

  !BWU variables (taking from Nyx test)
  real(rt), save :: p_0, u_x, u_y, u_z, rho_0, rhoe_0, frac, &
                    B_x, B_y, B_z, T_0
                   
  !real(rt), save :: center(3)
  integer, save :: idir

end module probdata_module
