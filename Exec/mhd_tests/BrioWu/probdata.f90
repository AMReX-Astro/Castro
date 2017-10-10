module probdata_module

  use network, only: nspec
  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: p_ambient, dens_ambient, exp_energy, temp_ambient, e_ambient
  real(rt), save :: xn_zone(nspec)
  real(rt), save :: r_init
  integer, save ::  nsub

  !BWU variables (taking from Nyx test)
  real(rt), save :: p_l, u_l, rho_l, p_r, u_r, rho_r, rhoe_l, rhoe_r, frac, &
                    B_x_l, B_x_r, B_y_l, B_y_r, B_z_l, B_z_r
  !real(rt), save :: center(3)

end module probdata_module
