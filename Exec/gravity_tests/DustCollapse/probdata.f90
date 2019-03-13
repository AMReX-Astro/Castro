module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: r_old, r_old_s
  real(rt), save :: rho_0, rho_ambient, r_0, p_0, T_0, T_ambient, smooth_delta
  real(rt), allocatable, save :: X_0(:)

  real(rt), save :: center_x, center_y, center_z

  logical,  save :: is_3d_fullstar

end module probdata_module
