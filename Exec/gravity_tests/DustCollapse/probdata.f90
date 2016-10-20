module probdata_module

  ! hold the state at the top of the initial model for the boundary
  ! conditions
  double precision, save :: r_old, r_old_s
  double precision, save :: rho_0, rho_ambient, r_0, p_0, T_0, T_ambient, smooth_delta
  double precision, allocatable, save :: X_0(:)
  
  double precision, save :: xmin, xmax, ymin, ymax, zmin, zmax

end module probdata_module
