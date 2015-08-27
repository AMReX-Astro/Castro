module probdata_module

  double precision, save :: pert_factor, dens_base, pres_base, y_pert_center
  double precision, save :: gravity, pert_width
  
  logical,          save :: do_isentropic

  integer,          save :: boundary_type

  double precision, save :: frac

end module probdata_module
