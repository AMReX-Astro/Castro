module probdata_module

  double precision, save :: pert_factor, dens_base, pres_base, y_pert_center
  double precision, save :: gravity, pert_width
  
  logical,          save :: do_isentropic

  integer,          save :: boundary_type

  double precision, save :: frac

  double precision, save :: INLET_RHO, INLET_EINT, INLET_TEMP, INLET_CS

  integer,          save :: npts_model
  double precision, save :: model_r(2560)
  double precision, save :: model_rho(2560)

end module probdata_module
