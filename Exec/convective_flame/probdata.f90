module probdata_module
  implicit none

  double precision, save :: pert_factor, dens_base, pres_base
  double precision, save :: x_pert_loc, pert_width
  double precision, save :: cutoff_density

  double precision, save :: thermal_conductivity

  logical,          save :: do_isentropic

  integer,          save :: boundary_type

  logical         , save :: zero_vels
end module probdata_module
