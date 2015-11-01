module probdata_module
  implicit none

  character(len=80), save :: model_name

  double precision, save :: pert_factor, x_pert_loc, pert_width
  double precision, save :: cutoff_density

  double precision, save :: thermal_conductivity

  logical         , save :: zero_vels

end module probdata_module
