module probdata_module

  use network, only : nspec

  implicit none

  character(len=80), save :: model_name

  double precision, save :: pert_factor, x_pert_loc, pert_width
  double precision, save :: cutoff_density

  double precision, save :: thermal_conductivity

  logical         , save :: zero_vels

  double precision, save :: rho_ambient, T_ambient, e_ambient, xn_ambient(nspec)

  double precision, save :: refine_cutoff_height
  
end module probdata_module
