module probdata_module

  use network, only : nspec

  use bl_fort_module, only : rt => c_real
  implicit none

  character(len=80), save :: model_name

  real(rt)        , save :: pert_factor, x_pert_loc, pert_width
  real(rt)        , save :: cutoff_density

  real(rt)        , save :: thermal_conductivity

  logical         , save :: zero_vels

  real(rt)        , save :: rho_ambient, T_ambient, e_ambient, xn_ambient(nspec)

  real(rt)        , save :: refine_cutoff_height
  
end module probdata_module
