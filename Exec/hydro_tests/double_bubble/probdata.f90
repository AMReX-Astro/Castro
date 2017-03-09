module probdata_module

  use network

  use amrex_fort_module, only : rt => amrex_real
  real(rt)        , save :: dens_base, pres_base
  real(rt)        , save :: pert_factor, y_pert_center, pert_width

  logical,          save :: do_isentropic

  integer,          save :: boundary_type

  real(rt)        , save :: xn_model(nspec)

  real(rt)        , save :: ymin, ymax

  real(rt)        , save :: left_bubble_x_center, right_bubble_x_center

  logical         , save :: single

end module probdata_module
