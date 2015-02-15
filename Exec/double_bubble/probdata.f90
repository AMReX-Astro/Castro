module probdata_module

  use network

  double precision, save :: dens_base, pres_base
  double precision, save :: pert_factor, y_pert_center, pert_width

  logical,          save :: do_isentropic

  integer,          save :: boundary_type

  double precision, save :: xn_model(nspec)

  double precision, save :: ymin, ymax

  double precision, save :: left_bubble_x_center, right_bubble_x_center

  logical         , save :: single

end module probdata_module
