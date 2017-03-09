module probdata_module

  use amrex_fort_module, only : rt => c_real
  real(rt)        , save :: pert_factor, dens_base, pres_base, y_pert_center
  real(rt)        , save :: pert_width
  
  logical,          save :: do_isentropic

  integer,          save :: boundary_type

  real(rt)        , save :: frac

end module probdata_module
