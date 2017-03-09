module probdata_module

  use amrex_fort_module, only : rt => amrex_real
  real(rt)        , save :: wref_l1, wref_l2

  ! rad shock parameters
  real(rt)        , save ::  rho_0, T_0, kappa_0, x_jump, R

  real(rt)        , save :: xmin,xmax,ymin,ymax,zmin,zmax
      
end module probdata_module
