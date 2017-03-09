module probdata_module

  ! rad shock parameters
  use amrex_fort_module, only : rt => amrex_real
  real(rt)        , save ::  rho_0, T_0, rhoe_0
  
  real(rt)        , save :: xmin,xmax,ymin,ymax,zmin,zmax
      
end module probdata_module
