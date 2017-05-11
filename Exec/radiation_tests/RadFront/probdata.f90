module probdata_module
  
  use bl_fort_module, only : rt => c_real
  real(rt)        , save :: rho_0, T_0, rhoe_0
  
  real(rt)        , save :: xmin,xmax,ymin,ymax,zmin,zmax
  
end module probdata_module
