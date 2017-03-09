module probdata_module
  
  use amrex_fort_module, only : rt => amrex_real
  real(rt)        , save :: x0, tau0, Q, Temp0, kapbar, epsilon, p0, p1
  real(rt)        , save :: t0, qn(0:1)

  real(rt)        , save :: xmin,xmax,ymin,ymax,zmin,zmax
  
end module probdata_module
