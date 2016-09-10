module probdata_module
  
  double precision, save :: rwind0, rwind1, rhowind1, Twind1, rbasefac

  double precision, save :: xmin,xmax,ymin,ymax,zmin,zmax

  integer, save :: npts_model
  integer, parameter :: npts_max = 1000
  double precision, save :: model_r(npts_max), model_rho(npts_max), &
       model_v(npts_max), model_T(npts_max), model_Ye(npts_max), model_Abar(npts_max)

  ! filter is used only when rho or time is below.
  double precision, save :: filter_rhomax, filter_timemax
  
end module probdata_module
