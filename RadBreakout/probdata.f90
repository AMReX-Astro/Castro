module probdata_module
  
  !     These determine the refinement criteria
  double precision, save :: denerr,  dengrad
  double precision, save :: velerr,  velgrad
  double precision, save :: presserr,pressgrad
  double precision, save :: temperr, tempgrad
  double precision, save :: raderr,radgrad
  integer         , save :: max_denerr_lev   ,max_dengrad_lev
  integer         , save :: max_velerr_lev   ,max_velgrad_lev
  integer         , save :: max_presserr_lev, max_pressgrad_lev
  integer         , save :: max_temperr_lev, max_tempgrad_lev
  integer         , save :: max_raderr_lev, max_radgrad_lev
  double precision, save :: rwind0, rwind1, rhowind1, Twind1, rbasefac

  double precision, save :: center(3)
  double precision, save :: xmin,xmax,ymin,ymax,zmin,zmax

  integer, save :: npts_model
  integer, parameter :: npts_max = 1000
  double precision, save :: model_r(npts_max), model_rho(npts_max), &
       model_v(npts_max), model_T(npts_max), model_Ye(npts_max), model_Abar(npts_max)

  ! filter is used only when rho or time is below.
  double precision, save :: filter_rhomax, filter_timemax
  
end module probdata_module
