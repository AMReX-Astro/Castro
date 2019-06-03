module probdata_module
  
  use amrex_fort_module, only : rt => amrex_real
  real(rt)        , save :: rwind0, rwind1, rhowind1, Twind1, rbasefac

  real(rt)        , save :: xmin,xmax,ymin,ymax,zmin,zmax

  integer, save :: npts_model
  integer, parameter :: npts_max = 1000
  real(rt)        , save :: model_r(npts_max), model_rho(npts_max), &
       model_v(npts_max), model_T(npts_max), model_Ye(npts_max), model_Abar(npts_max)

  ! filter is used only when rho or time is below.
  real(rt)        , save :: filter_rhomax, filter_timemax

  character (len=128), save :: model_file

end module probdata_module
