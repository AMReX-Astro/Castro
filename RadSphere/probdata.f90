module probdata_module

  ! These determine the refinement criteria
  double precision, save :: raderr,radgrad
  integer         , save :: max_raderr_lev, max_radgrad_lev

  ! needed for the generic tagging routine -- not used in this problem
  double precision, save :: denerr,dengrad
  double precision, save :: presserr,pressgrad
  double precision, save :: temperr, tempgrad
  double precision, save :: velerr  ,velgrad
  integer         , save :: max_denerr_lev,  max_dengrad_lev
  integer         , save :: max_presserr_lev,max_pressgrad_lev
  integer         , save :: max_temperr_lev, max_tempgrad_lev
  integer         , save :: max_velerr_lev  ,max_velgrad_lev

  ! rad shock parameters
  double precision, save ::  rho_0, T_0, rhoe_0
  
  double precision, save ::  center(3)
  double precision, save :: xmin,xmax,ymin,ymax,zmin,zmax
      
      
end module probdata_module
