module probdata_module
  use amrex_fort_module, only : rt => amrex_real
 
  ! These determine the refinement criteria
  real(rt), save :: raderr,radgrad
  integer         , save :: max_raderr_lev, max_radgrad_lev

  ! needed for the generic tagging routine -- not used in this problem
  real(rt), save :: denerr,dengrad
  real(rt), save :: presserr,pressgrad
  real(rt), save :: temperr, tempgrad
  real(rt), save :: velerr  ,velgrad
  integer         , save :: max_denerr_lev,  max_dengrad_lev
  integer         , save :: max_presserr_lev,max_pressgrad_lev
  integer         , save :: max_temperr_lev, max_tempgrad_lev
  integer         , save :: max_velerr_lev  ,max_velgrad_lev

  ! rad shock parameters
  real(rt), save ::  rho_0, T_0, rhoe_0, E_rad

  real(rt), save ::  center(3)
  real(rt), save :: xmin,xmax,ymin,ymax,zmin,zmax

end module probdata_module
