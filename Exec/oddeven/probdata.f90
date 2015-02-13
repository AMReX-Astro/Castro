module probdata_module

!     These determine the refinement criteria
      double precision, save ::    denerr,  dengrad
      double precision, save ::    velerr,  velgrad
      double precision, save ::  presserr,pressgrad
      double precision, save ::   temperr, tempgrad
      double precision, save ::    raderr,  radgrad
      integer         , save ::  max_denerr_lev   ,max_dengrad_lev
      integer         , save ::  max_velerr_lev   ,max_velgrad_lev
      integer         , save ::  max_presserr_lev, max_pressgrad_lev
      integer         , save ::  max_temperr_lev,  max_tempgrad_lev
      integer         , save ::  max_raderr_lev,   max_radgrad_lev

      double precision, save ::  center(3), xmin, xmax, ymin, ymax

! odd even variables
      double precision, save ::  p_ambient, dens_ambient, dens_pert_factor, vel_pert
      
end module probdata_module
