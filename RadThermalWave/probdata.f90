module probdata_module

      double precision Pi
      parameter (Pi=3.1415926535897932384d0)

      double precision, save :: rhocv, T0, Eexp, rexp
      
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

      double precision, save ::  center(3)
      double precision, save :: xmin,xmax,ymin,ymax,zmin,zmax

end module probdata_module
