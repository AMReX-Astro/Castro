module probdata_module

!     These determine the refinement criteria
      use amrex_fort_module, only : rt => amrex_real
      real(rt)        , save :: denerr,  dengrad
      real(rt)        , save :: velerr,  velgrad
      real(rt)        , save :: presserr,pressgrad
      real(rt)        , save :: temperr,tempgrad
      real(rt)        , save :: raderr,radgrad
      integer         , save :: max_denerr_lev   ,max_dengrad_lev
      integer         , save :: max_velerr_lev   ,max_velgrad_lev
      integer         , save :: max_presserr_lev, max_pressgrad_lev
      integer         , save :: max_temperr_lev, max_tempgrad_lev
      integer         , save :: max_raderr_lev, max_radgrad_lev

!     Sod variables
      real(rt)        , save ::  T_l, T_r, dens, frac, cfrac, w_T, center_T


!     These help specify which specific problem
      integer        , save ::  probtype,idir

      real(rt)        , save ::  center(3)
      real(rt)        , save :: xmin, xmax, ymin, ymax, zmin, zmax
      
      integer, save :: ihe4, ic12, io16
      real(rt)        , save, allocatable :: xn(:)

end module probdata_module
