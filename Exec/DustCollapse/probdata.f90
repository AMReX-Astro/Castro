module probdata_module

!     These determine the refinement criteria
      double precision, save ::   denerr,  dengrad
      double precision, save ::   velerr,  velgrad
      double precision, save :: presserr,pressgrad
      double precision, save ::   raderr,  radgrad
      double precision, save ::  temperr, tempgrad
      integer         , save :: max_denerr_lev   ,max_dengrad_lev
      integer         , save :: max_velerr_lev   ,max_velgrad_lev
      integer         , save :: max_presserr_lev, max_pressgrad_lev
      integer         , save :: max_raderr_lev,   max_radgrad_lev
      integer         , save :: max_temperr_lev,  max_tempgrad_lev

      integer         , save :: npts_model

      logical         , save :: is_3d_fullstar

      ! make the model name enter through the probin file
      character (len=80), save  :: model_name

      ! hold the state at the top of the initial model for the boundary
      ! conditions
      double precision, save :: r_old, r_old_s
      double precision, save :: rho_0, rho_ambient, r_0, p_0, T_0, T_ambient, smooth_delta
      double precision, allocatable, save :: X_0(:)

      double precision, save ::  center(3)

      double precision, save :: xmin, xmax, ymin, ymax, zmin, zmax

      !$omp threadprivate(r_old_s)
end module probdata_module
