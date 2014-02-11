module probdata_module

      double precision Pi
      parameter (Pi=3.1415926535897932384d0)

      double precision, save :: denerr,dengrad
      double precision, save :: presserr,pressgrad
      double precision, save :: temperr, tempgrad
      double precision, save :: velerr  ,velgrad
      double precision, save :: raderr  ,radgrad
      double precision, save :: enterr,  entgrad
      double precision, save :: yeerr,  yegrad
      double precision, save :: masserr
      integer         , save :: max_denerr_lev,  max_dengrad_lev
      integer         , save :: max_presserr_lev,max_pressgrad_lev
      integer         , save :: max_temperr_lev, max_tempgrad_lev
      integer         , save :: max_velerr_lev  ,max_velgrad_lev
      integer         , save :: max_raderr_lev  ,max_radgrad_lev
      integer         , save :: max_enterr_lev  ,max_entgrad_lev
      integer         , save :: max_yeerr_lev  ,max_yegrad_lev
      integer         , save :: max_masserr_lev
 
      double precision, save ::  center(3), corner(3)

      integer, save :: init_smoothing

      integer         , save :: npts_model, model_has_neut_data, ngr_model
      integer, save :: model_temp_in_K
      double precision, save, allocatable :: model_rad(:)
      double precision, save, allocatable :: model_rho(:)
      double precision, save, allocatable :: model_vel(:)
      double precision, save, allocatable :: model_pres(:)
      double precision, save, allocatable :: model_temp(:)
      double precision, save, allocatable :: model_energy(:)
      double precision, save, allocatable :: model_entropy(:)
      double precision, save, allocatable :: model_ye(:)
      double precision, save, allocatable :: model_neut(:,:)

      double precision, save ::  rho_bndry, T_bndry, Ye_bndry, v_bndry

end module probdata_module
