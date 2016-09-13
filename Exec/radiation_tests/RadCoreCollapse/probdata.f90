module probdata_module

      double precision Pi
      parameter (Pi=3.1415926535897932384d0)

      double precision, save :: yeerr,  yegrad
      double precision, save :: masserr
      integer         , save :: max_yeerr_lev  ,max_yegrad_lev
      integer         , save :: max_masserr_lev
 
      double precision, save :: corner(3)

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
