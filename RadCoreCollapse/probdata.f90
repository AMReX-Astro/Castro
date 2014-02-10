module probdata_module

      double precision Pi
      parameter (Pi=3.1415926535897932384d0)

      double precision, save :: rhoamb,pamb,unamb,ut1amb,ut2amb,eamb,Tamb
      double precision, save :: rhocld,denfact,pcld,uncld,ut1cld,ut2cld,ecld,Tcld
      double precision, save :: ashk,rhoshk,rmach,pshk,unshk,ut1shk,ut2shk,eshk,Tshk
      double precision, save :: shockpos,xcloud,ycloud,zcloud,ecent,radius,Tjunk
      double precision, save :: pertmag,ranampl(4:8,4:8,4:8),ranphse(4:8,4:8,4:8)
      double precision, save :: denerr,dengrad
      double precision, save :: presserr,pressgrad
      double precision, save :: temperr, tempgrad
      double precision, save :: velerr  ,velgrad
      double precision, save :: raderr  ,radgrad
      double precision, save :: enterr,  entgrad
      double precision, save :: yeerr,  yegrad
      double precision, save :: masserr
      double precision, save :: time_0
      integer         , save :: shockdir
      integer         , save :: max_denerr_lev,  max_dengrad_lev
      integer         , save :: max_presserr_lev,max_pressgrad_lev
      integer         , save :: max_temperr_lev, max_tempgrad_lev
      integer         , save :: max_velerr_lev  ,max_velgrad_lev
      integer         , save :: max_raderr_lev  ,max_radgrad_lev
      integer         , save :: max_enterr_lev  ,max_entgrad_lev
      integer         , save :: max_yeerr_lev  ,max_yegrad_lev
      integer         , save :: max_masserr_lev
      integer         , save :: do_thermal_wave, do_clouds, do_cloud_shock, do_light_front
      integer         , save :: do_gaussian_pulse, do_rad_sphere, do_divergence_test
      integer         , save :: direction
      integer         , save :: nmodes,probtype
 
      integer         , save :: do_sedov, do_mix
      double precision, save :: xmid,ymid,zmid,sed_rad
      double precision, save :: xmin,xmax,ymin,ymax,zmin,zmax
      integer         , save :: do_linearMGD_tp
      double precision, save :: x0, T0

      double precision, save ::  pref1, pref2

      double precision, save ::  center(3), corner(3)

      integer, save :: init_smoothing

! neutrino test problem
      integer         , save :: do_neutrino_test
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
