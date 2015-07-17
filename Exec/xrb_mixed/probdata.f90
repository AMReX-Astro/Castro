module probdata_module

  double precision, save :: H_min, cutoff_density

  character(len=80), save :: model_name

  ! Velocity perturbation
  logical         , save :: apply_vel_field

  double precision, save :: velpert_scale, velpert_amplitude, velpert_height_loc
  integer         , save :: num_vortices

  double precision, allocatable, save :: xloc_vortices(:)

  ! lower boundary
  logical         , save :: interp_BC, zero_vels

  ! sponge
  double precision, save :: sponge_center_density, sponge_start_factor, sponge_kappa
  
end module probdata_module
