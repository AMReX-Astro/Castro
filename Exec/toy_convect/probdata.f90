module probdata_module

  character(len=80), save :: model_name

  ! Velocity perturbation
  logical         , save :: apply_vel_field

  double precision, save :: velpert_scale, velpert_amplitude, velpert_height_loc
  integer         , save :: num_vortices

  double precision, allocatable, save :: xloc_vortices(:)

  ! lower boundary
  logical         , save :: interp_BC, zero_vels

end module probdata_module
