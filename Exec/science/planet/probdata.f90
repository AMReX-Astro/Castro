module probdata_module

  double precision, save :: cutoff_density

  character(len=80), save :: model_name

  ! Velocity perturbation
  logical         , save :: apply_vel_field, shear_vel_field

  double precision, save :: velpert_scale, velpert_amplitude, velpert_height_loc
  integer         , save :: num_vortices

  double precision, save :: shear_height_loc, shear_amplitude
  double precision, save :: shear_height, shear_width_x,shear_width_y
  double precision, allocatable, save :: xloc_vortices(:)

  ! Domain stuff
  double precision, save :: center(3)

  ! lower boundary
  logical         , save :: interp_BC, zero_vels

end module probdata_module
