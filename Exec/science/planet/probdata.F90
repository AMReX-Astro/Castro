module probdata_module

  use amrex_fort_module, only : rt => amrex_real
  real(rt)        , save :: cutoff_density

  character(len=80), save :: model_name

  ! Velocity perturbation
  logical         , save :: apply_vel_field, shear_vel_field

  real(rt)        , save :: velpert_scale, velpert_amplitude, velpert_height_loc
  integer         , save :: num_vortices

  real(rt)        , save :: shear_height_loc, shear_amplitude
  real(rt)        , save :: shear_height, shear_width_x,shear_width_y
  real(rt)        , allocatable, save :: xloc_vortices(:)

  ! Domain stuff
  real(rt)        , save :: center(3)

  ! lower boundary
  logical         , save :: interp_BC, zero_vels

end module probdata_module
