module probdata_module

  use bl_fort_module, only : rt => c_real
  real(rt)        , save :: H_min, cutoff_density

  character(len=80), save :: model_name

  ! Velocity perturbation
  logical         , save :: apply_vel_field

  real(rt)        , save :: velpert_scale, velpert_amplitude, velpert_height_loc
  integer         , save :: num_vortices

  real(rt)        , allocatable, save :: xloc_vortices(:)

  ! lower boundary
  logical         , save :: interp_BC, zero_vels

end module probdata_module
