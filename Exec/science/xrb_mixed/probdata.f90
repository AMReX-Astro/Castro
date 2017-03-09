module probdata_module

  use amrex_fort_module, only : rt => amrex_real
  real(rt)        , save :: H_min, cutoff_density

  character(len=80), save :: model_name

  ! Velocity perturbation
  logical         , save :: apply_vel_field

  real(rt)        , save :: velpert_scale, velpert_amplitude, velpert_height_loc
  integer         , save :: num_vortices

  real(rt)        , allocatable, save :: xloc_vortices(:)

end module probdata_module
