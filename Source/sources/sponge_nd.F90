module sponge_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), allocatable, save :: sponge_lower_factor, sponge_upper_factor
  real(rt), allocatable, save :: sponge_lower_radius, sponge_upper_radius
  real(rt), allocatable, save :: sponge_lower_density, sponge_upper_density
  real(rt), allocatable, save :: sponge_lower_pressure, sponge_upper_pressure
  real(rt), allocatable, save :: sponge_target_velocity(:)
  real(rt), allocatable, save :: sponge_timescale

#ifdef AMREX_USE_CUDA
  attributes(managed) :: sponge_lower_factor, sponge_upper_factor
  attributes(managed) :: sponge_lower_radius, sponge_upper_radius
  attributes(managed) :: sponge_lower_density, sponge_upper_density
  attributes(managed) :: sponge_lower_pressure, sponge_upper_pressure
  attributes(managed) :: sponge_target_velocity
  attributes(managed) :: sponge_timescale
#endif

  public

end module sponge_module
