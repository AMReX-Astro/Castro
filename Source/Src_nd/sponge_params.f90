module sponge_params_module

  implicit none

  double precision, save :: sponge_lower_radius, sponge_upper_radius
  double precision, save :: sponge_lower_density, sponge_upper_density
  double precision, save :: sponge_timescale

  !$omp threadprivate(sponge_lower_radius, sponge_upper_radius)
  !$omp threadprivate(sponge_lower_density, sponge_upper_density)
  !$omp threadprivate(sponge_timescale)
  
end module sponge_params_module
