module probdata_module

  use bl_fort_module, only : rt => c_real
  real(rt)        , save :: inner_radius = 0.75e0_rt
  real(rt)        , save :: outer_radius = 1.50e0_rt

  real(rt)        , save :: density_maximum_radius

  real(rt)        , save :: torus_width
  real(rt)        , save :: torus_center

  real(rt)        , save :: ambient_density = 1.0e-8_rt
  
end module probdata_module
