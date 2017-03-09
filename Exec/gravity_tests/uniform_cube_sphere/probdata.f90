module probdata_module

  use bl_fort_module, only : rt => c_real
  real(rt)        , save :: ambient_dens
  real(rt)        , save :: density, diameter
  real(rt)        , save :: problem

end module probdata_module
