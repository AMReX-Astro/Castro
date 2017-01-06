module probdata_module

  ! odd even variables
  use bl_fort_module, only : rt => c_real
  real(rt)        , save ::  p_ambient, dens_ambient, &
                             dens_pert_factor, vel_pert
      
end module probdata_module
