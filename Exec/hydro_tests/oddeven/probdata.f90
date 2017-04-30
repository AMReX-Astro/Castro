module probdata_module

  ! odd even variables
  use amrex_fort_module, only : rt => amrex_real
  real(rt)        , save ::  p_ambient, dens_ambient, &
                             dens_pert_factor, vel_pert
      
end module probdata_module
