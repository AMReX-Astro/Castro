module probdata_module

!     These determine the refinement criteria
      use amrex_fort_module, only : rt => amrex_real

      real(rt), save, allocatable :: heating_factor, g0, rho0, p0, gamma1

#ifdef AMREX_USE_CUDA
      attributes(managed) :: heating_factor, g0, rho0, p0, gamma1
#endif

end module probdata_module
