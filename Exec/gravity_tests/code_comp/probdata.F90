module probdata_module

!     These determine the refinement criteria
      use amrex_fort_module, only : rt => amrex_real

      real(rt), save :: heating_factor, g0, rho0, p0

end module probdata_module
