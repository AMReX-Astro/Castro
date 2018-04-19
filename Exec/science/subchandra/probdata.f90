module probdata_module

!     These determine the refinement criteria
      use amrex_fort_module, only : rt => amrex_real
      character (len=80), save  :: model_name

      real(rt), save :: R_pert
      real(rt), save :: pert_temp_factor
      real(rt), save :: pert_rad_factor

end module probdata_module
