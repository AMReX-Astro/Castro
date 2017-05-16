
! This module stores the runtime parameters that define the problem domain.  
! These parameter are initialized in set_problem_params().

module prob_rad_params_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none


  !radiation info
   real(rt)       ,save :: bcval_lo(3)
   real(rt)       ,save :: bcval_hi(3)

end module prob_rad_params_module
