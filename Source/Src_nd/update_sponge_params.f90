! This function is designed to permit you to update
! the sponge parameters as a time-dependent process.
! By default it does nothing, meaning that the probin
! parameters are used.

subroutine update_sponge_params(time) bind(C, name="update_sponge_params")

  use sponge_module

  use amrex_fort_module, only : rt => c_real
  implicit none
  
  real(rt)        , intent(in) :: time
  
end subroutine update_sponge_params
