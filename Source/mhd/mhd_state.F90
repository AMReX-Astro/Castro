module mhd_state_module

  use meth_params_module, only: NVAR
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: UMAGX = NVAR+1
  integer, parameter :: UMAGY = NVAR+2
  integer, parameter :: UMAGZ = NVAR+3

end module mhd_state_module
