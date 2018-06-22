
! This module contains information about the state of the Amr class in C++.

module amrinfo_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, allocatable :: amr_level

#ifdef CUDA
  attributes(managed) :: amr_level
#endif

  integer :: amr_iteration = 0
  integer :: amr_ncycle = 0

  real(rt)         :: amr_time = 0.0
  real(rt)         :: amr_dt = 0.0
  
end module amrinfo_module
