module amrinfo_module
  ! This module contains information about the state of the Amr class in C++.

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, allocatable :: amr_level
  integer, allocatable :: amr_iteration
  integer, allocatable :: amr_ncycle

  real(rt), allocatable         :: amr_time
  real(rt), allocatable         :: amr_dt

#ifdef AMREX_USE_CUDA
  attributes(managed) :: amr_level, amr_iteration, amr_ncycle
  attributes(managed) :: amr_time, amr_dt
#endif

end module amrinfo_module
