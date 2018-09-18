
! This module contains information about the state of the Amr class in C++.

module amrinfo_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, allocatable :: amr_level
  integer, allocatable :: amr_iteration
  integer, allocatable :: amr_ncycle

  real(rt), allocatable         :: amr_time
  real(rt), allocatable         :: amr_dt

end module amrinfo_module
