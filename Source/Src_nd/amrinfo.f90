
! This module contains information about the state of the Amr class in C++.

module amrinfo_module

  implicit none

  integer :: amr_level = 0
  integer :: amr_iteration = 0
  integer :: amr_ncycle = 0

  double precision :: amr_time = 0.0
  double precision :: amr_dt = 0.0
  
end module amrinfo_module
