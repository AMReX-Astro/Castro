module prob_params_module
    ! This module stores the runtime parameters that define the problem domain.
    ! These parameter are initialized in set_problem_params().

  use meth_params_module, only: UMX, UMZ
  use amrex_fort_module, only: rt => amrex_real
  use probdata_module, only: center ! for backwards compatibility

  implicit none

  ! geometry information
  integer,  allocatable, save :: coord_type
  real(rt), allocatable :: problo(:), probhi(:)

  ! dimension information
  integer, save, allocatable :: dim

end module prob_params_module
