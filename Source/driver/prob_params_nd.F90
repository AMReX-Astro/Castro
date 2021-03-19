module prob_params_module
    ! This module stores the runtime parameters that define the problem domain.
    ! These parameter are initialized in set_problem_params().

  use meth_params_module, only: UMX, UMZ
  use amrex_fort_module, only: rt => amrex_real
  use probdata_module, only: center ! for backwards compatibility

  implicit none

  ! boundary condition information
  integer, allocatable :: physbc_lo(:)
  integer, allocatable :: physbc_hi(:)

  ! geometry information
  integer,  allocatable, save :: coord_type
  real(rt), allocatable :: problo(:), probhi(:)

  ! dimension information
  integer, save, allocatable :: dim

  ! indices that we use for dimension agnostic routines
  ! to ensure we don't illegally access non-existent ghost cells
  ! the format is dg(1:dim) = 1, dg(dim+1:3) = 0
  integer, save, allocatable :: dg(:)

  !$acc declare create(physbc_lo, physbc_hi)
  !$acc declare create(dim)
  !$acc declare create(dg)
  !$acc declare create(coord_type)
  !$acc declare create(problo, probhi)

end module prob_params_module
