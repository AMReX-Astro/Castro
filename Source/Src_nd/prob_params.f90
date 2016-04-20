
! This module stores the runtime parameters that define the problem domain.  
! These parameter are initialized in set_problem_params().

module prob_params_module

  implicit none

  ! boundary condition information
  integer, save :: physbc_lo(3)
  integer, save :: physbc_hi(3)
  integer, save :: Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall

  ! geometry information
  integer         , save :: coord_type
  double precision, save :: center(3), problo(3), probhi(3)

  ! dimension information
  integer         , save :: dim

  ! indices that we use for dimension agnostic routines 
  ! to ensure we don't illegally access non-existent ghost cells
  ! the format is dg(1:dim) = 1, dg(dim+1:3) = 0
  integer         , save :: dg(3)

  ! grid information
  integer         , save              :: max_level
  double precision, save, allocatable :: dx_level(:,:)
  integer         , save, allocatable :: domlo_level(:,:)
  integer         , save, allocatable :: domhi_level(:,:)
  
end module prob_params_module
