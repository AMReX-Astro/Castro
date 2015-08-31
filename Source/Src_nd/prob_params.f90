
! This module stores the runtime parameters that define the problem domain.  
! These parameter are initialized in set_problem_params().

module prob_params_module

  implicit none

  ! boundary condition information
  integer :: physbc_lo(3)
  integer :: physbc_hi(3)
  integer :: Outflow, Symmetry, SlipWall, NoSlipWall

  ! geometry information
  integer          :: coord_type
  double precision :: center(3), problo(3), probhi(3)

  ! dimension information
  integer          :: dim

  ! indices that we use for dimension agnostic routines 
  ! to ensure we don't illegally access non-existent ghost cells
  integer          :: dg(3)
  
end module prob_params_module
