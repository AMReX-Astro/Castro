
! This module stores the runtime parameters that define the problem domain.  
! These parameter are initialized in set_problem_params().

module prob_params_module

  implicit none

  ! boundary condition information
  integer         , save, allocatable :: physbc_lo(:)
  integer         , save, allocatable :: physbc_hi(:)
  integer         , save :: Outflow, Symmetry, SlipWall, NoSlipWall

  ! geometry information
  integer         , save :: coord_type

  ! domain physical coordinates
  double precision, save :: xmin, xmax, ymin, ymax, zmin, zmax
  
end module prob_params_module
