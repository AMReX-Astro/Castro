
! This module stores the runtime parameters that define the problem domain.  
! These parameter are initialized in set_problem_params().

module prob_params_module

  implicit none
  integer         , save, allocatable :: physbc_lo(:)
  integer         , save, allocatable :: physbc_hi(:)
  integer         , save :: Outflow, Symmetry, SlipWall, NoSlipWall
  integer         , save :: coord_type
  double precision, save :: xmin, xmax, ymin, ymax, zmin, zmax
  
end module prob_params_module
