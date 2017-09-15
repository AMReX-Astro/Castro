
! This module stores the runtime parameters that define the problem domain.  
! These parameter are initialized in set_problem_params().

module prob_params_module

  use meth_params_module, only : UMX, UMZ

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! boundary condition information
  integer, save :: physbc_lo(3)
  integer, save :: physbc_hi(3)
  integer, save :: Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall

  ! geometry information
  integer         , save :: coord_type
  real(rt)        , save :: center(3), problo(3), probhi(3)

  ! dimension information
  integer         , save :: dim

  ! indices that we use for dimension agnostic routines 
  ! to ensure we don't illegally access non-existent ghost cells
  ! the format is dg(1:dim) = 1, dg(dim+1:3) = 0
  integer         , save :: dg(3)

  ! grid information
  integer         , save              :: max_level
  real(rt)        , save, allocatable :: dx_level(:,:)
  integer         , save, allocatable :: domlo_level(:,:)
  integer         , save, allocatable :: domhi_level(:,:)
  integer         , save, allocatable :: ref_ratio(:,:)
  integer         , save, allocatable :: n_error_buf(:)
  integer         , save, allocatable :: blocking_factor(:)

  integer, parameter :: MAX_MOM_INDEX = 5
  type momflux_t
     ! we want this to be able to use UMX, UMY, and UMZ to index here, but
     ! we can't use those to allocate, since they are not know until runtime.
     ! dynamic allocation might mess with GPUs, so we make this big enough
     ! to definitely contain UMX, UMY, and UMZ, and then check this when
     ! we fill it
     logical :: comp(MAX_MOM_INDEX)
  end type momflux_t

  ! one component for each coordinate direction flux
  type (momflux_t), save :: mom_flux_has_p(3)
  
end module prob_params_module
