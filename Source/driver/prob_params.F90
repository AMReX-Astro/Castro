module prob_params_module
    ! This module stores the runtime parameters that define the problem domain.
    ! These parameter are initialized in set_problem_params().

  use meth_params_module, only: UMX, UMZ
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  ! boundary condition information
  integer, allocatable :: physbc_lo(:)
  integer, allocatable :: physbc_hi(:)
  integer, allocatable :: Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall

  ! geometry information
  integer,  allocatable, save :: coord_type
  real(rt), allocatable :: center(:), problo(:), probhi(:)

  ! dimension information
  integer, save, allocatable :: dim

  ! indices that we use for dimension agnostic routines
  ! to ensure we don't illegally access non-existent ghost cells
  ! the format is dg(1:dim) = 1, dg(dim+1:3) = 0
  integer, save, allocatable :: dg(:)

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
  type (momflux_t), save, allocatable :: mom_flux_has_p(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: physbc_lo, physbc_hi
  attributes(managed) :: Interior, Inflow, Outflow, Symmetry, Slipwall, NoSlipWall
  attributes(managed) :: dim
  attributes(managed) :: dg
  attributes(managed) :: coord_type
  attributes(managed) :: center, problo, probhi
  attributes(managed) :: domlo_level, domhi_level, dx_level
  attributes(managed) :: ref_ratio, n_error_buf, blocking_factor
  attributes(managed) :: mom_flux_has_p
#endif

  !$acc declare create(physbc_lo, physbc_hi)
  !$acc declare create(Interior, Inflow, Outflow, Symmetry, Slipwall, NoSlipWall)
  !$acc declare create(dim)
  !$acc declare create(dg)
  !$acc declare create(coord_type)
  !$acc declare create(center, problo, probhi)
  !$acc declare create(domlo_level, domhi_level, dx_level)
  !$acc declare create(ref_ratio, n_error_buf, blocking_factor)
  !$acc declare create(mom_flux_has_p)

end module prob_params_module
