module ambient_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  ! This is a state vector that contains "ambient" material
  ! that will be used for filling material in regions that
  ! are not of interest, like at the edges of the domain.

  real(rt), allocatable :: ambient_state(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: ambient_state
#endif

end module ambient_module
