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

contains

  ! This routine returns the ambient state.
  ! This file is designed so that it can be overriden
  ! by specific problems that modify this routines with
  ! logic specific to that problem for ambient states
  ! that depend on position and/or time.

  subroutine get_ambient_state(ambient_state_out, loc, time)

    use state_sizes_module, only: NVAR

    implicit none

    real(rt), intent(inout) :: ambient_state_out(NVAR)
    real(rt), intent(in   ) :: loc(3), time

    !$gpu

    ambient_state_out(:) = ambient_state(:)

  end subroutine get_ambient_state

end module ambient_module
