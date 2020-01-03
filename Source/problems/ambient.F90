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

  subroutine get_ambient_eos(eos_state)
    ! Return an EOS state corresponding to the ambient_state.

    use eos_type_module, only: eos_t
    use meth_params_module, only: URHO, UTEMP, UEINT, UFS
    use network, only: nspec

    implicit none

    type(eos_t), intent(inout) :: eos_state

    !$gpu

    eos_state%rho = ambient_state(URHO)
    eos_state%T   = ambient_state(UTEMP)
    eos_state%e   = ambient_state(UEINT)
    eos_state%xn  = ambient_state(UFS:UFS+nspec-1) / ambient_state(URHO)

  end subroutine get_ambient_eos

end module ambient_module
