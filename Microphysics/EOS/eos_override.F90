module eos_override_module

  implicit none

  public eos_override

contains

  ! This is a user hook to override the details of the EOS state.

  AMREX_DEVICE subroutine eos_override(state)

    !$acc routine seq

    use eos_type_module, only: eos_t

    implicit none

    type (eos_t) :: state

  end subroutine eos_override

end module eos_override_module
