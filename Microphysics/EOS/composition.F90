module eos_composition_module

  use eos_type_module, only : eos_t
  use network, only: nspec, aion, zion
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

  ! Given a set of mass fractions, calculate quantities that depend
  ! on the composition like abar and zbar.

  subroutine composition(state)

    use amrex_constants_module, only: ONE
    use network, only: aion_inv, zion

    implicit none

    type (eos_t), intent(inout) :: state

    !$gpu

    ! Calculate abar, the mean nucleon number,
    ! zbar, the mean proton number,
    ! mu, the mean molecular weight,
    ! mu_e, the mean number of nucleons per electron, and
    ! y_e, the electron fraction.

    state % mu_e = ONE / (sum(state % xn(:) * zion(:) * aion_inv(:)))
    state % y_e = ONE / state % mu_e

    state % abar = ONE / (sum(state % xn(:) * aion_inv(:)))
    state % zbar = state % abar / state % mu_e

  end subroutine composition

end module eos_composition_module
