module burn_type_module

  use amrex_fort_module, only : rt => amrex_real
  use actual_network, only: nspec, nspec_evolve, naux

  implicit none

  ! A generic structure holding data necessary to do a nuclear burn.

  ! Set the number of independent variables -- this should be
  ! temperature, enuc + the number of species which participate
  ! in the evolution equations.

  integer, parameter :: neqs = 2 + nspec_evolve

  ! Indices of the temperature and energy variables in the work arrays.

  integer, parameter :: net_itemp = nspec_evolve + 1
  integer, parameter :: net_ienuc = nspec_evolve + 2

  type :: burn_t

    real(rt) :: rho
    real(rt) :: T
    real(rt) :: e
    real(rt) :: xn(nspec)
#if naux > 0
    real(rt) :: aux(naux)
#endif

    real(rt) :: cv
    real(rt) :: cp
    real(rt) :: y_e
    real(rt) :: eta
    real(rt) :: cs
    real(rt) :: dx
    real(rt) :: abar
    real(rt) :: zbar

    ! Last temperature we evaluated the EOS at
    real(rt) :: T_old

    ! Temperature derivatives of specific heat
    real(rt) :: dcvdT
    real(rt) :: dcpdT

    ! The following are the actual integration data.
    ! To avoid potential incompatibilities we won't
    ! include the integration vector y itself here.
    ! It can be reconstructed from all of the above
    ! data, particularly xn, e, and T.

    real(rt) :: ydot(neqs)
    real(rt) :: jac(neqs, neqs)

    ! Whether we are self-heating or not.

    logical          :: self_heat

    ! Zone index information.

    integer          :: i
    integer          :: j
    integer          :: k

    ! diagnostics
    integer :: n_rhs
    integer :: n_jac

    ! Integration time.

    real(rt) :: time

    ! Was the burn successful?

    logical :: success

  end type burn_t

contains

  ! Given an eos type, copy the data relevant to the burn type.

  subroutine eos_to_burn(eos_state, burn_state)

    !$acc routine seq

    use eos_type_module, only: eos_t

    implicit none

    type (eos_t)  :: eos_state
    type (burn_t) :: burn_state

    burn_state % rho  = eos_state % rho
    burn_state % T    = eos_state % T
    burn_state % e    = eos_state % e
    burn_state % xn   = eos_state % xn
#if naux > 0
    burn_state % aux  = eos_state % aux
#endif
    burn_state % cv   = eos_state % cv
    burn_state % cp   = eos_state % cp
    burn_state % y_e  = eos_state % y_e
    burn_state % eta  = eos_state % eta
    burn_state % cs   = eos_state % cs
    burn_state % abar = eos_state % abar
    burn_state % zbar = eos_state % zbar

  end subroutine eos_to_burn



  ! Given a burn type, copy the data relevant to the eos type.

  subroutine burn_to_eos(burn_state, eos_state)

    !$acc routine seq

    use eos_type_module, only: eos_t

    implicit none

    type (burn_t) :: burn_state
    type (eos_t)  :: eos_state

    eos_state % rho  = burn_state % rho
    eos_state % T    = burn_state % T
    eos_state % e    = burn_state % e
    eos_state % xn   = burn_state % xn
#if naux > 0
    eos_state % aux  = burn_state % aux
#endif
    eos_state % cv   = burn_state % cv
    eos_state % cp   = burn_state % cp
    eos_state % y_e  = burn_state % y_e
    eos_state % eta  = burn_state % eta
    eos_state % cs   = burn_state % cs
    eos_state % abar = burn_state % abar
    eos_state % zbar = burn_state % zbar

  end subroutine burn_to_eos


  subroutine normalize_abundances_burn(state)

    !$acc routine seq

    use amrex_constants_module, only: ONE
    use extern_probin_module, only: small_x

    implicit none

    type (burn_t), intent(inout) :: state

    state % xn(:) = max(small_x, min(ONE, state % xn(:)))
    state % xn(:) = state % xn(:) / sum(state % xn(:))

  end subroutine normalize_abundances_burn

end module burn_type_module
