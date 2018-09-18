module burn_type_module

  use actual_network, only: nspec, nspec_evolve, naux
  use amrex_fort_module, only : rt => amrex_real

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
    ! include the integration array y itself here.
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

  ! Implement a manual copy routine since CUDA Fortran doesn't
  ! (yet) support derived type copying on the device.
  subroutine copy_burn_t(to_state, from_state)

    implicit none

    type (burn_t), intent(in   ) :: from_state
    type (burn_t), intent(  out) :: to_state

    !$gpu

    to_state % rho = from_state % rho
    to_state % T   = from_state % T
    to_state % e   = from_state % e
    to_state % xn(1:nspec) = from_state % xn(1:nspec)

#if naux > 0
    to_state % aux(1:naux) = from_state % aux(1:naux)
#endif

    to_state % cv  = from_state % cv
    to_state % cp  = from_state % cp
    to_state % y_e = from_state % y_e
    to_state % eta = from_state % eta
    to_state % cs  = from_state % cs
    to_state % dx  = from_state % dx

    to_state % abar = from_state % abar
    to_state % zbar = from_state % zbar

    to_state % T_old = from_state % T_old

    to_state % dcvdT = from_state % dcvdT
    to_state % dcpdT = from_state % dcpdT

    to_state % ydot(1:neqs) = from_state % ydot(1:neqs)
    to_state % jac(1:neqs, 1:neqs) = from_state % jac(1:neqs, 1:neqs)

    to_state % self_heat = from_state % self_heat

    to_state % i = from_state % i
    to_state % j = from_state % j
    to_state % k = from_state % k

    to_state % n_rhs = from_state % n_rhs
    to_state % n_jac = from_state % n_jac

    to_state % time = from_state % time

    to_state % success = from_state % success

  end subroutine copy_burn_t

  ! Given an eos type, copy the data relevant to the burn type.

  subroutine eos_to_burn(eos_state, burn_state)

    !$acc routine seq

    use eos_type_module, only: eos_t

    implicit none

    type (eos_t)  :: eos_state
    type (burn_t) :: burn_state

    !$gpu

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

    !$gpu

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

    !$gpu

    state % xn(:) = max(small_x, min(ONE, state % xn(:)))
    state % xn(:) = state % xn(:) / sum(state % xn(:))

  end subroutine normalize_abundances_burn

end module burn_type_module
