module burn_type_module

  use network, only: nspec, naux
  use actual_burner_data, only: nrates
  use eos_module, only: eos_t

  implicit none

  ! A generic structure holding data necessary to do a nuclear burn.

  ! Initialize the main quantities to an unphysical number
  ! so that we know if the user forgot to initialize them
  ! when calling the burner.

  double precision, parameter :: init_num  = -1.0d200
  double precision, parameter :: init_test = -1.0d199

  type :: burn_t

    double precision :: rho              = init_num
    double precision :: T                = init_num
    double precision :: e                = init_num
    double precision :: h                = init_num
    double precision :: xn(nspec)        = init_num
    double precision :: aux(naux)        = init_num
    double precision :: cv               = init_num
    double precision :: cp               = init_num
    double precision :: dedX(nspec)      = init_num
    double precision :: dhdX(nspec)      = init_num
    double precision :: abar             = init_num
    double precision :: zbar             = init_num
    double precision :: rates(nrates)    = init_num
    double precision :: dratesdT(nrates) = init_num

    ! The following are time derivatives

    double precision :: rhodot           = init_num
    double precision :: Tdot             = init_num
    double precision :: xndot(nspec)     = init_num
    double precision :: auxdot(naux)     = init_num
    double precision :: edot             = init_num
    double precision :: hdot             = init_num

  end type burn_t

contains

  ! Given an eos type, copy the data relevant to the burn type.

  subroutine eos_to_burn(eos_state, burn_state)

    implicit none

    type (eos_t)  :: eos_state
    type (burn_t) :: burn_state

    burn_state % rho  = eos_state % rho
    burn_state % T    = eos_state % T
    burn_state % e    = eos_state % e
    burn_state % xn   = eos_state % xn
    burn_state % aux  = eos_state % aux
    burn_state % cv   = eos_state % cv
    burn_state % cp   = eos_state % cp
    burn_state % dedX = eos_state % dedX
    burn_state % dhdX = eos_state % dhdX
    burn_state % abar = eos_state % abar
    burn_state % zbar = eos_state % zbar

  end subroutine eos_to_burn



  ! Given a burn type, copy the data relevant to the eos type.

  subroutine burn_to_eos(burn_state, eos_state)

    use meth_params_module, only: allow_negative_energy

    implicit none

    type (burn_t) :: burn_state
    type (eos_t)  :: eos_state

    eos_state % rho  = burn_state % rho
    eos_state % T    = burn_state % T
    eos_state % e    = burn_state % e
    eos_state % xn   = burn_state % xn
    eos_state % aux  = burn_state % aux
    eos_state % cv   = burn_state % cv
    eos_state % cp   = burn_state % cp
    eos_state % dedX = burn_state % dedX
    eos_state % dhdX = burn_state % dhdX
    eos_state % abar = burn_state % abar
    eos_state % zbar = burn_state % zbar

    if (allow_negative_energy .eq. 0) eos_state % reset = .true.

  end subroutine burn_to_eos

end module burn_type_module
