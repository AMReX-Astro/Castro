module eos_type_module

  use bl_types
  use network
  use eos_data_module
  use mempool_module
  use bl_constants_module

  implicit none

  ! A generic structure holding thermodynamic quantities and their derivatives,
  ! plus some other quantities of interest.

  ! rho      -- mass density (g/cm**3)
  ! T        -- temperature (K)
  ! xn       -- the mass fractions of the individual isotopes
  ! p        -- the pressure (dyn/cm**2)
  ! h        -- the enthalpy (erg/g)
  ! e        -- the internal energy (erg/g)
  ! s        -- the entropy (erg/g/K)
  ! c_v      -- specific heat at constant volume
  ! c_p      -- specific heat at constant pressure
  ! ne       -- number density of electrons + positrons
  ! np       -- number density of positrons only
  ! eta      -- degeneracy parameter
  ! pele     -- electron pressure + positron pressure
  ! ppos     -- position pressure only
  ! mu       -- mean molecular weight
  ! mu_e     -- mean number of nucleons per electron
  ! y_e      -- electron fraction == 1 / mu_e
  ! dPdT     -- d pressure/ d temperature
  ! dPdr     -- d pressure/ d density
  ! dedT     -- d energy/ d temperature
  ! dedr     -- d energy/ d density
  ! dsdT     -- d entropy/ d temperature
  ! dsdr     -- d entropy/ d density
  ! dhdT     -- d enthalpy/ d temperature
  ! dhdr     -- d enthalpy/ d density
  ! dPdX     -- d pressure / d xmass
  ! dhdX     -- d enthalpy / d xmass at constant pressure
  ! gam1     -- first adiabatic index (d log P/ d log rho) |_s
  ! cs       -- sound speed
  ! abar     -- average atomic number ( sum_k {X_k} ) / ( sum_k {X_k/A_k} )
  ! zbar     -- average proton number ( sum_k {Z_k X_k/ A_k} ) / ( sum_k {X_k/A_k} )
  ! dpdA     -- d pressure/ d abar
  ! dpdZ     -- d pressure/ d zbar
  ! dedA     -- d energy/ d abar
  ! dedZ     -- d energy/ d zbar

  ! Initialize the main quantities to an unphysical number
  ! so that we know if the user forgot to initialize them
  ! when calling the EOS in a particular mode.

  double precision, parameter :: init_num  = -1.0d200
  double precision, parameter :: init_test = -1.0d199

  type :: eos_t

    double precision :: rho         = init_num
    double precision :: T           = init_num
    double precision :: p           = init_num
    double precision :: e           = init_num
    double precision :: h           = init_num
    double precision :: s           = init_num
    double precision :: dpdT        = init_num
    double precision :: dpdr        = init_num
    double precision :: dedT        = init_num
    double precision :: dedr        = init_num
    double precision :: dhdT        = init_num
    double precision :: dhdr        = init_num
    double precision :: dsdT        = init_num
    double precision :: dsdr        = init_num
    double precision :: dpde        = init_num
    double precision :: dpdr_e      = init_num

    double precision :: xn(nspec)   = init_num
    double precision :: aux(naux)   = init_num
    double precision :: cv          = init_num
    double precision :: cp          = init_num
    double precision :: xne         = init_num
    double precision :: xnp         = init_num
    double precision :: eta         = init_num
    double precision :: pele        = init_num
    double precision :: ppos        = init_num
    double precision :: mu          = init_num
    double precision :: mu_e        = init_num
    double precision :: y_e         = init_num
    double precision :: dedX(nspec) = init_num
    double precision :: dpdX(nspec) = init_num
    double precision :: dhdX(nspec) = init_num
    double precision :: gam1        = init_num
    double precision :: cs          = init_num

    double precision :: abar        = init_num
    double precision :: zbar        = init_num
    double precision :: dpdA        = init_num

    double precision :: dpdZ        = init_num
    double precision :: dedA        = init_num
    double precision :: dedZ        = init_num

    logical :: reset                = .false.
    logical :: check_small          = .true.

  end type eos_t

contains

  ! Given a set of mass fractions, calculate quantities that depend
  ! on the composition like abar and zbar.

  subroutine composition(state)

    use bl_constants_module
    use network

    implicit none

    type (eos_t), intent(inout) :: state

    ! Calculate abar, the mean nucleon number,
    ! zbar, the mean proton number,
    ! mu, the mean molecular weight,
    ! mu_e, the mean number of nucleons per electron, and
    ! y_e, the electron fraction.

    state % mu_e = ONE / (sum(state % xn(:) * zion(:) / aion(:)))
    state % y_e = ONE / state % mu_e

    state % abar = ONE / (sum(state % xn(:) / aion(:)))
    state % zbar = state % abar / state % mu_e

  end subroutine composition

  ! Compute thermodynamic derivatives with respect to xn(:)

  subroutine composition_derivatives(state)

    use bl_constants_module
    use network

    implicit none

    type (eos_t), intent(inout) :: state

    integer :: i

    state % dpdX(:) = state % dpdA * (state % abar/aion(:))   &
                                   * (aion(:) - state % abar) &
                    + state % dpdZ * (state % abar/aion(:))   &
                                   * (zion(:) - state % zbar)

    state % dEdX(:) = state % dedA * (state % abar/aion(:))   &
                                   * (aion(:) - state % abar) &
                    + state % dedZ * (state % abar/aion(:))   &
                                   * (zion(:) - state % zbar)

    if (state % dPdr > ZERO) then

       state % dhdX(:) = state % dedX(:) &
                       + (state % p / state % rho**2 - state % dedr) &
                       *  state % dPdX(:) / state % dPdr

    endif

  end subroutine composition_derivatives


  ! Normalize the mass fractions: they must be individually positive
  ! and less than one, and they must all sum to unity.

  subroutine normalize_abundances(state)

    use bl_constants_module
    use network

    implicit none

    type (eos_t), intent(inout) :: state

    state % xn(:) = max(smallx, min(ONE, state % xn(:)))

    state % xn(:) = state % xn(:) / sum(state % xn(:))

  end subroutine normalize_abundances

end module eos_type_module
