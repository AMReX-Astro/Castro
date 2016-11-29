module eos_type_module

  use bl_types, only: dp_t
  use network, only: nspec, naux

  implicit none

  integer, parameter :: eos_input_rt = 1  ! rho, T are inputs
  integer, parameter :: eos_input_rh = 2  ! rho, h are inputs
  integer, parameter :: eos_input_tp = 3  ! T, p are inputs
  integer, parameter :: eos_input_rp = 4  ! rho, p are inputs
  integer, parameter :: eos_input_re = 5  ! rho, e are inputs
  integer, parameter :: eos_input_ps = 6  ! p, s are inputs
  integer, parameter :: eos_input_ph = 7  ! p, h are inputs
  integer, parameter :: eos_input_th = 8  ! T, h are inputs

  ! these are used to allow for a generic interface to the 
  ! root finding
  integer, parameter :: itemp = 1
  integer, parameter :: idens = 2
  integer, parameter :: iener = 3
  integer, parameter :: ienth = 4
  integer, parameter :: ientr = 5
  integer, parameter :: ipres = 6

  ! error codes
  integer, parameter :: ierr_general         = 1
  integer, parameter :: ierr_input           = 2
  integer, parameter :: ierr_iter_conv       = 3
  integer, parameter :: ierr_neg_e           = 4
  integer, parameter :: ierr_neg_p           = 5
  integer, parameter :: ierr_neg_h           = 6
  integer, parameter :: ierr_neg_s           = 7
  integer, parameter :: ierr_iter_var        = 8
  integer, parameter :: ierr_init            = 9
  integer, parameter :: ierr_init_xn         = 10
  integer, parameter :: ierr_out_of_bounds   = 11
  integer, parameter :: ierr_not_implemented = 12

  ! Minimum and maximum thermodynamic quantities permitted by the EOS.

  real(dp_t), save :: mintemp = 1.d-200
  real(dp_t), save :: maxtemp = 1.d200
  real(dp_t), save :: mindens = 1.d-200
  real(dp_t), save :: maxdens = 1.d200
  real(dp_t), save :: minx    = 1.d-200
  real(dp_t), save :: maxx    = 1.d0 + 1.d-12
  real(dp_t), save :: minye   = 1.d-200
  real(dp_t), save :: maxye   = 1.d0 + 1.d-12
  real(dp_t), save :: mine    = 1.d-200
  real(dp_t), save :: maxe    = 1.d200
  real(dp_t), save :: minp    = 1.d-200
  real(dp_t), save :: maxp    = 1.d200
  real(dp_t), save :: mins    = 1.d-200
  real(dp_t), save :: maxs    = 1.d200
  real(dp_t), save :: minh    = 1.d-200
  real(dp_t), save :: maxh    = 1.d200

  !$acc declare &
  !$acc create(mintemp, maxtemp, mindens, maxdens, minx, maxx, minye, maxye) &
  !$acc create(mine, maxe, minp, maxp, mins, maxs, minh, maxh)

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

  type :: eos_t

    real(dp_t) :: rho
    real(dp_t) :: T
    real(dp_t) :: p
    real(dp_t) :: e
    real(dp_t) :: h
    real(dp_t) :: s
    real(dp_t) :: xn(nspec)
    real(dp_t) :: aux(naux)

    real(dp_t) :: dpdT
    real(dp_t) :: dpdr
    real(dp_t) :: dedT
    real(dp_t) :: dedr
    real(dp_t) :: dhdT
    real(dp_t) :: dhdr
    real(dp_t) :: dsdT
    real(dp_t) :: dsdr
    real(dp_t) :: dpde
    real(dp_t) :: dpdr_e

    real(dp_t) :: cv
    real(dp_t) :: cp
    real(dp_t) :: xne
    real(dp_t) :: xnp
    real(dp_t) :: eta
    real(dp_t) :: pele
    real(dp_t) :: ppos
    real(dp_t) :: mu
    real(dp_t) :: mu_e
    real(dp_t) :: y_e
    real(dp_t) :: dedX(nspec)
    real(dp_t) :: dpdX(nspec)
    real(dp_t) :: dhdX(nspec)
    real(dp_t) :: gam1
    real(dp_t) :: cs

    real(dp_t) :: abar
    real(dp_t) :: zbar
    real(dp_t) :: dpdA

    real(dp_t) :: dpdZ
    real(dp_t) :: dedA
    real(dp_t) :: dedZ

  end type eos_t

contains

  ! Given a set of mass fractions, calculate quantities that depend
  ! on the composition like abar and zbar.

  subroutine composition(state)

    !$acc routine seq

    use bl_constants_module, only: ONE
    use network, only: aion, aion_inv, zion

    implicit none

    type (eos_t), intent(inout) :: state

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

  ! Compute thermodynamic derivatives with respect to xn(:)

  subroutine composition_derivatives(state)

    !$acc routine seq

    use bl_constants_module, only: ZERO
    use network, only: aion, aion_inv, zion

    implicit none

    type (eos_t), intent(inout) :: state

    state % dpdX(:) = state % dpdA * (state % abar * aion_inv(:))   &
                                   * (aion(:) - state % abar) &
                    + state % dpdZ * (state % abar * aion_inv(:))   &
                                   * (zion(:) - state % zbar)

    state % dEdX(:) = state % dedA * (state % abar * aion_inv(:))   &
                                   * (aion(:) - state % abar) &
                    + state % dedZ * (state % abar * aion_inv(:))   &
                                   * (zion(:) - state % zbar)

    if (state % dPdr .ne. ZERO) then

       state % dhdX(:) = state % dedX(:) &
                       + (state % p / state % rho**2 - state % dedr) &
                       *  state % dPdX(:) / state % dPdr

    endif

  end subroutine composition_derivatives



  ! Normalize the mass fractions: they must be individually positive
  ! and less than one, and they must all sum to unity.

  subroutine normalize_abundances(state)

    !$acc routine seq

    use bl_constants_module, only: ONE
    use extern_probin_module, only: small_x

    implicit none

    type (eos_t), intent(inout) :: state

    state % xn = max(small_x, min(ONE, state % xn))

    state % xn = state % xn / sum(state % xn)

  end subroutine normalize_abundances



  ! Ensure that inputs are within reasonable limits.

  subroutine clean_state(state)

    !$acc routine seq

    implicit none

    type (eos_t), intent(inout) :: state

    state % T = min(maxtemp, max(mintemp, state % T))
    state % rho = min(maxdens, max(mindens, state % rho))

  end subroutine clean_state



  ! Print out details of the state.

  subroutine print_state(state)

    implicit none

    type (eos_t), intent(in) :: state

    print *, 'DENS = ', state % rho
    print *, 'TEMP = ', state % T
    print *, 'X    = ', state % xn
    print *, 'Y_E  = ', state % y_e

  end subroutine print_state



  subroutine eos_get_small_temp(small_temp_out)

    !$acc routine seq

    implicit none

    real(dp_t), intent(out) :: small_temp_out

    small_temp_out = mintemp

  end subroutine eos_get_small_temp



  subroutine eos_get_small_dens(small_dens_out)

    !$acc routine seq

    implicit none

    real(dp_t), intent(out) :: small_dens_out

    small_dens_out = mindens

  end subroutine eos_get_small_dens



  subroutine eos_get_max_temp(max_temp_out)

    !$acc routine seq

    implicit none

    real(dp_t), intent(out) :: max_temp_out

    max_temp_out = maxtemp

  end subroutine eos_get_max_temp



  subroutine eos_get_max_dens(max_dens_out)

    !$acc routine seq

    implicit none

    real(dp_t), intent(out) :: max_dens_out

    max_dens_out = maxdens

  end subroutine eos_get_max_dens

end module eos_type_module
