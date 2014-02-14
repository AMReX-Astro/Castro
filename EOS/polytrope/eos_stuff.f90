! This is the equation of state for a polytropic fluid:
! P = K rho^gamma
!
! The internal energy and pressure are related via a gamma law:
!
! P = (gamma - 1) rho e
!
! gamma and K are fixed quantities for the run, and must either be
! supplied by the user or selected from a list of available options.
! Currently, we have fully degenerate ionized gases (both relativistic
! and non-relativistic), where the pressure is supplied by electrons.
!
! We compute the mean molecular weight, mu, based on the mass fractions
! of the different species. We assume complete ionization, so:
!
!   1/mu = sum_k { (1 + Z_k) X_k / A_k }
!
! The mean number of electrons per ion is:
!
!   1/mu_e = sum_k { X_k Z_k / A_k }
!
! This is assumed to be constant for the degenerate gases.

module eos_module

  use bl_types
  use bl_space
  use bl_constants_module, only: M_PI, ZERO, ONE, FOUR3RD, FIVE3RD
  use network, only: nspec, aion, zion
  use eos_type_module

  implicit none

  integer, parameter, public :: eos_input_rt = 1   ! density, temperature are inputs
  integer, parameter, public :: eos_input_rh = 2   ! density, enthalpy are inputs
  integer, parameter, public :: eos_input_tp = 3   ! temperature, pressure are inputs
  integer, parameter, public :: eos_input_rp = 4   ! density, pressure are inputs
  integer, parameter, public :: eos_input_re = 5   ! density, internal energy are inputs
  integer, parameter, public :: eos_input_ps = 6   ! pressure, entropy are inputs
  integer, parameter, public :: eos_input_ph = 7   ! pressure, enthalpy are inputs
  integer, parameter, public :: eos_input_th = 8   ! temperature, enthalpy are inputs

  real(kind=dp_t), save, private :: smallt, smalld
  real(kind=dp_t), save, private :: gamma_const, K_const
  real(kind=dp_t), save, private :: mu_e
  integer        , save, private :: polytrope

  logical, save, private :: initialized = .false.

  private nspec, aion, zion

  public eos_init, eos_get_small_temp, eos_get_small_dens, eos

contains

  ! EOS initialization routine -- this is used by both MAESTRO and Castro
  subroutine eos_init(small_temp, small_dens)

    use extern_probin_module, only: polytrope_gamma, polytrope_K, polytrope_type, polytrope_mu_e
    use bl_error_module

    implicit none
 
    real(kind=dp_t), intent(in), optional :: small_temp
    real(kind=dp_t), intent(in), optional :: small_dens

    ! Available pre-defined polytrope options:

    ! 1: Non-relativistic, fully degenerate electron gas
    ! 2: Relativistic, fully degenerate electron gas 

    if (polytrope_type > 0) then
      mu_e = polytrope_mu_e

      polytrope = polytrope_type
      if (polytrope .eq. 1) then
        gamma_const = FIVE3RD
        K_const     = 9.9154d12 ! (3 / pi)^(2/3) * h^2 / (20 * m_e * m_p^(5/3))
        K_const     = K_const / mu_e**gamma_const
      elseif (polytrope .eq. 2) then
        gamma_const = FOUR3RD
        K_const     = 1.2316d15 ! (3 / pi)^(1/3) * h c / (8 * m_p^(4/3))
        K_const     = K_const / mu_e**gamma_const
      else
        call bl_error('EOS: Polytrope type currently not defined')
      endif
    elseif (polytrope_gamma .gt. 0.d0 .and. polytrope_K .gt. 0.d0) then
      gamma_const = polytrope_gamma
      K_const     = polytrope_K
      mu_e        = 2.0d0 ! This will not be used
    else
      call bl_error('EOS: Neither polytrope type nor both gamma and K are defined')
    endif
 
    ! small temperature and density parameters
    if (present(small_temp)) then
       smallt = small_temp
    else
       smallt = 0.d0
    endif
 
    if (present(small_dens)) then
       smalld = small_dens
    else
       smalld = 0.d0
    endif

    initialized = .true.
 
  end subroutine eos_init


  !---------------------------------------------------------------------------
  ! Castro interfaces 
  !---------------------------------------------------------------------------
  subroutine eos_get_small_temp(small_temp_out)
 
    real(kind=dp_t), intent(out) :: small_temp_out
 
    small_temp_out = smallt
 
  end subroutine eos_get_small_temp
 
  subroutine eos_get_small_dens(small_dens_out)
 
    real(kind=dp_t), intent(out) :: small_dens_out
 
    small_dens_out = smalld
 
  end subroutine eos_get_small_dens

  subroutine eos_get_polytrope_parameters(polytrope_out,gamma_out,K_out,mu_e_out)

    integer,         intent(out) :: polytrope_out
    real(kind=dp_t), intent(out) :: gamma_out, K_out, mu_e_out

    polytrope_out = polytrope
    gamma_out     = gamma_const
    K_out         = K_const
    mu_e_out      = mu_e

  end subroutine eos_get_polytrope_parameters

  subroutine eos_set_polytrope_parameters(polytrope_in,gamma_in,K_in,mu_e_in)

    integer,         intent(in) :: polytrope_in
    real(kind=dp_t), intent(in) :: gamma_in, K_in, mu_e_in

    polytrope   = polytrope_in
    gamma_const = gamma_in
    K_const     = K_in
    mu_e        = mu_e_in

  end subroutine eos_set_polytrope_parameters

  !---------------------------------------------------------------------------
  ! The main interface
  !---------------------------------------------------------------------------
  subroutine eos(input, state, do_eos_diag_in, pt_index)

    use bl_error_module
    use bl_constants_module
    use fundamental_constants_module, only: k_B, n_A, hbar

! dens     -- mass density (g/cc)
! temp     -- temperature (K) -- not well-defined for a polytropic fluid
! xmass    -- the mass fractions of the individual isotopes
! pres     -- the pressure (dyn/cm**2)
! enthalpy -- the enthalpy (erg/g)
! eint     -- the internal energy (erg/g)
! c_v      -- specific heat at constant volume
! c_p      -- specific heat at constant pressure
! ne       -- number density of electrons + positrons
! eta      -- degeneracy parameter
! pele     -- electron pressure + positron pressure
! dPdT     -- d pressure/ d temperature
! dPdR     -- d pressure/ d density
! dEdT     -- d energy/ d temperature
! dEdR     -- d energy/ d density
! dPdX     -- d pressure / d xmass(k)
! dhdX     -- d enthalpy / d xmass(k)  -- AT CONSTANT PRESSURE!!!
! gam1     -- first adiabatic index (d log P/ d log rho) |_s
! cs       -- sound speed -- sqrt(gam1 p /rho) 
! entropy  -- entropy (erg/g/K) -- not well-defined for a polytropic fluid
!
! input = 1 means dens, temp    , and xmass are inputs
!       = 2 means dens, enthalpy, and xmass are inputs
!       = 3 means temp, pres    , and xmass are inputs
!       = 4 means dens, pres    , and xmass are inputs
!       = 5 means dens, eint    , and xmass are inputs
!       = 6 means pres, entr    , and xmass are inputs
!
!
! derivatives wrt X_k:
!
!   For an ideal gas, the thermodynamic quantities only depend on composition
!   through the mean molecular weight, mu.
!
!   Using the chain rule:
!
!   dp/dX_k = dp/d(mu) d(mu)/dX_k
!

    implicit none

    logical do_eos_diag
    integer, intent(in) :: input

    real(kind=dp_t) :: dens, temp, enth, pres, eint, entr
    real(kind=dp_t) :: xmass(nspec)
    real(kind=dp_t) :: c_v, c_p
    real(kind=dp_t) :: ne, eta, pele
    real(kind=dp_t) :: dPdT, dPdr, dedT, dedr
    real(kind=dp_t) :: gam1, cs

    integer, optional, intent(in   ) :: pt_index(:)


    ! local variables
    real(kind=dp_t) :: ymass(nspec)    
    real(kind=dp_t) :: mu
    real(kind=dp_t) :: dmudX, sum_y

    ! get the mass of a nucleon from Avogadro's number.
    real(kind=dp_t), parameter :: m_nucleon = 1.d0/n_A

    integer :: k, n

    ! general sanity checks
    if (.not. initialized) call bl_error('EOS: not initialized')
      
    !-------------------------------------------------------------------------
    ! compute mu -- the mean molecular weight
    !-------------------------------------------------------------------------

    ! Assume completely ionized species.

    sum_y  = ZERO
          
    do n = 1, nspec
       ymass(n) = xmass(n)*(1.d0 + zion(n))/aion(n)
       sum_y = sum_y + ymass(n)
    enddo
          
    mu = ONE/sum_y

    ! Sanity check: make sure that the mu_e calculated from the species
    ! is equal to the input parameter. This only matters for the polytropic gases
    ! where we used mu_e to calculate K_const.

    if (polytrope .eq. 1 .or. polytrope .eq. 2) then

      sum_y  = ZERO
          
      do n = 1, nspec
         ymass(n) = xmass(n)*zion(n)/aion(n)
         sum_y = sum_y + ymass(n)
      enddo
    
      if (abs(mu_e - one/sum_y) .gt. 1.d-8) then
        print *, mu_e, sum_y
        call bl_error("Calculated mu_e is not equal to the input parameter.")
      endif

    endif

    dens = state % rho
    temp = state % T
    pres = state % pres
    enth = state % h
    eint = state % e
    entr = state % s

    select case (input)

    !-------------------------------------------------------------------------
    ! Now do the calculations. In every case,
    ! make sure we have pressure, density, energy, and enthalpy.
    ! Relevant equations:
    ! h   = e + p / rho = (p / rho) * gamma / (gamma - 1) = e * gamma
    ! p   = K * (rho ** gamma) = (gamma - 1) * rho * e
    ! rho = (p / K)**(1 / gamma)
    ! e   = h - p / rho = (p / rho) / (gamma - 1)         = h / gamma
    !-------------------------------------------------------------------------
    case (eos_input_rh)

       ! dens, enthalpy, and xmass are inputs

       ! Solve for the pressure and energy:

       pres = enth * dens * (gamma_const - ONE) / gamma_const
       eint = enth / gamma_const


    case (eos_input_rt)

       ! dens, temp, and xmass are inputs
          
       ! Solve for the pressure, energy and enthalpy:

       pres = K_const * dens**gamma_const
       enth = pres / dens * gamma_const / (gamma_const - ONE)
       eint = enth / gamma_const


    case (eos_input_tp)

       ! temp, pres, and xmass are inputs
          
       ! Solve for the density, energy and enthalpy:

       dens = (pres / K_const)**(ONE / gamma_const)
       enth = pres / dens * gamma_const / (gamma_const - ONE)
       eint = enth / gamma_const


    case (eos_input_rp)

       ! dens, pres, and xmass are inputs

       ! Solve for the enthalpy and energy:

       enth = (pres / dens) * gamma_const / (gamma_const - ONE)
       eint = enth / gamma_const


    case (eos_input_re)

       ! dens, energy, and xmass are inputs

       ! Solve for the pressure and enthalpy:

       pres = (gamma_const - one) * dens * eint


    case (eos_input_ps)
       
       ! pressure, entropy and xmass are inputs

       ! Solve for the density, energy and enthalpy:

       dens = (pres / K_const)**(ONE / gamma_const)
       enth = pres / dens * gamma_const / (gamma_const - ONE)
       eint = enth / gamma_const



    case (eos_input_ph)

       ! pressure, enthalpy and xmass are inputs

       ! Solve for the density and energy:

       eint = enth / gamma_const
       dens = (pres / K_const)**(ONE / gamma_const)



    case (eos_input_th)

       ! temperature, enthalpy and xmass are inputs

       ! Solve for the density, energy and pressure:

       eint = enth / gamma_const
       pres = (gamma_const - ONE) * dens * eint
       dens = (pres / K_const)**(ONE / gamma_const)



    case default

       call bl_error('EOS: invalid input.')

    end select

    !-------------------------------------------------------------------------
    ! Now we have all relevant quantities, regardless of the inputs.
    !-------------------------------------------------------------------------

    state % T   = temp
    state % rho = dens
    state % h   = enth
    state % s   = entr
    state % e   = eint
    state % p   = pres

    ! Compute the thermodynamic derivatives and specific heats 
    state % dPdT = ZERO
    state % dPdr = gamma_const * pres / dens
    state % dedT = ZERO
    state % dedr = pres / (dens * dens)
    state % dsdT = ZERO
    state % dsdr = ZERO
    state % dhdT = ZERO
    state % dhdr = state % dedr + (gamma_const - ONE) * pres / dens**2

    state % c_v = state % dedT
    state % c_p = gamma_const * state % c_v
 
    state % gam1 = gamma_const

    do n = 1, nspec

       ! the species only come into p and e (and therefore h)
       ! through mu, so first compute dmu/dX
       !
       ! NOTE: an extra, constant term appears in dmudx, which
       ! results from writing mu = sum { X_k} / sum {X_k / A_k}
       ! (for the neutral, analogous for the ionized).  The
       ! numerator is simply 1, but we can differentiate
       ! wrt it, giving the constant mu(k) term in dmudx.  Since
       ! dPdX only appears in a sum over species creation rate 
       ! (omegadot) and sum{omegadot} = 0, this term has no effect.
       ! If is added simply for completeness.

       dmudX =  (mu/aion(n))*(aion(n) - mu*(ONE + zion(n)))

       state % dPdX(n) = -(pres/mu)*dmudX
       state % dedX(n) = -(eint/mu)*dmudX
          
       ! dhdX is at constant pressure -- see paper III for details
       state % dhdX(n) = dedX(n) + &
            (pres/dens**2 - dedR)*dPdX(n)/dPdr
    enddo

    ! electron-specific stuff (really for the degenerate EOS)
    state % ne   = ZERO
    state % eta  = ZERO
    state % pele = ZERO

    ! sound speed
    state % cs = sqrt(gamma_const*pres/dens)

    return
  end subroutine eos_old

end module eos_module
