! This is a constant gamma equation of state, using an ideal gas.
!
! We compute the mean molecular weight, mu, based on the mass fractions
! of the different species.
!
!   NOTE: in the helmholtz EOS, we use Abar, which is the weighted
!   ion mass.  Abar does not include any free electron contribution.
!
! There are 2 ways to compute mu, depending on whether the gas is
! completely ionized or not.
!
! For a neutral gas, the mean molecular weight is:
!
!   1/mu = sum_k { X_k / A_k }
!
! For a completely ionized gas, the mean molecular weight is:
!
!   1/mu = sum_k { (1 + Z_k) X_k / A_k }
!
! The ratio of specific heats (gamma) is allowed to vary.  NOTE: the
! expression for entropy is only valid for an ideal MONATOMIC gas
! (gamma = 5/3).  

module eos_module

  use bl_types
  use bl_space
  use bl_constants_module
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
  real(kind=dp_t), save, public  :: gamma_const
  logical        , save, public  :: assume_neutral

  logical, save, private :: initialized = .false.

  private nspec, aion, zion

  public eos_init, eos_get_small_temp, eos_get_small_dens, eos

contains

  ! EOS initialization routine -- this is used by both MAESTRO and Castro
  subroutine eos_init(small_temp, small_dens)

    use extern_probin_module, only: eos_gamma, eos_assume_neutral

    implicit none
 
    real(kind=dp_t), intent(in), optional :: small_temp
    real(kind=dp_t), intent(in), optional :: small_dens
 
    ! constant ratio of specific heats
    if (eos_gamma .gt. 0.d0) then
       gamma_const = eos_gamma
    else
       gamma_const = 5.d0/3.d0
    end if
 
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

    assume_neutral = eos_assume_neutral

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

  !---------------------------------------------------------------------------
  ! The main interface -- this is used directly by MAESTRO
  !---------------------------------------------------------------------------
  subroutine eos(input, state, do_eos_diag_in, pt_index)

    use bl_error_module
    use bl_constants_module
    use fundamental_constants_module, only: k_B, n_A, hbar

! dens     -- mass density (g/cc)
! temp     -- temperature (K)
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
! cs       -- sound speed -- note that this is the non-relativistic one
!             (we compute it in this wrapper as sqrt(gam1 p /rho) instead
!             of taking the relativistic version from helmeos.
! entropy  -- entropy (erg/g/K)  NOTE: presently the entropy expression is 
!             valid only for an ideal MONATOMIC gas (gamma = 5/3).
!
! input = 1 means dens, temp    , and xmass are inputs, return enthalpy, eint
!       = 2 means dens, enthalpy, and xmass are inputs, return temp    , eint
!                (note, temp should be filled with an initial guess)
!       = 3 means temp, pres    , and xmass are inputs, return dens    , etc
!       = 4 means dens, pres    , and xmass are inputs, return temp    , etc
!       = 5 means dens, eint    , and xmass are inputs, return temp    , etc
!       = 6 means pres, entr    , and xmass are inputs, return temp    , etc
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

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state
    logical, optional, intent(in   ) :: do_eos_diag_in
    integer, optional, intent(in   ) :: pt_index(:)

    ! local variables
    real(kind=dp_t) :: ymass(nspec)
    real(kind=dp_t) :: mu
    real(kind=dp_t) :: dmudX, sum_y
    real(kind=dp_t) :: dens, temp

    ! get the mass of a nucleon from Avogadro's number.
    real(kind=dp_t), parameter :: m_nucleon = 1.d0/n_A

    integer :: k, n

    ! general sanity checks
    if (.not. initialized) call bl_error('EOS: not initialized')
      

    !-------------------------------------------------------------------------
    ! compute mu -- the mean molecular weight
    !-------------------------------------------------------------------------
    if (assume_neutral) then
       ! assume completely neutral atoms

       sum_y  = 0.d0
          
       do n = 1, nspec
          ymass(n) = xmass(n)/aion(n)
          sum_y = sum_y + ymass(n)
       enddo
          
       mu = 1.d0/sum_y

    else
       ! assume completely ionized species

       sum_y  = 0.d0
          
       do n = 1, nspec
          ymass(n) = xmass(n)*(1.d0 + zion(n))/aion(n)
          sum_y = sum_y + ymass(n)
       enddo
          
       mu = 1.d0/sum_y

    endif

    !-------------------------------------------------------------------------
    ! for all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs
    !-------------------------------------------------------------------------

    select case (input)

    case (eos_input_rt)

       ! We don't need to do anything here


    case (eos_input_rh)

       ! dens, enthalpy, and xmass are inputs

       ! Solve for the temperature:
       ! h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
       temp = (enthalpy*mu*m_nucleon/k_B)*(gamma_const - 1.0_dp_t)/gamma_const


    case (eos_input_tp)

       ! temp, pres, and xmass are inputs
          
       ! Solve for the density:
       ! p = rho k T / (mu m_nucleon)
       dens = pres*mu*m_nucleon/(k_B*temp)


    case (eos_input_rp)

       ! dens, pres, and xmass are inputs

       ! Solve for the temperature:
       ! p = rho k T / (mu m_nucleon)
       temp = pres*mu*m_nucleon/(k_B*dens)


    case (eos_input_re)

       ! dens, energy, and xmass are inputs

       ! Solve for the temperature
       ! e = k T / [(mu m_nucleon)*(gamma-1)]
       temp = eint*mu*m_nucleon*(gamma_const-1.0_dp_t)/k_B


    case (eos_input_ps)
       
       ! pressure and entropy are inputs

       ! Solve for the temperature
       ! Invert Sackur-Tetrode eqn (below) using 
       ! rho = p mu m_nucleon / (k T)
       temp = pres**(2.0_dp_t/5.0_dp_t) * &
            ( 2.0_dp_t*M_PI*hbar*hbar/(mu*m_nucleon) )**(3.0_dp_t/5.0_dp_t) * &
            dexp(2.0_dp_t*mu*m_nucleon*entropy/(5.0_dp_t*k_B) - 1.0_dp_t) / &
            k_B

       ! Solve for the density
       ! rho = p mu m_nucleon / (k T)
       dens = pres*mu*m_nucleon/(k_B*temp)


    case (eos_input_ph)

       call bl_error('EOS: This input is not currently supported in the gamma law EOS.')



    case (eos_input_th)

       call bl_error('EOS: This input is not currently supported in the gamma law EOS.')


    case default

       call bl_error('EOS: invalid input.')

    end select

    !-------------------------------------------------------------------------
    ! now we have the density and temperature (and mass fractions /
    ! mu), regardless of the inputs.
    !-------------------------------------------------------------------------

    state % T = temp
    state % rho = dens

    ! compute the pressure simply from the ideal gas law, and the
    ! specific internal energy using the gamma-law EOS relation
    state % p = dens*k_B*temp/(mu*m_nucleon)
    state % e = state % p/(gamma_const - 1.0_dp_t)/dens

    ! enthalpy is h = e + p/rho
    state % h = state % e + state % p / dens

    ! entropy (per gram) of an ideal monoatomic gas (the Sackur-Tetrode equation)
    ! NOTE: this expression is only valid for gamma = 5/3.
    state % s = (k_B/(mu*m_nucleon))*(2.5_dp_t + &
         log( ( (mu*m_nucleon)**2.5/dens )*(k_B*temp)**1.5_dp_t / (2.0_dp_t*M_PI*hbar*hbar)**1.5_dp_t ) )

    ! compute the thermodynamic derivatives and specific heats 
    state % dpdT = state % p / temp
    state % dpdr = state % p / dens
    state % dedT = state % e / temp
    state % dedr = 0.d0
    state % dsdT = 0.d0
    state % dsdr = 0.d0

    state % c_v = dedT
    state % c_p = gamma_const * c_v

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

       if (assume_neutral) then
          dmudX =  (mu/aion(n))*(aion(n) - mu)
       else
          dmudX =  (mu/aion(n))*(aion(n) - mu*(1.0_dp_t + zion(n)))
       endif

       state % dPdX(n) = -(pres/mu)*dmudX
       state % dedX(n) = -(eint/mu)*dmudX
          
       ! dhdX is at constant pressure -- see paper III for details
       state % dhdX(n) = dedX(n) + &
            (pres/dens**2 - dedR)*dPdX(n)/dPdr
    enddo

    ! electron-specific stuff (really for the degenerate EOS)
    state % ne = 0.d0
    state % eta = 0.d0
    state % pele = 0.d0

    ! sound speed
    state % cs = sqrt(gamma_const*pres/dens)

    return
  end subroutine eos

end module eos_module
