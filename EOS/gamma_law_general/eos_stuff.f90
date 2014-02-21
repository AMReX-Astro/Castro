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
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module

  implicit none

  double precision, save, public  :: gamma_const
  logical         , save, public  :: assume_neutral

  public eos_init, eos

contains

  ! EOS initialization routine
  subroutine eos_init(small_temp, small_dens)

    use extern_probin_module, only: eos_gamma, eos_assume_neutral

    implicit none
 
    double precision, intent(in), optional :: small_temp
    double precision, intent(in), optional :: small_dens
 
    ! constant ratio of specific heats
    if (eos_gamma .gt. 0.d0) then
       gamma_const = eos_gamma
    else
       gamma_const = FIVE3RD
    end if
 
    ! Small temperature and density parameters

    smallt = ZERO

    if (present(small_temp)) smallt = small_temp

    smalld = ZERO

    if (present(small_dens)) smalld = small_dens

    assume_neutral = eos_assume_neutral

    initialized = .true.
 
  end subroutine eos_init



  !---------------------------------------------------------------------------
  ! The main interface
  !---------------------------------------------------------------------------
  subroutine eos(input, state, do_eos_diag_in, pt_index)

    use fundamental_constants_module, only: k_B, n_A, hbar

    implicit none

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state
    logical, optional, intent(in   ) :: do_eos_diag_in
    integer, optional, intent(in   ) :: pt_index(:)

    ! Local variables
    double precision :: ymass(nspec)
    double precision :: mu
    double precision :: dmudX, sum_y
    double precision :: dens, temp

    ! Get the mass of a nucleon from Avogadro's number.
    double precision, parameter :: m_nucleon = ONE / n_A

    integer :: k, n

    ! general sanity checks
    if (.not. initialized) call bl_error('EOS: not initialized')

    do_eos_diag = .false.

    if (present(do_eos_diag_in)) do_eos_diag = do_eos_diag_in

    ! Get abar, zbar, mu.

    call composition(state, assume_neutral)

    !-------------------------------------------------------------------------
    ! For all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs.
    !-------------------------------------------------------------------------

    select case (input)

    case (eos_input_rt)

       ! dens, temp and xmass are inputs

       ! We don't need to do anything here


    case (eos_input_rh)

       ! dens, enthalpy, and xmass are inputs

       ! Solve for the temperature:
       ! h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
       temp = (enthalpy*mu*m_nucleon/k_B)*(gamma_const - ONE)/gamma_const


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
       temp = eint*mu*m_nucleon*(gamma_const-ONE)/k_B


    case (eos_input_ps)
       
       ! pressure entropy, and xmass are inputs

       ! Solve for the temperature
       ! Invert Sackur-Tetrode eqn (below) using 
       ! rho = p mu m_nucleon / (k T)
       temp = pres**(TWO/FIVE) * &
            ( TWO*M_PI*hbar*hbar/(mu*m_nucleon) )**(THREE/FIVE) * &
            dexp(TWO*mu*m_nucleon*entropy/(FIVE*k_B) - ONE) / k_B

       ! Solve for the density
       ! rho = p mu m_nucleon / (k T)
       dens = pres*mu*m_nucleon/(k_B*temp)



    case (eos_input_ph)

       ! pressure, enthalpy and xmass are inputs

       ! Solve for temperature and density
       dens = pres / enthalpy * gamma_const / (gamma_const - ONE)

       temp = pres*mu*m_nucleon/(k_B*dens)



    case (eos_input_th)

       ! temperature, enthalpy and xmass are inputs

       ! Solve for density
       dens = pres*mu*m_nucleon/(k_B*temp)



    case default

       call bl_error('EOS: invalid input.')

    end select

    !-------------------------------------------------------------------------
    ! Now we have the density and temperature (and mass fractions /
    ! mu), regardless of the inputs.
    !-------------------------------------------------------------------------

    state % T   = temp
    state % rho = dens

    ! Compute the pressure simply from the ideal gas law, and the
    ! specific internal energy using the gamma-law EOS relation.
    state % p = dens*k_B*temp/(mu*m_nucleon)
    state % e = state % p/(gamma_const - ONE)/dens

    ! enthalpy is h = e + p/rho
    state % h = state % e + state % p / dens

    ! entropy (per gram) of an ideal monoatomic gas (the Sackur-Tetrode equation)
    ! NOTE: this expression is only valid for gamma = 5/3.
    state % s = (k_B/(mu*m_nucleon))*(2.5_dp_t + &
         log( ( (mu*m_nucleon)**2.5/dens )*(k_B*temp)**1.5_dp_t / (TWO*M_PI*hbar*hbar)**1.5_dp_t ) )

    ! Compute the thermodynamic derivatives and specific heats 
    state % dpdT = state % p / temp
    state % dpdr = state % p / dens
    state % dedT = state % e / temp
    state % dedr = 0.d0
    state % dsdT = 0.d0
    state % dsdr = 0.d0
    state % dhdT = state % dedT + state % dpdT / dens
    state % dhdr = state % dedr + (gamma_const - ONE) * pres / dens**2

    state % cv = dedT
    state % cp = gamma_const * state % cv

    state % gam1 = gamma_const

    ! Get dpdX, dedX, dhdX.

    state % dpdA = - state % p / state % abar
    state % dedA = - state % e / state % abar

    if (assume_neutral) then
      state % dpdZ = ZERO
      state % dedZ = ZERO
    else
      state % dpdZ = state % p / (ONE + state % zbar)
      state % dedZ = state % e / (ONE + state % zbar)
    endif

    call composition_derivatives(state, assume_neutral)

    ! sound speed
    state % cs = sqrt(gamma_const*pres/dens)

    return
  end subroutine eos

end module eos_module
