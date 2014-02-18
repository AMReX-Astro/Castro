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
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module

  implicit none

  double precision, save, private :: gamma_const, K_const
  double precision, save, private :: mu_e
  integer         , save, private :: polytrope
  logical         , save, private :: assume_neutral

  public eos_init, eos

contains

  ! EOS initialization routine
  subroutine eos_init(small_temp, small_dens)

    use extern_probin_module, only: polytrope_gamma, polytrope_K, polytrope_type, polytrope_mu_e, eos_assume_neutral

    implicit none
 
    double precision, intent(in), optional :: small_temp
    double precision, intent(in), optional :: small_dens

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

    assume_neutral = eos_assume_neutral

    ! Small temperature and density parameters
 
    smallt = ZERO

    if (present(small_temp)) smallt = small_temp

    smalld = ZERO
 
    if (present(small_dens)) smalld = small_dens

    initialized = .true.
 
  end subroutine eos_init


  !---------------------------------------------------------------------------
  ! Castro interfaces 
  !---------------------------------------------------------------------------

  subroutine eos_get_polytrope_parameters(polytrope_out,gamma_out,K_out,mu_e_out)

    integer,          intent(out) :: polytrope_out
    double precision, intent(out) :: gamma_out, K_out, mu_e_out

    polytrope_out = polytrope
    gamma_out     = gamma_const
    K_out         = K_const
    mu_e_out      = mu_e

  end subroutine eos_get_polytrope_parameters

  subroutine eos_set_polytrope_parameters(polytrope_in,gamma_in,K_in,mu_e_in)

    integer,          intent(in) :: polytrope_in
    double precision, intent(in) :: gamma_in, K_in, mu_e_in

    polytrope   = polytrope_in
    gamma_const = gamma_in
    K_const     = K_in
    mu_e        = mu_e_in

  end subroutine eos_set_polytrope_parameters



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
    double precision :: dens, temp, enth, pres, eint, entr
    double precision :: dmudX

    ! get the mass of a nucleon from Avogadro's number.
    double precision, parameter :: m_nucleon = ONE / n_A

    integer :: k, n

    ! Make sure EOS is initialized before coming here.
    if (.not. initialized) call bl_error('EOS: not initialized')
    
    ! Calculate composition information
    call composition(state, assume_neutral)

    ! Sanity check: make sure that the mu_e calculated from the species
    ! is equal to the input parameter. This only matters for the polytropic gases
    ! where we used mu_e to calculate K_const.

    if (polytrope .eq. 1 .or. polytrope .eq. 2) then
    
      if (abs(state % mu_e - HALF) .gt. 1.d-8) then
        call bl_error("Calculated mu_e is not equal to the input parameter.")
      endif

    endif

    dens = state % rho
    temp = state % T
    pres = state % p
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

    state % cv = state % dedT
    state % cp = gamma_const * state % cv
 
    state % gam1 = gamma_const

    ! Compute dPdX, dedX, dhdX
    call composition_derivatives(state, assume_neutral)

    ! sound speed
    state % cs = sqrt(gamma_const*pres/dens)

    return
  end subroutine eos

end module eos_module
