! This is a constant gamma equation of state, using an ideal gas.
!
! The gas may either be completely ionized or completely neutral.
!
! The ratio of specific heats (gamma) is allowed to vary.  NOTE: the
! expression for entropy is only valid for an ideal MONATOMIC gas
! (gamma = 5/3).  

module actual_eos_module

  use bl_types
  use bl_space
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module

  implicit none

  character (len=64) :: eos_name = "gamma_law_general"  
  
  double precision  :: gamma_const

contains

  subroutine actual_eos_init

    use extern_probin_module, only: eos_gamma

    implicit none
 
    ! constant ratio of specific heats
    if (eos_gamma .gt. 0.d0) then
       gamma_const = eos_gamma
    else
       gamma_const = FIVE3RD
    end if
 
  end subroutine actual_eos_init



  subroutine actual_eos(input, state)

    use fundamental_constants_module, only: k_B, n_A, hbar

    implicit none

    integer,             intent(in   ) :: input
    type (eos_t_vector), intent(inout) :: state

    ! Get the mass of a nucleon from Avogadro's number.
    double precision, parameter :: m_nucleon = ONE / n_A
    double precision, parameter :: fac = ONE / (TWO*M_PI*hbar*hbar)**1.5d0

    integer :: j
    double precision :: Tinv, rhoinv

    !-------------------------------------------------------------------------
    ! For all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs.
    !-------------------------------------------------------------------------

    select case (input)

    case (eos_input_rt)

       ! dens, temp and xmass are inputs

       ! We don't need to do anything here
       continue

    case (eos_input_rh)

       ! dens, enthalpy, and xmass are inputs

       ! Solve for the temperature:
       ! h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
       do j = 1, state % N
          state % T(j)   = (state % h(j) * state % mu(j) * m_nucleon / k_B)*(gamma_const - ONE)/gamma_const
       end do

    case (eos_input_tp)

       ! temp, pres, and xmass are inputs
       
       ! Solve for the density:
       ! p = rho k T / (mu m_nucleon)
       do j = 1, state % N
          state % rho(j) =  state % p(j) * state % mu(j) * m_nucleon / (k_B * state % T(j))
       end do

    case (eos_input_rp)

       ! dens, pres, and xmass are inputs

       ! Solve for the temperature:
       ! p = rho k T / (mu m_nucleon)
       do j = 1, state % N
          state % T(j)   = state % p(j) * state % mu(j) * m_nucleon / (k_B * state % rho(j))
       end do

    case (eos_input_re)

       ! dens, energy, and xmass are inputs
       
       ! Solve for the temperature
       ! e = k T / [(mu m_nucleon)*(gamma-1)]
       do j = 1, state % N
          state % T(j)   = state % e(j) * state % mu(j) * m_nucleon * (gamma_const-ONE) / k_B
       end do

    case (eos_input_ps)

       ! pressure entropy, and xmass are inputs
       
       ! Solve for the temperature
       ! Invert Sackur-Tetrode eqn (below) using 
       ! rho = p mu m_nucleon / (k T)
       do j = 1, state % N
          state % T(j)   = state % p(j)**(TWO/FIVE) * &
               ( TWO*M_PI*hbar*hbar/(state % mu(j)*m_nucleon) )**(THREE/FIVE) * &
               dexp(TWO*state % mu(j)*m_nucleon*state % s(j)/(FIVE*k_B) - ONE) / k_B

          ! Solve for the density
          ! rho = p mu m_nucleon / (k T)
          state % rho(j) =  state % p(j) * state % mu(j) * m_nucleon / (k_B * state % T(j))
       end do

    case (eos_input_ph)

       ! pressure, enthalpy and xmass are inputs
       
       ! Solve for temperature and density
       do j = 1, state % N
          state % rho(j) =  state % p(j) / state % h(j) * gamma_const / (gamma_const - ONE)
          state % T(j)   = state % p(j) * state % mu(j) * m_nucleon / (k_B * state % rho(j))
       end do

    case (eos_input_th)

       ! temperature, enthalpy and xmass are inputs

       ! This system is underconstrained.
       
       call bl_error('EOS: eos_input_th is not a valid input for the gamma law EOS.')

    case default
       
       call bl_error('EOS: invalid input.')
       
    end select
    
    !-------------------------------------------------------------------------
    ! Now we have the density and temperature (and mass fractions /
    ! mu), regardless of the inputs.
    !-------------------------------------------------------------------------
    
    do j = 1, state % N
       
       Tinv = ONE / state % T(j)
       rhoinv = ONE / state % rho(j)

       ! Compute the pressure simply from the ideal gas law, and the
       ! specific internal energy using the gamma-law EOS relation.
       state % p(j) = state % rho(j) * state % T(j) * (k_B / (state % mu(j)*m_nucleon))
       state % e(j) = state % p(j)/(gamma_const - ONE) * rhoinv

       ! enthalpy is h = e + p/rho
       state % h(j) = state % e(j) + state % p(j) * rhoinv

       ! entropy (per gram) of an ideal monoatomic gas (the Sackur-Tetrode equation)
       ! NOTE: this expression is only valid for gamma = 5/3.
       state % s(j) = (k_B/(state % mu(j)*m_nucleon))*(2.5_dp_t + &
            log( ( (state % mu(j)*m_nucleon)**2.5_dp_t * rhoinv )*(k_B * state % T(j))**1.5_dp_t * fac ) )

       ! Compute the thermodynamic derivatives and specific heats 
       state % dpdT(j) = state % p(j) * Tinv
       state % dpdr(j) = state % p(j) * rhoinv
       state % dedT(j) = state % e(j) * Tinv
       state % dedr(j) = ZERO
       state % dsdT(j) = 1.5_dp_t * (k_B / (state % mu(j) * m_nucleon)) * Tinv
       state % dsdr(j) = - (k_B / (state % mu(j) * m_nucleon)) * rhoinv
       state % dhdT(j) = state % dedT(j) + state % dpdT(j) * rhoinv
       state % dhdr(j) = ZERO

       state % cv(j) = state % dedT(j)
       state % cp(j) = gamma_const * state % cv(j)

       state % gam1(j) = gamma_const

       state % dpdr_e(j) = state % dpdr(j) - state % dpdT(j) * state % dedr(j) * (ONE/state % dedT(j))
       state % dpde(j)   = state % dpdT(j) * (ONE/ state % dedT(j))

       ! sound speed
       state % cs(j) = sqrt(gamma_const * state % p(j) * rhoinv)

       state % dpdA(j) = - state % p(j) * (ONE/state % abar(j))
       state % dedA(j) = - state % e(j) * (ONE/state % abar(j))

       if (assume_neutral) then
         state % dpdZ(j) = ZERO
         state % dedZ(j) = ZERO
       else
         state % dpdZ(j) = state % p(j) * (ONE/(ONE + state % zbar(j)))
         state % dedZ(j) = state % e(j) * (ONE/(ONE + state % zbar(j)))
       endif

    enddo

  end subroutine actual_eos

end module actual_eos_module
