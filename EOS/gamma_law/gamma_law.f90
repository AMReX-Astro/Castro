! This is a constant gamma equation of state, using an ideal gas.
!
! This a simplified version of the more general eos_gamma_general.
!

module actual_eos_module

  use bl_types
  use bl_error_module
  use bl_constants_module
  use eos_type_module
  use eos_data_module

  implicit none

  character (len=64) :: eos_name = "gamma_law"  
  
  double precision, save :: gamma_const

  logical, save :: assume_neutral

contains

  subroutine actual_eos_init

    use extern_probin_module, only: eos_gamma, eos_assume_neutral

    implicit none
 
    ! constant ratio of specific heats
    if (eos_gamma .gt. 0.d0) then
       gamma_const = eos_gamma
    else
       gamma_const = FIVE3RD
    end if

    assume_neutral = eos_assume_neutral

  end subroutine actual_eos_init



  subroutine actual_eos(input, state)

    use fundamental_constants_module, only: k_B, n_A

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    double precision, parameter :: R = k_B*n_A

    double precision :: poverrho

    ! Calculate mu.
    
    if (assume_neutral) then
       state % mu = state % abar
    else
       state % mu = ONE / sum( (ONE + zion(:)) * state % xn(:) / aion(:) )
    endif    
    
    select case (input)

    case (eos_input_rt)

       ! dens, temp and xmass are inputs
       state % cv = R / (state % mu * (gamma_const-ONE)) 
       state % e = state % cv * state % T
       state % p = (gamma_const-ONE) * state % rho * state % e
       state % gam1 = gamma_const

    case (eos_input_rh)

       ! dens, enthalpy, and xmass are inputs

!       call bl_error('EOS: eos_input_rh is not supported in this EOS.')

    case (eos_input_tp)

       ! temp, pres, and xmass are inputs

!       call bl_error('EOS: eos_input_tp is not supported in this EOS.')
       
    case (eos_input_rp)

       ! dens, pres, and xmass are inputs

       poverrho = state % p / state % rho
       state % T = poverrho * state % mu * (ONE/R)
       state % e = poverrho * (ONE/(gamma_const-ONE))
       state % gam1 = gamma_const

    case (eos_input_re)

       ! dens, energy, and xmass are inputs

       poverrho = (gamma_const - ONE) * state % e

       state % p = poverrho * state % rho
       state % T = poverrho * state % mu * (ONE/R)
       state % gam1 = gamma_const
       
       ! sound speed
       state % cs = sqrt(gamma_const * poverrho)

       state % dpdr_e = poverrho
       state % dpde = (gamma_const-ONE) * state % rho

       ! Try to avoid the expensive log function.  Since we don't need entropy 
       ! in hydro solver, set it to an invalid but "nice" value for the plotfile.
       state % s = ONE  

    case (eos_input_ps)

       ! pressure entropy, and xmass are inputs

       call bl_error('EOS: eos_input_ps is not supported in this EOS.')
       
    case (eos_input_ph)

       ! pressure, enthalpy and xmass are inputs

       call bl_error('EOS: eos_input_ph is not supported in this EOS.')
       
    case (eos_input_th)

       ! temperature, enthalpy and xmass are inputs

       ! This system is underconstrained.
       
       call bl_error('EOS: eos_input_th is not a valid input for the gamma law EOS.')

    case default
       
       call bl_error('EOS: invalid input.')
       
    end select
    
  end subroutine actual_eos

end module actual_eos_module
