! This is the equation of state for zero-temperature white dwarf 
! matter composed of degenerate electrons:
! P = A ( x * (2x**2 - 3)(x**2 + 1)**1/2 + 3 sinh**-1(x) )
! 
! where rho = B x**3 and the constants are given by:
!
! A = pi m_e**4 c**5 / (3 h**3) = 6.0 x 10^22 dyne cm**-2
! B = 8 pi m_e**3 c**3 mu_e m_p  / (3 h**3) = 9.8 x 10^5 mu_e g cm**-3
!
! The equation of state comes from Chandrasekhar (1935), and the enthalpy
! is calculated by Hachisu (1986):
!
! h = (8A / B) (1 + x**2)**(1/2)
!
! The internal energy is calculated using the standard relation:
! 
! h = e + P / rho

module eos_module

  use bl_types
  use bl_space
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module

  implicit none

  double precision, private :: A, B, B2
  double precision, parameter, private :: iter_tol = 1.d-10
  integer,          parameter, private :: max_iter = 1000
  !$OMP THREADPRIVATE(B)

  private :: enthalpy, pressure
  public  :: eos_init, eos

contains

  ! EOS initialization routine
  subroutine eos_init(small_temp, small_dens)

    use fundamental_constants_module, only: m_e, m_p, c_light, hplanck

    implicit none
 
    double precision, intent(in), optional :: small_temp
    double precision, intent(in), optional :: small_dens

    ! Small temperature and density parameters
 
    smallt = ZERO

    if (present(small_temp)) smallt = small_temp

    smalld = ZERO
 
    if (present(small_dens)) smalld = small_dens

    A = M_PI * m_e**4 * c_light**5 / (THREE * hplanck**3)
    B2 = EIGHT * M_PI * m_e**3 * c_light**3 * m_p  / (THREE * hplanck**3)

    initialized = .true.
 
  end subroutine eos_init



  !---------------------------------------------------------------------------
  ! The main interface
  !---------------------------------------------------------------------------
  subroutine eos(input, state, do_eos_diag_in, pt_index)

    implicit none

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state
    logical, optional, intent(in   ) :: do_eos_diag_in
    integer, optional, intent(in   ) :: pt_index(:)

    ! Local variables
    double precision :: dens, temp, enth, pres, eint, entr
    double precision :: x, dxdr

    integer :: n

    ! Make sure EOS is initialized before coming here.
    if (.not. initialized) call bl_error('EOS: not initialized')

    ! Make sure that the composition was set properly.

    do n = 1, nspec
      if (state % xn(n) .lt. init_test) call bl_error("EOS: species abundances not set.")
    enddo
    
    ! Calculate composition information

    call composition(state, .false.)

    dens = state % rho
    temp = state % T
    pres = state % p
    enth = state % h
    eint = state % e
    entr = state % s

    B = B2 * state % mu_e

    select case (input)

    !-------------------------------------------------------------------------
    ! Now do the calculations. In every case,
    ! make sure we have pressure, density, energy, and enthalpy.
    ! Relevant equations:
    ! rho = B x**3
    ! p   = A ( x * (2x**2 - 3)(x**2 + 1)**1/2 + 3 sinh**-1(x) )
    ! h   = (8A / B) * (1 + x**2)**1/2
    ! e   = h - p / rho
    !-------------------------------------------------------------------------

    case (eos_input_rh)

       ! dens, enthalpy, and xmass are inputs

       if (state % rho .lt. init_test .or. state % h .lt. init_test) then
         call bl_error("EOS called with rho and enthalpy as inputs, but these were not initialized.")
       endif

       ! Solve for the pressure and energy:

       x = (dens / B)**THIRD
       pres = pressure(x)
       eint = enth - pres / dens


    case (eos_input_rt)

       ! dens, temp, and xmass are inputs

       if (state % rho .lt. init_test .or. state % T .lt. init_test) then
         call bl_error("EOS called with rho and T as inputs, but these were not initialized.")
       endif
          
       ! Solve for the pressure, energy and enthalpy:

       x = (dens / B)**THIRD
       pres = pressure(x)
       enth = enthalpy(x)
       eint = enth - pres / dens


    case (eos_input_tp)

       ! temp, pres, and xmass are inputs

       if (state % T .lt. init_test .or. state % p .lt. init_test) then
         call bl_error("EOS called with temp and pressure as inputs, but these were not initialized.")
       endif

       ! Solve for the density, energy and enthalpy:

       pres = state % p
       call pres_iter(pres, dens)

       x = (dens / B)**THIRD
       enth = enthalpy(x)
       eint = enth - pres / dens
       

    case (eos_input_rp)

       ! dens, pres, and xmass are inputs

       if (state % rho .lt. init_test .or. state % p .lt. init_test) then
         call bl_error("EOS called with rho and pressure as inputs, but these were not initialized.")
       endif

       ! Solve for the enthalpy and energy:

       x = (dens / B)**THIRD
       enth = enthalpy(x)
       eint = enth - pres / dens


    case (eos_input_re)

       ! dens, energy, and xmass are inputs

       if (state % rho .lt. init_test .or. state % e .lt. init_test) then
         call bl_error('EOS called with rho and e as inputs, but these were not initialized.')
       endif

       ! Solve for the pressure and enthalpy:

       x = (dens / B)**THIRD
       pres = pressure(x)
       enth = enthalpy(x)


    case (eos_input_ps)
       
       ! pressure, entropy and xmass are inputs

       if (state % p .lt. init_test .or. state % s .lt. init_test) then
         call bl_error("EOS called with pressure and entropy as inputs, but these were not initialized.")
       endif

       ! Solve for the density, energy and enthalpy:

       pres = state % p
       call pres_iter(pres, dens)

       x = (dens / B)**THIRD
       enth = enthalpy(x)
       eint = enth - pres / dens


    case (eos_input_ph)

       ! pressure, enthalpy and xmass are inputs

       if (state % p .lt. init_test .or. state % s .lt. init_test) then
         call bl_error("EOS called with pressure and enthalpy as inputs, but these were not initialized.")
       endif

       ! Solve for the density and energy:

       x = ( ( (B * enth) / (EIGHT * A) )**2 - ONE )**HALF
       dens = B * x**3
       eint = enth - pres / dens


    case (eos_input_th)

       ! temperature, enthalpy and xmass are inputs

       if (state % T .lt. init_test .or. state % h .lt. init_test) then
         call bl_error("EOS called with temperature and enthalpy as inputs, but these were not initialized.")
       endif

       ! Solve for the density, energy and pressure:

       x = ( ( (B * enth) / (EIGHT * A) )**2 - ONE )**HALF
       dens = B * x**3
       pres = pressure(x)
       eint = enth - pres / dens


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

    ! All temperature derivatives are zero since the gas is temperature-independent.

    state % dPdT = ZERO
    state % dhdT = ZERO
    state % dedT = ZERO
    state % dsdT = ZERO

    ! Density derivatives are computed using the chain rule, e.g. dPdr = dPdx * dxdr.

    x = (dens / B)**THIRD
    dxdr = THIRD * x / dens

    state % dPdr = dxdr * dpdx(x)
    state % dhdr = dxdr * dhdx(x)
    state % dedr = state % dhdr - state % dpdr / state % rho + state % p / (state % rho)**2
    state % dsdr = ZERO

    ! Heat capacities are zero: the gas properties don't change when the temperature changes.

    state % cv = ZERO
    state % cp = ZERO

    ! Adiabatic gamma_1 == d(log p) / d(log rho) |_s.
 
    state % gam1 = state % dpdr * (state % rho / state % p)

    ! Compute dPdX, dedX, dhdX.

    state % dpdA = - state % p / state % abar
    state % dpdZ =   state % p / (ONE + state % zbar)

    state % dedA = - state % e / state % abar
    state % dedZ =   state % e / (ONE + state % zbar)

    call composition_derivatives(state, .false.)

    ! Sound speed.

    state % cs = sqrt(state % dpdr)

  end subroutine eos



  double precision function pressure(x)

    implicit none

    double precision, intent(in)  :: x

    pressure = A * ( x * (TWO * x**2 - THREE) * (x**2 + ONE)**HALF + THREE * asinh(x) )

  end function pressure



  double precision function enthalpy(x)

    implicit none

    double precision, intent(in)  :: x

    enthalpy = (EIGHT * A / B) * (ONE + x**2)**HALF

  end function enthalpy



  double precision function dpdx(x)

    implicit none

    double precision, intent(in) :: x

    dpdx = A * ( (TWO * x**2 - THREE)*(x**2 + ONE)**HALF + &
                 x * (4*x) * (x**2 + ONE)**HALF + &
                 x**2 * (TWO * x**2 - THREE) * (x**2 + ONE)**(-HALF) + &
                 THREE * (x**2 + ONE)**(-HALF) )

  end function dpdx



  double precision function dhdx(x)

    implicit none

    double precision, intent(in) :: x

    dhdx = enthalpy(x) * (x / (x**2 + ONE))

  end function dhdx



  subroutine pres_iter(pres, dens)

    implicit none

    double precision, intent(inout) :: pres, dens

    double precision :: x, dx
    integer          :: iter

    ! Starting guess for the iteration.

    x = ONE

    ! We are solving the equation:
    ! f(x) = p_want - p(x) = 0.
    ! Then we can use Newton's method, with dfdx = -dpdx.
    ! We iterate until the density change is close enough to zero.

    do iter = 1, max_iter

       dx = (pres - pressure(x)) / dpdx(x)

       x  = x + dx

       if ( abs(dx) / x .lt. iter_tol ) then
          exit
       endif

    enddo

    if (iter .eq. max_iter) then
       call bl_error("EOS: pres_iter failed to converge.")
    endif

    dens = B * x**3

  end subroutine pres_iter



end module eos_module
