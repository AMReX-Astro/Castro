module rotation_frequency_module

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine get_omega(time, omega) bind(C, name="get_omega")

    use prob_params_module, only: coord_type
    use meth_params_module, only: rot_period, rot_period_dot, rot_axis
    use amrex_constants_module, only: ZERO, TWO, M_PI

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt), intent(in), value :: time
    real(rt), intent(inout) :: omega(3)

    real(rt) :: curr_period

    ! If rot_period is less than zero, that means rotation is disabled, and so we should effectively
    ! shut off the source term by setting omega = 0. Note that by default rot_axis == 3 for Cartesian
    ! coordinates and rot_axis == 2 for cylindrical coordinates.

    !$gpu

    omega(:) = ZERO

    if (coord_type == 0 .or. coord_type == 1) then

       if (rot_period > ZERO) then

          ! If we have a time rate of change of the rotational period,
          ! adjust it accordingly in the calculation of omega. We assume
          ! that the change has been linear and started at t == 0.

          curr_period = rot_period + rot_period_dot * time

          omega(rot_axis) = TWO * M_PI / curr_period

       endif

    else
#ifndef AMREX_USE_GPU
       call castro_error("Error:: rotation_nd.f90 :: invalid coord_type")
#endif
    endif

  end subroutine get_omega



  function get_domegadt(time) result(domegadt)

    use prob_params_module, only: coord_type
    use meth_params_module, only: rot_period, rot_period_dot
    use amrex_constants_module, only: ZERO, TWO, M_PI

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: time
    real(rt)         :: domegadt(3)

    real(rt)         :: curr_period, curr_omega(3)

    !$gpu

    domegadt(:) = ZERO

    if (coord_type == 0 .or. coord_type .eq. 1) then

       if (rot_period > ZERO) then

          ! Rate of change of the rotational frequency is given by
          ! d( ln(period) ) / dt = - d( ln(omega) ) / dt

          curr_period = rot_period + rot_period_dot * time
          call get_omega(time, curr_omega)

          domegadt = -curr_omega * (rot_period_dot / curr_period)

       endif

    else
#ifndef AMREX_USE_GPU
       call castro_error("Error:: rotation_nd.f90 :: unknown coord_type")
#endif
    endif

  end function get_domegadt

  subroutine set_rot_period(period) bind(C, name='set_rot_period')

    use meth_params_module, only: rot_period

    implicit none

    real(rt), intent(in) :: period

    rot_period = period

  end subroutine set_rot_period
  
end module rotation_frequency_module
