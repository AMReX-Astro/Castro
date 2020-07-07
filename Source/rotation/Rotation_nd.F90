module rotation_module

  use meth_params_module, only: rotation_include_centrifugal, rotation_include_coriolis

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

  subroutine inertial_to_rotational_velocity(idx, time, v, idir)
    ! Given a velocity vector in the inertial frame, transform it to a
    ! velocity vector in the rotating frame.


    use prob_params_module, only: center
    use castro_util_module, only: position ! function
    use math_module, only: cross_product ! function
    use rotation_frequency_module, only: get_omega
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: idx(3)
    real(rt)        , intent(in   ) :: time
    real(rt)        , intent(inout) :: v(3)
    integer, intent(in), optional :: idir

    real(rt)         :: loc(3), omega(3)

    !$gpu

    if (present(idir)) then
       if (idir .eq. 1) then
          loc = position(idx(1),idx(2),idx(3),ccx=.false.) - center
       else if (idir .eq. 2) then
          loc = position(idx(1),idx(2),idx(3),ccy=.false.) - center
       else if (idir .eq. 3) then
          loc = position(idx(1),idx(2),idx(3),ccz=.false.) - center
       else
#ifndef AMREX_USE_GPU
          call castro_error("Error: unknown direction in inertial_to_rotational_velocity.")
#endif
       endif
    else
       loc = position(idx(1),idx(2),idx(3)) - center
    endif

    call get_omega(omega)

    v = v - cross_product(omega, loc)

  end subroutine inertial_to_rotational_velocity

 


  function rotational_potential(r, time) result(phi)
    ! Construct rotational potential, phi_R = -1/2 | omega x r |**2
    !

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: state_in_rotating_frame, rotation_include_centrifugal
    use math_module, only: cross_product ! function
    use rotation_frequency_module, only: get_omega
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt)         :: r(3), time
    real(rt)         :: phi

    real(rt)         :: omega(3), omegacrossr(3)

    !$gpu

    if (state_in_rotating_frame .eq. 1) then

       call get_omega(omega)

       phi = ZERO

       if (rotation_include_centrifugal == 1) then

          omegacrossr = cross_product(omega, r)

          phi = phi - HALF * dot_product(omegacrossr,omegacrossr)

       endif

    else

       phi = ZERO

    endif

  end function rotational_potential




end module rotation_module
