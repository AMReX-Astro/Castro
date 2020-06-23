module rotation_module

  use meth_params_module, only: rotation_include_centrifugal, rotation_include_coriolis

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains


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
