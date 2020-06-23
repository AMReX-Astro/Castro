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



  subroutine ca_fill_rotational_potential(lo,hi,phi,phi_lo,phi_hi,dx,time) &
       bind(C, name="ca_fill_rotational_potential")
    !
    ! .. note::
    !    Binds to C function ``ca_fill_rotational_potential``

    use prob_params_module, only: problo, center
    use amrex_constants_module, only: HALF

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: phi_lo(3), phi_hi(3)

    real(rt)        , intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt)        , intent(in   ) :: dx(3)
    real(rt), value , intent(in   ) :: time

    integer          :: i, j, k
    real(rt)         :: r(3)

    !$gpu

    do k = lo(3), hi(3)
       r(3) = problo(3) + dx(3)*(dble(k)+HALF) - center(3)

       do j = lo(2), hi(2)
          r(2) = problo(2) + dx(2)*(dble(j)+HALF) - center(2)

          do i = lo(1), hi(1)
             r(1) = problo(1) + dx(1)*(dble(i)+HALF) - center(1)

             phi(i,j,k) = rotational_potential(r,time)

          enddo
       enddo
    enddo

  end subroutine ca_fill_rotational_potential


  subroutine ca_fill_rotational_psi(lo, hi, &
                                    psi, psi_lo, psi_hi, &
                                    dx, time) bind(C, name='ca_fill_rotational_psi')
    ! Construct psi, which is the distance-related part
    ! of the rotation law. See e.g. Hachisu 1986a, Equation 15.
    ! For rigid-body rotation, psi = -R^2 / 2, where R is the
    ! distance orthogonal to the rotation axis. There are also
    ! v-constant and j-constant rotation laws that we do not
    ! implement here. We will construct this as potential / omega**2,
    ! so that the rotational_potential routine uniquely determines
    ! the rotation law. For the other rotation laws, we would
    ! simply divide by v_0^2 or j_0^2 instead.

    use amrex_constants_module, only: HALF
    use prob_params_module, only: problo, center
    use rotation_frequency_module, only: get_omega

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: psi_lo(3), psi_hi(3)
    real(rt), intent(inout) :: psi(psi_lo(1):psi_hi(1),psi_lo(2):psi_hi(2),psi_lo(3):psi_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: time

    integer  :: i, j, k
    real(rt) :: r(3), omega(3)

    !$gpu

    call get_omega(omega)

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             psi(i,j,k) = rotational_potential(r, time) / sum(omega**2)

          enddo
       enddo
    enddo

  end subroutine ca_fill_rotational_psi

end module rotation_module
