module rotation_module

  use meth_params_module, only: rotation_include_centrifugal, rotation_include_coriolis, &
       rotation_include_domegadt

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
    use rotation_frequency_module, only: get_omega ! function

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

    omega = get_omega(time)

    v = v - cross_product(omega, loc)

  end subroutine inertial_to_rotational_velocity



  function rotational_acceleration(r, v, time, centrifugal, coriolis, domegadt) result(Sr)
    ! Given a position and velocity, calculate
    ! the rotational acceleration. This is the sum of:
    ! the Coriolis force (-2 omega x v),
    ! the centrifugal force (- omega x ( omega x r)),
    ! and a changing rotation rate (-d(omega)/dt x r).
    !

    use amrex_constants_module, only: ZERO, TWO
    use meth_params_module, only: state_in_rotating_frame
    use math_module, only: cross_product ! function
    use rotation_frequency_module, only: get_omega ! function
    use rotation_frequency_module, only: get_domegadt ! function

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: r(3), v(3), time
    real(rt)         :: Sr(3)

    real(rt)         :: omega(3), domega_dt(3), omegacrossr(3), omegacrossv(3)

    logical, optional :: centrifugal, coriolis, domegadt
    logical :: c1, c2, c3

    !$gpu

    omega = get_omega(time)

    if (state_in_rotating_frame .eq. 1) then

       ! Allow the various terms to be turned off.
       ! This is often used for diagnostic purposes,
       ! but there are genuine science cases for
       ! disabling certain terms in some cases (in
       ! particular, when obtaining a system in
       ! rotational equilibrium through a relaxation
       ! process involving damping or sponging, one
       ! may want to turn off the Coriolis force
       ! during the relaxation process, on the
       ! basis that the true equilibrium state
       ! will have zero velocity anyway).

       if (rotation_include_centrifugal == 1) then
          c1 = .true.
       else
          c1 = .false.
       endif

       if (present(centrifugal)) then
          if (.not. centrifugal) c1 = .false.
       endif

       if (rotation_include_coriolis == 1) then
          c2 = .true.
       else
          c2 = .false.
       endif

       if (present(coriolis)) then
          if (.not. coriolis) c2 = .false.
       endif

       if (rotation_include_domegadt == 1) then
          c3 = .true.
       else
          c3 = .false.
       endif

       if (present(domegadt)) then
          if (.not. domegadt) c3 = .false.
       endif

       domega_dt = get_domegadt(time)

       omegacrossr = cross_product(omega,r)
       omegacrossv = cross_product(omega,v)

       Sr = ZERO

       if (c1) then
          Sr = Sr - cross_product(omega, omegacrossr)
       endif

       if (c2) then
          Sr = Sr - TWO * omegacrossv
       endif

       if (c3) then
          Sr = Sr - cross_product(domega_dt, r)
       endif

    else

       ! The source term for the momenta when we're not measuring
       ! state variables in the rotating frame is not strictly the
       ! traditional Coriolis force, but we'll still allow it to
       ! be disabled with the same parameter.

       if (rotation_include_coriolis == 1) then
          c2 = .true.
       else
          c2 = .false.
       endif

       if (present(coriolis)) then
          if (.not. coriolis) c2 = .false.
       endif

       Sr = ZERO

       if (c2) then
          Sr = Sr - cross_product(omega, v)
       endif

    endif

  end function rotational_acceleration




  function rotational_potential(r, time) result(phi)
    ! Construct rotational potential, phi_R = -1/2 | omega x r |**2
    !

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: state_in_rotating_frame, rotation_include_centrifugal
    use math_module, only: cross_product ! function
    use rotation_frequency_module, only: get_omega ! function

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: r(3), time
    real(rt)         :: phi

    real(rt)         :: omega(3), omegacrossr(3)

    !$gpu

    if (state_in_rotating_frame .eq. 1) then

       omega = get_omega(time)

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


  subroutine ca_fill_rotational_acceleration(lo,hi,rot,rot_lo,rot_hi,state,state_lo,state_hi,dx,time) &
       bind(C, name="ca_fill_rotational_acceleration")
    !
    ! .. note::
    !    Binds to C function ``ca_fill_rotational_acceleration``

    use meth_params_module, only: NVAR, URHO, UMX, UMZ
    use prob_params_module, only: problo, center
    use amrex_constants_module, only: HALF

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: rot_lo(3), rot_hi(3)
    integer         , intent(in   ) :: state_lo(3), state_hi(3)

    real(rt)        , intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3),3)
    real(rt)        , intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    real(rt)        , intent(in   ) :: dx(3)
    real(rt), value , intent(in   ) :: time

    integer          :: i, j, k
    real(rt)         :: r(3), v(3)

    !$gpu

    do k = lo(3), hi(3)
       r(3) = problo(3) + dx(3)*(dble(k)+HALF) - center(3)

       do j = lo(2), hi(2)
          r(2) = problo(2) + dx(2)*(dble(j)+HALF) - center(2)

          do i = lo(1), hi(1)
             r(1) = problo(1) + dx(1)*(dble(i)+HALF) - center(1)

             v(:) = state(i,j,k,UMX:UMZ) / state(i,j,k,URHO)

             rot(i,j,k,:) = rotational_acceleration(r, v, time)

          enddo
       enddo
    enddo

  end subroutine ca_fill_rotational_acceleration



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

    omega = get_omega(time)

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
