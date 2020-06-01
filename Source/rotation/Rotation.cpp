void
Castro::rotational_acceleration(GpuArray<Real, 3>& r, GpuArray<Real, 3>& v,
                                GpuArray<Real, 3>& omega, GpuArray<Real, 3>& domega_dt.
                                const Real time, const bool coriolis, Real* Sr) {

  // Given a position and velocity, calculate
  // the rotational acceleration. This is the sum of:
  // the Coriolis force (-2 omega x v),
  // the centrifugal force (- omega x ( omega x r)),
  // and a changing rotation rate (-d(omega)/dt x r).

  Sr[0] = 0.0;
  Sr[1] = 0.0;
  Sr[2] = 0.0;

  if (state_in_rotating_frame == 1) {

    // Allow the various terms to be turned off.  This is often used
    // for diagnostic purposes, but there are genuine science cases
    // for disabling certain terms in some cases (in particular, when
    // obtaining a system in rotational equilibrium through a
    // relaxation process involving damping or sponging, one may want
    // to turn off the Coriolis force during the relaxation process,
    // on the basis that the true equilibrium state will have zero
    // velocity anyway).

    bool c1 = (rotation_include_centrifugal == 1) ? true : false;

    bool c2 = (rotation_include_coriolis == 1) ? true : false;

    if (! coriolis) {
      c2 = false;
    }

    bool c3 = (rotation_include_domegadt == 1) ? true : false;

    GpuArray<Real, 3> omega_cross_v;
    cross_product(omega, v, omega_cross_v);

    if (c1) {
      GpuArray<Real, 3> omega_cross_r;
      cross_product(omega, r, omega_cross_r);

      GpuArray<Real, 3> omega_cross_omega_cross_r;
      cross_product(omega, omega_cross_r, omega_cross_omega_cross_r);

      for (int idir = 0; idir < 3; idir++) {
        Sr[idir] -= omega_cross_omega_cross_r[idir];
      }
    }

    if (c2) {
      for (int idir = 0; idir < 3; idir++) {
        Sr[idir] -= 2.0_rt * omega_cross_v[idir];
      }
    }

    if (c3) {
      GpuArray<Real, 3> domega_dt_cross_r;
      cross_product(domega_dt, r, domega_dt_cross_r);

      for (int idir = 0; idir < 3; idir++) {
        Sr[idir] -= domega_dt_cross_r[idir];
      }
    }

  } else {

    // The source term for the momenta when we're not measuring state
    // variables in the rotating frame is not strictly the traditional
    // Coriolis force, but we'll still allow it to be disabled with
    // the same parameter.

    bool c2 = (rotation_include_coriolis == 1) ? true : false;

    if (! coriolis) {
      c2 = false;
    }

    if (c2) {
      GpuArray<Real, 3> omega_cross_v;
      cross_product(omega, v, omega_cross_v);

      for (int idir = 0; idir < 3; idir++) {
        Sr[idir] -= omega_cross_v[idir];
      }
    }

  }
}

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

       call get_omega(time, omega)

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

    !$gpu

    call get_omega(time, omega)

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
