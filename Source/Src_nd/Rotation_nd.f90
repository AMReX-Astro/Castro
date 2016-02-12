module rotation_module

  use math_module, only: cross_product
  use rotation_frequency_module, only: get_omega, get_domegadt

  implicit none

  private

  public rotational_acceleration, rotational_potential

contains

  ! Given a position and velocity, calculate 
  ! the rotational acceleration. This is the sum of:
  ! the Coriolis force (-2 omega x v),
  ! the centrifugal force (- omega x ( omega x r)),
  ! and a changing rotation rate (-d(omega)/dt x r).

  function rotational_acceleration(r, v, time) result(Sr)

    use bl_constants_module, only: ZERO, TWO

    implicit none

    double precision :: r(3), v(3), time
    double precision :: Sr(3)

    double precision :: omega(3), domega_dt(3), omegacrossr(3), omegacrossv(3)

    omega = get_omega(time)

    domega_dt = get_domegadt(time)

    omegacrossr = cross_product(omega,r)
    omegacrossv = cross_product(omega,v)

    Sr = -TWO * omegacrossv - cross_product(omega, omegacrossr) - cross_product(domega_dt, r)

  end function rotational_acceleration



  ! Construct rotational potential, phi_R = -1/2 | omega x r |**2

  function rotational_potential(r, time) result(phi)

    use bl_constants_module, only: HALF

    implicit none

    double precision :: r(3), time
    double precision :: phi

    double precision :: omega(3), omegacrossr(3)

    omega = get_omega(time)

    omegacrossr = cross_product(omega, r)

    phi = HALF * dot_product(omegacrossr,omegacrossr)

  end function rotational_potential



  subroutine ca_fill_rotational_potential(lo,hi,phi,phi_lo,phi_hi,dx,time) &
       bind(C, name="ca_fill_rotational_potential")

    use meth_params_module, only: NVAR, URHO, UMX, UMZ
    use prob_params_module, only: problo, center
    use bl_constants_module, only: HALF

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: phi_lo(3), phi_hi(3)

    double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    double precision, intent(in   ) :: dx(3), time

    integer          :: i, j, k
    double precision :: r(3)

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

    use meth_params_module, only: NVAR, URHO, UMX, UMZ
    use prob_params_module, only: problo, center
    use bl_constants_module, only: HALF

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: rot_lo(3), rot_hi(3)
    integer         , intent(in   ) :: state_lo(3), state_hi(3)

    double precision, intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3),3)
    double precision, intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    double precision, intent(in   ) :: dx(3), time

    integer          :: i, j, k
    double precision :: r(3)

    do k = lo(3), hi(3)
       r(3) = problo(3) + dx(3)*(dble(k)+HALF) - center(3)

       do j = lo(2), hi(2)
          r(2) = problo(2) + dx(2)*(dble(j)+HALF) - center(2)

          do i = lo(1), hi(1)
             r(1) = problo(1) + dx(1)*(dble(i)+HALF) - center(1)

             rot(i,j,k,:) = rotational_acceleration(r, state(i,j,k,UMX:UMZ) / state(i,j,k,URHO), time)

          enddo
       enddo
    enddo

  end subroutine ca_fill_rotational_acceleration

end module rotation_module
