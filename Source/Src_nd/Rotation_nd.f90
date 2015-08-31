module rotation_module

  implicit none

  private

  public cross_product, get_omega, rotational_acceleration, rotational_potential

contains

  function get_omega() result(omega)

    use prob_params_module, only: coord_type
    use meth_params_module, only: rot_period, rot_axis
    use bl_constants_module, only: ZERO, TWO, M_PI

    implicit none

    double precision :: omega(3)

    ! If rot_period is less than zero, that means rotation is disabled, and so we should effectively
    ! shut off the source term by setting omega = 0. Note that by default rot_axis == 3 for Cartesian
    ! coordinates and rot_axis == 2 for cylindrical coordinates.

    omega = ZERO

    if (coord_type == 0 .or. coord_type == 1) then

       if (rot_period > ZERO) then
          omega(rot_axis) = TWO * M_PI / rot_period
       endif

    else

       call bl_error("Error:: rotation_nd.f90 :: invalid coord_type")

    endif

  end function get_omega



  ! Given a position and velocity, calculate 
  ! the rotational acceleration.
  ! This is the Coriolis force  (-2 omega x v) 
  ! and the centrifugal force (- omega x ( omega x r)).
  ! For CASTRO we assume that the rotation rate is 
  ! unchanging so d(omega)/dt = 0.

  function rotational_acceleration(r, v) result(Sr)

    use bl_constants_module, only: ZERO, TWO

    implicit none

    double precision :: r(3), v(3)
    double precision :: Sr(3)

    double precision :: omega(3), domega_dt(3), omegacrossr(3), omegacrossv(3)

    omega = get_omega()

    domega_dt = ZERO

    omegacrossr = cross_product(omega,r)
    omegacrossv = cross_product(omega,v)

    Sr = -TWO * omegacrossv - cross_product(omega, omegacrossr) - cross_product(domega_dt, r)

  end function rotational_acceleration



  ! Construct rotational potential, phi_R = -1/2 | omega x r |**2

  function rotational_potential(r) result(phi)

    use bl_constants_module, only: HALF

    implicit none

    double precision :: r(3)
    double precision :: phi

    double precision :: omega(3), omegacrossr(3)

    omega = get_omega()

    omegacrossr = cross_product(omega, r)

    phi = - HALF * dot_product(omegacrossr,omegacrossr)

  end function rotational_potential



  ! Compute the standard cross-product of two three-vectors.

  function cross_product(A,B) result(C)

    implicit none

    double precision :: A(3), B(3)
    double precision :: C(3)

    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1) - A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)

  end function cross_product

end module rotation_module



  subroutine ca_fill_rotational_potential(lo,hi,phi,phi_lo,phi_hi,dx)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ
    use prob_params_module, only: problo, center
    use rotation_module, only: rotational_potential
    use bl_constants_module, only: HALF

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: phi_lo(3), phi_hi(3)

    double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    double precision, intent(in   ) :: dx(3)

    integer          :: i, j, k
    double precision :: r(3)

    do k = lo(3), hi(3)
       r(3) = problo(3) + dx(3)*(dble(k)+HALF) - center(3)       

       do j = lo(2), hi(2)
          r(2) = problo(2) + dx(2)*(dble(j)+HALF) - center(2)

          do i = lo(1), hi(1)
             r(1) = problo(1) + dx(1)*(dble(i)+HALF) - center(1)

             phi(i,j,k) = rotational_potential(r)
             
          enddo
       enddo
    enddo

  end subroutine ca_fill_rotational_potential



  subroutine ca_fill_rotational_acceleration(lo,hi,rot,rot_lo,rot_hi,state,state_lo,state_hi,dx)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ
    use prob_params_module, only: problo, center
    use rotation_module, only: rotational_acceleration
    use bl_constants_module, only: HALF

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: rot_lo(3), rot_hi(3)
    integer         , intent(in   ) :: state_lo(3), state_hi(3)

    double precision, intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3),3)
    double precision, intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    double precision, intent(in   ) :: dx(3)

    integer          :: i, j, k
    double precision :: r(3)

    do k = lo(3), hi(3)
       r(3) = problo(3) + dx(3)*(dble(k)+HALF) - center(3)       

       do j = lo(2), hi(2)
          r(2) = problo(2) + dx(2)*(dble(j)+HALF) - center(2)

          do i = lo(1), hi(1)
             r(1) = problo(1) + dx(1)*(dble(i)+HALF) - center(1)

             rot(i,j,k,:) = rotational_acceleration(r, state(i,j,k,UMX:UMZ) / state(i,j,k,URHO))
             
          enddo
       enddo
    enddo

  end subroutine ca_fill_rotational_acceleration
