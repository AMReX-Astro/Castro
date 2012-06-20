  subroutine ca_rotate(lo,hi, &
                       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                       rot_src,rot_src_l1,rot_src_l2,rot_src_l3,rot_src_h1,rot_src_h2,rot_src_h3, &
                       problo,dx,rot_freq)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN
    use probdata_module, only: center
    use prob_params_module, only: coord_type

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
    integer         , intent(in   ) :: rot_src_l1,rot_src_l2,rot_src_l3,rot_src_h1,rot_src_h2,rot_src_h3
    double precision, intent(in   ) :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
    double precision, intent(inout) :: rot_src(rot_src_l1:rot_src_h1,rot_src_l2:rot_src_h2,rot_src_l3:rot_src_h3,3+1)
    double precision, intent(in   ) :: problo(3), dx(3)
    double precision, intent(in   ) :: rot_freq

    integer          :: i,j,k
    double precision :: x,y,z,r(3)
    double precision :: v(3),omega(3)
    double precision :: dens
    double precision :: vdotr,omegadotr,omegadotv,omegacrossv(3),omega2

    double precision :: TWO_PI
    parameter (TWO_PI = 2.d0 * 3.1415926535d0)

    if (coord_type == 0) then
       omega = (/ 0.0d0, 0.0d0, TWO_PI*rot_freq /)
    else
       call bl_error("Error:: Rotate_2d.f90 :: unknown coord_type")
    endif

    omega2 = dot_product(omega,omega)
    
    do k = lo(3), hi(3)
       z = problo(3) + dx(3)*(float(k)+0.5d0) - center(3)
       do j = lo(2), hi(2)
          y = problo(2) + dx(2)*(float(j)+0.5d0) - center(2)
          do i = lo(1), hi(1)
             x = problo(1) + dx(1)*(float(i)+0.5d0) - center(1)

             r = (/ x, y, z /)

             dens = state(i,j,k,URHO)
             
             v = (/ state(i,j,k,UMX)/dens, &
                    state(i,j,k,UMY)/dens, &
                    state(i,j,k,UMZ)/dens /)

             omegacrossv = cross_product(omega,v)
             omegadotr   = dot_product(omega,r)
             omegadotv   = dot_product(omega,v)
             vdotr       = dot_product(v,r)

             rot_src(i,j,k,1:3) = -2.0d0 * dens * omegacrossv(:) - &
                  dens * (omegadotr * omega(:) - omega2 * r(:))
             rot_src(i,j,k,4) = -dens * omegadotv * omegadotr + &
                  dens * omega2 * vdotr

          enddo
       enddo
    enddo

    contains
      function cross_product(x,y) result(r)
        
        double precision :: x(3), y(3)
        double precision :: r(3)

        r(1) = x(2)*y(3) - x(3)*y(2)
        r(2) = x(3)*y(1) - x(1)*y(3)
        r(3) = x(1)*y(2) - x(2)*y(1)
      end function cross_product

    end subroutine ca_rotate
