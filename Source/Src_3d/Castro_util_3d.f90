module castro_util_3d_module

  implicit none

  public

contains
  
  subroutine ca_compute_avgstate(lo,hi,dx,dr,nc,&
                                 state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,radial_state, &
                                 vol,v_l1,v_l2,v_l3,v_h1,v_h2,v_h3,radial_vol, &
                                 problo,numpts_1d) &
                                 bind(C, name="ca_compute_avgstate")

    use meth_params_module, only : URHO, UMX, UMY, UMZ
    use prob_params_module, only : center
    use bl_constants_module

    implicit none

    integer          :: lo(3),hi(3),nc
    double precision :: dx(3),dr,problo(3)

    integer          :: numpts_1d
    double precision :: radial_state(nc,0:numpts_1d-1)
    double precision :: radial_vol(0:numpts_1d-1)

    integer          :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
    double precision :: state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,nc)

    integer          :: v_l1,v_l2,v_l3,v_h1,v_h2,v_h3
    double precision :: vol(v_l1:v_h1,v_l2:v_h2,v_l3:v_h3)

    integer          :: i,j,k,n,index
    double precision :: x,y,z,r
    double precision :: x_mom,y_mom,z_mom,radial_mom
    !
    ! Do not OMP this.
    !
    do k = lo(3), hi(3)
       z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
             r = sqrt(x**2 + y**2 + z**2)
             index = int(r/dr)
             if (index .gt. numpts_1d-1) then
                print *,'COMPUTE_AVGSTATE: INDEX TOO BIG ',index,' > ',numpts_1d-1
                print *,'AT (i,j,k) ',i,j,k
                print *,'R / DR ',r,dr
                call bl_error("Error:: Castro_3d.f90 :: ca_compute_avgstate")
             end if
             radial_state(URHO,index) = radial_state(URHO,index) &
                  + vol(i,j,k)*state(i,j,k,URHO)
             !
             ! Store the radial component of the momentum in the 
             ! UMX, UMY and UMZ components for now.
             !
             x_mom = state(i,j,k,UMX)
             y_mom = state(i,j,k,UMY)
             z_mom = state(i,j,k,UMZ)
             radial_mom = x_mom * (x/r) + y_mom * (y/r) + z_mom * (z/r)
             radial_state(UMX,index) = radial_state(UMX,index) + vol(i,j,k)*radial_mom
             radial_state(UMY,index) = radial_state(UMY,index) + vol(i,j,k)*radial_mom
             radial_state(UMZ,index) = radial_state(UMZ,index) + vol(i,j,k)*radial_mom

             do n = UMZ+1,nc
                radial_state(n,index) = radial_state(n,index) + vol(i,j,k)*state(i,j,k,n)
             end do
             radial_vol(index) = radial_vol(index) + vol(i,j,k)
          enddo
       enddo
    enddo

  end subroutine ca_compute_avgstate



  subroutine get_center(center_out) bind(C, name="get_center")

    use prob_params_module, only : center

    implicit none

    double precision, intent(inout) :: center_out(3)

    center_out(1:3) = center(1:3)

  end subroutine get_center


  
  subroutine set_center(center_in) bind(C, name="set_center")

    use prob_params_module, only : center

    implicit none

    double precision :: center_in(3)

    center(1:3) = center_in(1:3)

  end subroutine set_center



  subroutine find_center(data,new_center,icen,dx,problo) &
       bind(C, name="find_center")

    use bl_constants_module  

    implicit none

    double precision :: data(-1:1,-1:1,-1:1)
    double precision :: new_center(3)
    double precision :: dx(3),problo(3)
    double precision :: a,b,x,y,z,cen
    integer          :: icen(3)
    integer          :: i,j,k

    ! We do this to take care of precision issues
    cen = data(0,0,0)
    do k = -1,1
       do j = -1,1
          do i = -1,1
             data(i,j,k) = data(i,j,k) - cen 
          end do
       end do
    end do

    !       This puts the "center" at the cell center
    new_center(1) = problo(1) +  (icen(1)+HALF) * dx(1)
    new_center(2) = problo(2) +  (icen(2)+HALF) * dx(2)
    new_center(3) = problo(3) +  (icen(3)+HALF) * dx(3)

    ! Fit parabola y = a x^2  + b x + c through three points
    ! a = 1/2 ( y_1 + y_-1)
    ! b = 1/2 ( y_1 - y_-1)
    ! x_vertex = -b / 2a

    ! ... in x-direction
    a = HALF * (data(1,0,0) + data(-1,0,0)) - data(0,0,0)
    b = HALF * (data(1,0,0) - data(-1,0,0)) - data(0,0,0)
    x = -b / (TWO*a)
    new_center(1) = new_center(1) +  x*dx(1)

    ! ... in y-direction
    a = HALF * (data(0,1,0) + data(0,-1,0)) - data(0,0,0)
    b = HALF * (data(0,1,0) - data(0,-1,0)) - data(0,0,0)
    y = -b / (TWO*a)
    new_center(2) = new_center(2) +  y*dx(2)

    ! ... in z-direction
    a = HALF * (data(0,0,1) + data(0,0,-1)) - data(0,0,0)
    b = HALF * (data(0,0,1) - data(0,0,-1)) - data(0,0,0)
    z = -b / (TWO*a)
    new_center(3) = new_center(3) +  z*dx(3)

  end subroutine find_center

end module castro_util_3d_module
