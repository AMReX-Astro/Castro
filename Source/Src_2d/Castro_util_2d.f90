module castro_util_2d_module

  implicit none

  public

contains

  subroutine ca_check_initial_species(lo,hi, &
                                      state,state_l1,state_l2,state_h1,state_h2) &
                                      bind(C, name="ca_check_initial_species")

    use network           , only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(2), hi(2)
    integer          :: state_l1,state_l2,state_h1,state_h2
    double precision :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

    ! Local variables
    integer          :: i,j,n
    double precision :: sum

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          sum = ZERO
          do n = 1, nspec
             sum = sum + state(i,j,UFS+n-1)
          end do
          if (abs(state(i,j,URHO)-sum).gt. 1.d-8 * state(i,j,URHO)) then
             print *,'Sum of (rho X)_i vs rho at (i,j): ',i,j,sum,state(i,j,URHO)
             call bl_error("Error:: Failed check of initial species summing to 1")
          end if

       enddo
    enddo

  end subroutine ca_check_initial_species


  
  subroutine ca_compute_avgstate (lo,hi,dx,dr,nc,&
                                  state,s_l1,s_l2,s_h1,s_h2,radial_state, &
                                  vol,v_l1,v_l2,v_h1,v_h2,radial_vol, &
                                  problo,numpts_1d) &
                                  bind(C, name="ca_compute_avgstate")

    use meth_params_module, only: URHO, UMX, UMY
    use prob_params_module, only: center
    use bl_constants_module

    implicit none

    integer          :: lo(2),hi(2),nc
    double precision :: dx(2),dr,problo(2)

    integer          :: numpts_1d
    double precision :: radial_vol(0:numpts_1d-1)
    double precision :: radial_state(nc,0:numpts_1d-1)

    integer          :: s_l1,s_l2,s_h1,s_h2
    double precision :: state(s_l1:s_h1,s_l2:s_h2,nc)

    integer          :: v_l1,v_l2,v_h1,v_h2
    double precision :: vol(v_l1:v_h1,v_l2:v_h2)

    integer          :: i,j,n,index
    double precision :: x,y,r
    double precision :: x_mom,y_mom,radial_mom

    do j = lo(2), hi(2)
       y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
       do i = lo(1), hi(1)
          x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
          r = sqrt(x**2  + y**2)
          index = int(r/dr)
          if (index .gt. numpts_1d-1) then
             print *,'COMPUTE_AVGSTATE:INDEX TOO BIG ',index,' > ',numpts_1d-1
             print *,'AT (i,j) ',i,j
             print *,'R / DR IS ',r,dr
             call bl_error("Error:: Castro_2d.f90 :: ca_compute_avgstate")
          end if

          ! Store the radial component of the momentum in both the UMX and UMY components for now.
          x_mom = state(i,j,UMX)
          y_mom = state(i,j,UMY)
          radial_mom = x_mom * (x/r) + y_mom * (y/r)
          radial_state(UMX,index) = radial_state(UMX,index) + vol(i,j)*radial_mom
          radial_state(UMY,index) = radial_state(UMY,index) + vol(i,j)*radial_mom

          ! Store all the other variables as themselves
          radial_state(URHO,index) = radial_state(URHO,index) + vol(i,j)*state(i,j,URHO)
          do n = UMY+1,nc
             radial_state(n,index) = radial_state(n,index) + vol(i,j)*state(i,j,n)
          end do
          radial_vol(index) = radial_vol(index) + vol(i,j)
       enddo
    enddo

  end subroutine ca_compute_avgstate



  subroutine get_center(center_out) bind(C, name="get_center")

    use prob_params_module, only : center

    implicit none

    double precision :: center_out(2)

    center_out(1:2) = center(1:2)

  end subroutine get_center



  subroutine set_center(center_in) bind(C, name="set_center")

    use prob_params_module, only : center

    implicit none

    double precision :: center_in(2)

    center(1:2) = center_in(1:2)

  end subroutine set_center



  subroutine find_center(data,new_center,icen,dx,problo) &
       bind(C, name="find_center")

    use bl_constants_module

    implicit none

    double precision :: data(-1:1,-1:1)
    double precision :: new_center(2)
    double precision :: dx(2),problo(2)
    double precision :: a,b,x,y,cen
    integer          :: icen(2)
    integer          :: i,j

    ! We do this to take care of precision issues
    cen = data(0,0)
    do j = -1,1
       do i = -1,1
          data(i,j) = data(i,j) - cen 
       end do
    end do

    !       This puts the center at the cell center
    new_center(1) = problo(1) +  (icen(1)+HALF) * dx(1)
    new_center(2) = problo(2) +  (icen(2)+HALF) * dx(2)

    ! Fit parabola y = a x^2  + b x + c through three points
    ! a = 1/2 ( y_1 + y_-1)
    ! b = 1/2 ( y_1 - y_-1)
    ! x_vertex = -b / 2a

    ! ... in x-direction
    a = HALF * (data(1,0) + data(-1,0)) - data(0,0)
    b = HALF * (data(1,0) - data(-1,0)) - data(0,0)
    x = -b / (TWO*a)
    new_center(1) = new_center(1) +  x*dx(1)

    ! ... in y-direction
    a = HALF * (data(0,1) + data(0,-1)) - data(0,0)
    b = HALF * (data(0,1) - data(0,-1)) - data(0,0)
    y = -b / (TWO*a)
    new_center(2) = new_center(2) +  y*dx(2)

  end subroutine find_center

end module castro_util_2d_module

