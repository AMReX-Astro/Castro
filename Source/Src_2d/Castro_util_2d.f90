module castro_util_2d_module

  implicit none

  public

contains

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

end module castro_util_2d_module

