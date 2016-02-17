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

end module castro_util_3d_module
