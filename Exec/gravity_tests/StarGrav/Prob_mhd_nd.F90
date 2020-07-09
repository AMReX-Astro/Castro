
! :::
! ::: --------------------------------------------------------------------
! :::
subroutine ca_initmag(level, time, lo, hi, &
                      nbx, mag_x, bx_lo, bx_hi, &
                      nby, mag_y, by_lo, by_hi, &
                      nbz, mag_z, bz_lo, bz_hi, &
                      dx, xlo, xhi)

  use probdata_module
  use prob_params_module
  use model_parser_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer :: level, nbx, nby, nbz
  integer :: lo(3), hi(3)
  integer :: bx_lo(3), bx_hi(3)
  integer :: by_lo(3), by_hi(3)
  integer :: bz_lo(3), bz_hi(3)
  real(rt) :: xlo(3), xhi(3), time, dx(3)

  real(rt) :: mag_x(bx_lo(1):bx_hi(1), bx_lo(2):bx_hi(2), bx_lo(3):bx_hi(3), nbx)
  real(rt) :: mag_y(by_lo(1):by_hi(1), by_lo(2):by_hi(2), by_lo(3):by_hi(3), nby)
  real(rt) :: mag_z(bz_lo(1):bz_hi(1), bz_lo(2):bz_hi(2), bz_lo(3):bz_hi(3), nbz)

  real(rt) :: A_x(lo(1):hi(1), lo(2):hi(2)+1, lo(3):hi(3)+1)
  real(rt) :: A_y(lo(1):hi(1)+1, lo(2):hi(2), lo(3):hi(3)+1)
  real(rt) :: x, y, z, dist
  real(rt) :: theta, phi
  integer  :: i, j, k

  !get the vector potential A


  !get the A_x component
  do k = lo(3), hi(3)+1
     z = problo(3) + dx(3)*(dble(k)) - center(3)
     do j = lo(2), hi(2)+1
        y = problo(2) + dx(2)*(dble(j)) - center(2)
        do i = lo(1), hi(1)
           x = problo(1) + dx(1)*(dble(i)) - center(1)

           theta = atan2(sqrt(x**2+y**2),z)
           phi = atan2(y,x) 
 

           dist = sqrt(x**2 + y**2 + z**2)


           !inside the star radius treshold (R) the vector potential
           ! will be of the form A_x = constant*x 
           !                         = constant*rsin(theta)cos(phi)
           ! the constant is then obtained by matching B.C with
           ! the field outside at r = R  
           if (dist .le. 1.0d8 ) then
            
               A_x(i,j,k) = -m_0*dist*sin(theta)*sin(phi)/(1.0d24)

           else
           !outside the star rarius treshold
           ! x component of dipole vector potential        
               A_x(i,j,k) = -m_0*sin(theta)*sin(phi)/(dist*dist)
           endif

        end do
     end do
  end do

  !get the A_y component
  do k = lo(3), hi(3)+1
     z = problo(3) + dx(3)*(dble(k)) - center(3)
     do j = lo(2), hi(2)
        y = problo(2) + dx(2)*(dble(j)) - center(2)
        do i = lo(1), hi(1)+1
           x = problo(1) + dx(1)*(dble(i)) - center(1)

           theta = atan2(sqrt(x**2+y**2),z)
           phi = atan2(y,x) 

           dist = sqrt(x**2 + y**2 + z**2)

           !inside the star, doing the analogous as we did 
           ! for A_x
           if (dist .le. 1.0d8) then
               A_y(i,j,k) = m_0*dist*sin(theta)*cos(phi)/(1.0d24)
           else
           !outside the star, the y component of
           !a dipole vector potential  
               A_y(i,j,k) = m_0*sin(theta)*cos(phi)/(dist*dist)
           endif

        end do
     end do
  end do


 !since A is of the form A = (A_x, A_y, 0)
 ! B_x = -dA_x/dz
 ! B_y = dA_y/dz
 ! B_z = dA_y/dx - dA_x/dy 

 !Initialize the magnetic fields 
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           
           mag_x(i,j,k,1) = (A_y(i,j,k) - A_y(i,j,k+1))/dx(3)

        end do
     end do
  end do

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)

           mag_y(i,j,k,1) = (A_x(i,j,k+1) - A_x(i,j,k))/dx(3)

        end do
     end do
  end do

  do k = lo(3), hi(3)+1
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           mag_z(i,j,k,1) = (A_y(i+1,j,k) - A_y(i,j,k))/dx(1) &
                          - (A_x(i,j+1,k) - A_x(i,j,k))/dx(2) 

        end do
     end do
  end do


end subroutine ca_initmag

