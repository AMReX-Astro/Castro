   subroutine ca_initmag(level, time, lo, hi, &
                      nbx, mag_x, bx_lo, bx_hi, &
                      nby, mag_y, by_lo, by_hi, &
                      nbz, mag_z, bz_lo, bz_hi, &
                      dx, xlo, xhi)

   use probdata_module, only : B_x, B_y, B_z, m_0, center_P_initial, center_S_initial, & 
                               mass_P, mass_S, vel_P, vel_S
   use prob_params_module
   use amrex_fort_module, only : rt => amrex_real
   use castro_util_module, only: position
   use amrex_constants_module, only: HALF, ZERO
   use wdmerger_util_module
   use initial_model_module, only: model_P, model_S

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

   real(rt) :: theta, phi

   integer  :: i, j, k
   real(rt) :: loc(3)
   real(rt) :: x, y, z
   real(rt) :: dist_P, dist_S

   !Obtain the vector potential A
   !for each star we want to have
   !the magnetic field similarly
   !to what we do for StarGrav 
   
   !get the A_x component
   ! A_x, i, j-1/2, k-1/2
   do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
              loc = position(i,j,k, .false., .true., .true.)

              dist_P = sum((loc - center_P_initial)**2)**HALF
              dist_S = sum((loc - center_S_initial)**2)**HALF
             
              !vector potential inside the primary star    
              if (mass_P > ZERO .and. (dist_P < model_P % radius ) ) then
              
                  x = loc(1) - center_P_initial(1)
                  y = loc(2) - center_P_initial(2)
                  z = loc(3) - center_P_initial(3)
                  
                  theta = atan2(sqrt(x**2+y**2),z)
                  phi = atan2(y,x)

                  A_x(i,j,k) = -m_0*dist_P*sin(theta)*sin(phi)/(model_P %radius)**3
              
              !vector potential inside the secondary star
              else if (mass_S > ZERO .and. (dist_S < model_S % radius)) then
                  
                  x = loc(1) - center_S_initial(1)
                  y = loc(2) - center_S_initial(2)
                  z = loc(3) - center_S_initial(3)
                  
                  theta = atan2(sqrt(x**2+y**2),z)
                  phi = atan2(y,x)

                  A_x(i,j,k) = -m_0*dist_S*sin(theta)*sin(phi)/(model_S %radius)**3 

              
              else 
               !dipole outside the primary star
   
                  x = loc(1) - center_P_initial(1)
                  y = loc(2) - center_P_initial(2)
                  z = loc(3) - center_P_initial(3)
                 
                  theta = atan2(sqrt(x**2+y**2),z)
                  phi = atan2(y,x)

                  A_x(i,j,k) = -m_0*sin(theta)*sin(phi)/(dist_P*dist_P)

              
              ! add contribution of dipole outside the secondary star 
                  x = loc(1) - center_S_initial(1)
                  y = loc(2) - center_S_initial(2)
                  z = loc(3) - center_S_initial(3)
   
                  theta = atan2(sqrt(x**2+y**2),z)
                  phi = atan2(y,x)

                  A_x(i,j,k) = A_x(i,j,k) - m_0*sin(theta)*sin(phi)/(dist_S*dist_S)

              endif

            enddo
         enddo
     enddo

  !get the A_y component
  !A_y, i-1/2, j, k-1/2
   do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              loc = position(i,j,k, .true., .false., .true.)

              dist_P = sum((loc - center_P_initial)**2)**HALF
              dist_S = sum((loc - center_S_initial)**2)**HALF
            
              !vector potential inside the primary star  
              if (mass_P > ZERO .and. (dist_P < model_P % radius ) ) then 
                  
                  x = loc(1) - center_P_initial(1)
                  y = loc(2) - center_P_initial(2)
                  z = loc(3) - center_P_initial(3)
   
                  theta = atan2(sqrt(x**2+y**2),z)
                  phi = atan2(y,x)
    
                  A_y(i,j,k) = m_0*dist_P*sin(theta)*cos(phi)/(model_P % radius)**3

              !vector potentential inside the secondary star
              else if (mass_S > ZERO .and. (dist_S < model_S % radius)) then
                  
                  x = loc(1) - center_S_initial(1)
                  y = loc(2) - center_S_initial(2)
                  z = loc(3) - center_S_initial(3)
   
                  theta = atan2(sqrt(x**2+y**2),z)
                  phi = atan2(y,x)
    
                  A_y(i,j,k) = m_0*dist_S*sin(theta)*cos(phi)/(model_S % radius)**3

              else 

                  ! dipole outside the primary star
                  x = loc(1) - center_P_initial(1)
                  y = loc(2) - center_P_initial(2)
                  z = loc(3) - center_P_initial(3)
   
                  theta = atan2(sqrt(x**2+y**2),z)
                  phi = atan2(y,x)

                  A_y(i,j,k) = m_0*sin(theta)*cos(phi)/(dist_P*dist_P)
             
 
                  ! add contribution of dipole outside the secondary star 
                  x = loc(1) - center_S_initial(1)
                  y = loc(2) - center_S_initial(2)
                  z = loc(3) - center_S_initial(3)
   
                  theta = atan2(sqrt(x**2+y**2),z)
                  phi = atan2(y,x)
    
                  A_y(i,j,k) = A_y(i,j,k) + m_0*sin(theta)*cos(phi)/(dist_S*dist_S)   

              endif

              
            enddo
         enddo
     enddo



 ! A is of the form A = (A_x, A_y, 0)
 ! B_x = -dA_x/dz
 ! B_y = dA_y/dz
 ! B_z = dA_y/dx - dA_x/dy 

   !Initialize magnetic fields
   !B_x, i-1/2,j,k
   do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            mag_x(i,j,k,1) = (A_y(i,j,k) - A_y(i,j,k+1))/dx(3)
         enddo
      enddo
   enddo

  !B_y, i,j-1/2,k
   do k = lo(3), hi(3)
      do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)
            mag_y(i,j,k,1) = (A_x(i,j,k+1) - A_x(i,j,k))/dx(3)
         enddo
      enddo
   enddo

!B_z, i,j,k-1/2
   do k = lo(3), hi(3)+1
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            mag_z(i,j,k,1) = (A_y(i+1,j,k) - A_y(i,j,k))/dx(1) &
                           - (A_x(i,j+1,k) - A_x(i,j,k))/dx(2)

         enddo
      enddo
   enddo

   end subroutine ca_initmag
