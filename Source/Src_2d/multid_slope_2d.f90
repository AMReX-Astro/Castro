module multid_slope_module
  
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public multid_slope

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine multid_slope(q, flatn, q_lo, q_hi, &
                          dqx, dqy, qpd_lo, qpd_hi, &
                          dx, dy, &
                          ilo, jlo, ihi, jhi)

    use bl_constants_module        

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: ilo, jlo, ihi, jhi

    real(rt)        , intent(in) :: dx, dy

    real(rt)        , intent(in) ::     q( q_lo(1): q_hi(1), q_lo(2): q_hi(2))
    real(rt)        , intent(in) :: flatn( q_lo(1): q_hi(1), q_lo(2): q_hi(2))

    real(rt)        , intent(inout) :: dqx(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2))
    real(rt)        , intent(inout) :: dqy(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2))

    integer :: i, j, m, n
    
    real(rt)        , allocatable :: q_nd(:,:)

    integer, parameter :: ill = 1, ilh = 2, irl = 3, irh = 4
    real(rt)         :: ss(4), ss_temp(4), min_ss(4), max_ss(4), diff(4)
    real(rt)         :: A_x, A_y, A_xy
    real(rt)         :: sumdiff, sgndiff, redfac, redmax, div
    real(rt)        , parameter :: eps = 1.e-10_rt
    integer :: kdp
    integer, parameter :: niter = 3
    
  
    ! here we do a multidimensional reconstruction of the data, forming
    ! a bilinear polynomial: 
    !     q(x,y) = A_xy (x-x_i)(y-y_i) + A_x (x-x_i) + A_y (y-y_i) + A_avg
    !
    ! First find the values at the nodes.  Here q_nd(i,j) will refer to
    ! a_{i-1/2,j-1/2}

    allocate(q_nd(q_lo(1):q_hi(1),q_lo(2):q_hi(2)))

    do j = jlo-2, jhi+3
       do i = ilo-2, ihi+3
          
          q_nd(i,j) = &
             (         q(i-2,j-2) -  7.0e0_rt*(q(i-1,j-2) + q(i,j-2)) +       q(i+1,j-2) + &
              (-7.0e0_rt)*q(i-2,j-1) + 49.0e0_rt*(q(i-1,j-1) + q(i,j-1)) - 7.0e0_rt*q(i+1,j-1) + &
              (-7.0e0_rt)*q(i-2,j  ) + 49.0e0_rt*(q(i-1,j  ) + q(i,j  )) - 7.0e0_rt*q(i+1,j  ) + &
                       q(i-2,j+1) -  7.0e0_rt*(q(i-1,j+1) + q(i,j+1)) +       q(i+1,j+1))/ &
             144.0e0_rt

       enddo
    enddo

    ! now get the coefficients of the bilinear polynomial.  Here we loop
    ! ! over zone centers.
    do j = jlo-1, jhi+1
       do i = ilo-1, ihi+1
          
          ss(ill) = q_nd(i,j)     ! a_{i-1/2,j-1/2}
          ss(ilh) = q_nd(i,j+1)   ! a_{i-1/2,j+1/2}
          ss(irl) = q_nd(i+1,j)   ! a_{i+1/2,j-1/2}
          ss(irh) = q_nd(i+1,j+1) ! a_{i+1/2,j+1/2}
          
          A_x  = ((ss(irh) + ss(irl)) - (ss(ilh) + ss(ill)))/(2.0e0_rt*dx)
          A_y  = ((ss(ilh) + ss(irh)) - (ss(ill) + ss(irl)))/(2.0e0_rt*dy)
          A_xy = ((ss(irh) - ss(irl)) - (ss(ilh) - ss(ill)))/(dx*dy)
          
          ! check if we are within the cc-values
          ss_temp(ill) = q(i,j) - HALF*dx*A_x - HALF*dy*A_y + 0.25e0_rt*dx*dy*A_xy
          ss_temp(ilh) = q(i,j) - HALF*dx*A_x + HALF*dy*A_y - 0.25e0_rt*dx*dy*A_xy
          ss_temp(irl) = q(i,j) + HALF*dx*A_x - HALF*dy*A_y - 0.25e0_rt*dx*dy*A_xy
          ss_temp(irh) = q(i,j) + HALF*dx*A_x + HALF*dy*A_y + 0.25e0_rt*dx*dy*A_xy
          
          min_ss(ill) = min(q(i-1,j-1), q(i,j-1), q(i-1,j), q(i,j))
          max_ss(ill) = max(q(i-1,j-1), q(i,j-1), q(i-1,j), q(i,j))
          
          min_ss(ilh) = min(q(i-1,j+1), q(i,j+1), q(i-1,j), q(i,j))
          max_ss(ilh) = max(q(i-1,j+1), q(i,j+1), q(i-1,j), q(i,j))
          
          min_ss(irl) = min(q(i+1,j-1), q(i,j-1), q(i+1,j), q(i,j))
          max_ss(irl) = max(q(i+1,j-1), q(i,j-1), q(i+1,j), q(i,j))
          
          min_ss(irh) = min(q(i+1,j+1), q(i,j+1), q(i+1,j), q(i,j))
          max_ss(irh) = max(q(i+1,j+1), q(i,j+1), q(i+1,j), q(i,j))
          
          ! limit
          do m = 1, 4
             ss_temp(m) = max(min(ss_temp(m), max_ss(m)), min_ss(m))
          enddo

          do n = 1, niter
             sumdiff = (ss_temp(ill) + ss_temp(ilh) + &
                        ss_temp(irl) + ss_temp(irh)) - 4.0e0_rt*q(i,j)
             sgndiff = sign(1.e0_rt,sumdiff)

             do m = 1, 4
                diff(m) = (ss_temp(m) - q(i,j))*sgndiff
             enddo
             
             kdp = 0
             
             do m = 1, 4
                if (diff(m) > eps) kdp = kdp + 1
             enddo
             
             do m = 1, 4
                if (kdp < 1) then
                   div = 1.e0_rt
                else
                   div = dble(kdp)
                endif
                
                if (diff(m) > eps) then
                   redfac = sumdiff*sgndiff/div
                   kdp = kdp-1
                else
                   redfac = 0.0e0_rt
                endif
              
                if (sgndiff > 0.0e0_rt) then
                   redmax = ss_temp(m) - min_ss(m)
                else
                   redmax = max_ss(m) - ss_temp(m)
                endif
                
                redfac = min(redfac, redmax)
                sumdiff = sumdiff - redfac*sgndiff
                ss_temp(m) = ss_temp(m) - redfac*sgndiff
             enddo
          enddo
          
          ! construct the final slopes
          A_x  = ((ss_temp(irh) + ss_temp(irl)) - &
                  (ss_temp(ilh) + ss_temp(ill)))/(2.0e0_rt*dx)

          A_y  = ((ss_temp(ilh) + ss_temp(irh)) - &
                  (ss_temp(ill) + ss_temp(irl)))/(2.0e0_rt*dy)

          A_xy = ((ss_temp(irh) - ss_temp(irl)) - &
                  (ss_temp(ilh) - ss_temp(ill)))/(dx*dy)


          ! now construct the limited differences in x and y
          dqx(i,j) = flatn(i,j)*A_x*dx
          dqy(i,j) = flatn(i,j)*A_y*dy

       enddo
    enddo

    deallocate(q_nd)

  end subroutine multid_slope   

end module multid_slope_module
