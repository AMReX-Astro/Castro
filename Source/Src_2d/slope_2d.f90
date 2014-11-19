module slope_module
  
  implicit none

  private

  public uslope, pslope, multid_slope

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine uslope(q,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                    dq,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                    ilo1,ilo2,ihi1,ihi2,nv,idir)

    use bl_constants_module        

    implicit none

    integer ilo,ihi
    integer qd_l1,qd_l2,qd_h1,qd_h2
    integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer ilo1,ilo2,ihi1,ihi2,nv,idir
    
    double precision     q( qd_l1: qd_h1, qd_l2: qd_h2,nv)
    double precision flatn( qd_l1: qd_h1, qd_l2: qd_h2)
    double precision    dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,nv)
    
    ! local
    double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)
    
    integer i, j, n
    double precision dlft, drgt, dq1
    
    ilo = MIN(ilo1,ilo2)
    ihi = MAX(ihi1,ihi2)
    
    allocate (dsgn(ilo-2:ihi+2))
    allocate (dlim(ilo-2:ihi+2))
    allocate (  df(ilo-2:ihi+2))
    allocate (dcen(ilo-2:ihi+2))
    
    do n = 1, nv 
       
       if (idir .eq. 1) then
          
          ! Slopes in first coordinate direction
          do j = ilo2-1, ihi2+1
             ! First compute Fromm slopes
             do i = ilo1-2, ihi1+2
                dlft = q(i  ,j,n) - q(i-1,j,n)
                drgt = q(i+1,j,n) - q(i  ,j,n)
                
                dcen(i) = HALF*(dlft+drgt)
                dsgn(i) = sign(ONE, dcen(i))
                
                if (dlft*drgt .ge. ZERO) then
                   dlim(i) = TWO * min( abs(dlft), abs(drgt) )
                else
                   dlim(i) = ZERO
                endif
                df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
             enddo
             
             ! now limited fourth order slopes
             do i = ilo1-1, ihi1+1
                dq1 = FOUR3RD*dcen(i) - SIXTH*(df(i+1) + df(i-1))
                dq(i,j,n) = flatn(i,j)*dsgn(i)*min(dlim(i),abs(dq1))
             enddo
          enddo
          
       else

          ! Slopes in second coordinate direction
          do i = ilo1-1, ihi1+1
             ! First compute Fromm slopes
             do j = ilo2-2, ihi2+2
                dlft = q(i,j  ,n) - q(i,j-1,n)
                drgt = q(i,j+1,n) - q(i,j  ,n)
                
                dcen(j) = HALF*(dlft+drgt)
                dsgn(j) = sign( ONE, dcen(j) )
                
                if (dlft*drgt .ge. ZERO) then
                   dlim(j) = TWO * min( abs(dlft), abs(drgt) )
                else
                   dlim(j) = ZERO
                endif
                df(j) = dsgn(j)*min( dlim(j),abs(dcen(j)) )
             enddo
             
             ! now limited fourth order slopes
             do j = ilo2-1, ihi2+1
                dq1 = FOUR3RD*dcen(j) - SIXTH*(df(j+1) + df(j-1))
                dq(i,j,n) = flatn(i,j)*dsgn(j)*min(dlim(j),abs(dq1))
             enddo
          enddo
          
       endif
       
    enddo
    
    deallocate(dsgn,dlim,df,dcen)
    
  end subroutine uslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine pslope(p,rho,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                    dp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                    grav,gv_l1,gv_l2,gv_h1,gv_h2, &
                    ilo1,ilo2,ihi1,ihi2,dx,dy,idir)
    
    use bl_constants_module
    
    implicit none
    
    integer ilo,ihi
    integer qd_l1,qd_l2,qd_h1,qd_h2
    integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer gv_l1,gv_l2,gv_h1,gv_h2
    integer ilo1,ilo2,ihi1,ihi2,idir
    
    double precision, intent(in   ) ::      p( qd_l1: qd_h1, qd_l2: qd_h2)
    double precision, intent(in   ) ::    rho( qd_l1: qd_h1, qd_l2: qd_h2)
    double precision, intent(in   ) ::  flatn( qd_l1: qd_h1, qd_l2: qd_h2)
    double precision, intent(  out) ::     dp(qpd_l1:qpd_h1,qpd_l2:qpd_h2)
    double precision, intent(in   ) ::   grav( gv_l1: gv_h1, gv_l2: gv_h2,2)
    double precision, intent(in   ) ::  dx,dy
    
    ! local
    double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)
    
    integer i, j
    double precision dlft, drgt, dp1
    
    ilo = MIN(ilo1,ilo2)
    ihi = MAX(ihi1,ihi2)
    
    allocate (dsgn(ilo-2:ihi+2))
    allocate (dlim(ilo-2:ihi+2))
    allocate (  df(ilo-2:ihi+2))
    allocate (dcen(ilo-2:ihi+2))
    
    if (idir .eq. 1) then
       
       ! Slopes in first coordinate direction
       do j = ilo2-1, ihi2+1
          ! First compute Fromm slopes
          do i = ilo1-2, ihi1+2
             dlft = p(i  ,j) - p(i-1,j)
             drgt = p(i+1,j) - p(i  ,j)
             
             ! Subtract off (rho * grav) so as not to limit that part of the slope
             dlft = dlft - FOURTH * (rho(i,j)+rho(i-1,j))*(grav(i,j,1)+grav(i-1,j,1))*dx
             drgt = drgt - FOURTH * (rho(i,j)+rho(i+1,j))*(grav(i,j,1)+grav(i+1,j,1))*dx
             
             dcen(i) = HALF*(dlft+drgt)
             dsgn(i) = sign(ONE, dcen(i))
             
             if (dlft*drgt .ge. ZERO) then
                dlim(i) = TWO * min( abs(dlft), abs(drgt) )
             else
                dlim(i) = ZERO
             endif
             df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
          enddo
          
          ! now limited fourth order slopes
          do i = ilo1-1, ihi1+1
             dp1 = FOUR3RD*dcen(i) - SIXTH*(df(i+1) + df(i-1))
             dp(i,j) = flatn(i,j)*dsgn(i)*min(dlim(i),abs(dp1))
             dp(i,j) = dp(i,j) + rho(i,j)*grav(i,j,1)*dx
          enddo
       enddo
       
    else

       ! Slopes in second coordinate direction
       do i = ilo1-1, ihi1+1
          ! First compute Fromm slopes
          do j = ilo2-2, ihi2+2
             dlft = p(i,j  ) - p(i,j-1)
             drgt = p(i,j+1) - p(i,j  )
             
             ! Subtract off (rho * grav) so as not to limit that part of the slope
             dlft = dlft - FOURTH * (rho(i,j)+rho(i,j-1))*(grav(i,j,2)+grav(i,j-1,2))*dy
             drgt = drgt - FOURTH * (rho(i,j)+rho(i,j+1))*(grav(i,j,2)+grav(i,j+1,2))*dy
             
             dcen(j) = HALF*(dlft+drgt)
             dsgn(j) = sign( ONE, dcen(j) )
             
             if (dlft*drgt .ge. ZERO) then
                dlim(j) = TWO * min( abs(dlft), abs(drgt) )
             else
                dlim(j) = ZERO
             endif
             df(j) = dsgn(j)*min( dlim(j),abs(dcen(j)) )
          enddo
          
          ! now limited fourth order slopes
          do j = ilo2-1, ihi2+1
             dp1 = FOUR3RD*dcen(j) - SIXTH*(df(j+1) + df(j-1))
             dp(i,j) = flatn(i,j)*dsgn(j)*min(dlim(j),abs(dp1))
             dp(i,j) = dp(i,j) + rho(i,j)*grav(i,j,2)*dy
          enddo
       enddo
       
    endif
    
    deallocate(dsgn,dlim,df,dcen)
    
  end subroutine pslope


  subroutine multid_slope(q, flatn, &
                          qd_l1, qd_l2, qd_h1, qd_h2, &
                          dqx, dqy, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                          dx, dy, &
                          ilo, jlo, ihi, jhi)

    use bl_constants_module        

    implicit none

    integer, intent(in) :: qd_l1, qd_l2, qd_h1, qd_h2
    integer, intent(in) :: qpd_l1, qpd_l2, qpd_h1, qpd_h2
    integer, intent(in) :: ilo, jlo, ihi, jhi

    double precision, intent(in) :: dx, dy

    double precision, intent(in) ::     q( qd_l1: qd_h1, qd_l2: qd_h2)
    double precision, intent(in) :: flatn( qd_l1: qd_h1, qd_l2: qd_h2)

    double precision, intent(inout) :: dqx(qpd_l1:qpd_h1,qpd_l2:qpd_h2)
    double precision, intent(inout) :: dqy(qpd_l1:qpd_h1,qpd_l2:qpd_h2)

    integer :: i, j, m, n
    
    double precision, allocatable :: q_nd(:,:)

    integer, parameter :: ill = 1, ilh = 2, irl = 3, irh = 4
    double precision :: ss(4), ss_temp(4), min_ss(4), max_ss(4), diff(4)
    double precision :: A_x, A_y, A_xy
    double precision :: sumdiff, sgndiff, redfac, redmax, div
    double precision, parameter :: eps = 1.d-10
    integer :: kdp
    integer, parameter :: niter = 3
    
  
    ! here we do a multidimensional reconstruction of the data, forming
    ! a bilinear polynomial: 
    !     q(x,y) = A_xy (x-x_i)(y-y_i) + A_x (x-x_i) + A_y (y-y_i) + A_avg
    !
    ! First find the values at the nodes.  Here q_nd(i,j) will refer to
    ! a_{i-1/2,j-1/2}

    allocate(q_nd(qd_l1:qd_h1,qd_l2:qd_h2))

    do j = jlo-2, jhi+3
       do i = ilo-2, ihi+3
          
          q_nd(i,j) = &
             (         q(i-2,j-2) -  7.0d0*(q(i-1,j-2) + q(i,j-2)) +       q(i+1,j-2) + &
              (-7.0d0)*q(i-2,j-1) + 49.0d0*(q(i-1,j-1) + q(i,j-1)) - 7.0d0*q(i+1,j-1) + &
              (-7.0d0)*q(i-2,j  ) + 49.0d0*(q(i-1,j  ) + q(i,j  )) - 7.0d0*q(i+1,j  ) + &
                       q(i-2,j+1) -  7.0d0*(q(i-1,j+1) + q(i,j+1)) +       q(i+1,j+1))/ &
             144.0d0

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
          
          A_x  = ((ss(irh) + ss(irl)) - (ss(ilh) + ss(ill)))/(2.0d0*dx)
          A_y  = ((ss(ilh) + ss(irh)) - (ss(ill) + ss(irl)))/(2.0d0*dy)
          A_xy = ((ss(irh) - ss(irl)) - (ss(ilh) - ss(ill)))/(dx*dy)
          
          ! check if we are within the cc-values
          ss_temp(ill) = q(i,j) - HALF*dx*A_x - HALF*dy*A_y + 0.25d0*dx*dy*A_xy
          ss_temp(ilh) = q(i,j) - HALF*dx*A_x + HALF*dy*A_y - 0.25d0*dx*dy*A_xy
          ss_temp(irl) = q(i,j) + HALF*dx*A_x - HALF*dy*A_y - 0.25d0*dx*dy*A_xy
          ss_temp(irh) = q(i,j) + HALF*dx*A_x + HALF*dy*A_y + 0.25d0*dx*dy*A_xy
          
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
                        ss_temp(irl) + ss_temp(irh)) - 4.0d0*q(i,j)
             sgndiff = sign(1.d0,sumdiff)

             do m = 1, 4
                diff(m) = (ss_temp(m) - q(i,j))*sgndiff
             enddo
             
             kdp = 0
             
             do m = 1, 4
                if (diff(m) > eps) kdp = kdp + 1
             enddo
             
             do m = 1, 4
                if (kdp < 1) then
                   div = 1.d0
                else
                   div = dble(kdp)
                endif
                
                if (diff(m) > eps) then
                   redfac = sumdiff*sgndiff/div
                   kdp = kdp-1
                else
                   redfac = 0.0d0
                endif
              
                if (sgndiff > 0.0d0) then
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
                  (ss_temp(ilh) + ss_temp(ill)))/(2.0d0*dx)

          A_y  = ((ss_temp(ilh) + ss_temp(irh)) - &
                  (ss_temp(ill) + ss_temp(irl)))/(2.0d0*dy)

          A_xy = ((ss_temp(irh) - ss_temp(irl)) - &
                  (ss_temp(ilh) - ss_temp(ill)))/(dx*dy)


          ! now construct the limited differences in x and y
          dqx(i,j) = A_x*dx
          dqy(i,j) = A_y*dy

       enddo
    enddo

  end subroutine multid_slope   

end module slope_module
