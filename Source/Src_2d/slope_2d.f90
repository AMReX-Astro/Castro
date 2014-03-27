module slope_module
  
  implicit none

  private

  public uslope, pslope

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

end module slope_module
