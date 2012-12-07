module slope_module

  implicit none

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
      subroutine uslope(q,flatn,qd_l1,qd_h1,dq,qpd_l1,qpd_h1,ilo,ihi,nv)

      use meth_params_module, only : QPRES

      implicit none

      integer ilo, ihi, nv
      integer qd_l1,qd_h1,qpd_l1,qpd_h1
      double precision q(qd_l1:qd_h1, nv)
      double precision dq(qpd_l1:qpd_h1,nv)
      double precision flatn(qd_l1:qd_h1)

!     Local arrays 
      double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)

      integer i, n
      double precision dlft, drgt, slop, dq1
      double precision four3rd, sixth

      four3rd = 4.d0/3.d0
      sixth = 1.d0/6.d0

      allocate (dsgn(ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2))

      do n = 1, nv 
!       if (n .ne. QPRES) then

            ! first compute Fromm slopes
            do i = ilo-2, ihi+2 
                dlft = 2.d0*(q(i  ,n) - q(i-1,n))
                drgt = 2.d0*(q(i+1,n) - q(i  ,n))
                dcen(i) = .25d0 * (dlft+drgt)
                dsgn(i) = sign(1.d0, dcen(i))
                slop = min( abs(dlft), abs(drgt) )
!                dlim(i) = cvmgp( slop, 0.d0, dlft*drgt )
                if (dlft*drgt .ge. 0.d0) then
                   dlim(i) = slop
                else
                   dlim(i) = 0.d0
                endif
                df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
            enddo

            ! now limited fourth order slopes
            do i = ilo-1, ihi+1 
                dq1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
                dq(i,n) = flatn(i)* &
                            dsgn(i)*min(dlim(i),abs(dq1))
            enddo
!       end if
      enddo

      deallocate (dsgn,dlim,df,dcen)

      end subroutine uslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
      subroutine pslope(p,rho,flatn,qd_l1,qd_h1,dp,qpd_l1,qpd_h1,grav,gv_l1,gv_h1,ilo,ihi,dx)

      implicit none

      integer ilo, ihi
      integer  qd_l1, qd_h1
      integer qpd_l1,qpd_h1
      integer  gv_l1, gv_h1
      double precision, intent(in   ) ::      p( qd_l1: qd_h1)
      double precision, intent(in   ) ::    rho( qd_l1: qd_h1)
      double precision, intent(in   ) ::  flatn( qd_l1: qd_h1)
      double precision, intent(  out) ::     dp(qpd_l1:qpd_h1)
      double precision, intent(in   ) ::   grav( gv_l1: gv_h1)
      double precision, intent(in   ) ::  dx

!     Local arrays
      double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)

      integer i
      double precision dlft, drgt, dp1
      double precision four3rd, sixth

      four3rd = 4.d0/3.d0
      sixth = 1.d0/6.d0

      allocate (dsgn(ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2))

      ! first compute Fromm slopes
      do i = ilo-2, ihi+2 
          dlft = p(i  ) - p(i-1)
          drgt = p(i+1) - p(i  )

          ! Here we subtract off (rho * grav) so as not to limit that part of the slope
          dlft = dlft - 0.25d0 * (rho(i)+rho(i-1))*(grav(i)+grav(i-1))*dx
          drgt = drgt - 0.25d0 * (rho(i)+rho(i+1))*(grav(i)+grav(i+1))*dx
!         dlft = dlft - rho(i)*grav(i)*dx
!         drgt = drgt - rho(i)*grav(i)*dx

          dcen(i) = 0.5d0*(dlft+drgt)
          dsgn(i) = sign(1.d0, dcen(i))

          if (dlft*drgt .ge. 0.d0) then
             dlim(i) = 2.d0 * min( abs(dlft), abs(drgt) )
          else
             dlim(i) = 0.d0
          endif
          df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
      enddo

      if (ilo .eq. 0) then
        df(-1) = -df(0)
        df(-2) = -df(1)
      end if

      ! now limited fourth order slopes
      do i = ilo-1, ihi+1 
          dp1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
          dp(i) = flatn(i)* &
                      dsgn(i)*min(dlim(i),abs(dp1))
          dp(i) = dp(i) + rho(i)*grav(i)*dx
      enddo

      ! Here we are assuming a symmetry boundary condition
      if (ilo .eq. 0) dp(-1) = -dp(0)

      deallocate (dsgn,dlim,df,dcen)

      end subroutine pslope

end module slope_module
