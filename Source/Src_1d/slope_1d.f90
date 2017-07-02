module slope_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
      subroutine uslope(q, flatn, q_lo, q_hi, &
                        dq, qpd_lo, qpd_hi, &
                        ilo, ihi)

      use meth_params_module, only : NQ
      use bl_constants_module

      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer ilo, ihi, nv
      integer q_lo(3), q_hi(3)
      integer :: qpd_lo(3), qpd_hi(3)
      real(rt)         q(q_lo(1):q_hi(1),NQ)
      real(rt)         dq(qpd_lo(1):qpd_hi(1),NQ)
      real(rt)         flatn(q_lo(1):q_hi(1))

!     Local arrays 
      real(rt)        , allocatable::dsgn(:),dlim(:),df(:),dcen(:)

      integer i, n
      real(rt)         dlft, drgt, slop, dq1

      allocate (dsgn(ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2))

      do n = 1, nv 
!       if (n .ne. QPRES) then

            ! first compute Fromm slopes
            do i = ilo-2, ihi+2 
                dlft = TWO*(q(i  ,n) - q(i-1,n))
                drgt = TWO*(q(i+1,n) - q(i  ,n))
                dcen(i) = FOURTH * (dlft+drgt)
                dsgn(i) = sign(ONE, dcen(i))
                slop = min( abs(dlft), abs(drgt) )
!                dlim(i) = cvmgp( slop, ZERO, dlft*drgt )
                if (dlft*drgt .ge. ZERO) then
                   dlim(i) = slop
                else
                   dlim(i) = ZERO
                endif
                df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
            enddo

            ! now limited fourth order slopes
            do i = ilo-1, ihi+1 
                dq1 = FOUR3RD*dcen(i) - SIXTH*(df(i+1) + df(i-1))
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
      subroutine pslope(q, flatn, q_lo, q_hi, &
                        dq, qpd_lo, qpd_hi, &
                        src, src_lo, src_hi, &
                        ilo, ihi, dx)

      use bl_constants_module
      use meth_params_module, only: QU, QVAR, NQ, QPRES, QRHO
      
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer ilo, ihi
      integer  q_lo(3), q_hi(3)
      integer qpd_lo(3), qpd_hi(3)
      integer src_lo(3), src_hi(3)
      real(rt)        , intent(in   ) ::      q( q_lo(1):q_hi(1),NQ)
      real(rt)        , intent(in   ) ::  flatn( q_lo(1):q_hi(1))
      real(rt)        , intent(  out) ::     dq(qpd_lo(1):qpd_hi(1),NQ)
      real(rt)        , intent(in   ) ::    src(src_lo(1):src_hi(1),QVAR)
      real(rt)        , intent(in   ) ::  dx

!     Local arrays
      real(rt)        , allocatable :: dsgn(:), dlim(:), df(:), dcen(:)

      integer          :: i
      real(rt)         :: dlft, drgt, dp1

      allocate (dsgn(ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2))

      ! first compute Fromm slopes
      do i = ilo-2, ihi+2 
          dlft = q(i  ,QPRES) - q(i-1,QPRES)
          drgt = q(i+1,QPRES) - q(i  ,QPRES)

          ! Here we subtract off (rho * acceleration) so as not to limit that part of the slope
          dlft = dlft - FOURTH * (q(i,QRHO) + q(i-1,QRHO)) * (src(i,QU) + src(i-1,QU))*dx
          drgt = drgt - FOURTH * (q(i,QRHO) + q(i+1,QRHO)) * (src(i,QU) + src(i+1,QU))*dx

          dcen(i) = HALF*(dlft+drgt)
          dsgn(i) = sign(ONE, dcen(i))

          if (dlft*drgt .ge. ZERO) then
             dlim(i) = TWO * min( abs(dlft), abs(drgt) )
          else
             dlim(i) = ZERO
          endif
          df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
      enddo

      if (ilo .eq. 0) then
        df(-1) = -df(0)
        df(-2) = -df(1)
      end if

      ! now limited fourth order slopes
      do i = ilo-1, ihi+1 
          dp1 = FOUR3RD*dcen(i) - SIXTH*(df(i+1) + df(i-1))
          dq(i,QPRES) = flatn(i) * dsgn(i)*min(dlim(i),abs(dp1))
          dq(i,QPRES) = dq(i,QPRES) + q(i,QRHO)*src(i,QU)*dx
      enddo

      ! Here we are assuming a symmetry boundary condition
      if (ilo .eq. 0) dq(-1,QPRES) = -dq(0,QPRES)

      deallocate (dsgn,dlim,df,dcen)

      end subroutine pslope

end module slope_module
