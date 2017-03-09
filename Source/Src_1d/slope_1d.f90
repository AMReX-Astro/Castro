module slope_module

  use amrex_fort_module, only : rt => c_real
  implicit none

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
      subroutine uslope(q,flatn,qd_l1,qd_h1,dq,qpd_l1,qpd_h1,ilo,ihi,nv)

!      use meth_params_module, only : QPRES
      use bl_constants_module

      use amrex_fort_module, only : rt => c_real
      implicit none

      integer ilo, ihi, nv
      integer qd_l1,qd_h1,qpd_l1,qpd_h1
      real(rt)         q(qd_l1:qd_h1, nv)
      real(rt)         dq(qpd_l1:qpd_h1,nv)
      real(rt)         flatn(qd_l1:qd_h1)

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
      subroutine pslope(p,rho,flatn,qd_l1,qd_h1,dp,qpd_l1,qpd_h1, &
                        src,src_l1,src_h1,ilo,ihi,dx)

      use bl_constants_module
      use meth_params_module, only: QU, QVAR
      
      use amrex_fort_module, only : rt => c_real
      implicit none

      integer ilo, ihi
      integer  qd_l1, qd_h1
      integer qpd_l1,qpd_h1
      integer src_l1,src_h1
      real(rt)        , intent(in   ) ::      p( qd_l1: qd_h1)
      real(rt)        , intent(in   ) ::    rho( qd_l1: qd_h1)
      real(rt)        , intent(in   ) ::  flatn( qd_l1: qd_h1)
      real(rt)        , intent(  out) ::     dp(qpd_l1:qpd_h1)
      real(rt)        , intent(in   ) ::    src(src_l1:src_h1,QVAR)
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
          dlft = p(i  ) - p(i-1)
          drgt = p(i+1) - p(i  )

          ! Here we subtract off (rho * acceleration) so as not to limit that part of the slope
          dlft = dlft - FOURTH * (rho(i)+rho(i-1))*(src(i,QU)+src(i-1,QU))*dx
          drgt = drgt - FOURTH * (rho(i)+rho(i+1))*(src(i,QU)+src(i+1,QU))*dx
!         dlft = dlft - rho(i)*src(i,QU)*dx
!         drgt = drgt - rho(i)*src(i,QU)*dx

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
          dp(i) = flatn(i)* &
                      dsgn(i)*min(dlim(i),abs(dp1))
          dp(i) = dp(i) + rho(i)*src(i,QU)*dx
      enddo

      ! Here we are assuming a symmetry boundary condition
      if (ilo .eq. 0) dp(-1) = -dp(0)

      deallocate (dsgn,dlim,df,dcen)

      end subroutine pslope

end module slope_module
