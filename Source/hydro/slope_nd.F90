module slope_module

  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg

  implicit none

  private

  public uslope, pslope

contains

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine uslope(q, flatn, qd_lo, qd_hi, &
                    dqx, dqy, dqz, qpd_lo, qpd_hi, &
                    ilo1, ilo2, ihi1, ihi2, kc, k3d)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only: NQ, plm_iorder
    use bl_constants_module, only: ZERO, HALF, ONE, TWO, FOUR3RD, FOURTH, SIXTH

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qpd_lo(3),qpd_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2, kc, k3d

    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    real(rt), intent(inout) :: dqx(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: dqy(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: dqz(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    integer :: i, j, k, n

    real(rt) :: dlft, drgt, slop, dq1
    real(rt) :: dm, dp, dc, ds, sl, dl, dfm, dfp

    integer :: ilo, ihi

    real(rt), pointer::dsgn(:,:),dlim(:,:),df(:,:),dcen(:,:)

    ilo = min(ilo1, ilo2)
    ihi = max(ihi1, ihi2)

    call bl_allocate(dsgn, ilo-2, ihi+2, ilo-2*dg(2), ihi+2*dg(2))
    call bl_allocate(dlim, ilo-2, ihi+2, ilo-2*dg(2), ihi+2*dg(2))
    call bl_allocate(  df, ilo-2, ihi+2, ilo-2*dg(2), ihi+2*dg(2))
    call bl_allocate(dcen, ilo-2, ihi+2, ilo-2*dg(2), ihi+2*dg(2))

    if (plm_iorder == 1) then

       do n = 1, NQ
          do j = ilo2-dg(2), ihi2+dg(2)
             do i = ilo1-1, ihi1+1
                dqx(i,j,kc,n) = ZERO
#if (BL_SPACEDIM >= 2)
                dqy(i,j,kc,n) = ZERO
#if (BL_SPACEDIM == 3)
                dqz(i,j,kc,n) = ZERO
#endif
#endif
             enddo
          enddo
       enddo

    else

       do n = 1, NQ


          ! Compute slopes in first coordinate direction
          do j = ilo2-dg(2), ihi2+dg(2)

             ! First compute Fromm slopes
             do i = ilo1-2, ihi1+2
                dlft = TWO*(q(i ,j,k3d,n) - q(i-1,j,k3d,n))
                drgt = TWO*(q(i+1,j,k3d,n) - q(i ,j,k3d,n))
                dcen(i,j) = FOURTH * (dlft+drgt)
                dsgn(i,j) = sign(ONE, dcen(i,j))
                slop = min( abs(dlft), abs(drgt) )
                if (dlft*drgt .ge. ZERO) then
                   dlim(i,j) = slop
                else
                   dlim(i,j) = ZERO
                endif
                df(i,j) = dsgn(i,j)*min( dlim(i,j), abs(dcen(i,j)) )
             enddo

             ! Now compute limited fourth order slopes
             do i = ilo1-1, ihi1+1
                dq1       = FOUR3RD*dcen(i,j) - SIXTH*(df(i+1,j) + df(i-1,j))
                dqx(i,j,kc,n) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dq1))
             enddo

          enddo

#if (BL_SPACEDIM >= 2)
          ! Compute slopes in second coordinate direction
          do i = ilo1-1, ihi1+1
             ! First compute Fromm slopes for this column
             do j = ilo2-2*dg(2), ihi2+2*dg(2)
                dlft = TWO*(q(i,j ,k3d,n) - q(i,j-1,k3d,n))
                drgt = TWO*(q(i,j+1,k3d,n) - q(i,j ,k3d,n))
                dcen(i,j) = FOURTH * (dlft+drgt)
                dsgn(i,j) = sign( ONE, dcen(i,j) )
                slop = min( abs(dlft), abs(drgt) )
                if (dlft*drgt .ge. ZERO) then
                   dlim(i,j) = slop
                else
                   dlim(i,j) = ZERO
                endif
                df(i,j) = dsgn(i,j)*min( dlim(i,j),abs(dcen(i,j)) )
             enddo

             ! Now compute limited fourth order slopes
             do j = ilo2-dg(2), ihi2+dg(2)
                dq1 = FOUR3RD*dcen(i,j) - SIXTH*( df(i,j+1) + df(i,j-1) )
                dqy(i,j,kc,n) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dq1))
             enddo
          enddo
#endif

#if (BL_SPACEDIM == 3)
          ! Compute slopes in third coordinate direction
          do j = ilo2-1, ihi2+1
             do i = ilo1-1, ihi1+1

                ! Compute Fromm slope on slab below
                k = k3d-1
                dm = TWO*(q(i,j,k ,n) - q(i,j,k-1,n))
                dp = TWO*(q(i,j,k+1,n) - q(i,j,k, n))
                dc = FOURTH*(dm+dp)
                ds = sign( ONE, dc )
                sl = min( abs(dm), abs(dp) )
                if (dm*dp .ge. ZERO) then
                   dl = sl
                else
                   dl = ZERO
                endif
                dfm = ds*min(dl,abs(dc))

                ! Compute Fromm slope on slab above
                k = k3d+1
                dm = TWO*(q(i,j,k ,n) - q(i,j,k-1,n))
                dp = TWO*(q(i,j,k+1,n) - q(i,j,k, n))
                dc = FOURTH*(dm+dp)
                ds = sign( ONE, dc )
                sl = min( abs(dm), abs(dp) )
                if (dm*dp .ge. ZERO) then
                   dl = sl
                else
                   dl = ZERO
                endif
                dfp = ds*min(dl,abs(dc))

                ! Compute Fromm slope on current slab
                k = k3d
                dm = TWO*(q(i,j,k ,n) - q(i,j,k-1,n))
                dp = TWO*(q(i,j,k+1,n) - q(i,j,k, n))
                dc = FOURTH*(dm+dp)
                ds = sign( ONE, dc )
                sl = min( abs(dm), abs(dp) )
                if (dm*dp .ge. ZERO) then
                   dl = sl
                else
                   dl = ZERO
                endif

                ! Now compute limited fourth order slopes
                dq1 = FOUR3RD*dc - SIXTH*( dfp + dfm )
                dqz(i,j,kc,n) = flatn(i,j,k3d)*ds*min(dl,abs(dq1))
             enddo
          enddo
#endif

       enddo  ! component loop

    endif

    call bl_deallocate(dsgn)
    call bl_deallocate(dlim)
    call bl_deallocate(  df)
    call bl_deallocate(dcen)

  end subroutine uslope

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine pslope(q, flatn, q_lo, q_hi, &
                    dqx, dqy, dqz, qpd_lo, qpd_hi, &
                    src, src_lo, src_hi, &
                    ilo1, ilo2, ihi1, ihi2, kc, k3d, dx)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QRHO, QPRES, QU, QV, QW, NQ, QVAR, plm_iorder
    use bl_constants_module, only : ZERO, FOURTH, FOUR3RD, HALF, TWO, ONE, SIXTH

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qpd_lo(3),qpd_hi(3)
    integer, intent(in) :: src_lo(3),src_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2, kc, k3d

    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: flatn(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    real(rt), intent(inout) :: dqx(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: dqy(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: dqz(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    real(rt), intent(in) :: dx(3)

    integer :: i, j, k

    integer :: ilo, ihi

    real(rt) :: dlft, drgt, dp1
    real(rt) :: dm, dp, dc, dl, dfm, dfp, ds

    !     Local arrays
    real(rt), pointer :: dsgn(:,:), dlim(:,:), df(:,:), dcen(:,:)

    ilo = min(ilo1,ilo2)
    ihi = max(ihi1,ihi2)

    call bl_allocate(dsgn, ilo-2,ihi+2,ilo-2,ihi+2)
    call bl_allocate(dlim, ilo-2,ihi+2,ilo-2,ihi+2)
    call bl_allocate(  df, ilo-2,ihi+2,ilo-2,ihi+2)
    call bl_allocate(dcen, ilo-2,ihi+2,ilo-2,ihi+2)

    if (plm_iorder == 1) then

       do j = ilo2-dg(2), ihi2+dg(2)
          do i = ilo1-1, ihi1+1
             dqx(i,j,kc,QPRES) = ZERO
#if (BL_SPACEDIM >= 2)
             dqy(i,j,kc,QPRES) = ZERO
#if (BL_SPACEDIM == 3)
             dqz(i,j,kc,QPRES) = ZERO
#endif
#endif
          enddo
       enddo

    else
       ! Compute slopes in first coordinate direction
       do j = ilo2-dg(2), ihi2+dg(2)

          ! First compute Fromm slopes
          do i = ilo1-2, ihi1+2

             dlft = q(i  ,j,k3d,QPRES) - q(i-1,j,k3d,QPRES)
             drgt = q(i+1,j,k3d,QPRES) - q(i  ,j,k3d,QPRES)

             ! Subtract off (rho * acceleration) so as not to limit that part of the slope
             dlft = dlft - FOURTH * &
                  (q(i,j,k3d,QRHO) + q(i-1,j,k3d,QRHO))*(src(i,j,k3d,QU)+src(i-1,j,k3d,QU))*dx(1)
             drgt = drgt - FOURTH * &
                  (q(i,j,k3d,QRHO) + q(i+1,j,k3d,QRHO))*(src(i,j,k3d,QU)+src(i+1,j,k3d,QU))*dx(1)

             dcen(i,j) = HALF*(dlft+drgt)
             dsgn(i,j) = sign(ONE, dcen(i,j))
             if (dlft*drgt .ge. ZERO) then
                dlim(i,j) = TWO * min( abs(dlft), abs(drgt) )
             else
                dlim(i,j) = ZERO
             endif
             df(i,j) = dsgn(i,j)*min( dlim(i,j), abs(dcen(i,j)) )
          enddo

          ! Now limited fourth order slopes
          do i = ilo1-1, ihi1+1
             dp1         = FOUR3RD*dcen(i,j) - SIXTH*(df(i+1,j) + df(i-1,j))
             dqx(i,j,kc,QPRES) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dp1))
             dqx(i,j,kc,QPRES) = dqx(i,j,kc,QPRES) + q(i,j,k3d,QRHO)*src(i,j,k3d,QU)*dx(1)
          enddo
       enddo

#if (BL_SPACEDIM >= 2)
       ! Compute slopes in second coordinate direction
       do i = ilo1-1, ihi1+1

          ! First compute Fromm slopes
          do j = ilo2-2, ihi2+2
             dlft = q(i,j  ,k3d,QPRES) - q(i,j-1,k3d,QPRES)
             drgt = q(i,j+1,k3d,QPRES) - q(i,j  ,k3d,QPRES)

             ! Subtract off (rho * acceleration) so as not to limit that part of the slope
             dlft = dlft - FOURTH * &
                  (q(i,j,k3d,QRHO) + q(i,j-1,k3d,QRHO))*(src(i,j,k3d,QV)+src(i,j-1,k3d,QV))*dx(2)
             drgt = drgt - FOURTH * &
                  (q(i,j,k3d,QRHO) + q(i,j+1,k3d,QRHO))*(src(i,j,k3d,QV)+src(i,j+1,k3d,QV))*dx(2)

             dcen(i,j) = HALF*(dlft+drgt)
             dsgn(i,j) = sign( ONE, dcen(i,j) )
             if (dlft*drgt .ge. ZERO) then
                dlim(i,j) = TWO * min( abs(dlft), abs(drgt) )
             else
                dlim(i,j) = ZERO
             endif
             df(i,j) = dsgn(i,j)*min( dlim(i,j),abs(dcen(i,j)) )
          enddo

          ! Now limited fourth order slopes
          do j = ilo2-1, ihi2+1
             dp1 = FOUR3RD*dcen(i,j) - SIXTH*( df(i,j+1) + df(i,j-1) )
             dqy(i,j,kc,QPRES) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dp1))
             dqy(i,j,kc,QPRES) = dqy(i,j,kc,QPRES) + q(i,j,k3d,QRHO)*src(i,j,k3d,QV)*dx(2)
          enddo
       enddo
#endif

#if (BL_SPACEDIM == 3)
       ! Compute slopes in third coordinate direction
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! compute Fromm slopes on slab below
             k = k3d-1
             dm = q(i,j,k  ,QPRES) - q(i,j,k-1,QPRES)
             dp = q(i,j,k+1,QPRES) - q(i,j,k  ,QPRES)
             dm = dm - FOURTH * (q(i,j,k,QRHO) + q(i,j,k-1,QRHO))* &
                  (src(i,j,k,QW)+src(i,j,k-1,QW))*dx(3)
             dp = dp - FOURTH * (q(i,j,k,QRHO) + q(i,j,k+1,QRHO))* &
                  (src(i,j,k,QW)+src(i,j,k+1,QW))*dx(3)
             dc = HALF*(dm+dp)
             ds = sign( ONE, dc )
             if (dm*dp .ge. ZERO) then
                dl = TWO * min( abs(dm), abs(dp) )
             else
                dl = ZERO
             endif
             dfm = ds*min(dl,abs(dc))

             ! compute Fromm slopes on slab above
             k = k3d+1
             dm = q(i,j,k  ,QPRES) - q(i,j,k-1,QPRES)
             dp = q(i,j,k+1,QPRES) - q(i,j,k  ,QPRES)
             dm = dm - FOURTH * (q(i,j,k,QRHO) + q(i,j,k-1,QRHO))* &
                  (src(i,j,k,QW)+src(i,j,k-1,QW))*dx(3)
             dp = dp - FOURTH * (q(i,j,k,QRHO) + q(i,j,k+1,QRHO))* &
                  (src(i,j,k,QW)+src(i,j,k+1,QW))*dx(3)
             dc = HALF*(dm+dp)
             ds = sign( ONE, dc )
             if (dm*dp .ge. ZERO) then
                dl = TWO * min( abs(dm), abs(dp) )
             else
                dl = ZERO
             endif
             dfp = ds*min(dl,abs(dc))

             ! compute Fromm slopes on current slab
             k = k3d
             dm = q(i,j,k  ,QPRES) - q(i,j,k-1,QPRES)
             dp = q(i,j,k+1,QPRES) - q(i,j,k  ,QPRES)
             dm = dm - FOURTH * (q(i,j,k,QRHO) + q(i,j,k-1,QRHO))* &
                  (src(i,j,k,QW)+src(i,j,k-1,QW))*dx(3)
             dp = dp - FOURTH * (q(i,j,k,QRHO) + q(i,j,k+1,QRHO))* &
                  (src(i,j,k,QW)+src(i,j,k+1,QW))*dx(3)
             dc = HALF*(dm+dp)
             ds = sign( ONE, dc )
             if (dm*dp .ge. ZERO) then
                dl = TWO * min( abs(dm), abs(dp) )
             else
                dl = ZERO
             endif

             ! now limited fourth order slopes
             dp1 = FOUR3RD*dc - SIXTH*( dfp + dfm )
             dqz(i,j,kc,QPRES) = flatn(i,j,k3d)*ds*min(dl,abs(dp1))
             dqz(i,j,kc,QPRES) = dqz(i,j,kc,QPRES) + q(i,j,k3d,QRHO)*src(i,j,k3d,QW)*dx(3)
          enddo
       enddo
#endif

    endif

    call bl_deallocate(dsgn)
    call bl_deallocate(dlim)
    call bl_deallocate(  df)
    call bl_deallocate(dcen)

  end subroutine pslope

end module slope_module
