module slope_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public uslope, pslope

contains

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine uslope(q, qd_lo, qd_hi, n, &
                    flatn, f_lo, f_hi, &
                    dqx, &
#if AMREX_SPACEDIM >= 2
                    dqy, &
#endif
#if AMREX_SPACEDIM == 3
                    dqz, &
#endif
                    qpd_lo, qpd_hi, &
                    lo, hi)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only: NQ, plm_iorder
    use amrex_constants_module, only: ZERO, HALF, ONE, TWO, FOUR3RD, FOURTH, SIXTH

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: n

    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(inout) :: dqx(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: dqy(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: dqz(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
#endif
    integer, intent(in) :: lo(3), hi(3)

    integer :: i, j, k

    real(rt) :: dlft, drgt, dsgn, dcen, dlim, slop, dq1
    real(rt) :: dlftp1, drgtp1, dfp1
    real(rt) :: dlftm1, drgtm1, dfm1

    if (plm_iorder == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dqx(i,j,k,n) = ZERO
#if AMREX_SPACEDIM >= 2
                dqy(i,j,k,n) = ZERO
#endif
#if AMREX_SPACEDIM == 3
                dqz(i,j,k,n) = ZERO
#endif
             end do
          end do
       end do

    else

       ! Compute slopes in first coordinate direction
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! First compute Fromm slopes

                ! df at i+1
                dlftp1 = TWO*(q(i+1,j,k,n) - q(i,j,k,n))
                drgtp1 = TWO*(q(i+2,j,k,n) - q(i+1,j,k,n))
                dcen = FOURTH * (dlftp1 + drgtp1)
                dsgn = sign(ONE, dcen)
                slop = min(abs(dlftp1), abs(drgtp1))
                dlim = merge(slop, ZERO, dlftp1*drgtp1 >= ZERO)
                dfp1 = dsgn*min(dlim, abs(dcen))

                ! df at i-1
                dlftm1 = TWO*(q(i-1,j,k,n) - q(i-2,j,k,n))
                drgtm1 = TWO*(q(i,j,k,n) - q(i-1,j,k,n))
                dcen = FOURTH * (dlftm1 + drgtm1)
                dsgn = sign(ONE, dcen)
                slop = min(abs(dlftm1), abs(drgtm1))
                dlim = merge(slop, ZERO, dlftm1*drgtm1 >= ZERO)
                dfm1 = dsgn*min(dlim, abs(dcen))

                ! Now compute limited fourth order slopes at i
                dlft = drgtm1
                drgt = dlftp1
                dcen = FOURTH * (dlft + drgt)
                dsgn = sign(ONE, dcen)
                slop = min(abs(dlft), abs(drgt))
                dlim = merge(slop, ZERO, dlft*drgt >= ZERO)

                dq1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                dqx(i,j,k,n) = flatn(i,j,k)*dsgn*min(dlim, abs(dq1))

             end do
          end do
       end do

#if (AMREX_SPACEDIM >= 2)
       ! Compute slopes in second coordinate direction
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! First compute Fromm slopes

                ! df at j+1
                dlftp1 = TWO*(q(i,j+1,k,n) - q(i,j,k,n))
                drgtp1 = TWO*(q(i,j+2,k,n) - q(i,j+1,k,n))
                dcen = FOURTH * (dlftp1 + drgtp1)
                dsgn = sign(ONE, dcen)
                slop = min(abs(dlftp1), abs(drgtp1))
                dlim = merge(slop, ZERO, dlftp1*drgtp1 >= ZERO)
                dfp1 = dsgn*min(dlim, abs(dcen))

                ! df at j-1
                dlftm1 = TWO*(q(i,j-1,k,n) - q(i,j-2,k,n))
                drgtm1 = TWO*(q(i,j,k,n) - q(i,j-1,k,n))
                dcen = FOURTH * (dlftm1 + drgtm1)
                dsgn = sign(ONE, dcen)
                slop = min(abs(dlftm1), abs(drgtm1))
                dlim = merge(slop, ZERO, dlftm1*drgtm1 >= ZERO)
                dfm1 = dsgn*min(dlim, abs(dcen))

                ! Now compute limited fourth order slopes at j
                dlft = drgtm1
                drgt = dlftp1
                dcen = FOURTH * (dlft + drgt)
                dsgn = sign(ONE, dcen)
                slop = min(abs(dlft), abs(drgt))
                dlim = merge(slop, ZERO, dlft*drgt >= ZERO)

                dq1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                dqy(i,j,k,n) = flatn(i,j,k)*dsgn*min(dlim, abs(dq1))

             end do
          end do
       end do
#endif

#if (AMREX_SPACEDIM == 3)
       ! Compute slopes in third coordinate direction
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! First compute Fromm slopes

                ! df at k+1
                dlftp1 = TWO*(q(i,j,k+1,n) - q(i,j,k,n))
                drgtp1 = TWO*(q(i,j,k+2,n) - q(i,j,k+1,n))
                dcen = FOURTH * (dlftp1 + drgtp1)
                dsgn = sign(ONE, dcen)
                slop = min(abs(dlftp1), abs(drgtp1))
                dlim = merge(slop, ZERO, dlftp1*drgtp1 >= ZERO)
                dfp1 = dsgn*min(dlim, abs(dcen))

                ! df at k-1
                dlftm1 = TWO*(q(i,j,k-1,n) - q(i,j,k-2,n))
                drgtm1 = TWO*(q(i,j,k,n) - q(i,j,k-1,n))
                dcen = FOURTH * (dlftm1 + drgtm1)
                dsgn = sign(ONE, dcen)
                slop = min(abs(dlftm1), abs(drgtm1))
                dlim = merge(slop, ZERO, dlftm1*drgtm1 >= ZERO)
                dfm1 = dsgn*min(dlim, abs(dcen))

                ! Now compute limited fourth order slopes at k
                dlft = drgtm1
                drgt = dlftp1
                dcen = FOURTH * (dlft + drgt)
                dsgn = sign(ONE, dcen)
                slop = min(abs(dlft), abs(drgt))
                dlim = merge(slop, ZERO, dlft*drgt >= ZERO)

                dq1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                dqz(i,j,k,n) = flatn(i,j,k)*dsgn*min(dlim, abs(dq1))
             end do
          end do
       end do
#endif
    end if

  end subroutine uslope

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine pslope(q, q_lo, q_hi, &
                    flatn, f_lo, f_hi, &
                    dqx, &
#if AMREX_SPACEDIM >= 2
                    dqy, &
#endif
#if AMREX_SPACEDIM == 3
                    dqz, &
#endif
                    qpd_lo, qpd_hi, &
                    src, src_lo, src_hi, &
                    lo, hi, dx)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QRHO, QPRES, QU, QV, QW, NQ, QVAR, plm_iorder
    use amrex_constants_module, only : ZERO, FOURTH, FOUR3RD, HALF, TWO, ONE, SIXTH

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qpd_lo(3),qpd_hi(3)
    integer, intent(in) :: src_lo(3),src_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(inout) :: dqx(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: dqy(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: dqz(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
#endif
    real(rt), intent(in) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    real(rt), intent(in) :: dx(3)

    integer :: i, j, k

    real(rt) :: dlft, drgt, dp1
    real(rt) :: dm, dp, dc, dl, dfm, dfp

    !     Local arrays
    real(rt) :: dsgn, dlim, dcen
    real(rt) :: dfp1, dfm1, dlftp1, drgtp1, dlftm1, drgtm1

    if (plm_iorder == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dqx(i,j,k,QPRES) = ZERO
#if AMREX_SPACEDIM >= 2
                dqy(i,j,k,QPRES) = ZERO
#endif
#if AMREX_SPACEDIM == 3
                dqz(i,j,k,QPRES) = ZERO
#endif
             end do
          end do
       end do

    else
       ! Compute slopes in first coordinate direction
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! First compute Fromm slopes

                ! df at i+1
                dlftp1 = q(i+1,j,k,QPRES) - q(i,j,k,QPRES)
                drgtp1 = q(i+2,j,k,QPRES) - q(i+1,j,k,QPRES)

                ! Subtract off (rho * acceleration) so as not to limit that part of the slope
                dlftp1 = dlftp1 - FOURTH * &
                     (q(i+1,j,k,QRHO) + q(i,j,k,QRHO)) * (src(i+1,j,k,QU)+src(i,j,k,QU))*dx(1)
                drgtp1 = drgtp1 - FOURTH * &
                     (q(i+1,j,k,QRHO) + q(i+2,j,k,QRHO)) * (src(i+1,j,k,QU)+src(i+2,j,k,QU))*dx(1)

                dcen = HALF*(dlftp1 + drgtp1)
                dsgn = sign(ONE, dcen)
                dlim = merge(TWO * min(abs(dlftp1), abs(drgtp1)), ZERO, dlftp1*drgtp1 >= ZERO)
                dfp1 = dsgn*min(dlim, abs(dcen))

                ! df at i-1
                dlftm1 = q(i-1,j,k,QPRES) - q(i-2,j,k,QPRES)
                drgtm1 = q(i,j,k,QPRES) - q(i-1,j,k,QPRES)

                ! Subtract off (rho * acceleration) so as not to limit that part of the slope
                dlftm1 = dlftm1 - FOURTH * &
                     (q(i-1,j,k,QRHO) + q(i-2,j,k,QRHO)) * (src(i-1,j,k,QU)+src(i-2,j,k,QU))*dx(1)
                drgtm1 = drgtm1 - FOURTH * &
                     (q(i-1,j,k,QRHO) + q(i,j,k,QRHO)) * (src(i-1,j,k,QU)+src(i,j,k,QU))*dx(1)

                dcen = HALF*(dlftm1 + drgtm1)
                dsgn = sign(ONE, dcen)
                dlim = merge(TWO * min(abs(dlftm1), abs(drgtm1)), ZERO, dlftm1*drgtm1 >= ZERO)
                dfm1 = dsgn*min(dlim, abs(dcen))

                ! Now limited fourth order slopes at i
                dlft = drgtm1
                drgt = dlftp1
                dcen = HALF*(dlft + drgt)
                dsgn = sign(ONE, dcen)
                dlim = merge(TWO * min(abs(dlft), abs(drgt)), ZERO, dlft*drgt >= ZERO)

                dp1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                dqx(i,j,k,QPRES) = flatn(i,j,k)*dsgn*min(dlim, abs(dp1))
                dqx(i,j,k,QPRES) = dqx(i,j,k,QPRES) + q(i,j,k,QRHO)*src(i,j,k,QU)*dx(1)

             end do
          end do
       end do

#if AMREX_SPACEDIM >= 2
       ! Compute slopes in second coordinate direction
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! First compute Fromm slopes

                ! dt at j+1
                dlftp1 = q(i,j+1,k,QPRES) - q(i,j,k,QPRES)
                drgtp1 = q(i,j+2,k,QPRES) - q(i,j+1,k,QPRES)

                ! Subtract off (rho * acceleration) so as not to limit that part of the slope
                dlftp1 = dlftp1 - FOURTH * &
                     (q(i,j+1,k,QRHO) + q(i,j,k,QRHO)) * (src(i,j+1,k,QV)+src(i,j,k,QV))*dx(2)
                drgtp1 = drgtp1 - FOURTH * &
                     (q(i,j+1,k,QRHO) + q(i,j+2,k,QRHO)) * (src(i,j+1,k,QV)+src(i,j+2,k,QV))*dx(2)

                dcen = HALF*(dlftp1 + drgtp1)
                dsgn = sign(ONE, dcen)
                dlim = merge(TWO * min(abs(dlftp1), abs(drgtp1)), ZERO, dlftp1*drgtp1 >= ZERO)
                dfp1 = dsgn*min(dlim, abs(dcen))

                ! df at j-1
                dlftm1 = q(i,j-1,k,QPRES) - q(i,j-2,k,QPRES)
                drgtm1 = q(i,j,k,QPRES) - q(i,j-1,k,QPRES)

                ! Subtract off (rho * acceleration) so as not to limit that part of the slope
                dlftm1 = dlftm1 - FOURTH * &
                     (q(i,j-1,k,QRHO) + q(i,j-2,k,QRHO)) * (src(i,j-1,k,QV)+src(i,j-2,k,QV))*dx(2)
                drgtm1 = drgtm1 - FOURTH * &
                     (q(i,j-1,k,QRHO) + q(i,j,k,QRHO)) * (src(i,j-1,k,QV)+src(i,j,k,QV))*dx(2)

                dcen = HALF*(dlftm1 + drgtm1)
                dsgn = sign(ONE, dcen)
                dlim = merge(TWO * min(abs(dlftm1), abs(drgtm1)), ZERO, dlftm1*drgtm1 >= ZERO)
                dfm1 = dsgn*min(dlim, abs(dcen))

                ! Now limited fourth order slopes at j
                dlft = drgtm1
                drgt = dlftp1

                dcen = HALF*(dlft+drgt)
                dsgn = sign(ONE, dcen)
                dlim = merge(TWO * min(abs(dlft), abs(drgt)), ZERO, dlft*drgt >= ZERO)

                dp1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                dqy(i,j,k,QPRES) = flatn(i,j,k)*dsgn*min(dlim, abs(dp1))
                dqy(i,j,k,QPRES) = dqy(i,j,k,QPRES) + q(i,j,k,QRHO)*src(i,j,k,QV)*dx(2)
             end do
          end do
       end do
#endif

#if AMREX_SPACEDIM == 3
       ! Compute slopes in third coordinate direction
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! First compute Fromm slopes

                ! df at k+1
                dlftp1 = q(i,j,k+1,QPRES) - q(i,j,k,QPRES)
                drgtp1 = q(i,j,k+2,QPRES) - q(i,j,k+1,QPRES)

                ! Subtract off (rho * acceleration)
                dlftp1 = dlftp1 - FOURTH * &
                     (q(i,j,k+1,QRHO) + q(i,j,k,QRHO)) * (src(i,j,k+1,QW)+src(i,j,k,QW))*dx(3)
                drgtp1 = drgtp1 - FOURTH * &
                     (q(i,j,k+1,QRHO) + q(i,j,k+2,QRHO)) * (src(i,j,k+1,QW)+src(i,j,k+2,QW))*dx(3)

                dcen = HALF*(dlftp1 + drgtp1)
                dsgn = sign(ONE, dcen)
                dlim = merge(TWO * min(abs(dlftp1), abs(drgtp1)), ZERO, dlftp1*drgtp1 >= ZERO)
                dfp1 = dsgn*min(dlim, abs(dcen))

                ! df at k-1
                dlftm1 = q(i,j,k-1,QPRES) - q(i,j,k-2,QPRES)
                drgtm1 = q(i,j,k,QPRES) - q(i,j,k-1,QPRES)

                ! Subtract off (rho * acceleration)
                dlftm1 = dlftm1 - FOURTH * &
                     (q(i,j,k-1,QRHO) + q(i,j,k-2,QRHO)) * (src(i,j,k-1,QW)+src(i,j,k-2,QW))*dx(3)
                drgtm1 = drgtm1 - FOURTH * &
                     (q(i,j,k-1,QRHO) + q(i,j,k,QRHO)) * (src(i,j,k-1,QW)+src(i,j,k,QW))*dx(3)

                dcen = HALF*(dlftm1 + drgtm1)
                dsgn = sign(ONE, dcen)
                dlim = merge(TWO * min( abs(dlftm1), abs(drgtm1) ), ZERO, dlftm1*drgtm1 >= ZERO)
                dfm1 = dsgn*min(dlim, abs(dcen))

                ! now limited fourth order slopes at k
                dlft = drgtm1
                drgt = dlftp1

                dcen = HALF*(dlft+drgt)
                dsgn = sign(ONE, dcen)
                dlim = merge(TWO * min(abs(dlft), abs(drgt)), ZERO, dlft*drgt >= ZERO)

                dp1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                dqz(i,j,k,QPRES) = flatn(i,j,k)*dsgn*min(dlim, abs(dp1))
                dqz(i,j,k,QPRES) = dqz(i,j,k,QPRES) + q(i,j,k,QRHO)*src(i,j,k,QW)*dx(3)
             end do
          end do
       end do
#endif

    endif

  end subroutine pslope

end module slope_module
