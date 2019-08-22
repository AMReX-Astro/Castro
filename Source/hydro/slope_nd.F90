module slope_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine uslope(lo, hi, idir, &
                    q, qd_lo, qd_hi, n, &
                    flatn, f_lo, f_hi, &
                    dq, qpd_lo, qpd_hi, &
                    dx, domlo, domhi)

    use meth_params_module, only: NQ, plm_iorder, QU, QPRES, QRHO, &
                                  const_grav, plm_well_balanced
    use amrex_constants_module, only: ZERO, HALF, ONE, TWO, FOUR3RD, FOURTH, SIXTH
    use prob_params_module, only : Symmetry, physbc_lo, physbc_hi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: n

    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(inout) :: dq(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in), value :: idir
    real(rt), intent(in) :: dx(3)
    integer, intent(in) :: domlo(3), domhi(3)

    integer :: i, j, k

    real(rt) :: dlft, drgt, dsgn, dcen, dlim, slop, dq1
    real(rt) :: dlftp1, drgtp1, dfp1
    real(rt) :: dlftm1, drgtm1, dfm1

    real(rt) :: qm2, qm1, q0, qp1, qp2
    real(rt) :: pp1, p0, pm1

    !$gpu

    if (plm_iorder == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dq(i,j,k,n) = ZERO
             end do
          end do
       end do

    else

       if (idir == 1) then
          ! Compute slopes in first coordinate direction
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   if (plm_well_balanced == 1 .and. n == QPRES .and. idir == AMREX_SPACEDIM) then
                      ! we'll only do a second-order pressure slope,
                      ! but we'll follow the well-balanced scheme of
                      ! Kappeli.  Note at the moment we are assuming
                      ! constant gravity.
                      p0 = q(i,j,k,QPRES)
                      pp1 = q(i+1,j,k,QPRES) - (p0 + HALF*dx(1)*(q(i,j,k,QRHO) + q(i+1,j,k,QRHO))*const_grav)
                      pm1 = q(i-1,j,k,QPRES) - (p0 - HALF*dx(1)*(q(i,j,k,QRHO) + q(i-1,j,k,QRHO))*const_grav)

                      ! if (i == domlo(1) .and. physbc_lo(1) == Symmetry) then
                      !    pm1 = ZERO  ! HSE is perfectly satisfied
                      ! end if

                      ! if (i == domhi(1) .and. physbc_hi(1) == Symmetry) then
                      !    pp1 = ZERO
                      ! end if

                      dlft = TWO*(p0 - pm1)
                      drgt = TWO*(pp1 - p0)
                      dcen = FOURTH * (dlft + drgt)
                      dsgn = sign(ONE, dcen)
                      slop = min(abs(dlft), abs(drgt))
                      if (dlft*drgt >= ZERO) then
                         dlim = slop
                      else
                         dlim = ZERO
                      end if

                      dq(i,j,k,n) = flatn(i,j,k)*dsgn*min(dlim, abs(dcen))

                   else

                      ! the 4th order MC limiter

                      qm2 = q(i-2,j,k,n)
                      qm1 = q(i-1,j,k,n)
                      q0 = q(i,j,k,n)
                      qp1 = q(i+1,j,k,n)
                      qp2 = q(i+2,j,k,n)

                      ! special consideration for reflecting BCs -- see
                      ! Saltzmann p. 162 (but note that Saltzmann has a
                      ! sign error)
                      if (i == domlo(1) .and. n == QU .and. physbc_lo(1) == Symmetry) then
                         qm2 = -qp1
                         qm1 = -3.0_rt*q0 + qp1 - 0.125_rt*(qp2 + qp1)
                      end if

                      if (i == domhi(1) .and. n == QU .and. physbc_hi(1) == Symmetry) then
                         qp2 = -qm1
                         qp1 = -3.0_rt*q0 + qm1 - 0.125_rt*(qm2 + qm1)
                      end if

                      ! First compute Fromm slopes

                      ! df at i+1
                      dlftp1 = TWO*(qp1 - q0)
                      drgtp1 = TWO*(qp2 - qp1)
                      dcen = FOURTH * (dlftp1 + drgtp1)
                      dsgn = sign(ONE, dcen)
                      slop = min(abs(dlftp1), abs(drgtp1))
                      if (dlftp1*drgtp1 >= ZERO) then
                         dlim = slop
                      else
                         dlim = ZERO
                      end if
                      dfp1 = dsgn*min(dlim, abs(dcen))

                      ! df at i-1
                      dlftm1 = TWO*(qm1 - qm2)
                      drgtm1 = TWO*(q0 - qm1)
                      dcen = FOURTH * (dlftm1 + drgtm1)
                      dsgn = sign(ONE, dcen)
                      slop = min(abs(dlftm1), abs(drgtm1))
                      if (dlftm1*drgtm1 >= ZERO) then
                         dlim = slop
                      else
                         dlim = ZERO
                      end if
                      dfm1 = dsgn*min(dlim, abs(dcen))

                      ! Now compute limited fourth order slopes at i
                      dlft = drgtm1
                      drgt = dlftp1
                      dcen = FOURTH * (dlft + drgt)
                      dsgn = sign(ONE, dcen)
                      slop = min(abs(dlft), abs(drgt))
                      if (dlft*drgt >= ZERO) then
                         dlim = slop
                      else
                         dlim = ZERO
                      end if

                      dq1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                      dq(i,j,k,n) = flatn(i,j,k)*dsgn*min(dlim, abs(dq1))

                   end if

                end do
             end do
          end do

#if (AMREX_SPACEDIM >= 2)
       else if (idir == 2) then

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
                   if (dlftp1*drgtp1 >= ZERO) then
                      dlim = slop
                   else
                      dlim = ZERO
                   end if
                   dfp1 = dsgn*min(dlim, abs(dcen))

                   ! df at j-1
                   dlftm1 = TWO*(q(i,j-1,k,n) - q(i,j-2,k,n))
                   drgtm1 = TWO*(q(i,j,k,n) - q(i,j-1,k,n))
                   dcen = FOURTH * (dlftm1 + drgtm1)
                   dsgn = sign(ONE, dcen)
                   slop = min(abs(dlftm1), abs(drgtm1))
                   if (dlftm1*drgtm1 >= ZERO) then
                      dlim = slop
                   else
                      dlim = ZERO
                   end if
                   dfm1 = dsgn*min(dlim, abs(dcen))

                   ! Now compute limited fourth order slopes at j
                   dlft = drgtm1
                   drgt = dlftp1
                   dcen = FOURTH * (dlft + drgt)
                   dsgn = sign(ONE, dcen)
                   slop = min(abs(dlft), abs(drgt))
                   if (dlft*drgt >= ZERO) then
                      dlim = slop
                   else
                      dlim = ZERO
                   end if

                   dq1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                   dq(i,j,k,n) = flatn(i,j,k)*dsgn*min(dlim, abs(dq1))

                end do
             end do
          end do
#endif

#if (AMREX_SPACEDIM == 3)
       else

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
                   if (dlftp1*drgtp1 >= ZERO) then
                      dlim = slop
                   else
                      dlim = ZERO
                   end if
                   dfp1 = dsgn*min(dlim, abs(dcen))

                   ! df at k-1
                   dlftm1 = TWO*(q(i,j,k-1,n) - q(i,j,k-2,n))
                   drgtm1 = TWO*(q(i,j,k,n) - q(i,j,k-1,n))
                   dcen = FOURTH * (dlftm1 + drgtm1)
                   dsgn = sign(ONE, dcen)
                   slop = min(abs(dlftm1), abs(drgtm1))
                   if (dlftm1*drgtm1 >= ZERO) then
                      dlim = slop
                   else
                      dlim = ZERO
                   end if
                   dfm1 = dsgn*min(dlim, abs(dcen))

                   ! Now compute limited fourth order slopes at k
                   dlft = drgtm1
                   drgt = dlftp1
                   dcen = FOURTH * (dlft + drgt)
                   dsgn = sign(ONE, dcen)
                   slop = min(abs(dlft), abs(drgt))
                   if (dlft*drgt >= ZERO) then
                      dlim = slop
                   else
                      dlim = ZERO
                   end if

                   dq1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                   dq(i,j,k,n) = flatn(i,j,k)*dsgn*min(dlim, abs(dq1))
                end do
             end do
          end do
#endif
       end if

    end if

  end subroutine uslope

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine pslope(lo, hi, idir, &
                    q, q_lo, q_hi, &
                    flatn, f_lo, f_hi, &
                    dq, qpd_lo, qpd_hi, &
                    src, src_lo, src_hi, &
                    dx)

    use meth_params_module, only : QRHO, QPRES, QU, QV, QW, NQ, NQSRC, plm_iorder
    use amrex_constants_module, only : ZERO, FOURTH, FOUR3RD, HALF, TWO, ONE, SIXTH

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qpd_lo(3),qpd_hi(3)
    integer, intent(in) :: src_lo(3),src_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in), value :: idir

    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(inout) :: dq(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NQSRC)
    real(rt), intent(in) :: dx(3)

    integer :: i, j, k

    real(rt) :: dlft, drgt, dp1

    !     Local arrays
    real(rt) :: dsgn, dlim, dcen
    real(rt) :: dfp1, dfm1, dlftp1, drgtp1, dlftm1, drgtm1

    !$gpu

    if (plm_iorder == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dq(i,j,k,QPRES) = ZERO
             end do
          end do
       end do

    else
       ! Compute slopes in first coordinate direction
       if (idir == 1) then
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
                   if (dlftp1*drgtp1 >= ZERO) then
                      dlim = TWO * min(abs(dlftp1), abs(drgtp1))
                   else
                      dlim = ZERO
                   end if
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
                   if (dlftm1*drgtm1 >= ZERO) then
                      dlim = TWO * min(abs(dlftm1), abs(drgtm1))
                   else
                      dlim = ZERO
                   end if
                   dfm1 = dsgn*min(dlim, abs(dcen))

                   ! Now limited fourth order slopes at i
                   dlft = drgtm1
                   drgt = dlftp1
                   dcen = HALF*(dlft + drgt)
                   dsgn = sign(ONE, dcen)
                   if (dlft*drgt >= ZERO) then
                      dlim = TWO * min(abs(dlft), abs(drgt))
                   else
                      dlim = ZERO
                   end if

                   dp1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                   dq(i,j,k,QPRES) = flatn(i,j,k)*dsgn*min(dlim, abs(dp1))
                   dq(i,j,k,QPRES) = dq(i,j,k,QPRES) + q(i,j,k,QRHO)*src(i,j,k,QU)*dx(1)

                end do
             end do
          end do

#if AMREX_SPACEDIM >= 2
       else if (idir == 2) then

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
                   if (dlftp1*drgtp1 >= ZERO) then
                      dlim = TWO * min(abs(dlftp1), abs(drgtp1))
                   else
                      dlim = ZERO
                   end if
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
                   if (dlftm1*drgtm1 >= ZERO) then
                      dlim = TWO * min(abs(dlftm1), abs(drgtm1))
                   else
                      dlim = ZERO
                   end if
                   dfm1 = dsgn*min(dlim, abs(dcen))

                   ! Now limited fourth order slopes at j
                   dlft = drgtm1
                   drgt = dlftp1

                   dcen = HALF*(dlft+drgt)
                   dsgn = sign(ONE, dcen)
                   if (dlft*drgt >= ZERO) then
                      dlim = TWO * min(abs(dlft), abs(drgt))
                   else
                      dlim = ZERO
                   end if

                   dp1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                   dq(i,j,k,QPRES) = flatn(i,j,k)*dsgn*min(dlim, abs(dp1))
                   dq(i,j,k,QPRES) = dq(i,j,k,QPRES) + q(i,j,k,QRHO)*src(i,j,k,QV)*dx(2)
                end do
             end do
          end do
#endif

#if AMREX_SPACEDIM == 3
       else

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
                   if (dlftp1*drgtp1 >= ZERO) then
                      dlim = TWO * min(abs(dlftp1), abs(drgtp1))
                   else
                      dlim = ZERO
                   end if
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
                   if (dlftm1*drgtm1 >= ZERO) then
                      dlim = TWO * min(abs(dlftm1), abs(drgtm1))
                   else
                      dlim = ZERO
                   end if
                   dfm1 = dsgn*min(dlim, abs(dcen))

                   ! now limited fourth order slopes at k
                   dlft = drgtm1
                   drgt = dlftp1

                   dcen = HALF*(dlft+drgt)
                   dsgn = sign(ONE, dcen)
                   if (dlft*drgt >= ZERO) then
                      dlim = TWO * min(abs(dlft), abs(drgt))
                   else
                      dlim = ZERO
                   end if

                   dp1 = FOUR3RD*dcen - SIXTH*(dfp1 + dfm1)
                   dq(i,j,k,QPRES) = flatn(i,j,k)*dsgn*min(dlim, abs(dp1))
                   dq(i,j,k,QPRES) = dq(i,j,k,QPRES) + q(i,j,k,QRHO)*src(i,j,k,QW)*dx(3)
                end do
             end do
          end do
#endif
       end if
    end if

  end subroutine pslope

end module slope_module
