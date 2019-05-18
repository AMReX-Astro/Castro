module flatten_module

  use amrex_constants_module, only : ZERO
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

! :::
! ::: ------------------------------------------------------------------
! :::

#ifdef RADIATION
  subroutine ca_rad_flatten(lo, hi, &
                            q_core, qc_lo, qc_hi, &
                            q_rad, qr_lo, qr_hi, &
                            flatn, f_lo, f_hi, &
                            flatg, fg_lo, fg_hi) bind(C, name="ca_rad_flatten")

    use meth_params_module, only : QPRES, QPTOT, QU, QV, QW, flatten_pp_threshold, NQC, NQR

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: qc_lo(3), qc_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: fg_lo(3), fg_hi(3)

    real(rt), intent(in) :: q_core(qc_lo(1):qc_hi(1),qc_lo(2):qc_hi(2),qc_lo(3):qc_hi(3),NQC)
    real(rt), intent(in) :: q_rad(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),NQR)
    real(rt), intent(inout) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(inout) :: flatg(fg_lo(1):fg_hi(1),fg_lo(2):fg_hi(2),fg_lo(3):fg_hi(3))

    integer :: i, j, k

    call ca_uflatten(lo, hi, &
                     q_core, qc_lo, qc_hi, &
                     q_core, qc_lo, qc_hi, NQC, QPRES, &
                     flatn, f_lo, f_hi)

    call ca_uflatten(lo, hi, &
                     q_core, qc_lo, qc_hi, &
                     q_rad, qr_lo, qr_hi, NQR, QPTOT, &
                     flatg, fg_lo, fg_hi)


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             flatn(i,j,k) = flatn(i,j,k) * flatg(i,j,k)

             if (flatten_pp_threshold > ZERO) then
                if ( q_core(i-1,j,k,QU) + q_core(i,j-1,k,QV) + q_core(i,j,k-1,QW) > &
                     q_core(i+1,j,k,QU) + q_core(i,j+1,k,QV) + q_core(i,j,k+1,QW) ) then

                   if (q_core(i,j,k,QPRES) < flatten_pp_threshold * q_rad(i,j,k,QPTOT)) then
                      flatn(i,j,k) = ZERO
                   end if

                end if
             endif

          end do
       end do
    end do

  end subroutine ca_rad_flatten
#endif


  subroutine ca_uflatten(lo, hi, &
                         q_core, qc_lo, qc_hi, &
                         q_pres, qp_lo, qp_hi, ncomp, pres_comp, &
                         flatn, f_lo, f_hi) bind(c,name='ca_uflatten')

    ! here, q_pres is the primitive variable array that has the
    ! pressure we wish to consider.  It has ncomp components and the
    ! pressure is component pres_comp

    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQC, QU, QV, QW

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: qc_lo(3), qc_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: q_core(qc_lo(1):qc_hi(1),qc_lo(2):qc_hi(2),qc_lo(3):qc_hi(3),NQC)
    real(rt), intent(in   ) :: q_pres(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),ncomp)
    real(rt), intent(inout) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    integer, intent(in), value :: pres_comp, ncomp

    integer :: i, j, k, ishft

    real(rt), parameter :: small_pres = 1.e-200_rt

    real(rt) :: denom, zeta, tst, tmp, ftmp
    real(rt) :: dp, z, z2, chi, chi2

    ! Knobs for detection of strong shock
    real(rt), parameter :: shktst = 0.33e0_rt, zcut1 = 0.75e0_rt, zcut2 = 0.85e0_rt, dzcut = ONE/(zcut2-zcut1)

    !$gpu

    ! x-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             dp = q_pres(i+1,j,k,pres_comp) - q_pres(i-1,j,k,pres_comp)

             if (dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q_pres(i+2,j,k,pres_comp) - q_pres(i-2,j,k,pres_comp)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q_core(i-1,j,k,QU) - q_core(i+1,j,k,QU) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q_pres(i+1,j,k,pres_comp), q_pres(i-1,j,k,pres_comp))

             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q_pres(i+1-ishft,j,k,pres_comp) - q_pres(i-1-ishft,j,k,pres_comp)

             denom = max(small_pres, abs(q_pres(i+2-ishft,j,k,pres_comp) - q_pres(i-2-ishft,j,k,pres_comp)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q_core(i-1-ishft,j,k,QU) - q_core(i+1-ishft,j,k,QU) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q_pres(i+1-ishft,j,k,pres_comp), q_pres(i-1-ishft,j,k,pres_comp))

             if ((abs(dp)/tmp) .gt. shktst) then
                chi2 = tst
             else
                chi2 = ZERO
             endif

             flatn(i,j,k) = ONE - max(chi2 * z2, chi * z)

          end do
       end do
    end do

#if AMREX_SPACEDIM >= 2
    ! y-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             dp = q_pres(i,j+1,k,pres_comp) - q_pres(i,j-1,k,pres_comp)

             if (dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q_pres(i,j+2,k,pres_comp) - q_pres(i,j-2,k,pres_comp)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q_core(i,j-1,k,QV) - q_core(i,j+1,k,QV) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q_pres(i,j+1,k,pres_comp), q_pres(i,j-1,k,pres_comp))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q_pres(i,j+1-ishft,k,pres_comp) - q_pres(i,j-1-ishft,k,pres_comp)

             denom = max(small_pres, abs(q_pres(i,j+2-ishft,k,pres_comp) - q_pres(i,j-2-ishft,k,pres_comp)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q_core(i,j-1-ishft,k,QV) - q_core(i,j+1-ishft,k,QV) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q_pres(i,j+1-ishft,k,pres_comp), q_pres(i,j-1-ishft,k,pres_comp))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi2 = tst
             else
                chi2 = ZERO
             endif

             ftmp = ONE - max(chi2 * z2, chi * z)
             flatn(i,j,k) = min(flatn(i,j,k), ftmp)

          end do
       end do
    end do
#endif

#if AMREX_SPACEDIM == 3
    ! z-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             dp = q_pres(i,j,k+1,pres_comp) - q_pres(i,j,k-1,pres_comp)

             if(dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q_pres(i,j,k+2,pres_comp) - q_pres(i,j,k-2,pres_comp)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q_core(i,j,k-1,QW) - q_core(i,j,k+1,QW) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q_pres(i,j,k+1,pres_comp), q_pres(i,j,k-1,pres_comp))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q_pres(i,j,k+1-ishft,pres_comp) - q_pres(i,j,k-1-ishft,pres_comp)

             denom = max(small_pres, abs(q_pres(i,j,k+2-ishft,pres_comp) - q_pres(i,j,k-2-ishft,pres_comp)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q_core(i,j,k-1-ishft,QW) - q_core(i,j,k+1-ishft,QW) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q_pres(i,j,k+1-ishft,pres_comp), q_pres(i,j,k-1-ishft,pres_comp))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi2 = tst
             else
                chi2 = ZERO
             endif

             ftmp = ONE - max(chi2 * z2, chi * z)
             flatn(i,j,k) = min(flatn(i,j,k), ftmp)
          enddo
       enddo
    enddo
#endif

  end subroutine ca_uflatten

end module flatten_module
