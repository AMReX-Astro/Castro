module flatten_module

  implicit none

contains

  AMREX_DEVICE subroutine ca_uflaten(lo, hi, q, q_lo, q_hi, flatn, f_lo, f_hi) bind(c,name='ca_uflaten')

    use bl_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dg
    use meth_params_module, only: NQ, QU, QV, QW, QPRES

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(inout) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))

    integer :: i, j, k, ishft

    real(rt), parameter :: small_pres = 1.e-200_rt

    real(rt) :: denom, zeta, tst, tmp, ftmp
    real(rt) :: dp, z, z2, chi, chi2

    ! Knobs for detection of strong shock
    real(rt), parameter :: shktst = 0.33e0_rt, zcut1 = 0.75e0_rt, zcut2 = 0.85e0_rt, dzcut = ONE/(zcut2-zcut1)

    ! x-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             dp = q(i+1,j,k,QPRES) - q(i-1,j,k,QPRES)

             if (dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q(i+2,j,k,QPRES)-q(i-2,j,k,QPRES)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i-1,j,k,QU)-q(i+1,j,k,QU) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i+1,j,k,QPRES), q(i-1,j,k,QPRES))

             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q(i+1-ishft,j,k,QPRES) - q(i-1-ishft,j,k,QPRES)

             denom = max(small_pres, abs(q(i+2-ishft,j,k,QPRES)-q(i-2-ishft,j,k,QPRES)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i-1-ishft,j,k,QU)-q(i+1-ishft,j,k,QU) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i+1-ishft,j,k,QPRES), q(i-1-ishft,j,k,QPRES))

             if ((abs(dp)/tmp) .gt. shktst) then
                chi2 = tst
             else
                chi2 = ZERO
             endif

             flatn(i,j,k) = ONE - max(chi2 * z2, chi * z)

          end do
       end do
    end do

    ! y-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             dp = q(i,j+1,k,QPRES) - q(i,j-1,k,QPRES)

             if (dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q(i,j+2,k,QPRES)-q(i,j-2,k,QPRES)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j-1,k,QV)-q(i,j+1,k,QV) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j+1,k,QPRES), q(i,j-1,k,QPRES))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q(i,j+1-ishft,k,QPRES) - q(i,j-1-ishft,k,QPRES)

             denom = max(small_pres, abs(q(i,j+2-ishft,k,QPRES)-q(i,j-2-ishft,k,QPRES)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j-1-ishft,k,QV)-q(i,j+1-ishft,k,QV) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j+1-ishft,k,QPRES), q(i,j-1-ishft,k,QPRES))
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

    ! z-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             dp = q(i,j,k+1,QPRES) - q(i,j,k-1,QPRES)

             if(dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q(i,j,k+2,QPRES)-q(i,j,k-2,QPRES)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j,k-1,QW)-q(i,j,k+1,QW) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j,k+1,QPRES),q(i,j,k-1,QPRES))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q(i,j,k+1-ishft,QPRES) - q(i,j,k-1-ishft,QPRES)

             denom = max(small_pres, abs(q(i,j,k+2-ishft,QPRES)-q(i,j,k-2-ishft,QPRES)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j,k-1-ishft,QW)-q(i,j,k+1-ishft,QW) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j,k+1-ishft,QPRES),q(i,j,k-1-ishft,QPRES))
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

  end subroutine ca_uflaten

end module flatten_module
