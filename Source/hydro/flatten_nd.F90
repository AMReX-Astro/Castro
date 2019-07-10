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
                            q, q_lo, q_hi, &
                            flatn, f_lo, f_hi, &
                            flatg, fg_lo, fg_hi) bind(C, name="ca_rad_flatten")

    use meth_params_module, only : QPRES, QU, QV, QW, flatten_pp_threshold, QPTOT, NQ

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: fg_lo(3), fg_hi(3)

    real(rt)        , intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt)        , intent(inout) :: flatg(fg_lo(1):fg_hi(1),fg_lo(2):fg_hi(2),fg_lo(3):fg_hi(3))

    integer :: i, j, k

    call ca_uflatten(lo, hi, &
                     q, q_lo, q_hi, &
                     flatn, f_lo, f_hi, QPRES)

    call ca_uflatten(lo, hi, &
                     q, q_lo, q_hi, &
                     flatg, fg_lo, fg_hi, QPTOT)


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             flatn(i,j,k) = flatn(i,j,k) * flatg(i,j,k)

             if (flatten_pp_threshold > ZERO) then
                if ( q(i-1,j,k,QU) + q(i,j-1,k,QV) + q(i,j,k-1,QW) > &
                     q(i+1,j,k,QU) + q(i,j+1,k,QV) + q(i,j,k+1,QW) ) then

                   if (q(i,j,k,QPRES) < flatten_pp_threshold * q(i,j,k,QPTOT)) then
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
                         q, q_lo, q_hi, &
                         flatn, f_lo, f_hi, pres_comp) bind(c,name='ca_uflatten')

    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQ, QU, QV, QW

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(inout) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    integer, intent(in), value :: pres_comp

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

             dp = q(i+1,j,k,pres_comp) - q(i-1,j,k,pres_comp)

             if (dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q(i+2,j,k,pres_comp)-q(i-2,j,k,pres_comp)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i-1,j,k,QU)-q(i+1,j,k,QU) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i+1,j,k,pres_comp), q(i-1,j,k,pres_comp))

             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q(i+1-ishft,j,k,pres_comp) - q(i-1-ishft,j,k,pres_comp)

             denom = max(small_pres, abs(q(i+2-ishft,j,k,pres_comp)-q(i-2-ishft,j,k,pres_comp)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i-1-ishft,j,k,QU)-q(i+1-ishft,j,k,QU) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i+1-ishft,j,k,pres_comp), q(i-1-ishft,j,k,pres_comp))

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

             dp = q(i,j+1,k,pres_comp) - q(i,j-1,k,pres_comp)

             if (dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q(i,j+2,k,pres_comp)-q(i,j-2,k,pres_comp)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j-1,k,QV)-q(i,j+1,k,QV) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j+1,k,pres_comp), q(i,j-1,k,pres_comp))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q(i,j+1-ishft,k,pres_comp) - q(i,j-1-ishft,k,pres_comp)

             denom = max(small_pres, abs(q(i,j+2-ishft,k,pres_comp)-q(i,j-2-ishft,k,pres_comp)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j-1-ishft,k,QV)-q(i,j+1-ishft,k,QV) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j+1-ishft,k,pres_comp), q(i,j-1-ishft,k,pres_comp))
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

             dp = q(i,j,k+1,pres_comp) - q(i,j,k-1,pres_comp)

             if(dp .gt. ZERO) then
                ishft = 1
             else
                ishft = -1
             endif

             denom = max(small_pres, abs(q(i,j,k+2,pres_comp)-q(i,j,k-2,pres_comp)))
             zeta = abs(dp) / denom
             z = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j,k-1,QW)-q(i,j,k+1,QW) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j,k+1,pres_comp),q(i,j,k-1,pres_comp))
             if ((abs(dp)/tmp) .gt. shktst) then
                chi = tst
             else
                chi = ZERO
             endif

             dp = q(i,j,k+1-ishft,pres_comp) - q(i,j,k-1-ishft,pres_comp)

             denom = max(small_pres, abs(q(i,j,k+2-ishft,pres_comp)-q(i,j,k-2-ishft,pres_comp)))
             zeta = abs(dp) / denom
             z2 = min(ONE, max(ZERO, dzcut * (zeta - zcut1)))

             if (q(i,j,k-1-ishft,QW)-q(i,j,k+1-ishft,QW) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif

             tmp = min(q(i,j,k+1-ishft,pres_comp),q(i,j,k-1-ishft,pres_comp))
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
