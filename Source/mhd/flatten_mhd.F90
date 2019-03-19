module flatten_module_mhd

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_constants_module, only : ZERO

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public :: uflatten
#ifdef RADIATION
  public :: rad_flatten
#endif
#ifdef MHD
  public :: mhd_flatten 
#endif
contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine uflatten(lo, hi, q, flatn, q_lo, q_hi, ipres)

    ! here, ipres is the pressure variable we want to consider jumps on
    ! passing it in allows
    use meth_params_module, only : small_pres, QU, QV, QW, NQ
    use prob_params_module, only : dg
    use amrex_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: ipres

    real(rt)        , intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ+1)
    real(rt)        , intent(inout) :: flatn(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))

    integer :: i, j, k, ishft

    real(rt)         :: denom, zeta, tst, tmp, ftmp

    ! Local arrays
    real(rt)        , pointer :: dp(:,:,:), z(:,:,:), chi(:,:,:)

    ! Knobs for detection of strong shock
    real(rt)        , parameter :: shktst = 0.33e0_rt, zcut1 = 0.75e0_rt, zcut2 = 0.85e0_rt, dzcut = ONE/(zcut2-zcut1)

    call bl_allocate(dp ,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(z  ,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(chi,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)

    ! x-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          !dir$ ivdep
          do i = lo(1)-1*dg(1),hi(1)+1*dg(1)
             dp(i,j,k) = q(i+1*dg(1),j,k,ipres) - q(i-1*dg(1),j,k,ipres)
             denom = max(small_pres, abs(q(i+2*dg(1),j,k,ipres) - q(i-2*dg(1),j,k,ipres)))
             zeta = abs(dp(i,j,k))/denom
             z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (q(i-1*dg(1),j,k,QU) - q(i+1*dg(1),j,k,QU) >= ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(q(i+1*dg(1),j,k,ipres), q(i-1*dg(1),j,k,ipres))
             if ((abs(dp(i,j,k))/tmp) > shktst) then
                chi(i,j,k) = tst
             else
                chi(i,j,k) = ZERO
             endif
          enddo
          do i = lo(1), hi(1)
             if(dp(i,j,k) > ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             flatn(i,j,k) = ONE - &
                  max(chi(i-ishft*dg(1),j,k)*z(i-ishft*dg(1),j,k), &
                      chi(i,j,k)*z(i,j,k))
          enddo
       enddo
    enddo

    ! y-direction flattening coef
    do k = lo(3), hi(3)
       do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
          !dir$ ivdep
          do i = lo(1), hi(1)
             dp(i,j,k) = q(i,j+1*dg(2),k,ipres) - q(i,j-1*dg(2),k,ipres)
             denom = max(small_pres, abs(q(i,j+2*dg(2),k,ipres) - q(i,j-2*dg(2),k,ipres)))
             zeta = abs(dp(i,j,k))/denom
             z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (q(i,j-1*dg(2),k,QV) - q(i,j+1*dg(2),k,QV) >= ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(q(i,j+1*dg(2),k,ipres), q(i,j-1*dg(2),k,ipres))
             if ((abs(dp(i,j,k))/tmp) > shktst) then
                chi(i,j,k) = tst
             else
                chi(i,j,k) = ZERO
             endif
          enddo
       end do
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if(dp(i,j,k) > ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             ftmp = ONE - &
                  max(chi(i,j-ishft*dg(2),k)*z(i,j-ishft*dg(2),k), &
                      chi(i,j,k)*z(i,j,k))
             flatn(i,j,k) = min( flatn(i,j,k), ftmp )
          enddo
       enddo
    enddo

    ! z-direction flattening coef
    do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
       do j = lo(2), hi(2)
          !dir$ ivdep
          do i = lo(1), hi(1)
             dp(i,j,k) = q(i,j,k+1*dg(3),ipres) - q(i,j,k-1*dg(3),ipres)
             denom = max(small_pres, abs(q(i,j,k+2*dg(3),ipres) - q(i,j,k-2*dg(3),ipres)))
             zeta = abs(dp(i,j,k))/denom
             z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (q(i,j,k-1*dg(3),QW) - q(i,j,k+1*dg(3),QW) >= ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(q(i,j,k+1*dg(3),ipres), q(i,j,k-1*dg(3),ipres))
             if ((abs(dp(i,j,k))/tmp) > shktst) then
                chi(i,j,k) = tst
             else
                chi(i,j,k) = ZERO
             endif
          enddo
       enddo
    enddo
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if(dp(i,j,k) > ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             ftmp = ONE - &
                  max(chi(i,j,k-ishft*dg(3))*z(i,j,k-ishft*dg(3)), &
                      chi(i,j,k)*z(i,j,k))
             flatn(i,j,k) = min( flatn(i,j,k), ftmp )
          enddo
       enddo
    enddo

    call bl_deallocate(dp )
    call bl_deallocate(z  )
    call bl_deallocate(chi)

  end subroutine uflatten

#ifdef RADIATION
  subroutine rad_flatten(lo, hi, q, flatn, q_lo, q_hi)

    use meth_params_module, only : QPRES, QU, QV, QW, flatten_pp_threshold, QPTOT, NQ

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)

    real(rt)        , intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: flatn(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))

    integer :: i, j, k

    real(rt)        , pointer :: flatg(:,:,:)

    call uflatten(lo, hi, q, flatn, q_lo, q_hi, QPTOT)

    call bl_allocate(flatg, q_lo(1), q_hi(1), q_lo(2), q_hi(2), q_lo(3), q_hi(3))
    call uflatten(lo, hi, q, flatg, q_lo, q_hi, QPRES)

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

    call bl_deallocate(flatg)

  end subroutine rad_flatten
#endif

#ifdef MHD
  subroutine mhd_flatten(lo, hi, q, flatn, q_lo, q_hi)

    use meth_params_module, only : QPRES, QU, QV, QW, NQ, QMAGX, QMAGY, QMAGZ

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)

    real(rt)        , intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: flatn(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))

    integer :: i, j, k, QPTOT

    real(rt)        :: q2(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ+1)
    real(rt)        , pointer :: flatg(:,:,:)


    QPTOT = NQ+1
    
    q2(:,:,:,1:NQ) = q(:,:,:,:)
    q2(:,:,:,QPTOT) = q(:,:,:,QPRES) + 0.5d0 * ( q(:,:,:,QMAGX)**2 + q(:,:,:,QMAGY)**2 + q(:,:,:,QMAGZ)**2 )

    call uflatten(lo, hi, q2, flatn, q_lo, q_hi, QPTOT)

    call bl_allocate(flatg, q_lo(1), q_hi(1), q_lo(2), q_hi(2), q_lo(3), q_hi(3))
    call uflatten(lo, hi, q2, flatg, q_lo, q_hi, QPRES)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             flatn(i,j,k) = flatn(i,j,k) * flatg(i,j,k)

!             if (flatten_pp_threshold > ZERO) then
!                if ( q(i-1,j,k,QU) + q(i,j-1,k,QV) + q(i,j,k-1,QW) > &
!                     q(i+1,j,k,QU) + q(i,j+1,k,QV) + q(i,j,k+1,QW) ) then

!                   if (q(i,j,k,QPRES) < flatten_pp_threshold * q(i,j,k,QPTOT)) then
!                      flatn(i,j,k) = ZERO
!                   end if

!                end if
!             endif

          end do
       end do
    end do

    call bl_deallocate(flatg)

  end subroutine mhd_flatten
#endif



end module flatten_module_mhd
