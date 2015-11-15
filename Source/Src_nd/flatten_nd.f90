module flatten_module

  implicit none

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine uflaten(lo,hi,p,u,v,w,flatn,q_lo,q_hi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : iorder, small_pres
    use prob_params_module, only : dg
    use bl_constants_module

    implicit none

    integer :: lo(3), hi(3)
    integer :: q_lo(3), q_hi(3)
    
    double precision :: p(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision :: u(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision :: v(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision :: w(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision :: flatn(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))

    integer :: i, j, k, ishft

    double precision :: denom, zeta, tst, tmp, ftmp

    ! Local arrays
    double precision, pointer :: dp(:,:,:), z(:,:,:), chi(:,:,:)
    
    ! Knobs for detection of strong shock
    double precision, parameter :: shktst = 0.33d0, zcut1 = 0.75d0, zcut2 = 0.85d0, dzcut = ONE/(zcut2-zcut1)

    call bl_allocate(dp ,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(z  ,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
    call bl_allocate(chi,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)

    ! x-direction flattening coef
    do k = lo(3),hi(3)
       do j = lo(2),hi(2) 
          !dir$ ivdep
          do i = lo(1)-1*dg(1),hi(1)+1*dg(1)
             dp(i,j,k) = p(i+1*dg(1),j,k) - p(i-1*dg(1),j,k)
             denom = max(small_pres,abs(p(i+2*dg(2),j,k)-p(i-2*dg(2),j,k)))
             zeta = abs(dp(i,j,k))/denom
             z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (u(i-1*dg(1),j,k)-u(i+1*dg(1),j,k) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(p(i+1*dg(1),j,k),p(i-1*dg(1),j,k))
             if ((abs(dp(i,j,k))/tmp).gt.shktst) then
                chi(i,j,k) = tst
             else
                chi(i,j,k) = ZERO
             endif
          enddo
          do i = lo(1),hi(1)
             if(dp(i,j,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             flatn(i,j,k) = ONE - &
                  max(chi(i-ishft*dg(1),j,k)*z(i-ishft*dg(1),j,k),chi(i,j,k)*z(i,j,k))
          enddo
       enddo
    enddo

    ! y-direction flattening coef
    do k = lo(3),hi(3)
       do j = lo(2)-1*dg(2),hi(2)+1*dg(2)
          !dir$ ivdep
          do i = lo(1),hi(1)
             dp(i,j,k) = p(i,j+1*dg(2),k) - p(i,j-1*dg(2),k)
             denom = max(small_pres,abs(p(i,j+2*dg(2),k)-p(i,j-2*dg(2),k)))
             zeta = abs(dp(i,j,k))/denom
             z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (v(i,j-1*dg(2),k)-v(i,j+1*dg(2),k) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(p(i,j+1*dg(2),k),p(i,j-1*dg(2),k))
             if ((abs(dp(i,j,k))/tmp).gt.shktst) then
                chi(i,j,k) = tst
             else
                chi(i,j,k) = ZERO
             endif
          enddo
       end do
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if(dp(i,j,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             ftmp = ONE - &
                  max(chi(i,j-ishft*dg(2),k)*z(i,j-ishft*dg(2),k),chi(i,j,k)*z(i,j,k))
             flatn(i,j,k) = min( flatn(i,j,k), ftmp )
          enddo
       enddo
    enddo

    ! z-direction flattening coef
    do k = lo(3)-1*dg(3),hi(3)+1*dg(3)
       do j = lo(2),hi(2) 
          !dir$ ivdep
          do i = lo(1),hi(1)
             dp(i,j,k) = p(i,j,k+1*dg(3)) - p(i,j,k-1*dg(3))
             denom = max(small_pres,abs(p(i,j,k+2*dg(3))-p(i,j,k-2*dg(3))))
             zeta = abs(dp(i,j,k))/denom
             z(i,j,k) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
             if (w(i,j,k-1*dg(3))-w(i,j,k+1*dg(3)) .ge. ZERO) then
                tst = ONE
             else
                tst = ZERO
             endif
             tmp = min(p(i,j,k+1*dg(3)),p(i,j,k-1*dg(3)))
             if ((abs(dp(i,j,k))/tmp).gt.shktst) then
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
             if(dp(i,j,k).gt.ZERO)then
                ishft = 1
             else
                ishft = -1
             endif
             ftmp = ONE - &
                  max(chi(i,j,k-ishft*dg(3))*z(i,j,k-ishft*dg(3)),chi(i,j,k)*z(i,j,k))
             flatn(i,j,k) = min( flatn(i,j,k), ftmp )
          enddo
       enddo
    enddo
    
    call bl_deallocate(dp )
    call bl_deallocate(z  )
    call bl_deallocate(chi)

  end subroutine uflaten

end module flatten_module
