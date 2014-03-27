module flatten_module

  implicit none

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine uflaten(lo,hi,p,u,v,flatn, &
                     q_l1,q_l2,q_h1,q_h2)

    use meth_params_module, only : small_pres
    use bl_constants_module
    
    implicit none

    integer, intent(in) :: lo(2),hi(2)
    integer, intent(in) ::  q_l1,q_l2,q_h1,q_h2
    double precision, intent(in)    ::     p(q_l1:q_h1,q_l2:q_h2)
    double precision, intent(in)    ::     u(q_l1:q_h1,q_l2:q_h2)
    double precision, intent(in)    ::     v(q_l1:q_h1,q_l2:q_h2)
    double precision, intent(inout) :: flatn(q_l1:q_h1,q_l2:q_h2)
    
    double precision, allocatable :: dp(:), z(:), chi(:)

    ! Local arrays    
    integer i, j, ishft
    double precision dzcut
    double precision denom, zeta, tst, tmp, ftmp
    
    ! Knobs for detection of strong shock
    double precision, parameter :: shktst = 0.33d0
    double precision, parameter :: zcut1 = 0.75d0
    double precision, parameter :: zcut2 = 0.85d0

    dzcut = ONE/(zcut2-zcut1)
    
    ! x-direction flattening coef
    allocate(dp(q_l1:q_h1),z(q_l1:q_h1),chi(q_l1:q_h1))
    
    do j = lo(2),hi(2) 
       do i = lo(1)-1,hi(1)+1
          dp(i) = p(i+1,j) - p(i-1,j)
          denom = max(small_pres,abs(p(i+2,j)-p(i-2,j)))

          zeta = abs(dp(i))/denom
          z(i) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )

          if (u(i-1,j)-u(i+1,j) .ge. ZERO) then
             tst = ONE
          else
             tst = ZERO
          endif

          tmp = min(p(i+1,j),p(i-1,j))
          if ((abs(dp(i))/tmp).gt.shktst) then
             chi(i) = tst
          else
             chi(i) = ZERO
          endif
       enddo

       do i = lo(1),hi(1)
          if(dp(i).gt.ZERO)then
             ishft = 1
          else
             ishft = -1
          endif
          flatn(i,j) = ONE - &
               max(chi(i-ishft)*z(i-ishft),chi(i)*z(i))
       enddo
    enddo

    deallocate(dp,z,chi)

    ! y-direction flattening coef
    allocate(dp(q_l2:q_h2),z(q_l2:q_h2),chi(q_l2:q_h2))
    
    do i = lo(1),hi(1)
       do j = lo(2)-1,hi(2)+1
          dp(j) = p(i,j+1) - p(i,j-1)
          denom = max(small_pres,abs(p(i,j+2)-p(i,j-2)))

          zeta = abs(dp(j))/denom
          z(j) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )

          if (v(i,j-1)-v(i,j+1) .ge. ZERO) then
             tst = ONE
          else
             tst = ZERO
          endif

          tmp = min(p(i,j+1),p(i,j-1))
          if ((abs(dp(j))/tmp).gt.shktst) then
             chi(j) = tst
          else
             chi(j) = ZERO
          endif
       enddo

       do j = lo(2),hi(2)
          if(dp(j).gt.ZERO)then
             ishft = 1
          else
             ishft = -1
          endif
          ftmp = ONE - &
               max(chi(j-ishft)*z(j-ishft),chi(j)*z(j))

          ! merge the x- and y-directions into a single flattening
          ! coeff for the zone
          flatn(i,j) = min( flatn(i,j), ftmp )
       enddo
    enddo
    
    deallocate(dp,z,chi)
    
  end subroutine uflaten

end module flatten_module

