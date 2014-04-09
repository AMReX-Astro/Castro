module flatten_module

  implicit none

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine uflaten(lo,hi,p,u,flatn,q_l1,q_h1)

    use meth_params_module, only : small_pres
    use bl_constants_module

    implicit none

    integer lo(1),hi(1)
    integer q_l1,q_h1
    double precision p(q_l1:q_h1)
    double precision u(q_l1:q_h1)
    double precision flatn(q_l1:q_h1)

    ! Local arrays
    double precision, allocatable :: dp(:), z(:), chi(:)

    integer i, ishft
    double precision dzcut
    double precision denom, zeta, tst, tmp
    
    ! Knobs for detection of strong shock
    double precision, parameter :: shktst = 0.33d0
    double precision, parameter :: zcut1 = 0.75d0
    double precision, parameter :: zcut2 = 0.85d0
    
    allocate(dp(lo(1)-1:hi(1)+1),z(lo(1)-1:hi(1)+1),chi(lo(1)-1:hi(1)+1))

    dzcut = ONE/(zcut2-zcut1)

    
    ! x-direction flattening coef
    do i = lo(1)-1,hi(1)+1
       denom = max(small_pres,abs(p(i+2)-p(i-2)))
       dp(i) = p(i+1) - p(i-1)
       zeta = abs(dp(i))/denom
       z(i) = min( ONE, max( ZERO, dzcut*(zeta - zcut1) ) )
       if (u(i-1)-u(i+1) .ge. ZERO) then
          tst = ONE
       else
          tst = ZERO
       endif
       tmp = min(p(i+1),p(i-1))
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
       flatn(i) = ONE - &
            max(chi(i-ishft)*z(i-ishft),chi(i)*z(i))
    enddo
    
    deallocate(dp,z,chi)
    
  end subroutine uflaten

end module flatten_module
