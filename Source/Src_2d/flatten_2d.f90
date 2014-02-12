module flatten_module

  implicit none

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine uflaten(lo,hi,p,u,v,flatn, &
                     q_l1,q_l2,q_h1,q_h2)

    use meth_params_module, only : iorder, small_pres
    
    implicit none

    integer lo(2),hi(2)
    integer q_l1,q_l2,q_h1,q_h2
    double precision p(q_l1:q_h1,q_l2:q_h2)
    double precision u(q_l1:q_h1,q_l2:q_h2)
    double precision v(q_l1:q_h1,q_l2:q_h2)
    double precision flatn(q_l1:q_h1,q_l2:q_h2)
    
    ! Local arrays
    double precision, allocatable :: dp(:), z(:), chi(:)
    
    integer i, j, idx, ishft
    double precision shktst, zcut1, zcut2, dzcut
    double precision denom, zeta, tst, tmp, ftmp
    integer nx,ny,nmax
    
    ! Knobs for detection of strong shock
    data shktst /0.33d0/
    data zcut1 /0.75d0/
    data zcut2 /0.85d0/
    
    nx = hi(1)-lo(1)+3
    ny = hi(2)-lo(2)+3
    nmax = max(nx,ny)
    allocate(dp(0:nmax-1),z(0:nmax-1),chi(0:nmax-1))
    
    dzcut = 1.d0/(zcut2-zcut1)

    if (iorder .eq. 3) then
       do j = lo(2),hi(2) 
          do i = lo(1),hi(1) 
             flatn(i,j) = 1.d0
          enddo
       enddo
       return
    endif
    
    ! x-direction flattening coef
    do j = lo(2),hi(2) 
       do i = lo(1)-1,hi(1)+1
          idx = i-lo(1)+1
          dp(idx) = p(i+1,j) - p(i-1,j)
          denom = max(small_pres,abs(p(i+2,j)-p(i-2,j)))
          zeta = abs(dp(idx))/denom
          z(idx) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
          if (u(i-1,j)-u(i+1,j) .ge. 0.d0) then
             tst = 1.d0
          else
             tst = 0.d0
          endif
          tmp = min(p(i+1,j),p(i-1,j))
          if ((abs(dp(idx))/tmp).gt.shktst) then
             chi(idx) = tst
          else
             chi(idx) = 0.d0
          endif
       enddo
       do i = lo(1),hi(1)
          idx = i-lo(1)+1
          if(dp(idx).gt.0.d0)then
             ishft = 1
          else
             ishft = -1
          endif
          flatn(i,j) = 1.d0 - &
               max(chi(idx-ishft)*z(idx-ishft),chi(idx)*z(idx))
       enddo
    enddo
    
    ! y-direction flattening coef
    do i = lo(1),hi(1)
       do j = lo(2)-1,hi(2)+1
          idx = j-lo(2)+1
          dp(idx) = p(i,j+1) - p(i,j-1)
          denom = max(small_pres,abs(p(i,j+2)-p(i,j-2)))
          zeta = abs(dp(idx))/denom
          z(idx) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
          if (v(i,j-1)-v(i,j+1) .ge. 0.d0) then
             tst = 1.d0
          else
             tst = 0.d0
          endif
          tmp = min(p(i,j+1),p(i,j-1))
          if ((abs(dp(idx))/tmp).gt.shktst) then
             chi(idx) = tst
          else
             chi(idx) = 0.d0
          endif
       enddo
       do j = lo(2),hi(2)
          idx = j-lo(2)+1
          if(dp(idx).gt.0.d0)then
             ishft = 1
          else
             ishft = -1
          endif
          ftmp = 1.d0 - &
               max(chi(idx-ishft)*z(idx-ishft),chi(idx)*z(idx))
          flatn(i,j) = min( flatn(i,j), ftmp )
       enddo
    enddo
    
    deallocate(dp,z,chi)
    
  end subroutine uflaten

end module flatten_module

