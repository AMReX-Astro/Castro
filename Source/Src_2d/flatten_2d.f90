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
    
    integer i, j, ishft
    double precision dzcut
    double precision denom, zeta, tst, tmp, ftmp
    integer nx,ny,nmax
    
    ! Knobs for detection of strong shock
    double precision, parameter :: shktst = 0.33d0
    double precision, parameter :: zcut1 = 0.75d0
    double precision, parameter :: zcut2 = 0.85d0
    
    nx = hi(1)-lo(1)+3
    ny = hi(2)-lo(2)+3
    nmax = max(nx,ny)

    allocate(dp(q_l1:q_h1),z(q_l1:q_h1),chi(q_l1:q_h1))
    
    dzcut = 1.d0/(zcut2-zcut1)
    
    ! x-direction flattening coef
    do j = lo(2),hi(2) 
       do i = lo(1)-1,hi(1)+1
          dp(i) = p(i+1,j) - p(i-1,j)
          denom = max(small_pres,abs(p(i+2,j)-p(i-2,j)))

          zeta = abs(dp(i))/denom
          z(i) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )

          if (u(i-1,j)-u(i+1,j) .ge. 0.d0) then
             tst = 1.d0
          else
             tst = 0.d0
          endif

          tmp = min(p(i+1,j),p(i-1,j))
          if ((abs(dp(i))/tmp).gt.shktst) then
             chi(i) = tst
          else
             chi(i) = 0.d0
          endif
       enddo

       do i = lo(1),hi(1)
          if(dp(i).gt.0.d0)then
             ishft = 1
          else
             ishft = -1
          endif
          flatn(i,j) = 1.d0 - &
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
          z(j) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )

          if (v(i,j-1)-v(i,j+1) .ge. 0.d0) then
             tst = 1.d0
          else
             tst = 0.d0
          endif

          tmp = min(p(i,j+1),p(i,j-1))
          if ((abs(dp(j))/tmp).gt.shktst) then
             chi(j) = tst
          else
             chi(j) = 0.d0
          endif
       enddo

       do j = lo(2),hi(2)
          if(dp(j).gt.0.d0)then
             ishft = 1
          else
             ishft = -1
          endif
          ftmp = 1.d0 - &
               max(chi(j-ishft)*z(j-ishft),chi(j)*z(j))

          flatn(i,j) = min( flatn(i,j), ftmp )
       enddo
    enddo
    
    deallocate(dp,z,chi)
    
  end subroutine uflaten

end module flatten_module

