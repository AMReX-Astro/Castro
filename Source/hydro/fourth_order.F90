module fourth_order

  use amrex_constants_module, only : ZERO, TWO, ONE
  use prob_params_module, only : dg

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), parameter :: TWENTYFOURTH = ONE/24.0_rt

contains

  subroutine states(idir, &
                    a, a_lo, a_hi, ncomp, n, &
                    flatn, f_lo, f_hi, &
                    al, ar, ai_lo, ai_hi, &
                    lo, hi, &
                    domlo, domhi)

    use meth_params_module, only : limit_fourth_order
    use prob_params_module, only : Interior, physbc_lo, physbc_hi

    implicit none

    integer, intent(in) :: idir, n, ncomp
    integer, intent(in) :: a_lo(3), a_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: ai_lo(3), ai_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in) :: a(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3), ncomp)
    real(rt), intent(in) :: flatn( f_lo(1): f_hi(1), f_lo(2): f_hi(2), f_lo(3): f_hi(3))

    real(rt), intent(inout) :: al(ai_lo(1):ai_hi(1), ai_lo(2):ai_hi(2), ai_lo(3):ai_hi(3), ncomp)
    real(rt), intent(inout) :: ar(ai_lo(1):ai_hi(1), ai_lo(2):ai_hi(2), ai_lo(3):ai_hi(3), ncomp)
    integer, intent(in) :: domlo(3), domhi(3)

    ! local variables
    real(rt) :: a_int(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt) :: dafm(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt) :: dafp(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt) :: d2af(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt) :: d2ac(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt) :: d3a(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))

    real(rt), parameter :: C2 = 1.25_rt
    real(rt), parameter :: C3 = 0.1_rt

    integer :: i, j, k
    double precision :: rho, s

    double precision :: d2a_lim, d3a_min, d3a_max

    ! our convention here is that:
    !     al(i,j,k)   will be al_{i-1/2,j,k),
    !     al(i+1,j,k) will be al_{i+1/2,j,k)

    ! we need interface values on all faces of the domain
    if (idir == 1) then

       ! this loop is over interfaces
       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2)-dg(2), hi(2)+dg(2)
             do i = lo(1)-1, hi(1)+2

                ! interpolate to the edges -- this is a_{i-1/2}
                ! note for non-periodic physical boundaries, we use a special stencil

                if (i == domlo(1)+1 .and. physbc_lo(1) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the left physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i-1,j,k,n) + 13.0_rt*a(i,j,k,n) - &
                                                    5.0_rt*a(i+1,j,k,n) + a(i+2,j,k,n))

                else if (i == domlo(1) .and. physbc_lo(1) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! left physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k,n) - 23.0_rt*a(i+1,j,k,n) + &
                                                    13.0_rt*a(i+2,j,k,n) - 3.0_rt*a(i+3,j,k,n))

                else if (i == domhi(1) .and. physbc_hi(1) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the right physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k,n) + 13.0_rt*a(i-1,j,k,n) - &
                                                    5.0_rt*a(i-2,j,k,n) + a(i-3,j,k,n))

                else if (i == domhi(1)+1 .and. physbc_hi(1) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! right physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i-1,j,k,n) - 23.0_rt*a(i-2,j,k,n) + &
                                                    13.0_rt*a(i-3,j,k,n) - 3.0_rt*a(i-4,j,k,n))

                else
                   ! regular stencil
                   a_int(i,j,k) = (7.0_rt/12.0_rt)*(a(i-1,j,k,n) + a(i,j,k,n)) - &
                        (1.0_rt/12.0_rt)*(a(i-2,j,k,n) + a(i+1,j,k,n))
                end if

                al(i,j,k,n) = a_int(i,j,k)
                ar(i,j,k,n) = a_int(i,j,k)

             enddo
          enddo
       enddo

       ! the limiting loops are now over zones
       if (limit_fourth_order == 1) then
          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-dg(2), hi(2)+dg(2)
                do i = lo(1)-1, hi(1)+1
                   ! these live on cell-centers
                   dafm(i,j,k) = a(i,j,k,n) - a_int(i,j,k)
                   dafp(i,j,k) = a_int(i+1,j,k) - a(i,j,k,n)

                   ! these live on cell-centers
                   d2af(i,j,k) = 6.0_rt*(a_int(i,j,k) - 2.0_rt*a(i,j,k,n) + a_int(i+1,j,k))
                enddo
             enddo
          enddo

          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-dg(2), hi(2)+dg(2)
                do i = lo(1)-3, hi(1)+3
                   d2ac(i,j,k) = a(i-1,j,k,n) - 2.0_rt*a(i,j,k,n) + a(i+1,j,k,n)
                enddo
             enddo
          enddo

          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-dg(2), hi(2)+dg(2)
                do i = lo(1)-2, hi(1)+3
                   ! this lives on the interface
                   d3a(i,j,k) = d2ac(i,j,k) - d2ac(i-1,j,k)
                enddo
             enddo
          enddo

          ! this is a loop over cell centers, affecting
          ! i-1/2,R and i+1/2,L
          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-dg(2), hi(2)+dg(2)
                do i = lo(1)-1, hi(1)+1

                   ! limit? MC Eq. 24 and 25
                   if (dafm(i,j,k) * dafp(i,j,k) <= 0.0_rt .or. &
                        (a(i,j,k,n) - a(i-2,j,k,n))*(a(i+2,j,k,n) - a(i,j,k,n)) <= 0.0_rt) then

                      ! we are at an extrema

                      s = sign(1.0_rt, d2ac(i,j,k))
                      if ( s == sign(1.0_rt, d2ac(i-1,j,k)) .and. &
                           s == sign(1.0_rt, d2ac(i+1,j,k)) .and. &
                           s == sign(1.0_rt, d2af(i,j,k))) then
                         ! MC Eq. 26
                         d2a_lim = s*min(abs(d2af(i,j,k)), C2*abs(d2ac(i-1,j,k)), &
                              C2*abs(d2ac(i,j,k)), C2*abs(d2ac(i+1,j,k)))
                      else
                         d2a_lim = 0.0_rt
                      endif

                      if (abs(d2af(i,j,k)) <= 1.e-12*max(abs(a(i-2,j,k,n)), abs(a(i-1,j,k,n)), &
                           abs(a(i,j,k,n)), abs(a(i+1,j,k,n)), abs(a(i+2,j,k,n)))) then
                         rho = 0.0_rt
                      else
                         ! MC Eq. 27
                         rho = d2a_lim/d2af(i,j,k)
                      endif

                      if (rho < 1.0_rt - 1.e-12_rt) then
                         ! we may need to limit -- these quantities are at cell-centers
                         d3a_min = min(d3a(i-1,j,k), d3a(i,j,k), d3a(i+1,j,k), d3a(i+2,j,k))
                         d3a_max = max(d3a(i-1,j,k), d3a(i,j,k), d3a(i+1,j,k), d3a(i+2,j,k))

                         if (C3*max(abs(d3a_min), abs(d3a_max)) <= (d3a_max - d3a_min)) then
                            ! limit
                            if (dafm(i,j,k)*dafp(i,j,k) < 0.0_rt) then
                               ! Eqs. 29, 30
                               ar(i,j,k,n) = a(i,j,k,n) - rho*dafm(i,j,k)  ! note: typo in Eq 29
                               al(i+1,j,k,n) = a(i,j,k,n) + rho*dafp(i,j,k)
                            else if (abs(dafm(i,j,k)) >= 2.0_rt*abs(dafp(i,j,k))) then
                               ! Eq. 31
                               ar(i,j,k,n) = a(i,j,k,n) - 2.0_rt*(1.0_rt - rho)*dafp(i,j,k) - rho*dafm(i,j,k)
                            else if (abs(dafp(i,j,k)) >= 2.0_rt*abs(dafm(i,j,k))) then
                               ! Eq. 32
                               al(i+1,j,k,n) = a(i,j,k,n) + 2.0_rt*(1.0_rt - rho)*dafm(i,j,k) + rho*dafp(i,j,k)
                            endif

                         endif
                      endif

                   else
                      ! if Eqs. 24 or 25 didn't hold we still may need to limit
                      if (abs(dafm(i,j,k)) >= 2.0_rt*abs(dafp(i,j,k))) then
                         ar(i,j,k,n) = a(i,j,k,n) - 2.0_rt*dafp(i,j,k)
                      endif
                      if (abs(dafp(i,j,k)) >= 2.0_rt*abs(dafm(i,j,k))) then
                         al(i+1,j,k,n) = a(i,j,k,n) + 2.0_rt*dafm(i,j,k)
                      endif
                   endif

                   ! apply flattening
                   al(i+1,j,k,n) = flatn(i,j,k)*al(i+1,j,k,n) + (ONE - flatn(i,j,k))*a(i,j,k,n)
                   ar(i,j,k,n) = flatn(i,j,k)*ar(i,j,k,n) + (ONE - flatn(i,j,k))*a(i,j,k,n)
                enddo

             enddo
          enddo

       endif

    else if (idir == 2) then

       ! this loop is over interfaces
       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2)-1, hi(2)+2
             do i = lo(1)-1, hi(1)+1

                ! interpolate to the edges
                if (j == domlo(2)+1 .and. physbc_lo(2) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the left physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j-1,k,n) + 13.0_rt*a(i,j,k,n) - &
                                                    5.0_rt*a(i,j+1,k,n) + a(i,j+2,k,n))

                else if (j == domlo(2) .and. physbc_lo(2) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! left physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k,n) - 23.0_rt*a(i,j+1,k,n) + &
                                                    13.0_rt*a(i,j+2,k,n) - 3.0_rt*a(i,j+3,k,n))

                else if (j == domhi(2) .and. physbc_hi(2) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the right physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k,n) + 13.0_rt*a(i,j-1,k,n) - &
                                                    5.0_rt*a(i,j-2,k,n) + a(i,j-3,k,n))

                else if (j == domhi(2)+1 .and. physbc_hi(2) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! right physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j-1,k,n) - 23.0_rt*a(i,j-2,k,n) + &
                                                    13.0_rt*a(i,j-3,k,n) - 3.0_rt*a(i,j-4,k,n))

                else
                   ! regular stencil
                   a_int(i,j,k) = (7.0_rt/12.0_rt)*(a(i,j-1,k,n) + a(i,j,k,n)) - &
                        (1.0_rt/12.0_rt)*(a(i,j-2,k,n) + a(i,j+1,k,n))
                end if

                al(i,j,k,n) = a_int(i,j,k)
                ar(i,j,k,n) = a_int(i,j,k)

             enddo
          enddo
       enddo

       ! the limiting loops are now over zones
       if (limit_fourth_order == 1) then
          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1
                   ! these live on cell-centers
                   dafm(i,j,k) = a(i,j,k,n) - a_int(i,j,k)
                   dafp(i,j,k) = a_int(i,j+1,k) - a(i,j,k,n)

                   ! these live on cell-centers
                   d2af(i,j,k) = 6.0_rt*(a_int(i,j,k) - 2.0_rt*a(i,j,k,n) + a_int(i,j+1,k))
                enddo
             enddo
          enddo

          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-3, hi(2)+3
                do i = lo(1)-1, hi(1)+1
                   d2ac(i,j,k) = a(i,j-1,k,n) - 2.0_rt*a(i,j,k,n) + a(i,j+1,k,n)
                enddo
             enddo
          enddo

          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-2, hi(2)+3
                do i = lo(1)-1, hi(1)+1
                   ! this lives on the interface
                   d3a(i,j,k) = d2ac(i,j,k) - d2ac(i,j-1,k)
                enddo
             enddo
          enddo

          ! this is a loop over cell centers, affecting
          ! j-1/2,R and j+1/2,L
          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1

                   ! limit? MC Eq. 24 and 25
                   if (dafm(i,j,k) * dafp(i,j,k) <= 0.0_rt .or. &
                        (a(i,j,k,n) - a(i,j-2,k,n))*(a(i,j+2,k,n) - a(i,j,k,n)) <= 0.0_rt) then

                      ! we are at an extrema

                      s = sign(1.0_rt, d2ac(i,j,k))
                      if ( s == sign(1.0_rt, d2ac(i,j-1,k)) .and. &
                           s == sign(1.0_rt, d2ac(i,j+1,k)) .and. &
                           s == sign(1.0_rt, d2af(i,j,k))) then
                         ! MC Eq. 26
                         d2a_lim = s*min(abs(d2af(i,j,k)), C2*abs(d2ac(i,j-1,k)), &
                              C2*abs(d2ac(i,j,k)), C2*abs(d2ac(i,j+1,k)))
                      else
                         d2a_lim = 0.0_rt
                      endif

                      if (abs(d2af(i,j,k)) <= 1.e-12*max(abs(a(i,j-2,k,n)), abs(a(i,j-1,k,n)), &
                           abs(a(i,j,k,n)), abs(a(i,j+1,k,n)), abs(a(i,j+2,k,n)))) then
                         rho = 0.0_rt
                      else
                         ! MC Eq. 27
                         rho = d2a_lim/d2af(i,j,k)
                      endif

                      if (rho < 1.0_rt - 1.e-12_rt) then
                         ! we may need to limit -- these quantities are at cell-centers
                         d3a_min = min(d3a(i,j-1,k), d3a(i,j,k), d3a(i,j+1,k), d3a(i,j+2,k))
                         d3a_max = max(d3a(i,j-1,k), d3a(i,j,k), d3a(i,j+1,k), d3a(i,j+2,k))

                         if (C3*max(abs(d3a_min), abs(d3a_max)) <= (d3a_max - d3a_min)) then
                            ! limit
                            if (dafm(i,j,k)*dafp(i,j,k) < 0.0_rt) then
                               ! Eqs. 29, 30
                               ar(i,j,k,n) = a(i,j,k,n) - rho*dafm(i,j,k)  ! note: typo in Eq 29
                               al(i,j+1,k,n) = a(i,j,k,n) + rho*dafp(i,j,k)
                            else if (abs(dafm(i,j,k)) >= 2.0_rt*abs(dafp(i,j,k))) then
                               ! Eq. 31
                               ar(i,j,k,n) = a(i,j,k,n) - 2.0_rt*(1.0_rt - rho)*dafp(i,j,k) - rho*dafm(i,j,k)
                            else if (abs(dafp(i,j,k)) >= 2.0*abs(dafm(i,j,k))) then
                               ! Eq. 32
                               al(i,j+1,k,n) = a(i,j,k,n) + 2.0_rt*(1.0_rt - rho)*dafm(i,j,k) + rho*dafp(i,j,k)
                            endif

                         endif
                      endif

                   else
                      ! if Eqs. 24 or 25 didn't hold we still may need to limit
                      if (abs(dafm(i,j,k)) >= 2.0_rt*abs(dafp(i,j,k))) then
                         ar(i,j,k,n) = a(i,j,k,n) - 2.0_rt*dafp(i,j,k)
                      endif
                      if (abs(dafp(i,j,k)) >= 2.0_rt*abs(dafm(i,j,k))) then
                         al(i,j+1,k,n) = a(i,j,k,n) + 2.0_rt*dafm(i,j,k)
                      endif
                   endif

                   ! apply flattening
                   al(i,j+1,k,n) = flatn(i,j,k)*al(i,j+1,k,n) + (ONE - flatn(i,j,k))*a(i,j,k,n)
                   ar(i,j,k,n) = flatn(i,j,k)*ar(i,j,k,n) + (ONE - flatn(i,j,k))*a(i,j,k,n)

                enddo
             enddo
          enddo
       endif

    else if (idir == 3) then

       ! this loop is over interfaces
       do k = lo(3)-1, hi(3)+2
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

                ! interpolate to the edges
                if (k == domlo(3)+1 .and. physbc_lo(3) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the left physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k-1,n) + 13.0_rt*a(i,j,k,n) - &
                                                    5.0_rt*a(i,j,k+1,n) + a(i,j,k+2,n))

                else if (k == domlo(3) .and. physbc_lo(3) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! left physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k,n) - 23.0_rt*a(i,j,k+1,n) + &
                                                    13.0_rt*a(i,j,k+2,n) - 3.0_rt*a(i,j,k+3,n))

                else if (k == domhi(3) .and. physbc_hi(3) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the right physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k,n) + 13.0_rt*a(i,j,k-1,n) - &
                                                    5.0_rt*a(i,j,k-2,n) + a(i,j,k-3,n))

                else if (k == domhi(3)+1 .and. physbc_hi(3) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! right physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k-1,n) - 23.0_rt*a(i,j,k-2,n) + &
                                                    13.0_rt*a(i,j,k-3,n) - 3.0_rt*a(i,j,k-4,n))

                else
                   ! regular stencil
                   a_int(i,j,k) = (7.0_rt/12.0_rt)*(a(i,j,k-1,n) + a(i,j,k,n)) - &
                        (1.0_rt/12.0_rt)*(a(i,j,k-2,n) + a(i,j,k+1,n))
                end if

                al(i,j,k,n) = a_int(i,j,k)
                ar(i,j,k,n) = a_int(i,j,k)

             enddo
          enddo
       enddo

       ! the limiting loops are now over zones
       if (limit_fourth_order == 1) then

          do k = lo(3)-1, hi(3)+1
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1
                   ! these live on cell-centers
                   dafm(i,j,k) = a(i,j,k,n) - a_int(i,j,k)
                   dafp(i,j,k) = a_int(i,j,k+1) - a(i,j,k,n)

                   ! these live on cell-centers
                   d2af(i,j,k) = 6.0_rt*(a_int(i,j,k) - 2.0_rt*a(i,j,k,n) + a_int(i,j,k+1))
                enddo
             enddo
          enddo

          do k = lo(3)-3, hi(3)+3
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1
                   d2ac(i,j,k) = a(i,j,k-1,n) - 2.0_rt*a(i,j,k,n) + a(i,j,k+1,n)
                enddo
             enddo
          enddo

          do k = lo(3)-2, hi(3)+3
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1
                   ! this lives on the interface
                   d3a(i,j,k) = d2ac(i,j,k) - d2ac(i,j,k-1)
                enddo
             enddo
          enddo

          ! this is a loop over cell centers, affecting
          ! j-1/2,R and j+1/2,L
          do k = lo(3)-1, hi(3)+1
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1

                   ! limit? MC Eq. 24 and 25
                   if (dafm(i,j,k) * dafp(i,j,k) <= 0.0_rt .or. &
                        (a(i,j,k,n) - a(i,j,k-2,n))*(a(i,j,k+2,n) - a(i,j,k,n)) <= 0.0_rt) then

                      ! we are at an extrema

                      s = sign(1.0_rt, d2ac(i,j,k))
                      if ( s == sign(1.0_rt, d2ac(i,j,k-1)) .and. &
                           s == sign(1.0_rt, d2ac(i,j,k+1)) .and. &
                           s == sign(1.0_rt, d2af(i,j,k))) then
                         ! MC Eq. 26
                         d2a_lim = s*min(abs(d2af(i,j,k)), C2*abs(d2ac(i,j,k-1)), &
                              C2*abs(d2ac(i,j,k)), C2*abs(d2ac(i,j,k+1)))
                      else
                         d2a_lim = 0.0_rt
                      endif

                      if (abs(d2af(i,j,k)) <= 1.e-12*max(abs(a(i,j,k-2,n)), abs(a(i,j,k-1,n)), &
                           abs(a(i,j,k,n)), abs(a(i,j,k+1,n)), abs(a(i,j,k+2,n)))) then
                         rho = 0.0_rt
                      else
                         ! MC Eq. 27
                         rho = d2a_lim/d2af(i,j,k)
                      endif

                      if (rho < 1.0_rt - 1.e-12_rt) then
                         ! we may need to limit -- these quantities are at cell-centers
                         d3a_min = min(d3a(i,j,k-1), d3a(i,j,k), d3a(i,j,k+1), d3a(i,j,k+2))
                         d3a_max = max(d3a(i,j,k-1), d3a(i,j,k), d3a(i,j,k+1), d3a(i,j,k+2))

                         if (C3*max(abs(d3a_min), abs(d3a_max)) <= (d3a_max - d3a_min)) then
                            ! limit
                            if (dafm(i,j,k)*dafp(i,j,k) < 0.0_rt) then
                               ! Eqs. 29, 30
                               ar(i,j,k,n) = a(i,j,k,n) - rho*dafm(i,j,k)  ! note: typo in Eq 29
                               al(i,j,k+1,n) = a(i,j,k,n) + rho*dafp(i,j,k)
                            else if (abs(dafm(i,j,k)) >= 2.0_rt*abs(dafp(i,j,k))) then
                               ! Eq. 31
                               ar(i,j,k,n) = a(i,j,k,n) - 2.0_rt*(1.0_rt - rho)*dafp(i,j,k) - rho*dafm(i,j,k)
                            else if (abs(dafp(i,j,k)) >= 2.0*abs(dafm(i,j,k))) then
                               ! Eq. 32
                               al(i,j,k+1,n) = a(i,j,k,n) + 2.0_rt*(1.0_rt - rho)*dafm(i,j,k) + rho*dafp(i,j,k)
                            endif

                         endif
                      endif

                   else
                      ! if Eqs. 24 or 25 didn't hold we still may need to limit
                      if (abs(dafm(i,j,k)) >= 2.0_rt*abs(dafp(i,j,k))) then
                         ar(i,j,k,n) = a(i,j,k,n) - 2.0_rt*dafp(i,j,k)
                      endif
                      if (abs(dafp(i,j,k)) >= 2.0_rt*abs(dafm(i,j,k))) then
                         al(i,j,k+1,n) = a(i,j,k,n) + 2.0_rt*dafm(i,j,k)
                      endif
                   endif

                   ! apply flattening
                   al(i,j,k+1,n) = flatn(i,j,k)*al(i,j,k+1,n) + (ONE - flatn(i,j,k))*a(i,j,k,n)
                   ar(i,j,k,n) = flatn(i,j,k)*ar(i,j,k,n) + (ONE - flatn(i,j,k))*a(i,j,k,n)

                enddo
             enddo
          enddo
       endif
    endif

  end subroutine states


  ! Note: pretty much all of these routines below assume that dx(1) = dx(2) = dx(3)

  subroutine ca_make_cell_center(lo, hi, &
                                 U, U_lo, U_hi, nc, &
                                 U_cc, U_cc_lo, U_cc_hi, nc_cc, &
                                 domlo, domhi) &
                                 bind(C, name="ca_make_cell_center")
    ! Take a cell-average state U and a convert it to a cell-center
    ! state U_cc via U_cc = U - 1/24 L U

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3)
    integer, intent(in) :: U_cc_lo(3), U_cc_hi(3)
    integer, intent(in) :: nc, nc_cc
    real(rt), intent(in) :: U(U_lo(1):U_hi(1), U_lo(2):U_hi(2), U_lo(3):U_hi(3), nc)
    real(rt), intent(inout) :: U_cc(U_cc_lo(1):U_cc_hi(1), U_cc_lo(2):U_cc_hi(2), U_cc_lo(3):U_cc_hi(3), nc_cc)
    integer, intent(in) :: domlo(3), domhi(3)

    integer :: i, j, k, n
    real(rt) :: lapx, lapy, lapz

    if (U_lo(1) > lo(1)-1 .or. U_hi(1) < hi(1)+1 .or. &
        (AMREX_SPACEDIM >= 2 .and. (U_lo(2) > lo(2)-1 .or. U_hi(2) < hi(2)+1)) .or. &
        (AMREX_SPACEDIM == 3 .and. (U_lo(3) > lo(3)-1 .or. U_hi(3) < hi(3)+1))) then
       call bl_error("insufficient ghostcells in ca_make_cell_center")
    endif

    lapx = ZERO
    lapy = ZERO
    lapz = ZERO

    do n = 1, nc

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                lapx = U(i+1,j,k,n) - TWO*U(i,j,k,n) + U(i-1,j,k,n)
#if AMREX_SPACEDIM >= 2
                lapy = U(i,j+1,k,n) - TWO*U(i,j,k,n) + U(i,j-1,k,n)
#endif
#if AMREX_SPACEDIM == 3
                lapz = U(i,j,k+1,n) - TWO*U(i,j,k,n) + U(i,j,k-1,n)
#endif

                U_cc(i,j,k,n) = U(i,j,k,n) - TWENTYFOURTH * (lapx + lapy + lapz)

             end do
          end do
       end do
    end do

  end subroutine ca_make_cell_center

  subroutine ca_make_cell_center_in_place(lo, hi, &
                                          U, U_lo, U_hi, nc, &
                                          domlo, domhi) &
                                          bind(C, name="ca_make_cell_center_in_place")
    ! Take a cell-average state U and make it cell-centered in place
    ! via U <- U - 1/24 L U.  Note that this operation is not tile
    ! safe.

    use amrex_mempool_module, only : bl_allocate, bl_deallocate

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3)
    integer, intent(in) :: nc
    real(rt), intent(inout) :: U(U_lo(1):U_hi(1), U_lo(2):U_hi(2), U_lo(3):U_hi(3), nc)
    integer, intent(in) :: domlo(3), domhi(3)

    integer :: i, j, k, n

    real(rt) :: lapx, lapy, lapz
    real(rt), pointer :: lap(:,:,:)

    call bl_allocate(lap, lo, hi)

    lapx = ZERO
    lapy = ZERO
    lapz = ZERO

    do n = 1, nc

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                lapx = U(i+1,j,k,n) - TWO*U(i,j,k,n) + U(i-1,j,k,n)
#if AMREX_SPACEDIM >= 2
                lapy = U(i,j+1,k,n) - TWO*U(i,j,k,n) + U(i,j-1,k,n)
#endif
#if AMREX_SPACEDIM == 3
                lapz = U(i,j,k+1,n) - TWO*U(i,j,k,n) + U(i,j,k-1,n)
#endif

                lap(i,j,k) = lapx + lapy + lapz
             enddo
          enddo
       enddo

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                U(i,j,k,n) = U(i,j,k,n) - TWENTYFOURTH * lap(i,j,k)
             enddo
          enddo
       enddo
    enddo

    call bl_deallocate(lap)

  end subroutine ca_make_cell_center_in_place

  subroutine ca_compute_lap_term(lo, hi, &
                                 U, U_lo, U_hi, nc, &
                                 lap, lap_lo, lap_hi, ncomp, &
                                 domlo, domhi) &
                                 bind(C, name="ca_compute_lap_term")
    ! Computes the h**2/24 L U term that is used in correcting
    ! cell-center to averages (and back)

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3)
    integer, intent(in) :: lap_lo(3), lap_hi(3)
    integer, intent(in) :: nc
    real(rt), intent(in) :: U(U_lo(1):U_hi(1), U_lo(2):U_hi(2), U_lo(3):U_hi(3), nc)
    real(rt), intent(inout) :: lap(lap_lo(1):lap_hi(1), lap_lo(2):lap_hi(2), lap_lo(3):lap_hi(3))
    integer, intent(in) :: ncomp
    integer, intent(in) :: domlo(3), domhi(3)

    ! note: ncomp is C++ index (0-based)

    real(rt) :: lapx, lapy, lapz
    integer :: i, j, k, n

    lapx = ZERO
    lapy = ZERO
    lapz = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             lapx = U(i+1,j,k,ncomp+1) - TWO*U(i,j,k,ncomp+1) + U(i-1,j,k,ncomp+1)
#if AMREX_SPACEDIM >= 2
             lapy = U(i,j+1,k,ncomp+1) - TWO*U(i,j,k,ncomp+1) + U(i,j-1,k,ncomp+1)
#endif
#if AMREX_SPACEDIM == 3
             lapz = U(i,j,k+1,ncomp+1) - TWO*U(i,j,k,ncomp+1) + U(i,j,k-1,ncomp+1)
#endif

             lap(i,j,k) = lapx + lapy + lapz
          enddo
       enddo
    enddo

    lap(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = &
         TWENTYFOURTH * lap(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

  end subroutine ca_compute_lap_term

  subroutine ca_make_fourth_average(lo, hi, &
                                    q, q_lo, q_hi, nc, &
                                    q_bar, q_bar_lo, q_bar_hi, nc_bar, &
                                    domlo, domhi) &
                                    bind(C, name="ca_make_fourth_average")
    ! Take the cell-center state q and another state q_bar (e.g.,
    ! constructed from the cell-average U) and replace the cell-center
    ! q with a 4th-order accurate cell-average, q <- q + 1/24 L q_bar

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: q_bar_lo(3), q_bar_hi(3)
    integer, intent(in) :: nc, nc_bar
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), nc)
    real(rt), intent(in) :: q_bar(q_bar_lo(1):q_bar_hi(1), q_bar_lo(2):q_bar_hi(2), q_bar_lo(3):q_bar_hi(3), nc_bar)
    integer, intent(in) :: domlo(3), domhi(3)

    integer :: i, j, k, n
    real(rt) :: lapx, lapy, lapz

    lapx = ZERO
    lapy = ZERO
    lapz = ZERO

    do n = 1, nc
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                lapx = q_bar(i+1,j,k,n) - TWO*q_bar(i,j,k,n) + q_bar(i-1,j,k,n)
#if AMREX_SPACEDIM >= 2
                lapy = q_bar(i,j+1,k,n) - TWO*q_bar(i,j,k,n) + q_bar(i,j-1,k,n)
#endif
#if AMREX_SPACEDIM == 3
                lapz = q_bar(i,j,k+1,n) - TWO*q_bar(i,j,k,n) + q_bar(i,j,k-1,n)
#endif

                q(i,j,k,n) = q(i,j,k,n) + TWENTYFOURTH * (lapx + lapy + lapz)

             enddo
          enddo
       enddo
    enddo

  end subroutine ca_make_fourth_average

  subroutine ca_make_fourth_in_place(lo, hi, &
                                     q, q_lo, q_hi, nc, &
                                     domlo, domhi) &
                                     bind(C, name="ca_make_fourth_in_place")
    ! Take the cell-center q and makes it a cell-average q, in place
    ! (e.g. q is overwritten by its average), q <- q + 1/24 L q.
    ! Note: this routine is not tile safe.

    use amrex_mempool_module, only : bl_allocate, bl_deallocate

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: nc
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), nc)
    integer, intent(in) :: domlo(3), domhi(3)

    integer :: i, j, k, n
    real(rt) :: lapx, lapy, lapz
    real(rt), pointer :: lap(:,:,:)

    call bl_allocate(lap, lo, hi)

    lapx = ZERO
    lapy = ZERO
    lapz = ZERO

    do n = 1, nc

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                lapx = q(i+1,j,k,n) - TWO*q(i,j,k,n) + q(i-1,j,k,n)
#if AMREX_SPACEDIM >= 2
                lapy = q(i,j+1,k,n) - TWO*q(i,j,k,n) + q(i,j-1,k,n)
#endif
#if AMREX_SPACEDIM == 3
                lapz = q(i,j,k+1,n) - TWO*q(i,j,k,n) + q(i,j,k-1,n)
#endif

                lap(i,j,k) = lapx + lapy + lapz
             end do
          end do
       end do

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                q(i,j,k,n) = q(i,j,k,n) + TWENTYFOURTH * lap(i,j,k)
             end do
          end do
       end do
    end do

    call bl_deallocate(lap)

  end subroutine ca_make_fourth_in_place

  subroutine ca_make_fourth_in_place_n(lo, hi, &
                                       q, q_lo, q_hi, nc, ncomp, &
                                       domlo, domhi) &
                                       bind(C, name="ca_make_fourth_in_place_n")
    ! Take the cell-center q and makes it a cell-average q, in place
    ! (e.g. q is overwritten by its average), q <- q + 1/24 L q.
    ! Note: this routine is not tile safe.
    !
    ! This version operates on a single component.  Here ncomp is the
    ! component to update -- we expect this to come in 0-based from
    ! C++, so we add 1

    use amrex_mempool_module, only : bl_allocate, bl_deallocate

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: nc, ncomp
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), nc)
    integer, intent(in) :: domlo(3), domhi(3)

    integer :: i, j, k, n
    real(rt), pointer :: lap(:,:,:)
    real(rt) :: lapx, lapy, lapz

    call bl_allocate(lap, lo, hi)

    lapx = ZERO
    lapy = ZERO
    lapz = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             lapx = q(i+1,j,k,ncomp+1) - TWO*q(i,j,k,ncomp+1) + q(i-1,j,k,ncomp+1)
#if AMREX_SPACEDIM >= 2
             lapy = q(i,j+1,k,ncomp+1) - TWO*q(i,j,k,ncomp+1) + q(i,j-1,k,ncomp+1)
#endif
#if AMREX_SPACEDIM == 3
             lapz = q(i,j,k+1,ncomp+1) - TWO*q(i,j,k,ncomp+1) + q(i,j,k-1,ncomp+1)
#endif

             lap(i,j,k) = lapx + lapy + lapz
          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             q(i,j,k,ncomp+1) = q(i,j,k,ncomp+1) + TWENTYFOURTH * lap(i,j,k)
          enddo
       enddo
    enddo

    call bl_deallocate(lap)

  end subroutine ca_make_fourth_in_place_n


end module fourth_order
