module fourth_order

  use amrex_constants_module, only : ZERO, TWO, ONE
  use prob_params_module, only : dg
  use castro_error_module, only : castro_error
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), parameter :: TWENTYFOURTH = ONE/24.0_rt

contains

  subroutine fourth_interfaces(idir, ncomp, &
                               a, a_lo, a_hi, &
                               a_int,  ai_lo, ai_hi, &
                               lo, hi, &
                               domlo, domhi)

    ! this just computes the unlimited single-value interface state
    ! for the 4th order method.

    use meth_params_module, only : NQ
    use prob_params_module, only : Interior, Symmetry, Outflow, physbc_lo, physbc_hi

    implicit none

    integer, intent(in) :: idir, ncomp
    integer, intent(in) :: a_lo(3), a_hi(3)
    integer, intent(in) :: ai_lo(3), ai_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in) :: a(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3), NQ)
    real(rt), intent(inout) :: a_int(ai_lo(1):ai_hi(1), ai_lo(2):ai_hi(2), ai_lo(3):ai_hi(3))
    integer, intent(in) :: domlo(3), domhi(3)

    ! local variables
    integer :: i, j, k
    real(rt) :: ac

    ! our convention here is that:
    !     al(i,j,k)   will be al_{i-1/2,j,k),
    !     al(i+1,j,k) will be al_{i+1/2,j,k)

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
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i-1,j,k,ncomp) + 13.0_rt*a(i,j,k,ncomp) - &
                                                    5.0_rt*a(i+1,j,k,ncomp) + a(i+2,j,k,ncomp))

                else if (i == domlo(1) .and. physbc_lo(1) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! left physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k,ncomp) - 23.0_rt*a(i+1,j,k,ncomp) + &
                                                    13.0_rt*a(i+2,j,k,ncomp) - 3.0_rt*a(i+3,j,k,ncomp))

                else if (i == domhi(1) .and. physbc_hi(1) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the right physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k,ncomp) + 13.0_rt*a(i-1,j,k,ncomp) - &
                                                    5.0_rt*a(i-2,j,k,ncomp) + a(i-3,j,k,ncomp))

                else if (i == domhi(1)+1 .and. physbc_hi(1) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! right physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i-1,j,k,ncomp) - 23.0_rt*a(i-2,j,k,ncomp) + &
                                                    13.0_rt*a(i-3,j,k,ncomp) - 3.0_rt*a(i-4,j,k,ncomp))

                else if (i == domlo(1)-1 .and. physbc_lo(1) == Outflow) then
                   ! extrapolate to the domlo(1)-1 cell using a
                   ! conservative cubic polynomial averaged over
                   ! the cell
                   ac = 4.0_rt*a(domlo(1),j,k,ncomp) - 6.0_rt*a(domlo(1)+1,j,k,ncomp) + 4.0_rt*a(domlo(1)+2,j,k,ncomp) - a(domlo(1)+3,j,k,ncomp)

                   ! now use the 1-sided stencil from above with
                   ! this extrapolated value
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*ac - 23.0_rt*a(i+1,j,k,ncomp) + &
                                                    13.0_rt*a(i+2,j,k,ncomp) - 3.0_rt*a(i+3,j,k,ncomp))

                else if (i == domhi(1)+2 .and. physbc_hi(1) == Outflow) then
                   ! extrapolate to the domhi(1)+1 cell using a
                   ! conservative cubic polynomial averaged over
                   ! the cell
                   ac = 4.0_rt*a(domhi(1),j,k,ncomp) - 6.0_rt*a(domhi(1)-1,j,k,ncomp) + 4.0_rt*a(domhi(1)-2,j,k,ncomp) - a(domhi(1)-3,j,k,ncomp)

                   ! now use the 1-sided stencil from above with
                   ! this extrapolated value
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*ac - 23.0_rt*a(i-2,j,k,ncomp) + &
                                                    13.0_rt*a(i-3,j,k,ncomp) - 3.0_rt*a(i-4,j,k,ncomp))

                else
                   ! regular stencil
                   a_int(i,j,k) = (7.0_rt/12.0_rt)*(a(i-1,j,k,ncomp) + a(i,j,k,ncomp)) - &
                        (1.0_rt/12.0_rt)*(a(i-2,j,k,ncomp) + a(i+1,j,k,ncomp))
                end if

             end do
          end do
       end do

    else if (idir == 2) then

       ! this loop is over interfaces
       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2)-1, hi(2)+2
             do i = lo(1)-1, hi(1)+1

                ! interpolate to the edges
                if (j == domlo(2)+1 .and. physbc_lo(2) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the left physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j-1,k,ncomp) + 13.0_rt*a(i,j,k,ncomp) - &
                                                    5.0_rt*a(i,j+1,k,ncomp) + a(i,j+2,k,ncomp))

                else if (j == domlo(2) .and. physbc_lo(2) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! left physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k,ncomp) - 23.0_rt*a(i,j+1,k,ncomp) + &
                                                    13.0_rt*a(i,j+2,k,ncomp) - 3.0_rt*a(i,j+3,k,ncomp))

                else if (j == domhi(2) .and. physbc_hi(2) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the right physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k,ncomp) + 13.0_rt*a(i,j-1,k,ncomp) - &
                                                    5.0_rt*a(i,j-2,k,ncomp) + a(i,j-3,k,ncomp))

                else if (j == domhi(2)+1 .and. physbc_hi(2) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! right physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j-1,k,ncomp) - 23.0_rt*a(i,j-2,k,ncomp) + &
                                                    13.0_rt*a(i,j-3,k,ncomp) - 3.0_rt*a(i,j-4,k,ncomp))

                else if (j == domlo(2)-1 .and. physbc_lo(2) == Outflow) then
                   ! extrapolate to the domlo(1)-1 cell using a
                   ! conservative cubic polynomial averaged over
                   ! the cell
                   ac = 4.0_rt*a(i,domlo(2),k,ncomp) - 6.0_rt*a(i,domlo(2)+1,k,ncomp) + 4.0_rt*a(i,domlo(2)+2,k,ncomp) - a(i,domlo(2)+3,k,ncomp)

                   ! now use the 1-sided stencil from above with
                   ! this extrapolated value
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*ac - 23.0_rt*a(i,j+1,k,ncomp) + &
                                                    13.0_rt*a(i,j+2,k,ncomp) - 3.0_rt*a(i,j+3,k,ncomp))

                else if (j == domhi(2)+2 .and. physbc_hi(2) == Outflow) then
                   ! extrapolate to the domhi(1)+1 cell using a
                   ! conservative cubic polynomial averaged over
                   ! the cell
                   ac = 4.0_rt*a(i,domhi(2),k,ncomp) - 6.0_rt*a(i,domhi(2)-1,k,ncomp) + 4.0_rt*a(i,domhi(2)-2,k,ncomp) - a(i,domhi(2)-3,k,ncomp)

                   ! now use the 1-sided stencil from above with
                   ! this extrapolated value
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*ac - 23.0_rt*a(i,j-2,k,ncomp) + &
                                                    13.0_rt*a(i,j-3,k,ncomp) - 3.0_rt*a(i,j-4,k,ncomp))

                else
                   ! regular stencil
                   a_int(i,j,k) = (7.0_rt/12.0_rt)*(a(i,j-1,k,ncomp) + a(i,j,k,ncomp)) - &
                        (1.0_rt/12.0_rt)*(a(i,j-2,k,ncomp) + a(i,j+1,k,ncomp))
                end if

             end do
          end do
       end do

    else if (idir == 3) then

       ! this loop is over interfaces
       do k = lo(3)-1, hi(3)+2
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

                ! interpolate to the edges
                if (k == domlo(3)+1 .and. physbc_lo(3) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the left physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k-1,ncomp) + 13.0_rt*a(i,j,k,ncomp) - &
                                                    5.0_rt*a(i,j,k+1,ncomp) + a(i,j,k+2,ncomp))

                else if (k == domlo(3) .and. physbc_lo(3) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! left physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k,ncomp) - 23.0_rt*a(i,j,k+1,ncomp) + &
                                                    13.0_rt*a(i,j,k+2,ncomp) - 3.0_rt*a(i,j,k+3,ncomp))

                else if (k == domhi(3) .and. physbc_hi(3) /= Interior) then
                   ! use a stencil for the interface that is one zone
                   ! from the right physical boundary, MC Eq. 22
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k,ncomp) + 13.0_rt*a(i,j,k-1,ncomp) - &
                                                    5.0_rt*a(i,j,k-2,ncomp) + a(i,j,k-3,ncomp))

                else if (k == domhi(3)+1 .and. physbc_hi(3) /= Interior) then
                   ! use a stencil for when the interface is on the
                   ! right physical boundary MC Eq. 21
                   a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k-1,ncomp) - 23.0_rt*a(i,j,k-2,ncomp) + &
                                                    13.0_rt*a(i,j,k-3,ncomp) - 3.0_rt*a(i,j,k-4,ncomp))

                else
                   ! regular stencil
                   a_int(i,j,k) = (7.0_rt/12.0_rt)*(a(i,j,k-1,ncomp) + a(i,j,k,ncomp)) - &
                        (1.0_rt/12.0_rt)*(a(i,j,k-2,ncomp) + a(i,j,k+1,ncomp))
                end if

             end do
          end do
       end do

    end if

  end subroutine fourth_interfaces

  subroutine states(idir, ncomp, &
                    a, a_lo, a_hi, &
                    a_int, ai_lo, ai_hi, &
                    flatn, f_lo, f_hi, &
                    al, ar, as_lo, as_hi, &
                    lo, hi, &
                    domlo, domhi)

    use meth_params_module, only : NQ, QU, QV, QW, limit_fourth_order
    use prob_params_module, only : Interior, Symmetry, Outflow, physbc_lo, physbc_hi

    implicit none

    integer, intent(in) :: idir, ncomp
    integer, intent(in) :: a_lo(3), a_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: ai_lo(3), ai_hi(3)
    integer, intent(in) :: as_lo(3), as_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in) :: a(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3), NQ)
    real(rt), intent(in) :: a_int(ai_lo(1):ai_hi(1), ai_lo(2):ai_hi(2), ai_lo(3):ai_hi(3))
    real(rt), intent(in) :: flatn( f_lo(1): f_hi(1), f_lo(2): f_hi(2), f_lo(3): f_hi(3))

    real(rt), intent(inout) :: al(as_lo(1):as_hi(1), as_lo(2):as_hi(2), as_lo(3):as_hi(3), NQ)
    real(rt), intent(inout) :: ar(as_lo(1):as_hi(1), as_lo(2):as_hi(2), as_lo(3):as_hi(3), NQ)
    integer, intent(in) :: domlo(3), domhi(3)

    ! local variables
    real(rt) :: dafm, dafp, d2af
    real(rt) :: d2acm2, d2acm1, d2ac0, d2acp1, d2acp2
    real(rt) :: d3am1, d3a0, d3ap1, d3ap2
    real(rt) :: ac

    real(rt), parameter :: C2 = 1.25_rt
    real(rt), parameter :: C3 = 0.1_rt

    integer :: i, j, k
    real(rt) :: rho, s

    real(rt) :: d2a_lim, d3a_min, d3a_max

    ! our convention here is that:
    !     al(i,j,k)   will be al_{i-1/2,j,k),
    !     al(i+1,j,k) will be al_{i+1/2,j,k)

    ! we need interface values on all faces of the domain
    if (idir == 1) then

       if (limit_fourth_order == 0) then

          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-dg(2), hi(2)+dg(2)
                do i = lo(1)-1, hi(1)+1

                   al(i+1,j,k,ncomp) = a_int(i+1,j,k)
                   ar(i,j,k,ncomp) = a_int(i,j,k)

                end do
             end do
          end do

       else

          ! the limiting loops are now over zones

          ! this is a loop over cell centers, affecting
          ! i-1/2,R and i+1/2,L
          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-dg(2), hi(2)+dg(2)
                do i = lo(1)-1, hi(1)+1

                   ! these live on cell-centers
                   dafm = a(i,j,k,ncomp) - a_int(i,j,k)
                   dafp = a_int(i+1,j,k) - a(i,j,k,ncomp)

                   ! these live on cell-centers
                   d2af = 6.0_rt*(a_int(i,j,k) - 2.0_rt*a(i,j,k,ncomp) + a_int(i+1,j,k))

                   d2acm2 = a(i-3,j,k,ncomp) - 2.0_rt*a(i-2,j,k,ncomp) + a(i-1,j,k,ncomp)
                   d2acm1 = a(i-2,j,k,ncomp) - 2.0_rt*a(i-1,j,k,ncomp) + a(i,j,k,ncomp)
                   d2ac0  = a(i-1,j,k,ncomp) - 2.0_rt*a(i,j,k,ncomp) + a(i+1,j,k,ncomp)
                   d2acp1 = a(i,j,k,ncomp) - 2.0_rt*a(i+1,j,k,ncomp) + a(i+2,j,k,ncomp)
                   d2acp2 = a(i+1,j,k,ncomp) - 2.0_rt*a(i+2,j,k,ncomp) + a(i+3,j,k,ncomp)

                   ! limit? MC Eq. 24 and 25
                   if (dafm * dafp <= 0.0_rt .or. &
                        (a(i,j,k,ncomp) - a(i-2,j,k,ncomp))*(a(i+2,j,k,ncomp) - a(i,j,k,ncomp)) <= 0.0_rt) then

                      ! we are at an extrema

                      s = sign(1.0_rt, d2ac0)
                      if ( s == sign(1.0_rt, d2acm1) .and. &
                           s == sign(1.0_rt, d2acp1) .and. &
                           s == sign(1.0_rt, d2af)) then
                         ! MC Eq. 26
                         d2a_lim = s*min(abs(d2af), C2*abs(d2acm1), &
                                         C2*abs(d2ac0), C2*abs(d2acp1))
                      else
                         d2a_lim = 0.0_rt
                      end if

                      if (abs(d2af) <= 1.e-12_rt*max(abs(a(i-2,j,k,ncomp)), abs(a(i-1,j,k,ncomp)), &
                           abs(a(i,j,k,ncomp)), abs(a(i+1,j,k,ncomp)), abs(a(i+2,j,k,ncomp)))) then
                         rho = 0.0_rt
                      else
                         ! MC Eq. 27
                         rho = d2a_lim/d2af
                      end if

                      if (rho < 1.0_rt - 1.e-12_rt) then
                         ! we may need to limit -- these quantities are at cell-centers
                         d3am1 = d2acm1 - d2acm2
                         d3a0  = d2ac0 - d2acm1
                         d3ap1 = d2acp1 - d2ac0
                         d3ap2 = d2acp2 - d2acp1

                         d3a_min = min(d3am1, d3a0, d3ap1, d3ap2)
                         d3a_max = max(d3am1, d3a0, d3ap1, d3ap2)

                         if (C3*max(abs(d3a_min), abs(d3a_max)) <= (d3a_max - d3a_min)) then
                            ! limit
                            if (dafm*dafp < 0.0_rt) then
                               ! Eqs. 29, 30
                               ar(i,j,k,ncomp) = a(i,j,k,ncomp) - rho*dafm  ! note: typo in Eq 29
                               al(i+1,j,k,ncomp) = a(i,j,k,ncomp) + rho*dafp
                            else if (abs(dafm) >= 2.0_rt*abs(dafp)) then
                               ! Eq. 31
                               ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*(1.0_rt - rho)*dafp - rho*dafm
                            else if (abs(dafp) >= 2.0_rt*abs(dafm)) then
                               ! Eq. 32
                               al(i+1,j,k,ncomp) = a(i,j,k,ncomp) + 2.0_rt*(1.0_rt - rho)*dafm + rho*dafp
                            end if

                         end if
                      end if

                   else
                      ! if Eqs. 24 or 25 didn't hold we still may need to limit
                      if (abs(dafm) >= 2.0_rt*abs(dafp)) then
                         ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*dafp
                      end if
                      if (abs(dafp) >= 2.0_rt*abs(dafm)) then
                         al(i+1,j,k,ncomp) = a(i,j,k,ncomp) + 2.0_rt*dafm
                      end if
                   end if

                   ! apply flattening
                   al(i+1,j,k,ncomp) = flatn(i,j,k)*al(i+1,j,k,ncomp) + (ONE - flatn(i,j,k))*a(i,j,k,ncomp)
                   ar(i,j,k,ncomp) = flatn(i,j,k)*ar(i,j,k,ncomp) + (ONE - flatn(i,j,k))*a(i,j,k,ncomp)


                   ! now handle any physical boundaries here by modifying the interface values

                   if (i == domlo(1)) then

                      ! reset the left state at domlo(1) if needed -- it is outside the domain

                      if (physbc_lo(1) == Outflow) then
                         !al(domlo(1),j,k,:) = ar(domlo(1),j,k,:)
                         continue

                      else if (physbc_lo(1) == Symmetry) then
                         if (ncomp == QU) then
                            al(domlo(1),j,k,QU) = -ar(domlo(1),j,k,QU)
                         else
                            al(domlo(1),j,k,ncomp) = ar(domlo(1),j,k,ncomp)
                         end if

                      else if (physbc_lo(1) == Interior) then
                         ! we don't need to do anything here
                         continue

                      else
                         ! not supported
                         call castro_error("ERROR: boundary conditions not supported for -X in states")
                      end if

                   end if

                   if (i == domhi(1)+1) then

                      ! reset the right state at domhi(1)+1 if needed -- it is outside the domain

                      if (physbc_hi(1) == Outflow) then
                         !ar(domhi(1)+1,j,k,:) = al(domhi(1)+1,j,k,:)
                         continue

                      else if (physbc_hi(1) == Symmetry) then
                         if (ncomp == QU) then
                            ar(domhi(1)+1,j,k,QU) = -al(domhi(1)+1,j,k,QU)
                         else
                            ar(domhi(1)+1,j,k,ncomp) = al(domhi(1)+1,j,k,ncomp)
                         end if

                      else if (physbc_hi(1) == Interior) then
                         ! we don't need to do anything here
                         continue

                      else
                         ! not supported
                         call castro_error("ERROR: boundary conditions not supported for +X in states")
                      end if

                   end if

                end do
             end do
          end do

       end if


    else if (idir == 2) then

       if (limit_fourth_order == 0) then

          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1

                   al(i,j+1,k,ncomp) = a_int(i,j+1,k)
                   ar(i,j,k,ncomp) = a_int(i,j,k)

                end do
             end do
          end do

       else

          ! the limiting loops are now over zones

          ! this is a loop over cell centers, affecting
          ! j-1/2,R and j+1/2,L
          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1

                   ! these live on cell-centers
                   dafm = a(i,j,k,ncomp) - a_int(i,j,k)
                   dafp = a_int(i,j+1,k) - a(i,j,k,ncomp)

                   ! these live on cell-centers
                   d2af = 6.0_rt*(a_int(i,j,k) - 2.0_rt*a(i,j,k,ncomp) + a_int(i,j+1,k))

                   d2acm2 = a(i,j-3,k,ncomp) - 2.0_rt*a(i,j-2,k,ncomp) + a(i,j-1,k,ncomp)
                   d2acm1 = a(i,j-2,k,ncomp) - 2.0_rt*a(i,j-1,k,ncomp) + a(i,j,k,ncomp)
                   d2ac0  = a(i,j-1,k,ncomp) - 2.0_rt*a(i,j,k,ncomp) + a(i,j+1,k,ncomp)
                   d2acp1 = a(i,j,k,ncomp) - 2.0_rt*a(i,j+1,k,ncomp) + a(i,j+2,k,ncomp)
                   d2acp2 = a(i,j+1,k,ncomp) - 2.0_rt*a(i,j+2,k,ncomp) + a(i,j+3,k,ncomp)

                   ! limit? MC Eq. 24 and 25
                   if (dafm * dafp <= 0.0_rt .or. &
                        (a(i,j,k,ncomp) - a(i,j-2,k,ncomp))*(a(i,j+2,k,ncomp) - a(i,j,k,ncomp)) <= 0.0_rt) then

                      ! we are at an extrema

                      s = sign(1.0_rt, d2ac0)
                      if ( s == sign(1.0_rt, d2acm1) .and. &
                           s == sign(1.0_rt, d2acp1) .and. &
                           s == sign(1.0_rt, d2af)) then
                         ! MC Eq. 26
                         d2a_lim = s*min(abs(d2af), C2*abs(d2acm1), &
                                         C2*abs(d2ac0), C2*abs(d2acp1))
                      else
                         d2a_lim = 0.0_rt
                      end if

                      if (abs(d2af) <= 1.e-12_rt*max(abs(a(i,j-2,k,ncomp)), abs(a(i,j-1,k,ncomp)), &
                           abs(a(i,j,k,ncomp)), abs(a(i,j+1,k,ncomp)), abs(a(i,j+2,k,ncomp)))) then
                         rho = 0.0_rt
                      else
                         ! MC Eq. 27
                         rho = d2a_lim/d2af
                      end if

                      if (rho < 1.0_rt - 1.e-12_rt) then
                         ! we may need to limit -- these quantities are at cell-centers
                         d3am1 = d2acm1 - d2acm2
                         d3a0  = d2ac0 - d2acm1
                         d3ap1 = d2acp1 - d2ac0
                         d3ap2 = d2acp2 - d2acp1

                         d3a_min = min(d3am1, d3a0, d3ap1, d3ap2)
                         d3a_max = max(d3am1, d3a0, d3ap1, d3ap2)

                         if (C3*max(abs(d3a_min), abs(d3a_max)) <= (d3a_max - d3a_min)) then
                            ! limit
                            if (dafm*dafp < 0.0_rt) then
                               ! Eqs. 29, 30
                               ar(i,j,k,ncomp) = a(i,j,k,ncomp) - rho*dafm  ! note: typo in Eq 29
                               al(i,j+1,k,ncomp) = a(i,j,k,ncomp) + rho*dafp
                            else if (abs(dafm) >= 2.0_rt*abs(dafp)) then
                               ! Eq. 31
                               ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*(1.0_rt - rho)*dafp - rho*dafm
                            else if (abs(dafp) >= 2.0_rt*abs(dafm)) then
                               ! Eq. 32
                               al(i,j+1,k,ncomp) = a(i,j,k,ncomp) + 2.0_rt*(1.0_rt - rho)*dafm + rho*dafp
                            end if

                         end if
                      end if

                   else
                      ! if Eqs. 24 or 25 didn't hold we still may need to limit
                      if (abs(dafm) >= 2.0_rt*abs(dafp)) then
                         ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*dafp
                      end if
                      if (abs(dafp) >= 2.0_rt*abs(dafm)) then
                         al(i,j+1,k,ncomp) = a(i,j,k,ncomp) + 2.0_rt*dafm
                      end if
                   end if

                   ! apply flattening
                   al(i,j+1,k,ncomp) = flatn(i,j,k)*al(i,j+1,k,ncomp) + (ONE - flatn(i,j,k))*a(i,j,k,ncomp)
                   ar(i,j,k,ncomp) = flatn(i,j,k)*ar(i,j,k,ncomp) + (ONE - flatn(i,j,k))*a(i,j,k,ncomp)

                   ! now handle any physical boundaries here by modifying the interface values

                   if (j == domlo(2)) then
                      ! reset the left state at domlo(2) if needed -- it is outside the domain

                      if (physbc_lo(2) == Outflow) then
                         !al(i,domlo(2),k,:) = ar(i,domlo(2),k,:)
                         continue

                      else if (physbc_lo(2) == Symmetry) then
                         if (ncomp == QV) then
                            al(i,domlo(2),k,QV) = -ar(i,domlo(2),k,QV)
                         else
                            al(i,domlo(2),k,ncomp) = ar(i,domlo(2),k,ncomp)
                         end if

                      else if (physbc_lo(2) == Interior) then
                         ! we don't need to do anything here
                         continue

                      else
                         ! not supported
                         call castro_error("ERROR: boundary conditions not supported for -Y in states")
                      end if

                   end if

                   if (j == domhi(2)+1) then
                      ! reset the right state at domhi(2)+1 if needed -- it is outside the domain

                      if (physbc_hi(2) == Outflow) then
                         !ar(i,domhi(2)+1,k,:) = al(i,domhi(2)+1,k,:)
                         continue

                      else if (physbc_hi(2) == Symmetry) then
                         if (ncomp == QV) then
                            ar(i,domhi(2)+1,k,QV) = -al(i,domhi(2)+1,k,QV)
                         else
                            ar(i,domhi(2)+1,k,ncomp) = al(i,domhi(2)+1,k,ncomp)
                         end if

                      else if (physbc_hi(2) == Interior) then
                         ! we don't need to do anything here
                         continue

                      else
                         ! not supported
                         call castro_error("ERROR: boundary conditions not supported for +Y in states")
                      end if

                   end if

                end do
             end do
          end do

       end if

    else if (idir == 3) then

       if (limit_fourth_order == 0) then

          do k = lo(3)-1, hi(3)+1
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1

                   al(i,j,k+1,ncomp) = a_int(i,j,k+1)
                   ar(i,j,k,ncomp) = a_int(i,j,k)

                end do
             end do
          end do

       else

          ! this is a loop over cell centers, affecting
          ! j-1/2,R and j+1/2,L
          do k = lo(3)-1, hi(3)+1
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1

                   ! these live on cell-centers
                   dafm = a(i,j,k,ncomp) - a_int(i,j,k)
                   dafp = a_int(i,j,k+1) - a(i,j,k,ncomp)

                   ! these live on cell-centers
                   d2af = 6.0_rt*(a_int(i,j,k) - 2.0_rt*a(i,j,k,ncomp) + a_int(i,j,k+1))

                   d2acm2 = a(i,j,k-3,ncomp) - 2.0_rt*a(i,j,k-2,ncomp) + a(i,j,k-1,ncomp)
                   d2acm1 = a(i,j,k-2,ncomp) - 2.0_rt*a(i,j,k-1,ncomp) + a(i,j,k,ncomp)
                   d2ac0  = a(i,j,k-1,ncomp) - 2.0_rt*a(i,j,k,ncomp) + a(i,j,k+1,ncomp)
                   d2acp1 = a(i,j,k,ncomp) - 2.0_rt*a(i,j,k+1,ncomp) + a(i,j,k+2,ncomp)
                   d2acp2 = a(i,j,k+1,ncomp) - 2.0_rt*a(i,j,k+2,ncomp) + a(i,j,k+3,ncomp)

                   ! limit? MC Eq. 24 and 25
                   if (dafm * dafp <= 0.0_rt .or. &
                        (a(i,j,k,ncomp) - a(i,j,k-2,ncomp))*(a(i,j,k+2,ncomp) - a(i,j,k,ncomp)) <= 0.0_rt) then

                      ! we are at an extrema

                      s = sign(1.0_rt, d2ac0)
                      if ( s == sign(1.0_rt, d2acm1) .and. &
                           s == sign(1.0_rt, d2acp1) .and. &
                           s == sign(1.0_rt, d2af)) then
                         ! MC Eq. 26
                         d2a_lim = s*min(abs(d2af), C2*abs(d2acm1), &
                                         C2*abs(d2ac0), C2*abs(d2acp1))
                      else
                         d2a_lim = 0.0_rt
                      end if

                      if (abs(d2af) <= 1.e-12_rt*max(abs(a(i,j,k-2,ncomp)), abs(a(i,j,k-1,ncomp)), &
                           abs(a(i,j,k,ncomp)), abs(a(i,j,k+1,ncomp)), abs(a(i,j,k+2,ncomp)))) then
                         rho = 0.0_rt
                      else
                         ! MC Eq. 27
                         rho = d2a_lim/d2af
                      end if

                      if (rho < 1.0_rt - 1.e-12_rt) then
                         ! we may need to limit -- these quantities are at cell-centers
                         d3am1 = d2acm1 - d2acm2
                         d3a0  = d2ac0 - d2acm1
                         d3ap1 = d2acp1 - d2ac0
                         d3ap2 = d2acp2 - d2acp1

                         d3a_min = min(d3am1, d3a0, d3ap1, d3ap2)
                         d3a_max = max(d3am1, d3a0, d3ap1, d3ap2)

                         if (C3*max(abs(d3a_min), abs(d3a_max)) <= (d3a_max - d3a_min)) then
                            ! limit
                            if (dafm*dafp < 0.0_rt) then
                               ! Eqs. 29, 30
                               ar(i,j,k,ncomp) = a(i,j,k,ncomp) - rho*dafm  ! note: typo in Eq 29
                               al(i,j,k+1,ncomp) = a(i,j,k,ncomp) + rho*dafp
                            else if (abs(dafm) >= 2.0_rt*abs(dafp)) then
                               ! Eq. 31
                               ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*(1.0_rt - rho)*dafp - rho*dafm
                            else if (abs(dafp) >= 2.0_rt*abs(dafm)) then
                               ! Eq. 32
                               al(i,j,k+1,ncomp) = a(i,j,k,ncomp) + 2.0_rt*(1.0_rt - rho)*dafm + rho*dafp
                            end if

                         end if
                      end if

                   else
                      ! if Eqs. 24 or 25 didn't hold we still may need to limit
                      if (abs(dafm) >= 2.0_rt*abs(dafp)) then
                         ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*dafp
                      end if
                      if (abs(dafp) >= 2.0_rt*abs(dafm)) then
                         al(i,j,k+1,ncomp) = a(i,j,k,ncomp) + 2.0_rt*dafm
                      end if
                   end if

                   ! apply flattening
                   al(i,j,k+1,ncomp) = flatn(i,j,k)*al(i,j,k+1,ncomp) + (ONE - flatn(i,j,k))*a(i,j,k,ncomp)
                   ar(i,j,k,ncomp) = flatn(i,j,k)*ar(i,j,k,ncomp) + (ONE - flatn(i,j,k))*a(i,j,k,ncomp)


                   ! now handle any physical boundaries here by modifying the interface values

                   if (k == domlo(3)) then
                      ! reset the left state at domlo(3) if needed -- it is outside the domain

                      if (physbc_lo(3) == Outflow) then
                         !al(i,j,domlo(3),:) = ar(i,j,domlo(3),:)
                         continue

                      else if (physbc_lo(3) == Symmetry) then
                         if (ncomp == QW) then
                            al(i,j,domlo(3),QW) = -ar(i,j,domlo(3),QW)
                         else
                            al(i,j,domlo(3),ncomp) = ar(i,j,domlo(3),ncomp)
                         end if

                      else if (physbc_lo(3) == Interior) then
                         ! we don't need to do anything here
                         continue

                      else
                         ! not supported
                         call castro_error("ERROR: boundary conditions not supported for -Z in states")
                      end if

                   end if

                   if (k == domhi(3)+1) then
                      ! reset the right state at domhi(3)+1 if needed -- it is outside the domain

                      if (physbc_hi(3) == Outflow) then
                         !ar(i,j,domhi(3)+1,:) = al(i,j,domhi(3)+1,:)
                         continue

                      else if (physbc_hi(3) == Symmetry) then
                         if (ncomp == QW) then
                            ar(i,j,domhi(3)+1,QW) = -al(i,j,domhi(3)+1,QW)
                         else
                            ar(i,j,domhi(3)+1,ncomp) = al(i,j,domhi(3)+1,ncomp)
                         end if

                      else if (physbc_lo(3) == Interior) then
                         ! we don't need to do anything here
                         continue

                      else
                         ! not supported
                         call castro_error("ERROR: boundary conditions not supported for +Z in states")
                      end if

                   end if

                end do
             end do
          end do

       end if

    end if

  end subroutine states


  ! Note: pretty much all of these routines below assume that dx(1) = dx(2) = dx(3)
  pure function compute_laplacian(i, j, k, n, &
                                  a, a_lo, a_hi, nc, &
                                  domlo, domhi) result (lap)

    use prob_params_module, only : physbc_lo, physbc_hi, Interior
    implicit none

    integer, intent(in) :: i, j, k, n
    integer, intent(in) :: a_lo(3), a_hi(3)
    integer, intent(in) :: nc
    real(rt), intent(in) :: a(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3), nc)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt) :: lapx, lapy, lapz
    real(rt) :: lap

    lapx = ZERO
    lapy = ZERO
    lapz = ZERO

    ! we use 2nd-order accurate one-sided stencils at the physical
    ! boundaries note: this differs from the suggestion in MC2011 --
    ! they just suggest using the Laplacian from +1 off the interior.
    ! I like the one-sided better.

    if (i == domlo(1) .and. physbc_lo(1) /= Interior) then
       lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i+1,j,k,n) + 4.0_rt*a(i+2,j,k,n) - a(i+3,j,k,n)

    else if (i == domhi(1) .and. physbc_hi(1) /= Interior) then
       lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i-1,j,k,n) + 4.0_rt*a(i-2,j,k,n) - a(i-3,j,k,n)

    else
       lapx = a(i+1,j,k,n) - TWO*a(i,j,k,n) + a(i-1,j,k,n)
    end if

#if AMREX_SPACEDIM >= 2
    if (j == domlo(2) .and. physbc_lo(2) /= Interior) then
       lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j+1,k,n) + 4.0_rt*a(i,j+2,k,n) - a(i,j+3,k,n)

    else if (j == domhi(2) .and. physbc_hi(2) /= Interior) then
       lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j-1,k,n) + 4.0_rt*a(i,j-2,k,n) - a(i,j-3,k,n)

    else
       lapy = a(i,j+1,k,n) - TWO*a(i,j,k,n) + a(i,j-1,k,n)
    end if

#endif

#if AMREX_SPACEDIM == 3
    if (k == domlo(3) .and. physbc_lo(3) /= Interior) then
       lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k+1,n) + 4.0_rt*a(i,j,k+2,n) - a(i,j,k+3,n)

    else if (k == domhi(3) .and. physbc_hi(3) /= Interior) then
       lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k-1,n) + 4.0_rt*a(i,j,k-2,n) - a(i,j,k-3,n)

    else
       lapz = a(i,j,k+1,n) - TWO*a(i,j,k,n) + a(i,j,k-1,n)
    end if
#endif

    lap = lapx + lapy + lapz

  end function compute_laplacian

  pure function transx_laplacian(i, j, k, n, &
                                 a, a_lo, a_hi, nc, &
                                 domlo, domhi) result (lap)

    use prob_params_module, only : physbc_lo, physbc_hi, Interior
    implicit none

    integer, intent(in) :: i, j, k, n
    integer, intent(in) :: a_lo(3), a_hi(3)
    integer, intent(in) :: nc
    real(rt), intent(in) :: a(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3), nc)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt) :: lapx, lapy, lapz
    real(rt) :: lap

    lapy = ZERO
    lapz = ZERO

    ! we use 2nd-order accurate one-sided stencils at the physical boundaries

    if (j == domlo(2) .and. physbc_lo(2) /= Interior) then
       lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j+1,k,n) + 4.0_rt*a(i,j+2,k,n) - a(i,j+3,k,n)

    else if (j == domhi(2) .and. physbc_hi(2) /= Interior) then
       lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j-1,k,n) + 4.0_rt*a(i,j-2,k,n) - a(i,j-3,k,n)

    else
       lapy = a(i,j+1,k,n) - TWO*a(i,j,k,n) + a(i,j-1,k,n)
    end if

#if AMREX_SPACEDIM == 3
    if (k == domlo(3) .and. physbc_lo(3) /= Interior) then
       lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k+1,n) + 4.0_rt*a(i,j,k+2,n) - a(i,j,k+3,n)

    else if (k == domhi(3) .and. physbc_hi(3) /= Interior) then
       lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k-1,n) + 4.0_rt*a(i,j,k-2,n) - a(i,j,k-3,n)

    else
       lapz = a(i,j,k+1,n) - TWO*a(i,j,k,n) + a(i,j,k-1,n)
    end if
#endif

    lap = lapy + lapz

  end function transx_laplacian


  pure function transy_laplacian(i, j, k, n, &
                                 a, a_lo, a_hi, nc, &
                                 domlo, domhi) result (lap)

    use prob_params_module, only : physbc_lo, physbc_hi, Interior
    implicit none

    integer, intent(in) :: i, j, k, n
    integer, intent(in) :: a_lo(3), a_hi(3)
    integer, intent(in) :: nc
    real(rt), intent(in) :: a(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3), nc)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt) :: lapx, lapy, lapz
    real(rt) :: lap

    lapx = ZERO
    lapz = ZERO

    ! we use 2nd-order accurate one-sided stencils at the physical boundaries

    if (i == domlo(1) .and. physbc_lo(1) /= Interior) then
       lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i+1,j,k,n) + 4.0_rt*a(i+2,j,k,n) - a(i+3,j,k,n)

    else if (i == domhi(1) .and. physbc_hi(1) /= Interior) then
       lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i-1,j,k,n) + 4.0_rt*a(i-2,j,k,n) - a(i-3,j,k,n)

    else
       lapx = a(i+1,j,k,n) - TWO*a(i,j,k,n) + a(i-1,j,k,n)
    end if

#if AMREX_SPACEDIM == 3
    if (k == domlo(3) .and. physbc_lo(3) /= Interior) then
       lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k+1,n) + 4.0_rt*a(i,j,k+2,n) - a(i,j,k+3,n)

    else if (k == domhi(3) .and. physbc_hi(3) /= Interior) then
       lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k-1,n) + 4.0_rt*a(i,j,k-2,n) - a(i,j,k-3,n)

    else
       lapz = a(i,j,k+1,n) - TWO*a(i,j,k,n) + a(i,j,k-1,n)
    end if
#endif

    lap = lapx + lapz

  end function transy_laplacian


  pure function transz_laplacian(i, j, k, n, &
                                 a, a_lo, a_hi, nc, &
                                 domlo, domhi) result (lap)

    use prob_params_module, only : physbc_lo, physbc_hi, Interior
    implicit none

    integer, intent(in) :: i, j, k, n
    integer, intent(in) :: a_lo(3), a_hi(3)
    integer, intent(in) :: nc
    real(rt), intent(in) :: a(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3), nc)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt) :: lapx, lapy
    real(rt) :: lap

    lapx = ZERO
    lapy = ZERO

    ! we use 2nd-order accurate one-sided stencils at the physical boundaries

    if (i == domlo(1) .and. physbc_lo(1) /= Interior) then
       lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i+1,j,k,n) + 4.0_rt*a(i+2,j,k,n) - a(i+3,j,k,n)

    else if (i == domhi(1) .and. physbc_hi(1) /= Interior) then
       lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i-1,j,k,n) + 4.0_rt*a(i-2,j,k,n) - a(i-3,j,k,n)

    else
       lapx = a(i+1,j,k,n) - TWO*a(i,j,k,n) + a(i-1,j,k,n)
    end if

    if (j == domlo(2) .and. physbc_lo(2) /= Interior) then
       lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j+1,k,n) + 4.0_rt*a(i,j+2,k,n) - a(i,j+3,k,n)

    else if (j == domhi(2) .and. physbc_hi(2) /= Interior) then
       lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j-1,k,n) + 4.0_rt*a(i,j-2,k,n) - a(i,j-3,k,n)

    else
       lapy = a(i,j+1,k,n) - TWO*a(i,j,k,n) + a(i,j-1,k,n)
    end if

    lap = lapx + lapy

  end function transz_laplacian


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
    real(rt) :: lap

    if (U_lo(1) > lo(1)-1 .or. U_hi(1) < hi(1)+1 .or. &
        (AMREX_SPACEDIM >= 2 .and. (U_lo(2) > lo(2)-1 .or. U_hi(2) < hi(2)+1)) .or. &
        (AMREX_SPACEDIM == 3 .and. (U_lo(3) > lo(3)-1 .or. U_hi(3) < hi(3)+1))) then
       call bl_error("insufficient ghostcells in ca_make_cell_center")
    end if

    do n = 1, nc

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                lap = compute_laplacian(i, j, k, n, &
                                        U, U_lo, U_hi, nc, &
                                        domlo, domhi)

                U_cc(i,j,k,n) = U(i,j,k,n) - TWENTYFOURTH * lap

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

    real(rt), pointer :: lap(:,:,:)

    call bl_allocate(lap, lo, hi)

    do n = 1, nc

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                lap(i,j,k) = compute_laplacian(i, j, k, n, &
                                               U, U_lo, U_hi, nc, &
                                               domlo, domhi)

             end do
          end do
       end do

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                U(i,j,k,n) = U(i,j,k,n) - TWENTYFOURTH * lap(i,j,k)
             end do
          end do
       end do
    end do

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

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             lap(i,j,k) = compute_laplacian(i, j, k, ncomp+1, &
                                            U, U_lo, U_hi, nc, &
                                            domlo, domhi)

          end do
       end do
    end do

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
    real(rt) :: lap

    do n = 1, nc
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                lap = compute_laplacian(i, j, k, n, &
                                        q_bar, q_bar_lo, q_bar_hi, nc, &
                                        domlo, domhi)

                q(i,j,k,n) = q(i,j,k,n) + TWENTYFOURTH * lap

             end do
          end do
       end do
    end do

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
    real(rt), pointer :: lap(:,:,:)

    call bl_allocate(lap, lo, hi)

    do n = 1, nc

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                lap(i,j,k) = compute_laplacian(i, j, k, n, &
                                               q, q_lo, q_hi, nc, &
                                               domlo, domhi)

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

    integer :: i, j, k
    real(rt), pointer :: lap(:,:,:)

    call bl_allocate(lap, lo, hi)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             lap(i,j,k) = compute_laplacian(i, j, k, ncomp+1, &
                                            q, q_lo, q_hi, nc, &
                                            domlo, domhi)

          end do
       end do
    end do

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             q(i,j,k,ncomp+1) = q(i,j,k,ncomp+1) + TWENTYFOURTH * lap(i,j,k)
          end do
       end do
    end do

    call bl_deallocate(lap)

  end subroutine ca_make_fourth_in_place_n


end module fourth_order
