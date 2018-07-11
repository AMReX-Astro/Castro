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
                    lo, hi)

    use meth_params_module, only : limit_fourth_order

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

    ! local variables
    real(rt) :: a_int(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt) :: dafm(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt) :: dafp(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt) :: d2af(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt) :: d2ac(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt) :: d3a(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))

    real(rt), parameter :: C2 = 1.25d0
    real(rt), parameter :: C3 = 0.1d0

    integer :: i, j, k
    double precision :: rho, s

    double precision :: d2a_lim, d3a_min, d3a_max

    ! our convention here is that:
    !     al(i,j,k)   will be al_{i-1/2,j,k),
    !     al(i+1,j,k) will be al_{i+1/2,j,k)

    ! we need interface values on all faces of the domain
    if (idir == 1) then

       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2)-dg(2), hi(2)+dg(2)
             do i = lo(1)-2, hi(1)+3

                ! interpolate to the edges
                a_int(i,j,k) = (7.0d0/12.0d0)*(a(i-1,j,k,n) + a(i,j,k,n)) - &
                     (1.0d0/12.0d0)*(a(i-2,j,k,n) + a(i+1,j,k,n))

                al(i,j,k,n) = a_int(i,j,k)
                ar(i,j,k,n) = a_int(i,j,k)

             enddo
          enddo
       enddo

       if (limit_fourth_order == 1) then
          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-dg(2), hi(2)+dg(2)
                do i = lo(1)-2, hi(1)+3
                   ! these live on cell-centers
                   dafm(i,j,k) = a(i,j,k,n) - a_int(i,j,k)
                   dafp(i,j,k) = a_int(i+1,j,k) - a(i,j,k,n)
                   
                   ! these live on cell-centers
                   d2af(i,j,k) = 6.0d0*(a_int(i,j,k) - 2.0d0*a(i,j,k,n) + a_int(i+1,j,k))
                enddo
             enddo
          enddo

          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-dg(2), hi(2)+dg(2)
                do i = lo(1)-3, hi(1)+3
                   d2ac(i,j,k) = a(i-1,j,k,n) - 2.0d0*a(i,j,k,n) + a(i+1,j,k,n)
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

          ! this is a look over cell centers, affecting
          ! i-1/2,R and i+1/2,L
          do k = lo(3)-dg(3), hi(3)+dg(3)
             do j = lo(2)-dg(2), hi(2)+dg(2)
                do i = lo(1)-1, hi(1)+1
                   
                   ! limit? MC Eq. 24 and 25
                   if (dafm(i,j,k) * dafp(i,j,k) <= 0.0d0 .or. &
                        (a(i,j,k,n) - a(i-2,j,k,n))*(a(i+2,j,k,n) - a(i,j,k,n)) <= 0.0d0) then
                      
                      ! we are at an extrema
                      
                      s = sign(1.0d0, d2ac(i,j,k))
                      if ( s == sign(1.0d0, d2ac(i-1,j,k)) .and. &
                           s == sign(1.0d0, d2ac(i+1,j,k)) .and. &
                           s == sign(1.0d0, d2af(i,j,k))) then
                         ! MC Eq. 26
                         d2a_lim = s*min(abs(d2af(i,j,k)), C2*abs(d2ac(i-1,j,k)), &
                              C2*abs(d2ac(i,j,k)), C2*abs(d2ac(i+1,j,k)))
                      else
                         d2a_lim = 0.0d0
                      endif
                      
                      if (abs(d2af(i,j,k)) <= 1.e-12*max(abs(a(i-2,j,k,n)), abs(a(i-1,j,k,n)), &
                           abs(a(i,j,k,n)), abs(a(i+1,j,k,n)), abs(a(i+2,j,k,n)))) then
                         rho = 0.0d0
                      else
                         ! MC Eq. 27
                         rho = d2a_lim/d2af(i,j,k)
                      endif
                      
                      if (rho < 1.0d0 - 1.d-12) then
                         ! we may need to limit -- these quantities are at cell-centers
                         d3a_min = min(d3a(i-1,j,k), d3a(i,j,k), d3a(i+1,j,k), d3a(i+2,j,k))
                         d3a_max = max(d3a(i-1,j,k), d3a(i,j,k), d3a(i+1,j,k), d3a(i+2,j,k))
                         
                         if (C3*max(abs(d3a_min), abs(d3a_max)) <= (d3a_max - d3a_min)) then
                            ! limit
                            if (dafm(i,j,k)*dafp(i,j,k) < 0.0d0) then
                               ! Eqs. 29, 30
                               ar(i,j,k,n) = a(i,j,k,n) - rho*dafm(i,j,k)  ! note: typo in Eq 29
                               al(i+1,j,k,n) = a(i,j,k,n) + rho*dafp(i,j,k)
                            else if (abs(dafm(i,j,k)) >= 2.0d0*abs(dafp(i,j,k))) then
                               ! Eq. 31
                               ar(i,j,k,n) = a(i,j,k,n) - 2.0d0*(1.0d0 - rho)*dafp(i,j,k) - rho*dafm(i,j,k)
                            else if (abs(dafp(i,j,k)) >= 2.0d0*abs(dafm(i,j,k))) then
                               ! Eq. 32
                               al(i+1,j,k,n) = a(i,j,k,n) + 2.0d0*(1.0d0 - rho)*dafm(i,j,k) + rho*dafp(i,j,k)
                            endif
                            
                         endif
                      endif

                   else
                      ! if Eqs. 24 or 25 didn't hold we still may need to limit
                      if (abs(dafm(i,j,k)) >= 2.0d0*abs(dafp(i,j,k))) then
                         ar(i,j,k,n) = a(i,j,k,n) - 2.0d0*dafp(i,j,k)
                      endif
                      if (abs(dafp(i,j,k)) >= 2.0d0*abs(dafm(i,j,k))) then
                         al(i+1,j,k,n) = a(i,j,k,n) + 2.0d0*dafm(i,j,k)
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

       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2)-2, hi(2)+3
             do i = lo(1)-1, hi(1)+1

                ! interpolate to the edges
                a_int(i,j,k) = (7.0d0/12.0d0)*(a(i,j-1,k,n) + a(i,j,k,n)) - &
                     (1.0d0/12.0d0)*(a(i,j-2,k,n) + a(i,j+1,k,n))

                al(i,j,k,n) = a_int(i,j,k)
                ar(i,j,k,n) = a_int(i,j,k)

             enddo
          enddo
       enddo

       if (limit_fourth_order == 1) then
           do k = lo(3)-dg(3), hi(3)+dg(3)
              do j = lo(2)-2, hi(2)+3
                 do i = lo(1)-1, hi(1)+1
                    ! these live on cell-centers
                    dafm(i,j,k) = a(i,j,k,n) - a_int(i,j,k)
                    dafp(i,j,k) = a_int(i,j+1,k) - a(i,j,k,n)

                    ! these live on cell-centers
                    d2af(i,j,k) = 6.0d0*(a_int(i,j,k) - 2.0d0*a(i,j,k,n) + a_int(i,j+1,k))
                 enddo
              enddo
           enddo

           do k = lo(3)-dg(3), hi(3)+dg(3)
              do j = lo(2)-3, hi(2)+3
                 do i = lo(1)-1, hi(1)+1
                    d2ac(i,j,k) = a(i,j-1,k,n) - 2.0d0*a(i,j,k,n) + a(i,j+1,k,n)
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

           ! this is a look over cell centers, affecting
           ! j-1/2,R and j+1/2,L
           do k = lo(3)-dg(3), hi(3)+dg(3)
              do j = lo(2)-1, hi(2)+1
                 do i = lo(1)-1, hi(1)+1

                    ! limit? MC Eq. 24 and 25
                    if (dafm(i,j,k) * dafp(i,j,k) <= 0.0d0 .or. &
                         (a(i,j,k,n) - a(i,j-2,k,n))*(a(i,j+2,k,n) - a(i,j,k,n)) <= 0.0d0) then

                       ! we are at an extrema

                       s = sign(1.0d0, d2ac(i,j,k))
                       if ( s == sign(1.0d0, d2ac(i,j-1,k)) .and. &
                            s == sign(1.0d0, d2ac(i,j+1,k)) .and. &
                            s == sign(1.0d0, d2af(i,j,k))) then
                          ! MC Eq. 26
                          d2a_lim = s*min(abs(d2af(i,j,k)), C2*abs(d2ac(i,j-1,k)), &
                               C2*abs(d2ac(i,j,k)), C2*abs(d2ac(i,j+1,k)))
                       else
                          d2a_lim = 0.0d0
                       endif

                       if (abs(d2af(i,j,k)) <= 1.e-12*max(abs(a(i,j-2,k,n)), abs(a(i,j-1,k,n)), &
                            abs(a(i,j,k,n)), abs(a(i,j+1,k,n)), abs(a(i,j+2,k,n)))) then
                          rho = 0.0d0
                       else
                          ! MC Eq. 27
                          rho = d2a_lim/d2af(i,j,k)
                       endif

                       if (rho < 1.0d0 - 1.d-12) then
                          ! we may need to limit -- these quantities are at cell-centers
                          d3a_min = min(d3a(i,j-1,k), d3a(i,j,k), d3a(i,j+1,k), d3a(i,j+2,k))
                          d3a_max = max(d3a(i,j-1,k), d3a(i,j,k), d3a(i,j+1,k), d3a(i,j+2,k))

                          if (C3*max(abs(d3a_min), abs(d3a_max)) <= (d3a_max - d3a_min)) then
                             ! limit
                             if (dafm(i,j,k)*dafp(i,j,k) < 0.0d0) then
                                ! Eqs. 29, 30
                                ar(i,j,k,n) = a(i,j,k,n) - rho*dafm(i,j,k)  ! note: typo in Eq 29
                                al(i,j+1,k,n) = a(i,j,k,n) + rho*dafp(i,j,k)
                             else if (abs(dafm(i,j,k)) >= 2.0d0*abs(dafp(i,j,k))) then
                                ! Eq. 31
                                ar(i,j,k,n) = a(i,j,k,n) - 2.0d0*(1.0d0 - rho)*dafp(i,j,k) - rho*dafm(i,j,k)
                             else if (abs(dafp(i,j,k)) >= 2.0*abs(dafm(i,j,k))) then
                                ! Eq. 32
                                al(i,j+1,k,n) = a(i,j,k,n) + 2.0d0*(1.0d0 - rho)*dafm(i,j,k) + rho*dafp(i,j,k)
                             endif

                          endif
                       endif

                    else
                       ! if Eqs. 24 or 25 didn't hold we still may need to limit
                       if (abs(dafm(i,j,k)) >= 2.0d0*abs(dafp(i,j,k))) then
                          ar(i,j,k,n) = a(i,j,k,n) - 2.0d0*dafp(i,j,k)
                       endif
                       if (abs(dafp(i,j,k)) >= 2.0d0*abs(dafm(i,j,k))) then
                          al(i,j+1,k,n) = a(i,j,k,n) + 2.0d0*dafm(i,j,k)
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

       do k = lo(3)-2, hi(3)+3
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

                ! interpolate to the edges
                a_int(i,j,k) = (7.0d0/12.0d0)*(a(i,j,k-1,n) + a(i,j,k,n)) - &
                     (1.0d0/12.0d0)*(a(i,j,k-2,n) + a(i,j,k+1,n))

                al(i,j,k,n) = a_int(i,j,k)
                ar(i,j,k,n) = a_int(i,j,k)

             enddo
          enddo
       enddo

       if (limit_fourth_order == 1) then

           do k = lo(3)-2, hi(3)+3
              do j = lo(2)-1, hi(2)+1
                 do i = lo(1)-1, hi(1)+1
                    ! these live on cell-centers
                    dafm(i,j,k) = a(i,j,k,n) - a_int(i,j,k)
                    dafp(i,j,k) = a_int(i,j,k+1) - a(i,j,k,n)

                    ! these live on cell-centers
                    d2af(i,j,k) = 6.0d0*(a_int(i,j,k) - 2.0d0*a(i,j,k,n) + a_int(i,j,k+1))
                 enddo
              enddo
           enddo

           do k = lo(3)-3, hi(3)+3
              do j = lo(2)-1, hi(2)+1
                 do i = lo(1)-1, hi(1)+1
                    d2ac(i,j,k) = a(i,j,k-1,n) - 2.0d0*a(i,j,k,n) + a(i,j,k+1,n)
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

           ! this is a look over cell centers, affecting
           ! j-1/2,R and j+1/2,L
           do k = lo(3)-1, hi(3)+1
              do j = lo(2)-1, hi(2)+1
                 do i = lo(1)-1, hi(1)+1

                    ! limit? MC Eq. 24 and 25
                    if (dafm(i,j,k) * dafp(i,j,k) <= 0.0d0 .or. &
                         (a(i,j,k,n) - a(i,j,k-2,n))*(a(i,j,k+2,n) - a(i,j,k,n)) <= 0.0d0) then

                       ! we are at an extrema

                       s = sign(1.0d0, d2ac(i,j,k))
                       if ( s == sign(1.0d0, d2ac(i,j,k-1)) .and. &
                            s == sign(1.0d0, d2ac(i,j,k+1)) .and. &
                            s == sign(1.0d0, d2af(i,j,k))) then
                          ! MC Eq. 26
                          d2a_lim = s*min(abs(d2af(i,j,k)), C2*abs(d2ac(i,j,k-1)), &
                               C2*abs(d2ac(i,j,k)), C2*abs(d2ac(i,j,k+1)))
                       else
                          d2a_lim = 0.0d0
                       endif

                       if (abs(d2af(i,j,k)) <= 1.e-12*max(abs(a(i,j,k-2,n)), abs(a(i,j,k-1,n)), &
                            abs(a(i,j,k,n)), abs(a(i,j,k+1,n)), abs(a(i,j,k+2,n)))) then
                          rho = 0.0d0
                       else
                          ! MC Eq. 27
                          rho = d2a_lim/d2af(i,j,k)
                       endif

                       if (rho < 1.0d0 - 1.d-12) then
                          ! we may need to limit -- these quantities are at cell-centers
                          d3a_min = min(d3a(i,j,k-1), d3a(i,j,k), d3a(i,j,k+1), d3a(i,j,k+2))
                          d3a_max = max(d3a(i,j,k-1), d3a(i,j,k), d3a(i,j,k+1), d3a(i,j,k+2))

                          if (C3*max(abs(d3a_min), abs(d3a_max)) <= (d3a_max - d3a_min)) then
                             ! limit
                             if (dafm(i,j,k)*dafp(i,j,k) < 0.0d0) then
                                ! Eqs. 29, 30
                                ar(i,j,k,n) = a(i,j,k,n) - rho*dafm(i,j,k)  ! note: typo in Eq 29
                                al(i,j,k+1,n) = a(i,j,k,n) + rho*dafp(i,j,k)
                             else if (abs(dafm(i,j,k)) >= 2.0d0*abs(dafp(i,j,k))) then
                                ! Eq. 31
                                ar(i,j,k,n) = a(i,j,k,n) - 2.0d0*(1.0d0 - rho)*dafp(i,j,k) - rho*dafm(i,j,k)
                             else if (abs(dafp(i,j,k)) >= 2.0*abs(dafm(i,j,k))) then
                                ! Eq. 32
                                al(i,j,k+1,n) = a(i,j,k,n) + 2.0d0*(1.0d0 - rho)*dafm(i,j,k) + rho*dafp(i,j,k)
                             endif

                          endif
                       endif

                    else
                       ! if Eqs. 24 or 25 didn't hold we still may need to limit
                       if (abs(dafm(i,j,k)) >= 2.0d0*abs(dafp(i,j,k))) then
                          ar(i,j,k,n) = a(i,j,k,n) - 2.0d0*dafp(i,j,k)
                       endif
                       if (abs(dafp(i,j,k)) >= 2.0d0*abs(dafm(i,j,k))) then
                          al(i,j,k+1,n) = a(i,j,k,n) + 2.0d0*dafm(i,j,k)
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


  subroutine ca_make_cell_center(lo, hi, &
                                 U, U_lo, U_hi, nc, &
                                 U_cc, U_cc_lo, U_cc_hi, nc_cc) &
                                 bind(C, name="ca_make_cell_center")

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3)
    integer, intent(in) :: U_cc_lo(3), U_cc_hi(3)
    integer, intent(in) :: nc, nc_cc
    real(rt), intent(in) :: U(U_lo(1):U_hi(1), U_lo(2):U_hi(2), U_lo(3):U_hi(3), nc)
    real(rt), intent(inout) :: U_cc(U_cc_lo(1):U_cc_hi(1), U_cc_lo(2):U_cc_hi(2), U_cc_lo(3):U_cc_hi(3), nc_cc)

    integer :: i, j, k, n
    real(rt) :: lap

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             do n = 1, nc
                lap = U(i+1,j,k,n) - TWO*U(i,j,k,n) + U(i-1,j,k,n)
#if BL_SPACEDIM >= 2
                lap = lap + U(i,j+1,k,n) - TWO*U(i,j,k,n) + U(i,j-1,k,n)
#endif
#if BL_SPACEDIM == 3
                lap = lap + U(i,j,k+1,n) - TWO*U(i,j,k,n) + U(i,j,k-1,n)
#endif

                U_cc(i,j,k,n) = U(i,j,k,n) - TWENTYFOURTH * lap

             enddo
          enddo
       enddo
    enddo

  end subroutine ca_make_cell_center

  subroutine ca_make_cell_center_in_place(lo, hi, &
                                          U, U_lo, U_hi, nc) &
                                          bind(C, name="ca_make_cell_center_in_place")

    use amrex_mempool_module, only : bl_allocate, bl_deallocate

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3)
    integer, intent(in) :: nc
    real(rt), intent(inout) :: U(U_lo(1):U_hi(1), U_lo(2):U_hi(2), U_lo(3):U_hi(3), nc)

    integer :: i, j, k, n

    real(rt), pointer :: lap(:,:,:,:)

    call bl_allocate(lap, lo, hi, nc)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             do n = 1, nc
                lap(i,j,k,n) = U(i+1,j,k,n) - TWO*U(i,j,k,n) + U(i-1,j,k,n)
#if BL_SPACEDIM >= 2
                lap(i,j,k,n) = lap(i,j,k,n) + U(i,j+1,k,n) - TWO*U(i,j,k,n) + U(i,j-1,k,n)
#endif
#if BL_SPACEDIM == 3
                lap(i,j,k,n) = lap(i,j,k,n) + U(i,j,k+1,n) - TWO*U(i,j,k,n) + U(i,j,k-1,n)
#endif
             enddo
          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             do n = 1, nc
                U(i,j,k,n) = U(i,j,k,n) - TWENTYFOURTH * lap(i,j,k,n)
             enddo
          enddo
       enddo
    enddo

    call bl_deallocate(lap)

  end subroutine ca_make_cell_center_in_place

  subroutine ca_make_fourth_average(lo, hi, &
                                    q, q_lo, q_hi, nc, &
                                    q_bar, q_bar_lo, q_bar_hi, nc_bar) &
                                    bind(C, name="ca_make_fourth_average")

    ! this takes the cell-center q and the q constructed from the
    ! cell-average U (q_bar) and replaces the cell-center q with a
    ! proper 4th-order accurate cell-average

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: q_bar_lo(3), q_bar_hi(3)
    integer, intent(in) :: nc, nc_bar
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), nc)
    real(rt), intent(in) :: q_bar(q_bar_lo(1):q_bar_hi(1), q_bar_lo(2):q_bar_hi(2), q_bar_lo(3):q_bar_hi(3), nc_bar)

    integer :: i, j, k, n
    real(rt) :: lap

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             do n = 1, nc
                lap = q_bar(i+1,j,k,n) - TWO*q_bar(i,j,k,n) + q_bar(i-1,j,k,n)
#if BL_SPACEDIM >= 2
                lap = lap + q_bar(i,j+1,k,n) - TWO*q_bar(i,j,k,n) + q_bar(i,j-1,k,n)
#endif
#if BL_SPACEDIM == 3
                lap = lap + q_bar(i,j,k+1,n) - TWO*q_bar(i,j,k,n) + q_bar(i,j,k-1,n)
#endif

                q(i,j,k,n) = q(i,j,k,n) + TWENTYFOURTH * lap

             enddo
          enddo
       enddo
    enddo

  end subroutine ca_make_fourth_average

  subroutine ca_make_fourth_in_place(lo, hi, &
                                     q, q_lo, q_hi, nc) &
                                     bind(C, name="ca_make_fourth_in_place")

    use amrex_mempool_module, only : bl_allocate, bl_deallocate

    ! this takes the cell-center q and makes it a cell-average q, in
    ! place (e.g. q is overwritten by its average).  Note: this
    ! routine is not tile safe.

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: nc
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), nc)

    integer :: i, j, k, n
    real(rt), pointer :: lap(:,:,:,:)

    call bl_allocate(lap, lo, hi, nc)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             do n = 1, nc
                lap(i,j,k,n) = q(i+1,j,k,n) - TWO*q(i,j,k,n) + q(i-1,j,k,n)
#if BL_SPACEDIM >= 2
                lap(i,j,k,n) = lap(i,j,k,n) + q(i,j+1,k,n) - TWO*q(i,j,k,n) + q(i,j-1,k,n)
#endif
#if BL_SPACEDIM == 3
                lap(i,j,k,n) = lap(i,j,k,n) + q(i,j,k+1,n) - TWO*q(i,j,k,n) + q(i,j,k-1,n)
#endif
             enddo
          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             do n = 1, nc
                q(i,j,k,n) = q(i,j,k,n) + TWENTYFOURTH * lap(i,j,k,n)

             enddo
          enddo
       enddo
    enddo

    call bl_deallocate(lap)

  end subroutine ca_make_fourth_in_place


end module fourth_order

