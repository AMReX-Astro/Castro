module mhd_advection_module

  implicit none

contains

subroutine check_for_mhd_cfl_violation(lo, hi, &
                                       q, q_lo, q_hi, &
                                       cx, cx_lo, cx_hi, &
                                       cy, cy_lo, cy_hi, &
                                       cz, cz_lo, cz_hi, &
                                       courno, dx, dt) bind(C, name="check_for_mhd_cfl_violation")

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use amrex_error_module
  use meth_params_module, only : NQ, QRHO, QU, QV, QW, QC, NQAUX,&
       QREINT, QPRES, QFA, QFS, QTEMP, &
       QMAGX,  QMAGY, QMAGZ

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: cx_lo(3), cx_hi(3)
  integer, intent(in) :: cy_lo(3), cy_hi(3)
  integer, intent(in) :: cz_lo(3), cz_hi(3)

  real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
  real(rt), intent(in) :: cx(cx_lo(1):cx_hi(1), cx_lo(2):cx_hi(2), cx_lo(3):cx_hi(3))
  real(rt), intent(in) :: cy(cy_lo(1):cy_hi(1), cy_lo(2):cy_hi(2), cy_lo(3):cx_hi(3))
  real(rt), intent(in) :: cz(cz_lo(1):cz_hi(1), cz_lo(2):cz_hi(2), cz_lo(3):cx_hi(3))

  real(rt), intent(in) :: dx(AMREX_SPACEDIM)
  real(rt), intent(out) :: courno
  real(rt), intent(in), value :: dt

  integer :: i, j, k

  real(rt) :: courx, coury, courz, courmx, courmy, courmz
  real(rt) :: dtdx, dtdy, dtdz

  ! Compute running max of Courant number over grids
  courmx = courno
  courmy = courno
  courmz = courno

  dtdx = dt / dx(1)
  dtdy = dt / dx(2)
  dtdz = dt / dx(3)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           courx = ( cx(i,j,k)+abs(q(i,j,k,QU)) ) * dtdx
           coury = ( cy(i,j,k)+abs(q(i,j,k,QV)) ) * dtdy
           courz = ( cz(i,j,k)+abs(q(i,j,k,QW)) ) * dtdz

           courmx = max( courmx, courx )
           courmy = max( courmy, coury )
           courmz = max( courmz, courz )

           if (courx .gt. ONE) then
              !
              ! A critical region since we usually can't write from threads.
              !
              print *,'   '
              print *,'>>> ... (u+c) * a * dt / dx > 1 ', courx
              print *,'>>> ... at cell (i,j,k)   : ',i,j,k
              print *,'>>> ... u, c                ',q(i,j,k,QU), cx(i,j,k)
              print *,'>>> ... B                   ',q(i,j,k,QMAGX:QMAGZ)
              print *,'>>> ... Internal e          ',q(i,j,k,QREINT)
              print *,'>>> ... pressure            ',q(i,j,k,QPRES)
              print *,'>>> ... density             ',q(i,j,k,QRHO)
              call amrex_error("Error:: Castro_advection_mhd_3d.f90 :: CFL violation in x-dir in ctoprim")
           end if

           if (coury .gt. ONE) then
              !
              ! A critical region since we usually can't write from threads.
              !
              print *,'   '
              print *,'>>> ... (v+c) * a * dt / dx > 1 ', coury
              print *,'>>> ... at cell (i,j,k)   : ',i,j,k
              print *,'>>> ... v, c                ',q(i,j,k,QV), cy(i,j,k)
              print *,'>>> ... B                   ',q(i,j,k,QMAGX:QMAGZ)
              print *,'>>> ... pressure            ',q(i,j,k,QPRES)
              print *,'>>> ... density             ',q(i,j,k,QRHO)
              call amrex_error("Error:: Castro_advection_mhd_3d.f90 :: CFL violation in y-dir in ctoprim")
           end if

           if (courz .gt. ONE) then
              !
              ! A critical region since we usually can't write from threads.
              !
              print *,'   '
              print *,'>>> ... (w+c) * a * dt / dx > 1 ', courz
              print *,'>>> ... at cell (i,j,k)   : ',i,j,k
              print *,'>>> ... w, c                ',q(i,j,k,QW), cz(i,j,k)
              print *,'>>> ... B                   ',q(i,j,k,QMAGX:QMAGZ)
              print *,'>>> ... pressure            ',q(i,j,k,QPRES)
              print *,'>>> ... density             ',q(i,j,k,QRHO)
              call amrex_error("Error:: Castro_advection_mhd_3d.f90 :: CFL violation in z-dir in ctoprim")
           end if

        enddo
     enddo
  enddo
  courno = max( courmx, courmy, courmz )

end subroutine check_for_mhd_cfl_violation


end module mhd_advection_module
