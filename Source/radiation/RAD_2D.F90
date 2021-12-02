#include <AMReX_LO_BCTYPES.H>
#include <AMReX_ArrayLim.H>

module rad_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , parameter :: tiny = 1.e-50_rt
  real(rt)        , parameter :: BIGKR = 1.e25_rt

contains

subroutine sphc(r, s, &
                DIMS(reg), dx) bind(C, name="sphc")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(reg_l2:reg_h2)
  real(rt)         :: dx(2)
  real(rt)         :: h1, h2, d1, d2
  integer :: i, j
  h1 = 0.5e0_rt * dx(1)
  h2 = 0.5e0_rt * dx(2)
  d1 = 1.e0_rt / (3.e0_rt * dx(1))
  d2 = 1.e0_rt / dx(2)
  do i = reg_l1, reg_h1
     r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
  enddo
  do j = reg_l2, reg_h2
     s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
  enddo
end subroutine sphc

subroutine sphe(r, s, n, &
                DIMS(reg), dx) bind(C, name="sphe")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(reg_l2:reg_h2)
  integer :: n
  real(rt)         :: dx(2)
  real(rt)         :: h1, h2, d1, d2
  integer :: i, j
  if (n == 0) then
     do i = reg_l1, reg_h1
        r(i) = r(i)**2
     enddo
     h2 = 0.5e0_rt * dx(2)
     d2 = 1.e0_rt / dx(2)
     do j = reg_l2, reg_h2
        s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
     enddo
  else
     h1 = 0.5e0_rt * dx(1)
     d1 = 1.e0_rt / (3.e0_rt * dx(1))
     do i = reg_l1, reg_h1
        r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
     enddo
     do j = reg_l2, reg_h2
        s(j) = sin(s(j))
     enddo
  endif
end subroutine sphe

subroutine rfface(fine, &
                  DIMS(fbox), &
                  crse, &
                  DIMS(cbox), &
                  idim, irat) bind(C, name="rfface")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(fbox)
  integer :: DIMDEC(cbox)
  real(rt)         :: fine(DIMV(fbox))
  real(rt)         :: crse(DIMV(cbox))
  integer :: idim, irat(0:1)
  integer :: i, j
  if (idim == 0) then
     do j = fbox_l2, fbox_h2
        fine(fbox_l1,j) = crse(cbox_l1, j/irat(1)) / irat(1)
     enddo
  else
     do i = fbox_l1, fbox_h1
        fine(i,fbox_l2) = crse(i/irat(0), cbox_l2) / irat(0)
     enddo
  endif
end subroutine rfface


subroutine bextrp(f, fbox_l1, fbox_l2, fbox_h1, fbox_h2, &
                  reg_l1, reg_l2, reg_h1, reg_h2) bind(C, name="bextrp")

  use amrex_fort_module, only : rt => amrex_real
  integer :: fbox_l1, fbox_l2, fbox_h1, fbox_h2
  integer ::  reg_l1,  reg_l2,  reg_h1,  reg_h2
  real(rt)         :: f(fbox_l1:fbox_h1,fbox_l2:fbox_h2)
  integer :: i, j

  !     i direction first:
  do j = reg_l2, reg_h2
     i = reg_l1
     f(i-1,j) = 2.e0_rt * f(i,j) - f(i+1,j)
     i = reg_h1
     f(i+1,j) = 2.e0_rt * f(i,j) - f(i-1,j)
  enddo

  !     j direction second, including corners:
  do i = reg_l1 - 1, reg_h1 + 1
     j = reg_l2
     f(i,j-1) = 2.e0_rt * f(i,j) - f(i,j+1)
     j = reg_h2
     f(i,j+1) = 2.e0_rt * f(i,j) - f(i,j-1)
  enddo

  !  corner results are the same whichever direction we extrapolate first
end subroutine bextrp

end module rad_module
