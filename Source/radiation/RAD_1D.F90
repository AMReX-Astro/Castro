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
  implicit none
  integer :: DIMDEC(reg)
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(1)
  real(rt)         :: dx(1)
  real(rt)         :: h1, d1
  integer :: i
  h1 = 0.5e0_rt * dx(1)
  d1 = 1.e0_rt / (3.e0_rt * dx(1))
  do i = reg_l1, reg_h1
     r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
  enddo
end subroutine sphc

subroutine sphe(r, s, n, &
                DIMS(reg), dx) bind(C, name="sphe")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(reg)
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(1)
  integer :: n
  real(rt)         :: dx(1)
  integer :: i
  do i = reg_l1, reg_h1
     r(i) = r(i)**2
  enddo
end subroutine sphe

subroutine rfface(fine, &
                  DIMS(fbox), &
                  crse, &
                  DIMS(cbox), &
                  idim, irat) bind(C, name="rfface")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(fbox)
  integer :: DIMDEC(cbox)
  real(rt)         :: fine(DIMV(fbox))
  real(rt)         :: crse(DIMV(cbox))
  integer :: idim, irat(0:0)
  fine(fbox_l1) = crse(cbox_l1)
end subroutine rfface

subroutine bextrp(f, fbox_l1, fbox_h1, reg_l1, reg_h1) bind(C, name="bextrp")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: fbox_l1, fbox_h1
  integer ::  reg_l1,  reg_h1
  real(rt)         :: f(fbox_l1:fbox_h1)
  integer :: i

  !     i direction first:
  i = reg_l1
  f(i-1) = 2.e0_rt * f(i) - f(i+1)
  i = reg_h1
  f(i+1) = 2.e0_rt * f(i) - f(i-1)
end subroutine bextrp

end module rad_module
