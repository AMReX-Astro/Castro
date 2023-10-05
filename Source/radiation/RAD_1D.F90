#include <AMReX_LO_BCTYPES.H>
#include <AMReX_ArrayLim.H>

module rad_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , parameter :: tiny = 1.e-50_rt
  real(rt)        , parameter :: BIGKR = 1.e25_rt

contains

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
