#include <AMReX_LO_BCTYPES.H>
#include <AMReX_ArrayLim.H>

module rad_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , parameter :: tiny = 1.e-50_rt
  real(rt)        , parameter :: BIGKR = 1.e25_rt

contains

! The following routines implement metric terms in 2D and are included
! in the 3D source only to enable the code to link.

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
  integer :: idim, irat(0:2)
  integer :: i, j, k
  integer :: rfac
  if (idim == 0) then
     rfac = irat(1) * irat(2)
     do k = fbox_l3, fbox_h3
        do j = fbox_l2, fbox_h2
           fine(fbox_l1,j,k) = crse(cbox_l1, j/irat(1), k/irat(2)) / rfac
        enddo
     enddo
  else if (idim == 1) then
     rfac = irat(0) * irat(2)
     do k = fbox_l3, fbox_h3
        do i = fbox_l1, fbox_h1
           fine(i,fbox_l2,k) = crse(i/irat(0), cbox_l2, k/irat(2)) / rfac
        enddo
     enddo
  else
     rfac = irat(0) * irat(1)
     do j = fbox_l2, fbox_h2
        do i = fbox_l1, fbox_h1
           fine(i,j,fbox_l3) = crse(i/irat(0), j/irat(1), cbox_l3) / rfac
        enddo
     enddo
  endif
end subroutine rfface

end module rad_module

