#include "AMReX_LO_BCTYPES.H"
#include "AMReX_ArrayLim.H"

module rad_module

  use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS, UFX, NVAR

  use rad_util_module, only : FLDlambda

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , parameter :: tiny = 1.e-50_rt
  real(rt)        , parameter :: BIGKR = 1.e25_rt

contains

! The following routines implement metric terms in 2D and are included
! in the 3D source only to enable the code to link.

subroutine multrs(d, &
                  DIMS(dbox), &
                  DIMS(reg), &
                  r, s) bind(C, name="multrs")
  
  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(dbox)
  integer :: DIMDEC(reg)
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(reg_l2:reg_h2)
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           d(i,j,k) = d(i,j,k) * r(i) * s(j)
        enddo
     enddo
  enddo
end subroutine multrs

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

subroutine lacoef(a, &
                  DIMS(abox), &
                  DIMS(reg), &
                  fkp, eta, etainv, r, s, c, dt, theta) bind(C, name="lacoef")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  real(rt)         :: a(DIMV(abox))
  real(rt)         :: fkp(DIMV(abox))
  real(rt)         :: eta(DIMV(abox))
  real(rt)         :: etainv(DIMV(abox))
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(reg_l2:reg_h2)
  real(rt)         :: c, dt, theta
  integer :: i, j, k
  real(rt)         :: dtm
  dtm = 1.e0_rt / dt
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           a(i,j,k) = r(i) * s(j) * &
                (fkp(i,j,k) * etainv(i,j,k) * c + dtm) / &
                (1.e0_rt - (1.e0_rt - theta) * eta(i,j,k))
        enddo
     enddo
  enddo
end subroutine lacoef


subroutine bclim(b, &
                 lambda, DIMS(bbox), &
                 DIMS(reg), &
                 n, kappar, DIMS(kbox), &
                 r, s, c, dx) bind(C, name="bclim")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: lambda(DIMV(bbox))
  real(rt)         :: kappar(DIMV(kbox))
  real(rt)         :: r(reg_l1:reg_h1+1)
  real(rt)         :: s(reg_l2:reg_h2+1)
  real(rt)         :: c, dx(3)
  real(rt)         :: kavg
  integer :: i, j, k
  real(rt)         :: kap
  if (n == 0) then
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1 + 1
              kap = kavg(kappar(i-1,j,k), kappar(i,j,k), dx(1), -1)
              b(i,j,k) = r(i) * s(j) * c * lambda(i,j,k) / kap
           enddo
        enddo
     enddo
  else if (n == 1) then
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2 + 1
           do i = reg_l1, reg_h1
              kap = kavg(kappar(i,j-1,k), kappar(i,j,k), dx(2), -1)
              b(i,j,k) = r(i) * s(j) * c * lambda(i,j,k) / kap
           enddo
        enddo
     enddo
  else
     do k = reg_l3, reg_h3 + 1
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              kap = kavg(kappar(i,j,k-1), kappar(i,j,k), dx(3), -1)
              b(i,j,k) = r(i) * s(j) * c * lambda(i,j,k) / kap
           enddo
        enddo
     enddo
  endif
end subroutine bclim

subroutine flxlim(lambda, &
                  DIMS(rbox), & 
                  DIMS(reg), limiter) bind(C, name="flxlim")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: limiter
  real(rt)         :: lambda(DIMV(rbox))
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           lambda(i,j,k) = FLDlambda(lambda(i,j,k),limiter)
        enddo
     enddo
  enddo
end subroutine flxlim

subroutine eddfac(efact, &
                  DIMS(rbox), &
                  DIMS(reg), limiter, n) bind(C, name="eddfac")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: n, limiter
  real(rt)         :: efact(DIMV(rbox))
  integer :: i, j, k
  real(rt)         :: r, lambda
  integer :: dir(0:2)
  dir(0) = 0
  dir(1) = 0
  dir(2) = 0
  dir(n) = 1
  do k = reg_l3, reg_h3 + dir(2)
     do j = reg_l2, reg_h2 + dir(1)
        do i = reg_l1, reg_h1 + dir(0)
           r = efact(i,j,k)
           lambda = FLDlambda(r,limiter)
           efact(i,j,k) = lambda + (lambda * r)**2
        enddo
     enddo
  enddo
end subroutine eddfac

subroutine scgrd1(r, &
                  DIMS(rbox), &
                  DIMS(reg), &
                  n, kappar, DIMS(kbox), &
                  er, dx) bind(C, name="scgrd1")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real(rt)         :: r(DIMV(rbox))
  real(rt)         :: kappar(DIMV(kbox))
  real(rt)         :: er(DIMV(kbox))
  real(rt)         :: dx(3)
  real(rt)         :: kavg
  integer :: i, j, k
  real(rt)         :: kap
  if (n == 0) then
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           !     x derivatives, gradient assembly:
           do i = reg_l1, reg_h1 + 1
              r(i,j,k) = abs(er(i,j,k) - er(i-1,j,k)) / dx(1)
           enddo
           i = reg_l1
           if (er(i-1,j,k) == -1.e0_rt) then
              r(i,j,k) = abs(er(i+1,j,k) - er(i,j,k)) / dx(1)
           endif
           i = reg_h1 + 1
           if (er(i,j,k) == -1.e0_rt) then
              r(i,j,k) = abs(er(i-1,j,k) - er(i-2,j,k)) / dx(1)
           endif
           !     construct coefficients:
           do i = reg_l1, reg_h1 + 1
              kap = kavg(kappar(i-1,j,k), kappar(i,j,k), dx(1), -1)
              r(i,j,k) = r(i,j,k) / &
                   (kap * max(er(i-1,j,k), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  else if (n == 1) then
     do k = reg_l3, reg_h3
        !     y derivatives, gradient assembly:
        do j = reg_l2, reg_h2 + 1
           do i = reg_l1, reg_h1
              r(i,j,k) = abs(er(i,j,k) - er(i,j-1,k)) / dx(2)
           enddo
        enddo
        do i = reg_l1, reg_h1
           j = reg_l2
           if (er(i,j-1,k) == -1.e0_rt) then
              r(i,j,k) = abs(er(i,j+1,k) - er(i,j,k)) / dx(2)
           endif
           j = reg_h2 + 1
           if (er(i,j,k) == -1.e0_rt) then
              r(i,j,k) = abs(er(i,j-1,k) - er(i,j-2,k)) / dx(2)
           endif
        enddo
        ! construct coefficients:
        do j = reg_l2, reg_h2 + 1
           do i = reg_l1, reg_h1
              kap = kavg(kappar(i,j-1,k), kappar(i,j,k), dx(2), -1)
              r(i,j,k) = r(i,j,k) / &
                   (kap * max(er(i,j-1,k), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  else
     !     z derivatives, gradient assembly:
     do k = reg_l3, reg_h3 +1
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              r(i,j,k) = abs(er(i,j,k-1) - er(i,j,k)) / dx(3)
           enddo
        enddo
     enddo
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           k = reg_l3
           if (er(i,j,k-1) == -1.e0_rt) then
              r(i,j,k) = abs(er(i,j,k+1) - er(i,j,k)) / dx(3)
           endif
           k = reg_h3 + 1
           if (er(i,j,k) == -1.e0_rt) then
              r(i,j,k) = abs(er(i,j,k-2) - er(i,j,k-1)) / dx(3)
           endif
        enddo
     enddo
     ! construct coefficients:
     do k = reg_l3, reg_h3 + 1
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              kap = kavg(kappar(i,j,k-1), kappar(i,j,k), dx(3), -1)
              r(i,j,k) = r(i,j,k) / &
                   (kap * max(er(i,j,k-1), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  endif
end subroutine scgrd1

subroutine scgrd2(r, &
                  DIMS(rbox), &
                  DIMS(reg), &
                  n, kappar, DIMS(kbox), er, &
                  DIMS(dbox), da, db, dx) bind(C, name="scgrd2")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: DIMDEC(dbox)
  integer :: n
  real(rt)         :: r(DIMV(rbox))
  real(rt)         :: kappar(DIMV(kbox))
  real(rt)         :: er(DIMV(kbox))
  real(rt)         :: da(DIMV(dbox))
  real(rt)         :: db(DIMV(dbox))
  real(rt)         :: dx(3)
  real(rt)         :: kavg
  integer :: i, j, k
  real(rt)         :: kap
  if (n == 0) then
     !     y & z derivatives:
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1 - 1, reg_h1 + 1
              da(i,j,k) = er(i,j+1,k  ) - er(i,j-1,k  )
              db(i,j,k) = er(i,j  ,k+1) - er(i,j  ,k-1)
           enddo
        enddo
     enddo

     !     check y derivatives
     do k = reg_l3, reg_h3
        do i = reg_l1 - 1, reg_h1 + 1
           j = reg_l2
           if (er(i,j-1,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i,j+1,k) - er(i,j,k))
           endif
           j = reg_h2
           if (er(i,j+1,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j-1,k))
           endif
        enddo
        do j = reg_l2, reg_h2
           i = reg_l1 - 1
           !     (check at j-1 and j+1.  if value at j is bad it will not be used at all.)
           if (er(i,j-1,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i,j+1,k) - er(i,j,k))
           else if (er(i,j+1,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j-1,k))
           endif
           i = reg_h1 + 1
           if (er(i,j-1,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i,j+1,k) - er(i,j,k))
           else if (er(i,j+1,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j-1,k))
           endif
        enddo
     enddo

     !     check z derivatives
     do j = reg_l2, reg_h2
        do i = reg_l1-1, reg_h1 + 1
           k = reg_l3
           if (er(i,j,k-1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k+1) - er(i,j,k))
           endif
           k = reg_h3
           if (er(i,j,k+1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j,k-1))
           endif
        enddo
     enddo
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           i = reg_l1 - 1
           !     (check at k-1 and k+1.  if value at k is bad it will not be used at all.)
           if (er(i,j,k-1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k+1) - er(i,j,k))
           else if (er(i,j,k+1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j,k-1))
           endif
           i = reg_h1 + 1
           if (er(i,j,k-1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k+1) - er(i,j,k))
           else if (er(i,j,k+1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j,k-1))
           endif
        enddo
     enddo
     !     x derivatives
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1 + 1
              r(i,j,k) = ((er(i,j,k) - er(i-1,j,k)) / dx(1)) ** 2 &
                   + ((da(i-1,j,k) + da(i,j,k)) / &
                   (4.e0_rt * dx(2))) ** 2 &
                   + ((db(i-1,j,k) + db(i,j,k)) / &
                   (4.e0_rt * dx(3))) ** 2
           enddo
           i = reg_l1
           if (er(i-1,j,k) == -1.e0_rt) then
              r(i,j,k) = ((er(i+1,j,k) - er(i,j,k)) / dx(1)) ** 2 &
                   + (da(i,j,k) / (2.e0_rt * dx(2))) ** 2 &
                   + (db(i,j,k) / (2.e0_rt * dx(3))) ** 2
           endif
           i = reg_h1 + 1
           if (er(i,j,k) == -1.e0_rt) then
              r(i,j,k) = ((er(i-1,j,k) - er(i-2,j,k)) / dx(1)) ** 2 &
                   + (da(i-1,j,k) / (2.e0_rt * dx(2))) ** 2 &
                   + (db(i-1,j,k) / (2.e0_rt * dx(3))) ** 2
           endif
           !     construct scaled gradient
           do i = reg_l1, reg_h1 + 1
              kap = kavg(kappar(i-1,j,k), kappar(i,j,k), dx(1), -1)
              r(i,j,k) = sqrt(r(i,j,k)) / &
                   (kap * max(er(i-1,j,k), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  else if (n == 1) then
     !     x & z derivatives:
     do k = reg_l3, reg_h3
        do j = reg_l2-1, reg_h2+1
           do i = reg_l1, reg_h1
              da(i,j,k) = er(i+1,j,k) - er(i-1,j,k)
              db(i,j,k) = er(i,j,k+1) - er(i,j,k-1)
           enddo
        enddo
     enddo

     !     check x derivatives
     do k = reg_l3, reg_h3
        do j = reg_l2-1, reg_h2+1
           i = reg_l1
           if (er(i-1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i+1,j,k) - er(i,j,k))
           endif
           i = reg_h1
           if (er(i+1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i-1,j,k))
           endif
        enddo
        do i = reg_l1, reg_h1
           j = reg_l2-1
           !     (check at i-1 and i+1.  if value at i is bad it will not be used at all.)
           if (er(i-1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i+1,j,k) - er(i,j,k))
           else if (er(i+1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i-1,j,k))
           endif
           j = reg_h2+1
           if (er(i-1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i+1,j,k) - er(i,j,k))
           else if (er(i+1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i-1,j,k))
           endif
        enddo
     enddo

     !     check z derivatives
     do j = reg_l2-1, reg_h2+1
        do i = reg_l1, reg_h1
           k = reg_l3
           if (er(i,j,k-1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k+1) - er(i,j,k))
           endif
           k = reg_h3
           if (er(i,j,k+1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j,k-1))
           endif
        enddo
     enddo
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           j = reg_l2-1
           !     (check at k-1 and k+1.  if value at k is bad it will not be used at all.)
           if (er(i,j,k-1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k+1) - er(i,j,k))
           else if (er(i,j,k+1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j,k-1))
           endif
           j = reg_h2+1
           if (er(i,j,k-1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k+1) - er(i,j,k))
           else if (er(i,j,k+1) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j,k-1))
           endif
        enddo
     enddo

     !     y derivatives
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2 + 1
           do i = reg_l1, reg_h1
              r(i,j,k) = ((er(i,j,k) - er(i,j-1,k)) / dx(2)) ** 2 &
                   + ((da(i,j-1,k) + da(i,j,k)) / &
                   (4.e0_rt * dx(1))) ** 2 &
                   + ((db(i,j-1,k) + db(i,j,k)) / &
                   (4.e0_rt * dx(3))) ** 2
           enddo
        enddo
        do i = reg_l1, reg_h1
           j = reg_l2
           if (er(i,j-1,k) == -1.e0_rt) then
              r(i,j,k) = ((er(i,j+1,k) - er(i,j,k)) / dx(2)) ** 2 &
                   + (da(i,j,k) / (2.e0_rt * dx(1))) ** 2 &
                   + (db(i,j,k) / (2.e0_rt * dx(3))) ** 2
           endif
           j = reg_h2 + 1
           if (er(i,j,k) == -1.e0_rt) then
              r(i,j,k) = ((er(i,j-1,k) - er(i,j-2,k)) / dx(2)) ** 2 &
                   + (da(i,j-1,k) / (2.e0_rt * dx(1))) ** 2 &
                   + (db(i,j-1,k) / (2.e0_rt * dx(3))) ** 2
           endif
        enddo
        !     construct scaled gradient
        do j = reg_l2, reg_h2 + 1
           do i = reg_l1, reg_h1
              kap = kavg(kappar(i,j-1,k), kappar(i,j,k), dx(2), -1)
              r(i,j,k) = sqrt(r(i,j,k)) / &
                   (kap * max(er(i,j-1,k), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  else
     !     x & y derivatives
     do k = reg_l3-1, reg_h3+1
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              da(i,j,k) = er(i+1,j,k) - er(i-1,j,k)
              db(i,j,k) = er(i,j+1,k) - er(i,j-1,k)
           enddo
        enddo
     enddo

     !     check x derivatives
     do k = reg_l3-1, reg_h3+1
        do j = reg_l2, reg_h2
           i = reg_l1
           if (er(i-1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i+1,j,k) - er(i,j,k))
           endif
           i = reg_h1
           if (er(i+1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.0e0_rt * (er(i,j,k) - er(i-1,j,k))
           endif
        enddo
     enddo
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           k = reg_l3 - 1
           !     (check at i-1 and i+1.  if value at i is bad it will not be used at all.)
           if (er(i-1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i+1,j,k) - er(i,j,k))
           else if (er(i+1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.0e0_rt * (er(i,j,k) - er(i-1,j,k))
           endif
           k = reg_h3 + 1
           if (er(i-1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.e0_rt * (er(i+1,j,k) - er(i,j,k))
           else if (er(i+1,j,k) == -1.e0_rt) then
              da(i,j,k) = 2.0e0_rt * (er(i,j,k) - er(i-1,j,k))
           endif
        enddo
     enddo

     !     check y derivative
     do k = reg_l3-1, reg_h3+1
        do i = reg_l1, reg_h1
           j = reg_l2
           if (er(i,j-1,k) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j+1,k) - er(i,j,k))
           endif
           j = reg_h2
           if (er(i,j+1,k) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j-1,k))
           endif
        enddo
     enddo
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           k = reg_l3 - 1
           !     (check at j-1 and j+1.  if value at j is bad it will not be used at all.)
           if (er(i,j-1,k) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j+1,k) - er(i,j,k))
           else if (er(i,j+1,k) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j-1,k))
           endif
           k = reg_h3 + 1
           if (er(i,j-1,k) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j+1,k) - er(i,j,k))
           else if (er(i,j+1,k) == -1.e0_rt) then
              db(i,j,k) = 2.e0_rt * (er(i,j,k) - er(i,j-1,k))
           endif
        enddo
     enddo

     ! z derivatives
     do k = reg_l3, reg_h3+1
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              r(i,j,k) = ((er(i,j,k) - er(i,j,k-1)) / dx(3)) ** 2 &
                   + ((da(i,j,k) + da(i,j,k-1)) / &
                   (4.e0_rt * dx(1))) ** 2 &
                   + ((db(i,j,k) + db(i,j,k-1)) / &
                   (4.e0_rt * dx(2))) ** 2
           enddo
        enddo
     enddo

     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           k = reg_l3
           if (er(i,j,k-1) == -1.e0_rt) then
              r(i,j,k) = ((er(i,j,k+1) - er(i,j,k)) / dx(3)) ** 2 &
                   + (da(i,j,k) / (2.e0_rt * dx(1))) ** 2 &
                   + (db(i,j,k) / (2.e0_rt * dx(2))) ** 2
           endif
           k = reg_h3+1
           if (er(i,j,k) == -1.e0_rt) then
              r(i,j,k) = ((er(i,j,k-1) - er(i,j,k-2)) / dx(3)) ** 2 &
                   + (da(i,j,k-1) / (2.e0_rt * dx(1))) ** 2 &
                   + (db(i,j,k-1) / (2.e0_rt * dx(2))) ** 2
           endif
        enddo
     enddo

     !     gradient assembly
     do k = reg_l3, reg_h3+1
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              kap = kavg(kappar(i,j,k), kappar(i,j,k-1), dx(3), -1)
              r(i,j,k) = sqrt(r(i,j,k)) / &
                   (kap * max(er(i,j,k-1), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  endif
end subroutine scgrd2

subroutine scgrd3(r, &
                  DIMS(rbox), &
                  DIMS(reg), &
                  n, kappar, DIMS(kbox), er, &
                  DIMS(dbox), da, db, dx) bind(C, name="scgrd3")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: DIMDEC(dbox)
  integer :: n
  real(rt)         :: r(DIMV(rbox))
  real(rt)         :: kappar(DIMV(kbox))
  real(rt)         :: er(DIMV(kbox))
  real(rt)         :: da(DIMV(dbox))
  real(rt)         :: db(DIMV(dbox))
  real(rt)         :: dx(3)
  real(rt)         :: kavg
  integer :: i
  real(rt)         :: kap
  print *, "scgrd3 not implemented in 3d"
  stop
end subroutine scgrd3

subroutine lrhs(rhs, &
                DIMS(rbox), &
                DIMS(reg), &
                temp, fkp, eta, etainv, frhoem, frhoes, dfo, &
                ero, DIMS(ebox), edot, &
                r, s, dt, sigma, c, theta) bind(C, name="lrhs")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(reg)
  real(rt)         :: rhs(DIMV(rbox))
  real(rt)         :: temp(DIMV(rbox))
  real(rt)         :: fkp(DIMV(rbox))
  real(rt)         :: eta(DIMV(rbox))
  real(rt)         :: etainv(DIMV(rbox))
  real(rt)         :: frhoem(DIMV(rbox))
  real(rt)         :: frhoes(DIMV(rbox))
  real(rt)         :: dfo(DIMV(rbox))
  real(rt)         :: ero(DIMV(ebox))
  real(rt)         :: edot(DIMV(rbox))
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(reg_l2:reg_h2)
  real(rt)         :: dt, sigma, c, theta
  integer :: i, j, k
  real(rt)         :: dtm, ek, bs, es, ekt
  dtm = 1.e0_rt / dt
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           ek = fkp(i,j,k) * eta(i,j,k)
           bs = etainv(i,j,k) * &
                &               4.e0_rt * sigma * fkp(i,j,k) * temp(i,j,k)**4
           es = eta(i,j,k) * (frhoem(i,j,k) - frhoes(i,j,k))
           ekt = (1.e0_rt - theta) * eta(i,j,k)
           rhs(i,j,k) = (rhs(i,j,k) + r(i) * s(j) * &
                (bs + dtm * (ero(i,j,k) + es) + &
                ek * c * edot(i,j,k) - &
                ekt * dfo(i,j,k))) / &
                (1.e0_rt - ekt)
        enddo
     enddo
  enddo
end subroutine lrhs

subroutine anatw2(test, &
                  DIMS(reg), &
                  temp, p, xf, Tc, dx, xlo, lo) bind(C, name="anatw2")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  real(rt)         :: test(DIMV(reg), 0:1)
  real(rt)         :: temp(DIMV(reg))
  real(rt)         :: p, xf, Tc, dx(3), xlo(3)
  integer :: lo(3)
  integer :: i, j, k
  real(rt)         :: x, y, z, r2
  do k = reg_l3, reg_h3
     z = xlo(3) + dx(3) * ((k-lo(3)) + 0.5e0_rt)
     do j = reg_l2, reg_h2
        y = xlo(2) + dx(2) * ((j-lo(2)) + 0.5e0_rt)
        do i = reg_l1, reg_h1
           x  = xlo(1) + dx(1) * ((i-lo(1)) + 0.5e0_rt)
           r2 = x*x + y*y + z*z
           test(i,j,k,0) = Tc * max((1.e0_rt-r2/xf**2), 0.e0_rt)**(1.e0_rt/p)
           test(i,j,k,1) = temp(i,j,k) - test(i,j,k,0)
        enddo
     enddo
  enddo
end subroutine anatw2

! temp contains frhoe on input:

subroutine gtemp(DIMS(reg), &
                 temp, DIMS(tb), &
                 const, em, en, &
                 state, DIMS(sb)) bind(C, name="gtemp")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(sb)
  real(rt)         :: temp(DIMV(tb))
  real(rt)         :: const(0:1), em(0:1), en(0:1)
  real(rt)         :: state(DIMV(sb), NVAR)
  real(rt)         :: alpha, teff, ex, frhoal
  integer :: i, j, k
  if (en(0) >= 1.e0_rt) then
     print *, "Bad exponent for cv calculation"
     stop
  endif
  ex = 1.e0_rt / (1.e0_rt - en(0))
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (em(0) == 0.e0_rt) then
              alpha = const(0)
           else
              alpha = const(0) * state(i,j,k, URHO) ** em(0)
           endif
           frhoal = state(i,j,k, URHO) * alpha + tiny
           if (en(0) == 0.e0_rt) then
              temp(i,j,k) = temp(i,j,k) / frhoal
           else
              teff = max(temp(i,j,k), tiny)
              temp(i,j,k) = ((1.e0_rt - en(0)) * teff / frhoal) ** ex
           endif
        enddo
     enddo
  enddo
end subroutine gtemp

! temp contains temp on input:

subroutine gcv(DIMS(reg), &
               cv,DIMS(cbox), &
               temp, DIMS(tbox), &
               const, em, en, tf, &
               state, DIMS(sbox)) bind(C, name="gcv")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(cbox)
  integer :: DIMDEC(tbox)
  integer :: DIMDEC(sbox)
  real(rt)         :: cv(DIMV(cbox))
  real(rt)         :: temp(DIMV(tbox))
  real(rt)         :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real(rt)         :: state(DIMV(sbox), NVAR)
  real(rt)         :: alpha, teff, frhoal
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (em(0) == 0.e0_rt) then
              alpha = const(0)
           else
              alpha = const(0) * state(i,j,k, URHO) ** em(0)
           endif
           frhoal = state(i,j,k, URHO) * alpha + tiny
           if (en(0) == 0.e0_rt) then
              cv(i,j,k) = alpha
           else
              teff = max(temp(i,j,k), tiny)
              teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
              cv(i,j,k) = alpha * teff ** (-en(0))
           endif
        enddo
     enddo
  enddo
end subroutine gcv

! exch contains temp on input:

subroutine cexch(DIMS(reg), &
                 exch, DIMS(xbox), &
                 er  , DIMS(ebox), &
                 fkp , DIMS(kbox), &
                 sigma, c) bind(C, name="cexch")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(xbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(kbox)
  real(rt)         :: exch(DIMV(xbox))
  real(rt)         :: er  (DIMV(ebox))
  real(rt)         :: fkp (DIMV(kbox))
  real(rt)         :: sigma, c
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           exch(i,j,k) = fkp(i,j,k) * &
                (4.e0_rt * sigma * exch(i,j,k)**4 &
                - c * er(i,j,k))
        enddo
     enddo
  enddo
end subroutine cexch

subroutine ceta2(DIMS(reg), &
                 eta, etainv, DIMS(etab), &
                 frho, DIMS(sb), &
                 temp, DIMS(tb), &
                 cv, DIMS(cb), &
                 fkp, DIMS(fb), &
                 er, DIMS(ebox), &
                 dtemp, dtime, sigma, c, underr, lagpla) bind(C, name="ceta2")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(etab)
  integer :: DIMDEC(sb)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(cb)
  integer :: DIMDEC(fb)
  integer :: DIMDEC(ebox)
  real(rt)         :: eta(DIMV(etab))
  real(rt)         :: etainv(DIMV(etab))
  real(rt)         :: frho(DIMV(sb))
  real(rt)         :: temp(DIMV(tb))
  real(rt)         :: cv(DIMV(cb))
  real(rt)         :: fkp(DIMV(fb))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: dtemp, dtime, sigma, c, underr
  integer :: lagpla
  real(rt)         :: d, frc, fac0, fac1, fac2
  integer :: i, j, k
  fac1 = 16.e0_rt * sigma * dtime
  if (lagpla == 0) then
     fac0 = 0.25e0_rt * fac1 / dtemp
     fac2 = dtime * c / dtemp
  endif
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (lagpla /= 0) then
              ! assume eta and fkp are the same
              d = fac1 * fkp(i,j,k) * temp(i,j,k) ** 3
           else
              d = fac0 * (eta(i,j,k) * (temp(i,j,k) + dtemp) ** 4 - &
                   fkp(i,j,k) * (temp(i,j,k)        ) ** 4) - &
                   fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
              ! alternate form, sometimes worse, sometimes better:
              !                  d = fac1 * fkp(i,j,k) * temp(i,j,k) ** 3 +
              !     @                fac0 * (eta(i,j,k) - fkp(i,j,k)) * temp(i,j,k) ** 4 -
              !     @                fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
              ! another alternate form (much worse):
              !                  d = fac1 * fkp(i,j,k) * (temp(i,j,k) + dtemp) ** 3 +
              !     @                fac0 * (eta(i,j,k) - fkp(i,j,k))
              !     @                     * (temp(i,j,k) + dtemp) ** 4 -
              !     @                fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
           endif
           frc = frho(i,j,k) * cv(i,j,k) + tiny
           eta(i,j,k) = d / (d + frc)
           etainv(i,j,k) = underr * frc / (d + frc)
           eta(i,j,k) = 1.e0_rt - etainv(i,j,k)
           !               eta(i,j,k) = 1.e0_rt - underr * (1.e0_rt - eta(i,j,k))
        enddo
     enddo
  enddo
end subroutine ceta2

subroutine ceup(DIMS(reg), relres, absres, &
                frhoes, DIMS(grd), &
                frhoem, eta, etainv, dfo, dfn, exch, &
                dt, theta) bind(C, name="ceup")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(grd)
  real(rt)         :: frhoes(DIMV(grd))
  real(rt)         :: frhoem(DIMV(grd))
  real(rt)         :: eta(DIMV(grd))
  real(rt)         :: etainv(DIMV(grd))
  real(rt)         :: dfo(DIMV(grd))
  real(rt)         :: dfn(DIMV(grd))
  real(rt)         :: exch(DIMV(grd))
  real(rt)         :: dt, theta, relres, absres
  real(rt)         :: tmp, chg, tot
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           chg = 0.e0_rt
           tot = 0.e0_rt
           tmp = eta(i,j,k) * frhoes(i,j,k) + &
                etainv(i,j,k) * &
                (frhoem(i,j,k) - &
                dt * ((1.e0_rt - theta) * &
                (dfo(i,j,k) - dfn(i,j,k)) + &
                exch(i,j,k)))
           chg = abs(tmp - frhoes(i,j,k))
           tot = abs(frhoes(i,j,k))
           frhoes(i,j,k) = tmp
           absres = max(absres, chg)
           relres = max(relres, chg / (tot + tiny))
        enddo
     enddo
  enddo
end subroutine ceup

subroutine ceupdterm(DIMS(reg), relres, absres, &
                     frhoes, DIMS(grd), &
                     frhoem, eta, etainv, dfo, dfn, exch, dterm, &
                     dt, theta) bind(C, name="ceupdterm")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(grd)
  real(rt)         :: frhoes(DIMV(grd))
  real(rt)         :: frhoem(DIMV(grd))
  real(rt)         :: eta(DIMV(grd))
  real(rt)         :: etainv(DIMV(grd))
  real(rt)         :: dfo(DIMV(grd))
  real(rt)         :: dfn(DIMV(grd))
  real(rt)         :: exch(DIMV(grd))
  real(rt)         :: dterm(DIMV(grd))
  real(rt)         :: dt, theta, relres, absres
  real(rt)         :: tmp, chg, tot
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           chg = 0.e0_rt
           tot = 0.e0_rt
           tmp = eta(i,j,k) * frhoes(i,j,k) + &
                etainv(i,j,k) * &
                (frhoem(i,j,k) - &
                dt * ((1.e0_rt - theta) * &
                (dfo(i,j,k) - dfn(i,j,k)) + &
                exch(i,j,k))) &
                + dt * dterm(i,j,k)
           chg = abs(tmp - frhoes(i,j,k))
           tot = abs(frhoes(i,j,k))
           frhoes(i,j,k) = tmp
           absres = max(absres, chg)
           relres = max(relres, chg / (tot + tiny))
        enddo
     enddo
  enddo
end subroutine ceupdterm

! nonconservative form based on delta B
subroutine nceup(DIMS(reg), relres, absres, &
                 frhoes, DIMS(grd), &
                 frhoem, eta, etainv, &
                 er, DIMS(ebox), &
                 dfo, dfn, temp, fkp, cv, &
                 state, DIMS(sb), &
                 sigma, c, dt, theta) bind(C, name="nceup")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(grd)
  integer :: DIMDEC(sb)
  integer :: DIMDEC(ebox)
  real(rt)         :: frhoes(DIMV(grd))
  real(rt)         :: frhoem(DIMV(grd))
  real(rt)         :: eta(DIMV(grd))
  real(rt)         :: etainv(DIMV(grd))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: dfo(DIMV(grd))
  real(rt)         :: dfn(DIMV(grd))
  real(rt)         :: temp(DIMV(grd))
  real(rt)         :: fkp(DIMV(grd))
  real(rt)         :: cv(DIMV(reg))
  real(rt)         :: state(DIMV(sb), NVAR)
  real(rt)         :: sigma, c, dt, theta, relres, absres
  real(rt)         :: tmp, chg, tot, exch, b, db, dbdt, frhocv
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           chg = 0.e0_rt
           tot = 0.e0_rt
           frhocv = state(i,j,k, URHO) * cv(i,j,k)
           dbdt = 16.e0_rt * sigma * temp(i,j,k)**3
           b = 4.e0_rt * sigma * temp(i,j,k)**4
           exch = fkp(i,j,k) * (b - c * er(i,j,k))
           tmp = eta(i,j,k) * frhoes(i,j,k) + etainv(i,j,k) * &
                (frhoem(i,j,k) - &
                dt * ((1.e0_rt - theta) * &
                (dfo(i,j,k) - dfn(i,j,k)) + &
                exch))
#if 1
           if (frhocv > tiny .AND. tmp > frhoes(i,j,k)) then
              db = (tmp - frhoes(i,j,k)) * dbdt / frhocv
              if (b + db <= 0.e0_rt) then
                 print *, i, j, k, b, db, b+db
              endif
              tmp = ((b + db) / (4.e0_rt * sigma))**0.25e0_rt
              tmp = frhoes(i,j,k) + frhocv * (tmp - temp(i,j,k))
           endif
#endif
           chg = abs(tmp - frhoes(i,j,k))
           tot = abs(frhoes(i,j,k))
           frhoes(i,j,k) = tmp
           absres = max(absres, chg)
           relres = max(relres, chg / (tot + tiny))
        enddo
     enddo
  enddo
end subroutine nceup

subroutine cetot(DIMS(reg), &
                 state, DIMS(sb), &
                 frhoe, DIMS(fb)) bind(C, name="cetot")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(sb)
  integer :: DIMDEC(fb)
  real(rt)         :: state(DIMV(sb), NVAR)
  real(rt)         :: frhoe(DIMV(fb))
  real(rt)         :: kin
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           !               kin = 0.5e0_rt * (state(i,j,k,XMOM)   ** 2 +
           !     @                        state(i,j,k,XMOM+1) ** 2 +
           !     @                        state(i,j,k,XMOM+2) ** 2) /
           !     @                       state(i,j,k,DEN)
           kin = state(i,j,k, UEDEN) - state(i,j,k, UEINT)
           state(i,j,k, UEINT) = frhoe(i,j,k)
           state(i,j,k, UEDEN) = frhoe(i,j,k) + kin
        enddo
     enddo
  enddo
end subroutine cetot

subroutine fkpn(DIMS(reg), &
                fkp, DIMS(fb), &
                const, em, en, &
                ep, nu, tf, &
                temp, DIMS(tb), &
                state, DIMS(sb)) bind(C, name="fkpn")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(fb)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(sb)
  real(rt)         :: fkp(DIMV(fb))
  real(rt)         :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real(rt)         :: ep(0:1), nu
  real(rt)         :: temp(DIMV(tb))
  real(rt)         :: state(DIMV(sb), NVAR)
  real(rt)         :: teff
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           teff = max(temp(i,j,k), tiny)
           teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
           fkp(i,j,k) = const(0) * &
                (state(i,j,k, URHO) ** em(0)) * &
                (teff ** (-en(0))) * &
                (nu ** (ep(0)))
        enddo
     enddo
  enddo
end subroutine fkpn

subroutine rosse1(DIMS(reg), &
                  kappar, DIMS(kbox), &
                  const, em, en, &
                  ep, nu, &
                  tf, kfloor, &
                  temp, DIMS(tb), &
                  state, DIMS(sb)) bind(C, name="rosse1")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(sb)
  real(rt)         :: kappar(DIMV(kbox))
  real(rt)         :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real(rt)         :: ep(0:1), nu
  real(rt)         :: temp(DIMV(tb))
  real(rt)         :: state(DIMV(sb), NVAR)
  real(rt)         :: kfloor
  real(rt)         :: kf, teff
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           teff = max(temp(i,j,k), tiny)
           teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
           kf = const(0) * &
                (state(i,j,k, URHO) ** em(0)) * &
                (teff ** (-en(0))) * &
                (nu ** (ep(0)))
           kappar(i,j,k) = max(kf, kfloor)
        enddo
     enddo
  enddo
end subroutine rosse1

subroutine rosse1s(DIMS(reg), &
                   kappar, DIMS(kbox), &
                   const, em, en, &
                   ep, &
                   sconst, sem, sen, &
                   sep, &
                   nu, &
                   tf, kfloor, &
                   temp, DIMS(tb), &
                   state, DIMS(sb)) bind(C, name="rosse1s")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(sb)
  real(rt)         :: kappar(DIMV(kbox))
  real(rt)         :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real(rt)         :: ep(0:1), nu
  real(rt)         :: sconst(0:1), sem(0:1), sen(0:1), sep(0:1)
  real(rt)         :: temp(DIMV(tb))
  real(rt)         :: state(DIMV(sb), NVAR)
  real(rt)         :: kfloor
  real(rt)         :: kf, teff, sct
  integer :: i, j, k
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           teff = max(temp(i,j,k), tiny)
           teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
           kf = const(0) * &
                (state(i,j,k, URHO) ** em(0)) * &
                (teff ** (-en(0))) * &
                (nu ** (ep(0)))
           sct = sconst(0) * &
                (state(i,j,k, URHO) ** sem(0)) * &
                (teff ** (-sen(0))) * &
                (nu ** (sep(0)))
           kappar(i,j,k) = max(kf+sct, kfloor)
        enddo
     enddo
  enddo
end subroutine rosse1s

subroutine nfloor(dest, &
                  DIMS(dbox), &
                  DIMS(reg), &
                  nflr, flr, nvar) bind(C, name="nfloor")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(dbox)
  integer :: DIMDEC(reg)
  integer :: nvar, nflr
  real(rt)         :: dest(DIMV(dbox), 0:nvar-1)
  real(rt)         :: flr
  integer :: i, j, k, n
  nflr = 0
  do n = 0, nvar-1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              if (dest(i,j,k,n) < flr) then
                 dest(i,j,k,n) = flr
                 nflr = nflr + 1
              endif
           enddo
        enddo
     enddo
  enddo
end subroutine nfloor

! *********************************
! ** BEGIN MGFLD routines        **
! *********************************

subroutine lacoefmgfld(a, &
                       DIMS(abox), &
                       DIMS(reg), &
                       kappa, &
                       DIMS(kbox), &
                       r, s, &
                       dt, c) bind(C, name="lacoefmgfld")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)

  real(rt)         :: a(DIMV(abox))
  real(rt)         :: kappa(DIMV(kbox))
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(reg_l2:reg_h2)
  real(rt)         :: dt, c

  integer :: i, j, k

  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           a(i,j,k) = c*kappa(i,j,k) + 1.e0_rt/dt
        enddo
     enddo
  enddo

end subroutine lacoefmgfld

! *********************************
! ** END MGFLD routines          **
! *********************************

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


subroutine bextrp( &
                  f, fbox_l1, fbox_l2, fbox_l3, fbox_h1, fbox_h2, fbox_h3, &
                  reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3) bind(C, name="bextrp")

  use amrex_fort_module, only : rt => amrex_real
  integer :: fbox_l1, fbox_l2, fbox_l3, fbox_h1, fbox_h2, fbox_h3
  integer ::  reg_l1,  reg_l2,  reg_l3,  reg_h1,  reg_h2,  reg_h3
  real(rt)         :: f(fbox_l1:fbox_h1,fbox_l2:fbox_h2,fbox_l3:fbox_h3)
  integer :: i, j, k

  !     i direction first:
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        i = reg_l1
        f(i-1,j,k) = 2.e0_rt * f(i,j,k) - f(i+1,j,k)
        i = reg_h1
        f(i+1,j,k) = 2.e0_rt * f(i,j,k) - f(i-1,j,k)
     enddo
  enddo

  !     j direction second, including edges:
  do k = reg_l3, reg_h3
     do i = reg_l1 - 1, reg_h1 + 1
        j = reg_l2
        f(i,j-1,k) = 2.e0_rt * f(i,j,k) - f(i,j+1,k)
        j = reg_h2
        f(i,j+1,k) = 2.e0_rt * f(i,j,k) - f(i,j-1,k)
     enddo
  enddo

  !     k direction third, including corners:
  do j = reg_l2 - 1, reg_h2 + 1
     do i = reg_l1 - 1, reg_h1 + 1
        k = reg_l3
        f(i,j,k-1) = 2.e0_rt * f(i,j,k) - f(i,j,k+1)
        k = reg_h3
        f(i,j,k+1) = 2.e0_rt * f(i,j,k) - f(i,j,k-1)
     enddo
  enddo

  !   corner results are the same whichever direction we extrapolate first
end subroutine bextrp


subroutine lbcoefna(bcoef, &
                    bcgrp, bboxl0, bboxl1, bboxl2, bboxh0, bboxh1, bboxh2, &
                    reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3, &
                    spec, sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2, &
                    idim) bind(C, name="lbcoefna")

  use amrex_fort_module, only : rt => amrex_real
  integer :: idim
  integer ::  reg_l1,  reg_l2,  reg_l3,  reg_h1,  reg_h2,  reg_h3
  integer :: bboxl0, bboxl1, bboxl2, bboxh0, bboxh1, bboxh2
  integer :: sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2
  real(rt)         :: bcoef(bboxl0:bboxh0,bboxl1:bboxh1,bboxl2:bboxh2)
  real(rt)         :: bcgrp(bboxl0:bboxh0,bboxl1:bboxh1,bboxl2:bboxh2)
  real(rt)         :: spec(sboxl0:sboxh0,sboxl1:sboxh1,sboxl2:sboxh2)
  integer :: i, j, k
  if (idim == 0) then
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              bcoef(i,j,k) = bcoef(i,j,k) &
                   + 0.5e0_rt * (spec(i-1,j,k) + spec(i,j,k)) * bcgrp(i,j,k)
           enddo
        enddo
     enddo
  else if (idim == 1) then
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              bcoef(i,j,k) = bcoef(i,j,k) &
                   + 0.5e0_rt * (spec(i,j-1,k) + spec(i,j,k)) * bcgrp(i,j,k)
           enddo
        enddo
     enddo
  else
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              bcoef(i,j,k) = bcoef(i,j,k) &
                   + 0.5e0_rt * (spec(i,j,k-1) + spec(i,j,k)) * bcgrp(i,j,k)
           enddo
        enddo
     enddo
  endif

end subroutine lbcoefna


subroutine ljupna( &
                  jnew, jboxl0, jboxl1, jboxl2, jboxh0, jboxh1, jboxh2, &
                  reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3, &
                  spec, sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2, &
                  accel, aboxl0, aboxl1, aboxl2, aboxh0, aboxh1, aboxh2, &
                  nTotal) bind(C, name="ljupna")

  use amrex_fort_module, only : rt => amrex_real
  integer :: nTotal
  integer ::  reg_l1,  reg_l2,  reg_l3,  reg_h1,  reg_h2,  reg_h3
  integer :: jboxl0, jboxl1, jboxl2, jboxh0, jboxh1, jboxh2
  integer :: sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2
  integer :: aboxl0, aboxl1, aboxl2, aboxh0, aboxh1, aboxh2
  real(rt)         :: jnew(jboxl0:jboxh0,jboxl1:jboxh1,jboxl2:jboxh2,0:nTotal-1)
  real(rt)         :: spec(sboxl0:sboxh0,sboxl1:sboxh1,sboxl2:sboxh2,0:nTotal-1)
  real(rt)         :: accel(aboxl0:aboxh0,aboxl1:aboxh1,aboxl2:aboxh2)

  integer :: i, j, k, n
  do n = 0, nTotal - 1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              jnew(i,j,k,n) = jnew(i,j,k,n) + spec(i,j,k,n) * accel(i,j,k)
           enddo
        enddo
     enddo
  enddo

end subroutine ljupna

end module rad_module

