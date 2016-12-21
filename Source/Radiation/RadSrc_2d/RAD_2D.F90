#include "AMReX_LO_BCTYPES.H"
#include "AMReX_ArrayLim.H"

module rad_module

  use bl_types

  use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS, UFX, NVAR

  use rad_util_module, only : FLDlambda

  implicit none

  double precision, parameter :: tiny = 1.e-50_dp_t
  double precision, parameter :: BIGKR = 1.e25_dp_t

contains

subroutine multrs(d, &
                  DIMS(dbox), &
                  DIMS(reg), &
                  r, s) bind(C, name="multrs")
  implicit none
  integer :: DIMDEC(dbox)
  integer :: DIMDEC(reg)
  real*8 :: d(DIMV(dbox))
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(reg_l2:reg_h2)
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        d(i,j) = d(i,j) * r(i) * s(j)
     enddo
  enddo
end subroutine multrs

subroutine sphc(r, s, &
                DIMS(reg), dx) bind(C, name="sphc")

  integer :: DIMDEC(reg)
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(reg_l2:reg_h2)
  real*8 :: dx(2)
  real*8 :: h1, h2, d1, d2
  integer :: i, j
  h1 = 0.5d0 * dx(1)
  h2 = 0.5d0 * dx(2)
  d1 = 1.d0 / (3.d0 * dx(1))
  d2 = 1.d0 / dx(2)
  do i = reg_l1, reg_h1
     r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
  enddo
  do j = reg_l2, reg_h2
     s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
  enddo
end subroutine sphc

subroutine sphe(r, s, n, &
                DIMS(reg), dx) bind(C, name="sphe")

  integer :: DIMDEC(reg)
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(reg_l2:reg_h2)
  integer :: n
  real*8 :: dx(2)
  real*8 :: h1, h2, d1, d2
  integer :: i, j
  if (n == 0) then
     do i = reg_l1, reg_h1
        r(i) = r(i)**2
     enddo
     h2 = 0.5d0 * dx(2)
     d2 = 1.d0 / dx(2)
     do j = reg_l2, reg_h2
        s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
     enddo
  else
     h1 = 0.5d0 * dx(1)
     d1 = 1.d0 / (3.d0 * dx(1))
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

  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  real*8 :: a(DIMV(abox))
  real*8 :: fkp(DIMV(abox))
  real*8 :: eta(DIMV(abox))
  real*8 :: etainv(DIMV(abox))
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(reg_l2:reg_h2)
  real*8 :: c, dt, theta
  integer :: i, j
  real*8 :: dtm
  dtm = 1.d0 / dt
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        a(i,j) = r(i) * s(j) * &
             (fkp(i,j) * etainv(i,j) * c + dtm) / &
             (1.d0 - (1.d0 - theta) * eta(i,j))
     enddo
  enddo
end subroutine lacoef

subroutine bclim(b, &
                 lambda, DIMS(bbox), &
                 DIMS(reg), &
                 n, kappar, DIMS(kbox), &
                 r, s, c, dx) bind(C, name="bclim")

  integer :: DIMDEC(bbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real*8 :: b(DIMV(bbox))
  real*8 :: lambda(DIMV(bbox))
  real*8 :: kappar(DIMV(kbox))
  real*8 :: r(reg_l1:reg_h1+1)
  real*8 :: s(reg_l2:reg_h2+1)
  real*8 :: c, dx(2)
  real*8 :: kavg
  integer :: i, j
  real*8 :: kap
  if (n == 0) then
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1 + 1
           kap = kavg(kappar(i-1,j), kappar(i,j), dx(1), -1)
           b(i,j) = r(i) * s(j) * c * lambda(i,j) / kap
        enddo
     enddo
  else
     do j = reg_l2, reg_h2 + 1
        do i = reg_l1, reg_h1
           kap = kavg(kappar(i,j-1), kappar(i,j), dx(2), -1)
           b(i,j) = r(i) * s(j) * c * lambda(i,j) / kap
        enddo
     enddo
  endif
end subroutine bclim

subroutine flxlim(lambda, &
                  DIMS(rbox), &
                  DIMS(reg), limiter) bind(C, name="flxlim")

  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: limiter
  real*8 :: lambda(DIMV(rbox))
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        lambda(i,j) = FLDlambda(lambda(i,j),limiter)
     enddo
  enddo
end subroutine flxlim

subroutine eddfac(efact, &
                  DIMS(rbox), &
                  DIMS(reg), limiter, n) bind(C, name="eddfac")

  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: n, limiter
  real*8 :: efact(DIMV(rbox))
  integer :: i, j
  real*8 :: r, lambda
  if (n == 0) then
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1 + 1
           r = efact(i,j)
           lambda = FLDlambda(r,limiter)
           efact(i,j) = lambda + (lambda * r)**2
        enddo
     enddo
  else
     do j = reg_l2, reg_h2 + 1
        do i = reg_l1, reg_h1
           r = efact(i,j)
           lambda = FLDlambda(r,limiter)
           efact(i,j) = lambda + (lambda * r)**2
        enddo
     enddo
  endif
end subroutine eddfac

subroutine scgrd1(r, &
                  DIMS(rbox), &
                  DIMS(reg), &
                  n, kappar, DIMS(kbox), &
                  er, dx) bind(C, name="scgrd1")

  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real*8 :: r(DIMV(rbox))
  real*8 :: kappar(DIMV(kbox))
  real*8 :: er(DIMV(kbox))
  real*8 :: dx(2)
  real*8 :: kavg
  integer :: i, j
  real*8 :: kap
  if (n == 0) then
     do j = reg_l2, reg_h2
        ! x derivatives, gradient assembly:
        do i = reg_l1, reg_h1 + 1
           r(i,j) = abs(er(i,j) - er(i-1,j)) / dx(1)
        enddo
        i = reg_l1
        if (er(i-1,j) == -1.d0) then
           r(i,j) = abs(er(i+1,j) - er(i,j)) / dx(1)
        endif
        i = reg_h1 + 1
        if (er(i,j) == -1.d0) then
           r(i,j) = abs(er(i-1,j) - er(i-2,j)) / dx(1)
        endif
        ! construct coefficients:
        do i = reg_l1, reg_h1 + 1
           kap = kavg(kappar(i-1,j), kappar(i,j), dx(1), -1)
           r(i,j) = r(i,j) / &
                (kap * max(er(i-1,j), er(i,j), tiny))
        enddo
     enddo
  else
     ! y derivatives, gradient assembly:
     do j = reg_l2, reg_h2 + 1
        do i = reg_l1, reg_h1
           r(i,j) = abs(er(i,j) - er(i,j-1)) / dx(2)
        enddo
     enddo
     do i = reg_l1, reg_h1
        j = reg_l2
        if (er(i,j-1) == -1.d0) then
           r(i,j) = abs(er(i,j+1) - er(i,j)) / dx(2)
        endif
        j = reg_h2 + 1
        if (er(i,j) == -1.d0) then
           r(i,j) = abs(er(i,j-1) - er(i,j-2)) / dx(2)
        endif
     enddo
     ! construct coefficients:
     do j = reg_l2, reg_h2 + 1
        do i = reg_l1, reg_h1
           kap = kavg(kappar(i,j-1), kappar(i,j), dx(2), -1)
           r(i,j) = r(i,j) / &
                (kap * max(er(i,j-1), er(i,j), tiny))
        enddo
     enddo
  endif
end subroutine scgrd1

subroutine scgrd2(r, &
                  DIMS(rbox), &
                  DIMS(reg), &
                  n, kappar, DIMS(kbox), er, &
                  DIMS(dbox), d, dx) bind(C, name="scgrd2")

  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: DIMDEC(dbox)
  integer :: n
  real*8 :: r(DIMV(rbox))
  real*8 :: kappar(DIMV(kbox))
  real*8 :: er(DIMV(kbox))
  real*8 :: d(DIMV(dbox))
  real*8 :: dx(2)
  real*8 :: kavg
  integer :: i, j
  real*8 :: kap
  if (n == 0) then
     ! y derivatives:
     do j = reg_l2, reg_h2
        do i = reg_l1 - 1, reg_h1 + 1
           d(i,j) = er(i,j+1) - er(i,j-1)
        enddo
     enddo
     ! check y derivatives at reg_l2 - 1 and reg_h2 + 1:
     do i = reg_l1 - 1, reg_h1 + 1
        j = reg_l2
        if (er(i,j-1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
        endif
        j = reg_h2
        if (er(i,j+1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
        endif
     enddo
     do j = reg_l2, reg_h2
        ! check y derivatives at reg_l1 - 1 and reg_h1 + 1:
        ! (check at j-1 and j+1.  if value at j is bad it will not be used at all.)
        i = reg_l1 - 1
        if (er(i,j-1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
        else if (er(i,j+1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
        endif
        i = reg_h1 + 1
        if (er(i,j-1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
        else if (er(i,j+1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
        endif
        ! x derivatives, gradient assembly:
        do i = reg_l1, reg_h1 + 1
           r(i,j) = ((er(i,j) - er(i-1,j)) / dx(1)) ** 2 + &
                ((d(i-1,j) + d(i,j)) / &
                (4.d0 * dx(2))) ** 2
        enddo
        i = reg_l1
        if (er(i-1,j) == -1.d0) then
           r(i,j) = ((er(i+1,j) - er(i,j)) / dx(1)) ** 2 + &
                (d(i,j) / &
                (2.d0 * dx(2))) ** 2
        endif
        i = reg_h1 + 1
        if (er(i,j) == -1.d0) then
           r(i,j) = ((er(i-1,j) - er(i-2,j)) / dx(1)) ** 2 + &
                (d(i-1,j) / &
                (2.d0 * dx(2))) ** 2
        endif
        ! construct scaled gradient:
        do i = reg_l1, reg_h1 + 1
           kap = kavg(kappar(i-1,j), kappar(i,j), dx(1), -1)
           r(i,j) = sqrt(r(i,j)) / &
                (kap * max(er(i-1,j), er(i,j), tiny))
        enddo
     enddo
  else
     ! x derivatives:
     do j = reg_l2 - 1, reg_h2 + 1
        do i = reg_l1, reg_h1
           d(i,j) = er(i+1,j) - er(i-1,j)
        enddo
        ! check x derivatives at reg_l1 - 1 and reg_h1 + 1:
        i = reg_l1
        if (er(i-1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
        endif
        i = reg_h1
        if (er(i+1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
        endif
     enddo
     do i = reg_l1, reg_h1
        ! check x derivatives at reg_l2 - 1 and reg_h2 + 1:
        ! (check at i-1 and i+1.  if value at i is bad it will not be used at all.)
        j = reg_l2 - 1
        if (er(i-1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
        else if (er(i+1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
        endif
        j = reg_h2 + 1
        if (er(i-1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
        else if (er(i+1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
        endif
     enddo
     ! y derivatives, gradient assembly:
     do j = reg_l2, reg_h2 + 1
        do i = reg_l1, reg_h1
           r(i,j) = ((er(i,j) - er(i,j-1)) / dx(2)) ** 2 + &
                ((d(i,j-1) + d(i,j)) / &
                (4.d0 * dx(1))) ** 2
        enddo
     enddo
     do i = reg_l1, reg_h1
        j = reg_l2
        if (er(i,j-1) == -1.d0) then
           r(i,j) = ((er(i,j+1) - er(i,j)) / dx(2)) ** 2 + &
                (d(i,j) / &
                (2.d0 * dx(1))) ** 2
        endif
        j = reg_h2 + 1
        if (er(i,j) == -1.d0) then
           r(i,j) = ((er(i,j-1) - er(i,j-2)) / dx(2)) ** 2 + &
                (d(i,j-1) / &
                (2.d0 * dx(1))) ** 2
        endif
     enddo
     ! construct coefficients:
     do j = reg_l2, reg_h2 + 1
        do i = reg_l1, reg_h1
           kap = kavg(kappar(i,j-1), kappar(i,j), dx(2), -1)
           r(i,j) = sqrt(r(i,j)) / &
                (kap * max(er(i,j-1), er(i,j), tiny))
        enddo
     enddo
  endif
end subroutine scgrd2

subroutine scgrd3(r, &
                  DIMS(rbox), &
                  DIMS(reg), &
                  n, kappar, DIMS(kbox), er, &
                  DIMS(dbox), d, dx) bind(C, name="scgrd3")

  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: DIMDEC(dbox)
  integer :: n
  real*8 :: r(DIMV(rbox))
  real*8 :: kappar(DIMV(kbox))
  real*8 :: er(DIMV(kbox))
  real*8 :: d(DIMV(dbox))
  real*8 :: dx(2)
  real*8 :: kavg
  integer :: i, j
  real*8 :: kap
  if (n == 0) then
     ! y derivatives:
     do j = reg_l2, reg_h2
        do i = reg_l1 - 1, reg_h1 + 1
           d(i,j) = er(i,j+1) - er(i,j-1)
        enddo
     enddo
     ! check y derivatives at reg_l2 - 1 and reg_h2 + 1:
     do i = reg_l1 - 1, reg_h1 + 1
        j = reg_l2
        if (er(i,j-1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
        endif
        j = reg_h2
        if (er(i,j+1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
        endif
     enddo
     do j = reg_l2, reg_h2
        ! check y derivatives at reg_l1 - 1 and reg_h1 + 1:
        ! (check at j-1 and j+1.  if value at j is bad it will not be used at all.)
        i = reg_l1 - 1
        if (er(i,j-1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
        else if (er(i,j+1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
        endif
        i = reg_h1 + 1
        if (er(i,j-1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j+1) - er(i,j))
        else if (er(i,j+1) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i,j-1))
        endif
        ! x derivatives, gradient assembly:
        do i = reg_l1, reg_h1 + 1
           r(i,j) = ((er(i,j) - er(i-1,j)) / dx(1)) ** 2 + &
                ((d(i-1,j) + d(i,j)) / &
                (4.d0 * dx(2))) ** 2
        enddo
        i = reg_l1
        if (er(i-1,j) == -1.d0) then
           r(i,j) = ((er(i+1,j) - er(i,j)) / dx(1)) ** 2 + &
                (d(i,j) / &
                (2.d0 * dx(2))) ** 2
        endif
        i = reg_h1 + 1
        if (er(i,j) == -1.d0) then
           r(i,j) = ((er(i-1,j) - er(i-2,j)) / dx(1)) ** 2 + &
                (d(i-1,j) / &
                (2.d0 * dx(2))) ** 2
        endif
        ! construct scaled gradient:
        do i = reg_l1, reg_h1 + 1
           kap = kavg(kappar(i-1,j), kappar(i,j), dx(1), -1)
           r(i,j) = sqrt(r(i,j)) / &
                (kap * max(er(i-1,j), er(i,j), &
                er(i-1,j-1), er(i,j-1), &
                er(i-1,j+1), er(i,j+1), tiny))
        enddo
     enddo
  else
     ! x derivatives:
     do j = reg_l2 - 1, reg_h2 + 1
        do i = reg_l1, reg_h1
           d(i,j) = er(i+1,j) - er(i-1,j)
        enddo
        ! check x derivatives at reg_l1 - 1 and reg_h1 + 1:
        i = reg_l1
        if (er(i-1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
        endif
        i = reg_h1
        if (er(i+1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
        endif
     enddo
     do i = reg_l1, reg_h1
        ! check x derivatives at reg_l2 - 1 and reg_h2 + 1:
        ! (check at i-1 and i+1.  if value at i is bad it will not be used at all.)
        j = reg_l2 - 1
        if (er(i-1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
        else if (er(i+1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
        endif
        j = reg_h2 + 1
        if (er(i-1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i+1,j) - er(i,j))
        else if (er(i+1,j) == -1.d0) then
           d(i,j) = 2.d0 * (er(i,j) - er(i-1,j))
        endif
     enddo
     ! y derivatives, gradient assembly:
     do j = reg_l2, reg_h2 + 1
        do i = reg_l1, reg_h1
           r(i,j) = ((er(i,j) - er(i,j-1)) / dx(2)) ** 2 + &
                ((d(i,j-1) + d(i,j)) / &
                (4.d0 * dx(1))) ** 2
        enddo
     enddo
     do i = reg_l1, reg_h1
        j = reg_l2
        if (er(i,j-1) == -1.d0) then
           r(i,j) = ((er(i,j+1) - er(i,j)) / dx(2)) ** 2 + &
                (d(i,j) / &
                (2.d0 * dx(1))) ** 2
        endif
        j = reg_h2 + 1
        if (er(i,j) == -1.d0) then
           r(i,j) = ((er(i,j-1) - er(i,j-2)) / dx(2)) ** 2 + &
                (d(i,j-1) / &
                (2.d0 * dx(1))) ** 2
        endif
     enddo
     ! construct coefficients:
     do j = reg_l2, reg_h2 + 1
        do i = reg_l1, reg_h1
           kap = kavg(kappar(i,j-1), kappar(i,j), dx(2), -1)
           r(i,j) = sqrt(r(i,j)) / &
                (kap * max(er(i,j-1), er(i,j), &
                er(i-1,j-1), er(i-1,j), &
                er(i+1,j-1), er(i+1,j), tiny))
        enddo
     enddo
  endif
end subroutine scgrd3

subroutine lrhs(rhs, &
                DIMS(rbox), &
                DIMS(reg), &
                temp, fkp, eta, etainv, frhoem, frhoes, dfo, &
                ero, DIMS(ebox), edot, &
                r, s, dt, sigma, c, theta) bind(C, name="lrhs")

  integer :: DIMDEC(rbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(reg)
  real*8 :: rhs(DIMV(rbox))
  real*8 :: temp(DIMV(rbox))
  real*8 :: fkp(DIMV(rbox))
  real*8 :: eta(DIMV(rbox))
  real*8 :: etainv(DIMV(rbox))
  real*8 :: frhoem(DIMV(rbox))
  real*8 :: frhoes(DIMV(rbox))
  real*8 :: dfo(DIMV(rbox))
  real*8 :: ero(DIMV(ebox))
  real*8 :: edot(DIMV(rbox))
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(reg_l2:reg_h2)
  real*8 :: dt, sigma, c, theta
  integer :: i, j
  real*8 :: dtm, ek, bs, es, ekt
  dtm = 1.d0 / dt
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        ek = fkp(i,j) * eta(i,j)
        bs = etainv(i,j) * &
             &            4.d0 * sigma * fkp(i,j) * temp(i,j)**4
        es = eta(i,j) * (frhoem(i,j) - frhoes(i,j))
        ekt = (1.d0 - theta) * eta(i,j)
        rhs(i,j) = (rhs(i,j) + r(i) * s(j) * &
             (bs + dtm * (ero(i,j) + es) + &
             ek * c * edot(i,j) - &
             ekt * dfo(i,j))) / &
             (1.d0 - ekt)
     enddo
  enddo
end subroutine lrhs

subroutine anatw2(test, &
                  DIMS(reg), &
                  temp, p, xf, Tc, dx, xlo, lo) bind(C, name="anatw2")

  integer :: DIMDEC(reg)
  real*8 :: test(DIMV(reg), 0:1)
  real*8 :: temp(DIMV(reg))
  real*8 :: p, xf, Tc, dx(2), xlo(2)
  integer :: lo(2)
  integer :: i, j
  real*8 :: x, y, r2
  do j = reg_l2, reg_h2
     y = xlo(2) + dx(2) * ((j-lo(2)) + 0.5d0)
     do i = reg_l1, reg_h1
        x  = xlo(1) + dx(1) * ((i-lo(1)) + 0.5d0)
        r2 = x*x + y*y
        test(i,j,0) = Tc * max((1.d0-r2/xf**2), 0.d0)**(1.d0/p)
        test(i,j,1) = temp(i,j) - test(i,j,0)
     enddo
  enddo
end subroutine anatw2

subroutine cfrhoe(DIMS(reg), &
                  frhoe, &
                  DIMS(fb), &
                  state, &
                  DIMS(sb)) bind(C, name="cfrhoe")

  integer :: DIMDEC(reg)
  integer :: DIMDEC(fb)
  integer :: DIMDEC(sb)
  real*8 :: frhoe(DIMV(fb))
  real*8 :: state(DIMV(sb), NVAR)
  !      real*8 kin
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        !            kin = 0.5d0 * (state(i,j,XMOM)   ** 2 +
        !     @                     state(i,j,XMOM+1) ** 2) /
        !     @                    state(i,j,DEN)
        !            frhoe(i,j) = state(i,j,EDEN) - kin
        frhoe(i,j) = state(i,j, UEINT)
     enddo
  enddo
end subroutine cfrhoe

! temp contains frhoe on input:

subroutine gtemp(DIMS(reg), &
                 temp, DIMS(tb), &
                 const, em, en, &
                 state, DIMS(sb)) bind(C, name="gtemp")

  integer :: DIMDEC(reg)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(sb)
  real*8 :: temp(DIMV(tb))
  real*8 :: const(0:1), em(0:1), en(0:1)
  real*8 :: state(DIMV(sb),  NVAR)
  real*8 :: alpha, teff, ex, frhoal
  integer :: i, j
  if (en(0) >= 1.d0) then
     print *, "Bad exponent for cv calculation"
     stop
  endif
  ex = 1.d0 / (1.d0 - en(0))
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        if (em(0) == 0.d0) then
           alpha = const(0)
        else
           alpha = const(0) * state(i,j, URHO) ** em(0)
        endif
        frhoal = state(i,j, URHO) * alpha + tiny
        if (en(0) == 0.d0) then
           temp(i,j) = temp(i,j) / frhoal
        else
           teff = max(temp(i,j), tiny)
           temp(i,j) = ((1.d0 - en(0)) * teff / frhoal) ** ex
        endif
     enddo
  enddo
end subroutine gtemp

! temp contains temp on input:

subroutine gcv(DIMS(reg), &
     cv,DIMS(cbox), &
     temp, DIMS(tbox), &
     const, em, en, tf, &
     state, DIMS(sbox)) bind(C, name="gcv")

  integer :: DIMDEC(reg)
  integer :: DIMDEC(cbox)
  integer :: DIMDEC(tbox)
  integer :: DIMDEC(sbox)
  real*8 :: cv(DIMV(cbox))
  real*8 :: temp(DIMV(tbox))
  real*8 :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real*8 :: state(DIMV(sbox), NVAR)
  real*8 :: alpha, teff, frhoal
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        if (em(0) == 0.d0) then
           alpha = const(0)
        else
           alpha = const(0) * state(i,j, URHO) ** em(0)
        endif
        frhoal = state(i,j, URHO) * alpha + tiny
        if (en(0) == 0.d0) then
           cv(i,j) = alpha
        else
           teff = max(temp(i,j), tiny)
           teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
           cv(i,j) = alpha * teff ** (-en(0))
        endif
     enddo
  enddo
end subroutine gcv

! exch contains temp on input:

subroutine cexch( DIMS(reg), &
     exch, DIMS(xbox), &
     er  , DIMS(ebox), &
     fkp , DIMS(kbox), &
     sigma, c) bind(C, name="cexch")

  integer :: DIMDEC(reg)
  integer :: DIMDEC(xbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(kbox)
  real*8 :: exch(DIMV(xbox))
  real*8 :: er  (DIMV(ebox))
  real*8 :: fkp (DIMV(kbox))
  real*8 :: sigma, c
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        exch(i,j) = fkp(i,j) * &
             (4.d0 * sigma * exch(i,j)**4 &
             - c * er(i,j))
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

  integer :: DIMDEC(reg)
  integer :: DIMDEC(etab)
  integer :: DIMDEC(sb)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(cb)
  integer :: DIMDEC(fb)
  integer :: DIMDEC(ebox)
  real*8 :: eta(DIMV(etab))
  real*8 :: etainv(DIMV(etab))
  real*8 :: frho(DIMV(sb))
  real*8 :: temp(DIMV(tb))
  real*8 :: cv(DIMV(cb))
  real*8 :: fkp(DIMV(fb))
  real*8 :: er(DIMV(ebox))
  real*8 :: dtemp, dtime, sigma, c, underr
  integer :: lagpla
  real*8 :: d, frc, fac0, fac1, fac2
  integer :: i, j
  fac1 = 16.d0 * sigma * dtime
  if (lagpla == 0) then
     fac0 = 0.25d0 * fac1 / dtemp
     fac2 = dtime * c / dtemp
  endif
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        if (lagpla /= 0) then
           ! assume eta and fkp are the same
           d = fac1 * fkp(i,j) * temp(i,j) ** 3
        else
           d = fac0 * (eta(i,j) * (temp(i,j) + dtemp) ** 4 - &
                fkp(i,j) * (temp(i,j)        ) ** 4) - &
                fac2 * (eta(i,j) - fkp(i,j)) * er(i,j)
           ! alternate form, sometimes worse, sometimes better:
           !               d = fac1 * fkp(i,j) * temp(i,j) ** 3 +
           !     @             fac0 * (eta(i,j) - fkp(i,j)) * temp(i,j) ** 4 -
           !     @             fac2 * (eta(i,j) - fkp(i,j)) * er(i,j)
           ! analytic derivatives for specific test problem:
           !               d = (1.2d+6 * sigma * temp(i,j) ** 2 +
           !     @              1.d+5 * c * er(i,j) * (temp(i,j) + tiny) ** (-2)) * dtime
           ! another alternate form (much worse):
           !               d = fac1 * fkp(i,j) * (temp(i,j) + dtemp) ** 3 +
           !     @             fac0 * (eta(i,j) - fkp(i,j)) * (temp(i,j) + dtemp) ** 4 -
           !     @             fac2 * (eta(i,j) - fkp(i,j)) * er(i,j)
        endif
        frc = frho(i,j) * cv(i,j) + tiny
        eta(i,j) = d / (d + frc)
        etainv(i,j) = underr * frc / (d + frc)
        eta(i,j) = 1.d0 - etainv(i,j)
        !            eta(i,j) = 1.d0 - underr * (1.d0 - eta(i,j))
     enddo
  enddo
end subroutine ceta2

subroutine ceup(DIMS(reg), relres, absres, &
                frhoes, DIMS(grd), &
                frhoem, eta, etainv, dfo, dfn, exch, &
                dt, theta) bind(C, name="ceup")

  integer :: DIMDEC(reg)
  integer :: DIMDEC(grd)
  real*8 :: frhoes(DIMV(grd))
  real*8 :: frhoem(DIMV(grd))
  real*8 :: eta(DIMV(grd))
  real*8 :: etainv(DIMV(grd))
  real*8 :: dfo(DIMV(grd))
  real*8 :: dfn(DIMV(grd))
  real*8 :: exch(DIMV(grd))
  real*8 :: dt, theta, relres, absres
  real*8 :: tmp, chg, tot
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        chg = 0.d0
        tot = 0.d0
        tmp = eta(i,j) * frhoes(i,j) + &
             etainv(i,j) * &
             (frhoem(i,j) - &
             dt * ((1.d0 - theta) * &
             (dfo(i,j) - dfn(i,j)) + &
             exch(i,j)))
        chg = abs(tmp - frhoes(i,j))
        tot = abs(frhoes(i,j))
        frhoes(i,j) = tmp
        absres = max(absres, chg)
        relres = max(relres, chg / (tot + tiny))
     enddo
  enddo
end subroutine ceup

subroutine ceupdterm( DIMS(reg), relres, absres, &
     frhoes, DIMS(grd), &
     frhoem, eta, etainv, dfo, dfn, exch, dterm, &
     dt, theta) bind(C, name="ceupdterm")

  integer :: DIMDEC(reg)
  integer :: DIMDEC(grd)
  real*8 :: frhoes(DIMV(grd))
  real*8 :: frhoem(DIMV(grd))
  real*8 :: eta(DIMV(grd))
  real*8 :: etainv(DIMV(grd))
  real*8 :: dfo(DIMV(grd))
  real*8 :: dfn(DIMV(grd))
  real*8 :: exch(DIMV(grd))
  real*8 :: dterm(DIMV(grd))
  real*8 :: dt, theta, relres, absres
  real*8 :: tmp, chg, tot
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        chg = 0.d0
        tot = 0.d0
        tmp = eta(i,j) * frhoes(i,j) + &
             etainv(i,j) * &
             (frhoem(i,j) - &
             dt * ((1.d0 - theta) * &
             (dfo(i,j) - dfn(i,j)) + &
             exch(i,j))) &
             + dt * dterm(i,j)
        chg = abs(tmp - frhoes(i,j))
        tot = abs(frhoes(i,j))
        frhoes(i,j) = tmp
        absres = max(absres, chg)
        relres = max(relres, chg / (tot + tiny))
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

  integer :: DIMDEC(reg)
  integer :: DIMDEC(grd)
  integer :: DIMDEC(sb)
  integer :: DIMDEC(ebox)
  real*8 :: frhoes(DIMV(grd))
  real*8 :: frhoem(DIMV(grd))
  real*8 :: eta(DIMV(grd))
  real*8 :: etainv(DIMV(grd))
  real*8 :: er(DIMV(ebox))
  real*8 :: dfo(DIMV(grd))
  real*8 :: dfn(DIMV(grd))
  real*8 :: temp(DIMV(grd))
  real*8 :: fkp(DIMV(grd))
  real*8 :: cv(DIMV(reg))
  real*8 :: state(DIMV(sb), NVAR)
  real*8 :: sigma, c, dt, theta, relres, absres
  real*8 :: tmp, chg, tot, exch, b, db, dbdt, frhocv
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        chg = 0.d0
        tot = 0.d0
        frhocv = state(i,j, URHO) * cv(i,j)
        dbdt = 16.d0 * sigma * temp(i,j)**3
        b = 4.d0 * sigma * temp(i,j)**4
        exch = fkp(i,j) * (b - c * er(i,j))
        tmp = eta(i,j) * frhoes(i,j) + etainv(i,j) * &
             (frhoem(i,j) - &
             dt * ((1.d0 - theta) * &
             (dfo(i,j) - dfn(i,j)) + &
             exch))
#if 1
        if (frhocv > tiny .AND. tmp > frhoes(i,j)) then
           db = (tmp - frhoes(i,j)) * dbdt / frhocv
           if (b + db <= 0.d0) then
              print *, i, j, b, db, b+db
           endif
           tmp = ((b + db) / (4.d0 * sigma))**0.25d0
           tmp = frhoes(i,j) + frhocv * (tmp - temp(i,j))
        endif
#endif
        chg = abs(tmp - frhoes(i,j))
        tot = abs(frhoes(i,j))
        frhoes(i,j) = tmp
        absres = max(absres, chg)
        relres = max(relres, chg / (tot + tiny))
     enddo
  enddo
end subroutine nceup

subroutine cetot(DIMS(reg), &
                 state, DIMS(sb), &
                 frhoe, DIMS(fb)) bind(C, name="cetot")

  integer :: DIMDEC(reg)
  integer :: DIMDEC(sb)
  integer :: DIMDEC(fb)
  real*8 :: state(DIMV(sb), NVAR)
  real*8 :: frhoe(DIMV(fb))
  real*8 :: kin
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        !            kin = 0.5d0 * (state(i,j,XMOM)   ** 2 +
        !     @                     state(i,j,XMOM+1) ** 2) /
        !     @                    state(i,j,DEN)
        kin = state(i,j, UEDEN) - state(i,j, UEINT)
        state(i,j, UEINT) = frhoe(i,j)
        state(i,j, UEDEN) = frhoe(i,j) + kin
     enddo
  enddo
end subroutine cetot

subroutine fkpn( DIMS(reg), &
     fkp, DIMS(fb), &
     const, em, en, &
     ep, nu, tf, &
     temp, DIMS(tb), &
     state, DIMS(sb)) bind(C, name="fkpn")

  integer :: DIMDEC(reg)
  integer :: DIMDEC(fb)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(sb)
  real*8 :: fkp(DIMV(fb))
  real*8 :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real*8 :: ep(0:1), nu
  real*8 :: temp(DIMV(tb))
  real*8 :: state(DIMV(sb), NVAR)
  real*8 :: teff
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        teff = max(temp(i,j), tiny)
        teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
        fkp(i,j) = const(0) * &
             (state(i,j, URHO) ** em(0)) * &
             (teff ** (-en(0))) * &
             (nu ** (ep(0)))
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

  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(sb)
  real*8 :: kappar(DIMV(kbox))
  real*8 :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real*8 :: ep(0:1), nu
  real*8 :: temp(DIMV(tb))
  real*8 :: state(DIMV(sb), NVAR)
  real*8 :: kfloor
  real*8 :: kf, teff
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        teff = max(temp(i,j), tiny)
        teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
        kf = const(0) * &
             (state(i,j, URHO) ** em(0)) * &
             (teff ** (-en(0))) * &
             (nu ** (ep(0)))
        kappar(i,j) = max(kf, kfloor)
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

  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(sb)
  real*8 :: kappar(DIMV(kbox))
  real*8 :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real*8 :: ep(0:1), nu
  real*8 :: sconst(0:1), sem(0:1), sen(0:1), sep(0:1)
  real*8 :: temp(DIMV(tb))
  real*8 :: state(DIMV(sb), NVAR)
  real*8 :: kfloor
  real*8 :: kf, teff, sct
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        teff = max(temp(i,j), tiny)
        teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
        kf = const(0) * &
             (state(i,j, URHO) ** em(0)) * &
             (teff ** (-en(0))) * &
             (nu ** (ep(0)))
        sct = sconst(0) * &
             (state(i,j, URHO) ** sem(0)) * &
             (teff ** (-sen(0))) * &
             (nu ** (sep(0)))
        kappar(i,j) = max(kf+sct, kfloor)
     enddo
  enddo
end subroutine rosse1s

subroutine nfloor(dest, &
                  DIMS(dbox), &
                  DIMS(reg), &
                  nflr, flr, nvar) bind(C, name="nfloor")

  integer :: DIMDEC(dbox)
  integer :: DIMDEC(reg)
  integer :: nvar, nflr
  real*8 :: dest(DIMV(dbox), 0:nvar-1)
  real*8 :: flr
  integer :: i, j, n
  nflr = 0
  do n = 0, nvar-1
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (dest(i,j,n) < flr) then
              dest(i,j,n) = flr
              nflr = nflr + 1
           endif
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

  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)

  real*8 :: a(DIMV(abox))
  real*8 :: kappa(DIMV(kbox))
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(reg_l2:reg_h2)
  real*8 :: dt, c

  integer :: i, j

  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1

        a(i,j) = c*kappa(i,j) + 1.d0/dt
        a(i,j) = r(i) * s(j) * a(i,j)

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

  integer :: DIMDEC(fbox)
  integer :: DIMDEC(cbox)
  real*8 :: fine(DIMV(fbox))
  real*8 :: crse(DIMV(cbox))
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

  integer :: fbox_l1, fbox_l2, fbox_h1, fbox_h2
  integer ::  reg_l1,  reg_l2,  reg_h1,  reg_h2
  real*8 :: f(fbox_l1:fbox_h1,fbox_l2:fbox_h2)
  integer :: i, j

  !     i direction first:
  do j = reg_l2, reg_h2
     i = reg_l1
     f(i-1,j) = 2.d0 * f(i,j) - f(i+1,j)
     i = reg_h1
     f(i+1,j) = 2.d0 * f(i,j) - f(i-1,j)
  enddo

  !     j direction second, including corners:
  do i = reg_l1 - 1, reg_h1 + 1
     j = reg_l2
     f(i,j-1) = 2.d0 * f(i,j) - f(i,j+1)
     j = reg_h2
     f(i,j+1) = 2.d0 * f(i,j) - f(i,j-1)
  enddo

  !  corner results are the same whichever direction we extrapolate first
end subroutine bextrp


subroutine lbcoefna(bcoef, &
                    bcgrp, bboxl0, bboxl1, bboxh0, bboxh1, &
                    reg_l1, reg_l2, reg_h1, reg_h2, &
                    spec, sboxl0, sboxl1, sboxh0, sboxh1, &
                    idim) bind(C, name="lbcoefna")

  implicit none
  integer :: idim
  integer ::  reg_l1,  reg_l2,  reg_h1,  reg_h2
  integer :: bboxl0, bboxl1, bboxh0, bboxh1
  integer :: sboxl0, sboxl1, sboxh0, sboxh1
  real*8 :: bcoef(bboxl0:bboxh0,bboxl1:bboxh1)
  real*8 :: bcgrp(bboxl0:bboxh0,bboxl1:bboxh1)
  real*8 :: spec(sboxl0:sboxh0,sboxl1:sboxh1)
  integer :: i, j
  if (idim == 0) then
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           bcoef(i,j) = bcoef(i,j) &
                + 0.5d0 * (spec(i-1,j) + spec(i,j)) * bcgrp(i,j)
        enddo
     enddo
  else
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           bcoef(i,j) = bcoef(i,j) &
                + 0.5d0 * (spec(i,j-1) + spec(i,j)) * bcgrp(i,j)
        enddo
     enddo
  endif

end subroutine lbcoefna


subroutine ljupna(jnew, jboxl0, jboxl1, jboxh0, jboxh1, &
                  reg_l1, reg_l2, reg_h1, reg_h2, &
                  spec, sboxl0, sboxl1, sboxh0, sboxh1, &
                  accel, aboxl0, aboxl1, aboxh0, aboxh1, &
                  nTotal) bind(C, name="ljupna")

  integer :: nTotal
  integer ::  reg_l1,  reg_l2,  reg_h1, reg_h2
  integer :: jboxl0, jboxl1, jboxh0, jboxh1
  integer :: sboxl0, sboxl1, sboxh0, sboxh1
  integer :: aboxl0, aboxl1, aboxh0, aboxh1
  real*8 :: jnew(jboxl0:jboxh0,jboxl1:jboxh1,0:nTotal-1)
  real*8 :: spec(sboxl0:sboxh0,sboxl1:sboxh1,0:nTotal-1)
  real*8 :: accel(aboxl0:aboxh0,aboxl1:aboxh1)

  integer :: i, j, n

  do n = 0, nTotal - 1
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           jnew(i,j,n) = jnew(i,j,n) + spec(i,j,n) * accel(i,j)
        enddo
     enddo
  enddo

end subroutine ljupna

end module rad_module
