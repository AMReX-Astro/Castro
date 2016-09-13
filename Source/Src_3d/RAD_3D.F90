#include "LO_BCTYPES.H"

#define dims(a) a/**/l0, a/**/l1, a/**/l2, a/**/h0, a/**/h1, a/**/h2
#define dimdec(a) dims(a)
#define dimv(a) a/**/l0:a/**/h0, a/**/l1:a/**/h1, a/**/l2:a/**/h2

#define tiny 1.d-50
! define big  1.d+50

#define BIGKR 1.d25
! define BIGKP 1.d0

#define REAL_T real*8

#include "FLD_limiter.F"

subroutine setlayout(Den_in, Xmom_in, Eden_in, Eint_in, Temp_in, &
                     FirstSpec, FirstAux, NUM_STATE)
  implicit none
  INCLUDE 'StateLayout.H'
  integer :: Den_in, Xmom_in, Eden_in, Eint_in, Temp_in
  integer :: FirstSpec, FirstAux, NUM_STATE
  DEN = Den_in
  XMOM = Xmom_in
  YMOM = Xmom_in+1
  ZMOM = Xmom_in+2
  EDEN = Eden_in
  EINT = Eint_in
  ITEMP = Temp_in
  IFS = FirstSpec
  IFX = FirstAux
  SVSIZE = NUM_STATE
end subroutine setlayout

! The following routines implement metric terms in 2D and are included
! in the 3D source only to enable the code to link.

subroutine multrs(d, &
                  dims(dbox), &
                  dims(reg), &
                  r, s)
  implicit none
  integer :: dimdec(dbox)
  integer :: dimdec(reg)
  real*8 :: d(dimv(dbox))
  real*8 :: r(regl0:regh0)
  real*8 :: s(regl1:regh1)
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           d(i,j,k) = d(i,j,k) * r(i) * s(j)
        enddo
     enddo
  enddo
end subroutine multrs

subroutine sphc(r, s, &
                dims(reg), dx)
  implicit none
  integer :: dimdec(reg)
  real*8 :: r(regl0:regh0)
  real*8 :: s(regl1:regh1)
  real*8 :: dx(2)
  real*8 :: h1, h2, d1, d2
  integer :: i, j
  h1 = 0.5d0 * dx(1)
  h2 = 0.5d0 * dx(2)
  d1 = 1.d0 / (3.d0 * dx(1))
  d2 = 1.d0 / dx(2)
  do i = regl0, regh0
     r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
  enddo
  do j = regl1, regh1
     s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
  enddo
end subroutine sphc

subroutine sphe(r, s, n, &
                dims(reg), dx)
  implicit none
  integer :: dimdec(reg)
  real*8 :: r(regl0:regh0)
  real*8 :: s(regl1:regh1)
  integer :: n
  real*8 :: dx(2)
  real*8 :: h1, h2, d1, d2
  integer :: i, j
  if (n == 0) then
     do i = regl0, regh0
        r(i) = r(i)**2
     enddo
     h2 = 0.5d0 * dx(2)
     d2 = 1.d0 / dx(2)
     do j = regl1, regh1
        s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
     enddo
  else
     h1 = 0.5d0 * dx(1)
     d1 = 1.d0 / (3.d0 * dx(1))
     do i = regl0, regh0
        r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
     enddo
     do j = regl1, regh1
        s(j) = sin(s(j))
     enddo
  endif
end subroutine sphe

subroutine lacoef(a, &
                  dims(abox), &
                  dims(reg), &
                  fkp, eta, etainv, r, s, c, dt, theta)
  implicit none
  integer :: dimdec(abox)
  integer :: dimdec(reg)
  real*8 :: a(dimv(abox))
  real*8 :: fkp(dimv(abox))
  real*8 :: eta(dimv(abox))
  real*8 :: etainv(dimv(abox))
  real*8 :: r(regl0:regh0)
  real*8 :: s(regl1:regh1)
  real*8 :: c, dt, theta
  integer :: i, j, k
  real*8 :: dtm
  dtm = 1.d0 / dt
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           a(i,j,k) = r(i) * s(j) * &
                (fkp(i,j,k) * etainv(i,j,k) * c + dtm) / &
                (1.d0 - (1.d0 - theta) * eta(i,j,k))
        enddo
     enddo
  enddo
end subroutine lacoef

! arithmetic average, geometrically correct(?) but underestimates surface flux

#define KAVG0(a,b) (0.5d0 * (a + b + tiny))

! define KAVG(a,b,d) KAVG0(a,b)

! harmonic average, overestimates surface flux

#define KAVG1(a,b) ((2.d0 * a * b) / (a + b + tiny) + tiny)

! define KAVG(a,b,d) KAVG1(a,b)

! chooses arithmetic where optically thin, harmonic where optically thick,
! surface flux approximation at a thick/thin boundary

#define KAVG2(a,b,d) min(KAVG0(a,b), max(KAVG1(a,b), 4.d0 / (3.d0 * d)))

! define KAVG(a,b,d) KAVG2(a,b,d)

#define KAVG(a,b,d) kavg(a,b,d,-1)

real*8 function kavg(a, b, d, iopt)
  implicit none
  real*8 :: a, b, d
  integer :: iopt
  integer, save :: opt=100
  if (iopt >= 0) then
     opt = iopt
     if (opt > 2) then
        print *, "Fortran KAVG: invalid averaging option"
     endif
     return
  endif
  if (opt == 0) then
     kavg = KAVG0(a,b)
  else if (opt == 1) then
     kavg = KAVG1(a,b)
  else
     kavg = KAVG2(a,b,d)
  endif
end function kavg

subroutine bclim(b, &
                 lambda, dims(bbox), &
                 dims(reg), &
                 n, kappar, dims(kbox), &
                 r, s, c, dx)
  implicit none
  integer :: dimdec(bbox)
  integer :: dimdec(reg)
  integer :: dimdec(kbox)
  integer :: n
  real*8 :: b(dimv(bbox))
  real*8 :: lambda(dimv(bbox))
  real*8 :: kappar(dimv(kbox))
  real*8 :: r(regl0:regh0+1)
  real*8 :: s(regl1:regh1+1)
  real*8 :: c, dx(3)
  real*8 :: kavg
  integer :: i, j, k
  real*8 :: kap
  if (n == 0) then
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0 + 1
              kap = KAVG(kappar(i-1,j,k), kappar(i,j,k), dx(1))
              b(i,j,k) = r(i) * s(j) * c * lambda(i,j,k) / kap
           enddo
        enddo
     enddo
  else if (n == 1) then
     do k = regl2, regh2
        do j = regl1, regh1 + 1
           do i = regl0, regh0
              kap = KAVG(kappar(i,j-1,k), kappar(i,j,k), dx(2))
              b(i,j,k) = r(i) * s(j) * c * lambda(i,j,k) / kap
           enddo
        enddo
     enddo
  else
     do k = regl2, regh2 + 1
        do j = regl1, regh1
           do i = regl0, regh0
              kap = KAVG(kappar(i,j,k-1), kappar(i,j,k), dx(3))
              b(i,j,k) = r(i) * s(j) * c * lambda(i,j,k) / kap
           enddo
        enddo
     enddo
  endif
end subroutine bclim

subroutine flxlim(lambda, &
                  dims(rbox), &
                  dims(reg), limiter)
  implicit none
  integer :: dimdec(rbox)
  integer :: dimdec(reg)
  integer :: limiter
  real*8 :: lambda(dimv(rbox)), FLDlambda
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           lambda(i,j,k) = FLDlambda(lambda(i,j,k),limiter)
        enddo
     enddo
  enddo
end subroutine flxlim

subroutine eddfac(efact, &
                  dims(rbox), &
                  dims(reg), limiter, n)
  implicit none
  integer :: dimdec(rbox)
  integer :: dimdec(reg)
  integer :: n, limiter
  real*8 :: efact(dimv(rbox))
  integer :: i, j, k
  real*8 :: r, lambda, FLDlambda
  integer :: dir(0:2)
  dir(0) = 0
  dir(1) = 0
  dir(2) = 0
  dir(n) = 1
  do k = regl2, regh2 + dir(2)
     do j = regl1, regh1 + dir(1)
        do i = regl0, regh0 + dir(0)
           r = efact(i,j,k)
           lambda = FLDlambda(r,limiter)
           efact(i,j,k) = lambda + (lambda * r)**2
        enddo
     enddo
  enddo
end subroutine eddfac

subroutine scgrd1(r, &
                  dims(rbox), &
                  dims(reg), &
                  n, kappar, dims(kbox), &
                  er, dx)
  implicit none
  integer :: dimdec(rbox)
  integer :: dimdec(reg)
  integer :: dimdec(kbox)
  integer :: n
  real*8 :: r(dimv(rbox))
  real*8 :: kappar(dimv(kbox))
  real*8 :: er(dimv(kbox))
  real*8 :: dx(3)
  real*8 :: kavg
  integer :: i, j, k
  real*8 :: kap
  if (n == 0) then
     do k = regl2, regh2
        do j = regl1, regh1
           !     x derivatives, gradient assembly:
           do i = regl0, regh0 + 1
              r(i,j,k) = abs(er(i,j,k) - er(i-1,j,k)) / dx(1)
           enddo
           i = regl0
           if (er(i-1,j,k) == -1.d0) then
              r(i,j,k) = abs(er(i+1,j,k) - er(i,j,k)) / dx(1)
           endif
           i = regh0 + 1
           if (er(i,j,k) == -1.d0) then
              r(i,j,k) = abs(er(i-1,j,k) - er(i-2,j,k)) / dx(1)
           endif
           !     construct coefficients:
           do i = regl0, regh0 + 1
              kap = KAVG(kappar(i-1,j,k), kappar(i,j,k), dx(1))
              r(i,j,k) = r(i,j,k) / &
                   (kap * max(er(i-1,j,k), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  else if (n == 1) then
     do k = regl2, regh2
        !     y derivatives, gradient assembly:
        do j = regl1, regh1 + 1
           do i = regl0, regh0
              r(i,j,k) = abs(er(i,j,k) - er(i,j-1,k)) / dx(2)
           enddo
        enddo
        do i = regl0, regh0
           j = regl1
           if (er(i,j-1,k) == -1.d0) then
              r(i,j,k) = abs(er(i,j+1,k) - er(i,j,k)) / dx(2)
           endif
           j = regh1 + 1
           if (er(i,j,k) == -1.d0) then
              r(i,j,k) = abs(er(i,j-1,k) - er(i,j-2,k)) / dx(2)
           endif
        enddo
        ! construct coefficients:
        do j = regl1, regh1 + 1
           do i = regl0, regh0
              kap = KAVG(kappar(i,j-1,k), kappar(i,j,k), dx(2))
              r(i,j,k) = r(i,j,k) / &
                   (kap * max(er(i,j-1,k), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  else
     !     z derivatives, gradient assembly:
     do k = regl2, regh2 +1
        do j = regl1, regh1
           do i = regl0, regh0
              r(i,j,k) = abs(er(i,j,k-1) - er(i,j,k)) / dx(3)
           enddo
        enddo
     enddo
     do j = regl1, regh1
        do i = regl0, regh0
           k = regl2
           if (er(i,j,k-1) == -1.d0) then
              r(i,j,k) = abs(er(i,j,k+1) - er(i,j,k)) / dx(3)
           endif
           k = regh2 + 1
           if (er(i,j,k) == -1.d0) then
              r(i,j,k) = abs(er(i,j,k-2) - er(i,j,k-1)) / dx(3)
           endif
        enddo
     enddo
     ! construct coefficients:
     do k = regl2, regh2 + 1
        do j = regl1, regh1
           do i = regl0, regh0
              kap = KAVG(kappar(i,j,k-1), kappar(i,j,k), dx(3))
              r(i,j,k) = r(i,j,k) / &
                   (kap * max(er(i,j,k-1), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  endif
end subroutine scgrd1

subroutine scgrd2(r, &
                  dims(rbox), &
                  dims(reg), &
                  n, kappar, dims(kbox), er, &
                  dims(dbox), da, db, dx)
  implicit none
  integer :: dimdec(rbox)
  integer :: dimdec(reg)
  integer :: dimdec(kbox)
  integer :: dimdec(dbox)
  integer :: n
  real*8 :: r(dimv(rbox))
  real*8 :: kappar(dimv(kbox))
  real*8 :: er(dimv(kbox))
  real*8 :: da(dimv(dbox))
  real*8 :: db(dimv(dbox))
  real*8 :: dx(3)
  real*8 :: kavg
  integer :: i, j, k
  real*8 :: kap
  if (n == 0) then
     !     y & z derivatives:
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0 - 1, regh0 + 1
              da(i,j,k) = er(i,j+1,k  ) - er(i,j-1,k  )
              db(i,j,k) = er(i,j  ,k+1) - er(i,j  ,k-1)
           enddo
        enddo
     enddo

     !     check y derivatives
     do k = regl2, regh2
        do i = regl0 - 1, regh0 + 1
           j = regl1
           if (er(i,j-1,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
           endif
           j = regh1
           if (er(i,j+1,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
           endif
        enddo
        do j = regl1, regh1
           i = regl0 - 1
           !     (check at j-1 and j+1.  if value at j is bad it will not be used at all.)
           if (er(i,j-1,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
           else if (er(i,j+1,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
           endif
           i = regh0 + 1
           if (er(i,j-1,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
           else if (er(i,j+1,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
           endif
        enddo
     enddo

     !     check z derivatives
     do j = regl1, regh1
        do i = regl0-1, regh0 + 1
           k = regl2
           if (er(i,j,k-1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
           endif
           k = regh2
           if (er(i,j,k+1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
           endif
        enddo
     enddo
     do k = regl2, regh2
        do j = regl1, regh1
           i = regl0 - 1
           !     (check at k-1 and k+1.  if value at k is bad it will not be used at all.)
           if (er(i,j,k-1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
           else if (er(i,j,k+1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
           endif
           i = regh0 + 1
           if (er(i,j,k-1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
           else if (er(i,j,k+1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
           endif
        enddo
     enddo
     !     x derivatives
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0 + 1
              r(i,j,k) = ((er(i,j,k) - er(i-1,j,k)) / dx(1)) ** 2 &
                   + ((da(i-1,j,k) + da(i,j,k)) / &
                   (4.d0 * dx(2))) ** 2 &
                   + ((db(i-1,j,k) + db(i,j,k)) / &
                   (4.d0 * dx(3))) ** 2
           enddo
           i = regl0
           if (er(i-1,j,k) == -1.d0) then
              r(i,j,k) = ((er(i+1,j,k) - er(i,j,k)) / dx(1)) ** 2 &
                   + (da(i,j,k) / (2.d0 * dx(2))) ** 2 &
                   + (db(i,j,k) / (2.d0 * dx(3))) ** 2
           endif
           i = regh0 + 1
           if (er(i,j,k) == -1.d0) then
              r(i,j,k) = ((er(i-1,j,k) - er(i-2,j,k)) / dx(1)) ** 2 &
                   + (da(i-1,j,k) / (2.d0 * dx(2))) ** 2 &
                   + (db(i-1,j,k) / (2.d0 * dx(3))) ** 2
           endif
           !     construct scaled gradient
           do i = regl0, regh0 + 1
              kap = KAVG(kappar(i-1,j,k), kappar(i,j,k), dx(1))
              r(i,j,k) = sqrt(r(i,j,k)) / &
                   (kap * max(er(i-1,j,k), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  else if (n == 1) then
     !     x & z derivatives:
     do k = regl2, regh2
        do j = regl1-1, regh1+1
           do i = regl0, regh0
              da(i,j,k) = er(i+1,j,k) - er(i-1,j,k)
              db(i,j,k) = er(i,j,k+1) - er(i,j,k-1)
           enddo
        enddo
     enddo

     !     check x derivatives
     do k = regl2, regh2
        do j = regl1-1, regh1+1
           i = regl0
           if (er(i-1,j,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
           endif
           i = regh0
           if (er(i+1,j,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i,j,k) - er(i-1,j,k))
           endif
        enddo
        do i = regl0, regh0
           j = regl1-1
           !     (check at i-1 and i+1.  if value at i is bad it will not be used at all.)
           if (er(i-1,j,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
           else if (er(i+1,j,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i,j,k) - er(i-1,j,k))
           endif
           j = regh1+1
           if (er(i-1,j,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
           else if (er(i+1,j,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i,j,k) - er(i-1,j,k))
           endif
        enddo
     enddo

     !     check z derivatives
     do j = regl1-1, regh1+1
        do i = regl0, regh0
           k = regl2
           if (er(i,j,k-1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
           endif
           k = regh2
           if (er(i,j,k+1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
           endif
        enddo
     enddo
     do k = regl2, regh2
        do i = regl0, regh0
           j = regl1-1
           !     (check at k-1 and k+1.  if value at k is bad it will not be used at all.)
           if (er(i,j,k-1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
           else if (er(i,j,k+1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
           endif
           j = regh1+1
           if (er(i,j,k-1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
           else if (er(i,j,k+1) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
           endif
        enddo
     enddo

     !     y derivatives
     do k = regl2, regh2
        do j = regl1, regh1 + 1
           do i = regl0, regh0
              r(i,j,k) = ((er(i,j,k) - er(i,j-1,k)) / dx(2)) ** 2 &
                   + ((da(i,j-1,k) + da(i,j,k)) / &
                   (4.d0 * dx(1))) ** 2 &
                   + ((db(i,j-1,k) + db(i,j,k)) / &
                   (4.d0 * dx(3))) ** 2
           enddo
        enddo
        do i = regl0, regh0
           j = regl1
           if (er(i,j-1,k) == -1.d0) then
              r(i,j,k) = ((er(i,j+1,k) - er(i,j,k)) / dx(2)) ** 2 &
                   + (da(i,j,k) / (2.d0 * dx(1))) ** 2 &
                   + (db(i,j,k) / (2.d0 * dx(3))) ** 2
           endif
           j = regh1 + 1
           if (er(i,j,k) == -1.d0) then
              r(i,j,k) = ((er(i,j-1,k) - er(i,j-2,k)) / dx(2)) ** 2 &
                   + (da(i,j-1,k) / (2.d0 * dx(1))) ** 2 &
                   + (db(i,j-1,k) / (2.d0 * dx(3))) ** 2
           endif
        enddo
        !     construct scaled gradient
        do j = regl1, regh1 + 1
           do i = regl0, regh0
              kap = KAVG(kappar(i,j-1,k), kappar(i,j,k), dx(2))
              r(i,j,k) = sqrt(r(i,j,k)) / &
                   (kap * max(er(i,j-1,k), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  else
     !     x & y derivatives
     do k = regl2-1, regh2+1
        do j = regl1, regh1
           do i = regl0, regh0
              da(i,j,k) = er(i+1,j,k) - er(i-1,j,k)
              db(i,j,k) = er(i,j+1,k) - er(i,j-1,k)
           enddo
        enddo
     enddo

     !     check x derivatives
     do k = regl2-1, regh2+1
        do j = regl1, regh1
           i = regl0
           if (er(i-1,j,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
           endif
           i = regh0
           if (er(i+1,j,k) == -1.d0) then
              da(i,j,k) = 2.0d0 * (er(i,j,k) - er(i-1,j,k))
           endif
        enddo
     enddo
     do j = regl1, regh1
        do i = regl0, regh0
           k = regl2 - 1
           !     (check at i-1 and i+1.  if value at i is bad it will not be used at all.)
           if (er(i-1,j,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
           else if (er(i+1,j,k) == -1.d0) then
              da(i,j,k) = 2.0d0 * (er(i,j,k) - er(i-1,j,k))
           endif
           k = regh2 + 1
           if (er(i-1,j,k) == -1.d0) then
              da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
           else if (er(i+1,j,k) == -1.d0) then
              da(i,j,k) = 2.0d0 * (er(i,j,k) - er(i-1,j,k))
           endif
        enddo
     enddo

     !     check y derivative
     do k = regl2-1, regh2+1
        do i = regl0, regh0
           j = regl1
           if (er(i,j-1,k) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
           endif
           j = regh1
           if (er(i,j+1,k) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
           endif
        enddo
     enddo
     do j = regl1, regh1
        do i = regl0, regh0
           k = regl2 - 1
           !     (check at j-1 and j+1.  if value at j is bad it will not be used at all.)
           if (er(i,j-1,k) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
           else if (er(i,j+1,k) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
           endif
           k = regh2 + 1
           if (er(i,j-1,k) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
           else if (er(i,j+1,k) == -1.d0) then
              db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
           endif
        enddo
     enddo

     ! z derivatives
     do k = regl2, regh2+1
        do j = regl1, regh1
           do i = regl0, regh0
              r(i,j,k) = ((er(i,j,k) - er(i,j,k-1)) / dx(3)) ** 2 &
                   + ((da(i,j,k) + da(i,j,k-1)) / &
                   (4.d0 * dx(1))) ** 2 &
                   + ((db(i,j,k) + db(i,j,k-1)) / &
                   (4.d0 * dx(2))) ** 2
           enddo
        enddo
     enddo

     do j = regl1, regh1
        do i = regl0, regh0
           k = regl2
           if (er(i,j,k-1) == -1.d0) then
              r(i,j,k) = ((er(i,j,k+1) - er(i,j,k)) / dx(3)) ** 2 &
                   + (da(i,j,k) / (2.d0 * dx(1))) ** 2 &
                   + (db(i,j,k) / (2.d0 * dx(2))) ** 2
           endif
           k = regh2+1
           if (er(i,j,k) == -1.d0) then
              r(i,j,k) = ((er(i,j,k-1) - er(i,j,k-2)) / dx(3)) ** 2 &
                   + (da(i,j,k-1) / (2.d0 * dx(1))) ** 2 &
                   + (db(i,j,k-1) / (2.d0 * dx(2))) ** 2
           endif
        enddo
     enddo

     !     gradient assembly
     do k = regl2, regh2+1
        do j = regl1, regh1
           do i = regl0, regh0
              kap = KAVG(kappar(i,j,k), kappar(i,j,k-1), dx(3))
              r(i,j,k) = sqrt(r(i,j,k)) / &
                   (kap * max(er(i,j,k-1), er(i,j,k), tiny))
           enddo
        enddo
     enddo
  endif
end subroutine scgrd2

subroutine scgrd3(r, &
                  dims(rbox), &
                  dims(reg), &
                  n, kappar, dims(kbox), er, &
                  dims(dbox), da, db, dx)
  implicit none
  integer :: dimdec(rbox)
  integer :: dimdec(reg)
  integer :: dimdec(kbox)
  integer :: dimdec(dbox)
  integer :: n
  real*8 :: r(dimv(rbox))
  real*8 :: kappar(dimv(kbox))
  real*8 :: er(dimv(kbox))
  real*8 :: da(dimv(dbox))
  real*8 :: db(dimv(dbox))
  real*8 :: dx(3)
  real*8 :: kavg
  integer :: i
  real*8 :: kap
  print *, "scgrd3 not implemented in 3d"
  stop
end subroutine scgrd3

subroutine lrhs(rhs, &
                dims(rbox), &
                dims(reg), &
                temp, fkp, eta, etainv, frhoem, frhoes, dfo, &
                ero, dims(ebox), edot, &
                r, s, dt, sigma, c, theta)
  implicit none
  integer :: dimdec(rbox)
  integer :: dimdec(ebox)
  integer :: dimdec(reg)
  real*8 :: rhs(dimv(rbox))
  real*8 :: temp(dimv(rbox))
  real*8 :: fkp(dimv(rbox))
  real*8 :: eta(dimv(rbox))
  real*8 :: etainv(dimv(rbox))
  real*8 :: frhoem(dimv(rbox))
  real*8 :: frhoes(dimv(rbox))
  real*8 :: dfo(dimv(rbox))
  real*8 :: ero(dimv(ebox))
  real*8 :: edot(dimv(rbox))
  real*8 :: r(regl0:regh0)
  real*8 :: s(regl1:regh1)
  real*8 :: dt, sigma, c, theta
  integer :: i, j, k
  real*8 :: dtm, ek, bs, es, ekt
  dtm = 1.d0 / dt
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           ek = fkp(i,j,k) * eta(i,j,k)
           bs = etainv(i,j,k) * &
                &               4.d0 * sigma * fkp(i,j,k) * temp(i,j,k)**4
           es = eta(i,j,k) * (frhoem(i,j,k) - frhoes(i,j,k))
           ekt = (1.d0 - theta) * eta(i,j,k)
           rhs(i,j,k) = (rhs(i,j,k) + r(i) * s(j) * &
                (bs + dtm * (ero(i,j,k) + es) + &
                ek * c * edot(i,j,k) - &
                ekt * dfo(i,j,k))) / &
                (1.d0 - ekt)
        enddo
     enddo
  enddo
end subroutine lrhs

subroutine anatw2(test, &
                  dims(reg), &
                  temp, p, xf, Tc, dx, xlo, lo)
  implicit none
  integer :: dimdec(reg)
  real*8 :: test(dimv(reg), 0:1)
  real*8 :: temp(dimv(reg))
  real*8 :: p, xf, Tc, dx(3), xlo(3)
  integer :: lo(3)
  integer :: i, j, k
  real*8 :: x, y, z, r2
  do k = regl2, regh2
     z = xlo(3) + dx(3) * ((k-lo(3)) + 0.5d0)
     do j = regl1, regh1
        y = xlo(2) + dx(2) * ((j-lo(2)) + 0.5d0)
        do i = regl0, regh0
           x  = xlo(1) + dx(1) * ((i-lo(1)) + 0.5d0)
           r2 = x*x + y*y + z*z
           test(i,j,k,0) = Tc * max((1.d0-r2/xf**2), 0.d0)**(1.d0/p)
           test(i,j,k,1) = temp(i,j,k) - test(i,j,k,0)
        enddo
     enddo
  enddo
end subroutine anatw2

subroutine cfrhoe(dims(reg), &
                  frhoe, &
                  dims(fb), &
                  state, &
                  dims(sb))
  implicit none
  INCLUDE 'StateLayout.H'
  integer :: dimdec(reg)
  integer :: dimdec(fb)
  integer :: dimdec(sb)
  real*8 :: frhoe(dimv(fb))
  real*8 :: state(dimv(sb), 0:SVSIZE-1)
  !      real*8 kin
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           !               kin = 0.5d0 * (state(i,j,k,XMOM)   ** 2 +
           !     @                        state(i,j,k,XMOM+1) ** 2 +
           !     @                        state(i,j,k,XMOM+2) ** 2) /
           !     @                       state(i,j,k,DEN)
           !               frhoe(i,j,k) = state(i,j,k,EDEN) - kin
           frhoe(i,j,k) = state(i,j,k,EINT)
        enddo
     enddo
  enddo
end subroutine cfrhoe

! temp contains frhoe on input:

subroutine gtemp(dims(reg), &
                 temp, dims(tb), &
                 const, em, en, &
                 state, dims(sb))
  implicit none
  INCLUDE 'StateLayout.H'
  integer :: dimdec(reg)
  integer :: dimdec(tb)
  integer :: dimdec(sb)
  real*8 :: temp(dimv(tb))
  real*8 :: const(0:1), em(0:1), en(0:1)
  real*8 :: state(dimv(sb), 0:SVSIZE-1)
  real*8 :: alpha, teff, ex, frhoal
  integer :: i, j, k
  if (en(0) >= 1.d0) then
     print *, "Bad exponent for cv calculation"
     stop
  endif
  ex = 1.d0 / (1.d0 - en(0))
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           if (em(0) == 0.d0) then
              alpha = const(0)
           else
              alpha = const(0) * state(i,j,k,DEN) ** em(0)
           endif
           frhoal = state(i,j,k,DEN) * alpha + tiny
           if (en(0) == 0.d0) then
              temp(i,j,k) = temp(i,j,k) / frhoal
           else
              teff = max(temp(i,j,k), tiny)
              temp(i,j,k) = ((1.d0 - en(0)) * teff / frhoal) ** ex
           endif
        enddo
     enddo
  enddo
end subroutine gtemp

! temp contains temp on input:

subroutine gcv(dims(reg), &
               cv,dims(cbox), &
               temp, dims(tbox), &
               const, em, en, tf, &
               state, dims(sbox))
  implicit none
  INCLUDE 'StateLayout.H'
  integer :: dimdec(reg)
  integer :: dimdec(cbox)
  integer :: dimdec(tbox)
  integer :: dimdec(sbox)
  real*8 :: cv(dimv(cbox))
  real*8 :: temp(dimv(tbox))
  real*8 :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real*8 :: state(dimv(sbox), 0:SVSIZE-1)
  real*8 :: alpha, teff, frhoal
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           if (em(0) == 0.d0) then
              alpha = const(0)
           else
              alpha = const(0) * state(i,j,k,DEN) ** em(0)
           endif
           frhoal = state(i,j,k,DEN) * alpha + tiny
           if (en(0) == 0.d0) then
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

subroutine cexch(dims(reg), &
                 exch, dims(xbox), &
                 er  , dims(ebox), &
                 fkp , dims(kbox), &
                 sigma, c)
  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(xbox)
  integer :: dimdec(ebox)
  integer :: dimdec(kbox)
  real*8 :: exch(dimv(xbox))
  real*8 :: er  (dimv(ebox))
  real*8 :: fkp (dimv(kbox))
  real*8 :: sigma, c
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           exch(i,j,k) = fkp(i,j,k) * &
                (4.d0 * sigma * exch(i,j,k)**4 &
                - c * er(i,j,k))
        enddo
     enddo
  enddo
end subroutine cexch

subroutine ceta2(dims(reg), &
                 eta, etainv, dims(etab), &
                 frho, dims(sb), &
                 temp, dims(tb), &
                 cv, dims(cb), &
                 fkp, dims(fb), &
                 er, dims(ebox), &
                 dtemp, dtime, sigma, c, underr, lagpla)
  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(etab)
  integer :: dimdec(sb)
  integer :: dimdec(tb)
  integer :: dimdec(cb)
  integer :: dimdec(fb)
  integer :: dimdec(ebox)
  real*8 :: eta(dimv(etab))
  real*8 :: etainv(dimv(etab))
  real*8 :: frho(dimv(sb))
  real*8 :: temp(dimv(tb))
  real*8 :: cv(dimv(cb))
  real*8 :: fkp(dimv(fb))
  real*8 :: er(dimv(ebox))
  real*8 :: dtemp, dtime, sigma, c, underr
  integer :: lagpla
  real*8 :: d, frc, fac0, fac1, fac2
  integer :: i, j, k
  fac1 = 16.d0 * sigma * dtime
  if (lagpla == 0) then
     fac0 = 0.25d0 * fac1 / dtemp
     fac2 = dtime * c / dtemp
  endif
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
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
           eta(i,j,k) = 1.d0 - etainv(i,j,k)
           !               eta(i,j,k) = 1.d0 - underr * (1.d0 - eta(i,j,k))
        enddo
     enddo
  enddo
end subroutine ceta2

subroutine ceup(dims(reg), relres, absres, &
                frhoes, dims(grd), &
                frhoem, eta, etainv, dfo, dfn, exch, &
                dt, theta)
  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(grd)
  real*8 :: frhoes(dimv(grd))
  real*8 :: frhoem(dimv(grd))
  real*8 :: eta(dimv(grd))
  real*8 :: etainv(dimv(grd))
  real*8 :: dfo(dimv(grd))
  real*8 :: dfn(dimv(grd))
  real*8 :: exch(dimv(grd))
  real*8 :: dt, theta, relres, absres
  real*8 :: tmp, chg, tot
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           chg = 0.d0
           tot = 0.d0
           tmp = eta(i,j,k) * frhoes(i,j,k) + &
                etainv(i,j,k) * &
                (frhoem(i,j,k) - &
                dt * ((1.d0 - theta) * &
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

subroutine ceupdterm(dims(reg), relres, absres, &
                     frhoes, dims(grd), &
                     frhoem, eta, etainv, dfo, dfn, exch, dterm, &
                     dt, theta)
  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(grd)
  real*8 :: frhoes(dimv(grd))
  real*8 :: frhoem(dimv(grd))
  real*8 :: eta(dimv(grd))
  real*8 :: etainv(dimv(grd))
  real*8 :: dfo(dimv(grd))
  real*8 :: dfn(dimv(grd))
  real*8 :: exch(dimv(grd))
  real*8 :: dterm(dimv(grd))
  real*8 :: dt, theta, relres, absres
  real*8 :: tmp, chg, tot
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           chg = 0.d0
           tot = 0.d0
           tmp = eta(i,j,k) * frhoes(i,j,k) + &
                etainv(i,j,k) * &
                (frhoem(i,j,k) - &
                dt * ((1.d0 - theta) * &
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
subroutine nceup(dims(reg), relres, absres, &
                 frhoes, dims(grd), &
                 frhoem, eta, etainv, &
                 er, dims(ebox), &
                 dfo, dfn, temp, fkp, cv, &
                 state, dims(sb), &
                 sigma, c, dt, theta)
  implicit none
  INCLUDE 'StateLayout.H'
  integer :: dimdec(reg)
  integer :: dimdec(grd)
  integer :: dimdec(sb)
  integer :: dimdec(ebox)
  real*8 :: frhoes(dimv(grd))
  real*8 :: frhoem(dimv(grd))
  real*8 :: eta(dimv(grd))
  real*8 :: etainv(dimv(grd))
  real*8 :: er(dimv(ebox))
  real*8 :: dfo(dimv(grd))
  real*8 :: dfn(dimv(grd))
  real*8 :: temp(dimv(grd))
  real*8 :: fkp(dimv(grd))
  real*8 :: cv(dimv(reg))
  real*8 :: state(dimv(sb), 0:SVSIZE-1)
  real*8 :: sigma, c, dt, theta, relres, absres
  real*8 :: tmp, chg, tot, exch, b, db, dbdt, frhocv
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           chg = 0.d0
           tot = 0.d0
           frhocv = state(i,j,k,DEN) * cv(i,j,k)
           dbdt = 16.d0 * sigma * temp(i,j,k)**3
           b = 4.d0 * sigma * temp(i,j,k)**4
           exch = fkp(i,j,k) * (b - c * er(i,j,k))
           tmp = eta(i,j,k) * frhoes(i,j,k) + etainv(i,j,k) * &
                (frhoem(i,j,k) - &
                dt * ((1.d0 - theta) * &
                (dfo(i,j,k) - dfn(i,j,k)) + &
                exch))
#if 1
           if (frhocv > tiny .AND. tmp > frhoes(i,j,k)) then
              db = (tmp - frhoes(i,j,k)) * dbdt / frhocv
              if (b + db <= 0.d0) then
                 print *, i, j, k, b, db, b+db
              endif
              tmp = ((b + db) / (4.d0 * sigma))**0.25d0
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

subroutine cetot(dims(reg), &
                 state, dims(sb), &
                 frhoe, dims(fb))
  implicit none
  INCLUDE 'StateLayout.H'
  integer :: dimdec(reg)
  integer :: dimdec(sb)
  integer :: dimdec(fb)
  real*8 :: state(dimv(sb), 0:SVSIZE-1)
  real*8 :: frhoe(dimv(fb))
  real*8 :: kin
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           !               kin = 0.5d0 * (state(i,j,k,XMOM)   ** 2 +
           !     @                        state(i,j,k,XMOM+1) ** 2 +
           !     @                        state(i,j,k,XMOM+2) ** 2) /
           !     @                       state(i,j,k,DEN)
           kin = state(i,j,k,EDEN) - state(i,j,k,EINT)
           state(i,j,k,EINT) = frhoe(i,j,k)
           state(i,j,k,EDEN) = frhoe(i,j,k) + kin
        enddo
     enddo
  enddo
end subroutine cetot

subroutine fkpn(dims(reg), &
                fkp, dims(fb), &
                const, em, en, &
                ep, nu, tf, &
                temp, dims(tb), &
                state, dims(sb))
  implicit none
  INCLUDE 'StateLayout.H'
  integer :: dimdec(reg)
  integer :: dimdec(fb)
  integer :: dimdec(tb)
  integer :: dimdec(sb)
  real*8 :: fkp(dimv(fb))
  real*8 :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real*8 :: ep(0:1), nu
  real*8 :: temp(dimv(tb))
  real*8 :: state(dimv(sb), 0:SVSIZE-1)
  real*8 :: teff
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           teff = max(temp(i,j,k), tiny)
           teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
           fkp(i,j,k) = const(0) * &
                (state(i,j,k,DEN) ** em(0)) * &
                (teff ** (-en(0))) * &
                (nu ** (ep(0)))
        enddo
     enddo
  enddo
end subroutine fkpn

subroutine rosse1(dims(reg), &
                  kappar, dims(kbox), &
                  const, em, en, &
                  ep, nu, &
                  tf, kfloor, &
                  temp, dims(tb), &
                  state, dims(sb))
  implicit none
  INCLUDE 'StateLayout.H'
  integer :: dimdec(reg)
  integer :: dimdec(kbox)
  integer :: dimdec(tb)
  integer :: dimdec(sb)
  real*8 :: kappar(dimv(kbox))
  real*8 :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real*8 :: ep(0:1), nu
  real*8 :: temp(dimv(tb))
  real*8 :: state(dimv(sb), 0:SVSIZE-1)
  real*8 :: kfloor
  real*8 :: kf, teff
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           teff = max(temp(i,j,k), tiny)
           teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
           kf = const(0) * &
                (state(i,j,k,DEN) ** em(0)) * &
                (teff ** (-en(0))) * &
                (nu ** (ep(0)))
           kappar(i,j,k) = max(kf, kfloor)
        enddo
     enddo
  enddo
end subroutine rosse1

subroutine rosse1s(dims(reg), &
                   kappar, dims(kbox), &
                   const, em, en, &
                   ep, &
                   sconst, sem, sen, &
                   sep, &
                   nu, &
                   tf, kfloor, &
                   temp, dims(tb), &
                   state, dims(sb))
  implicit none
  INCLUDE 'StateLayout.H'
  integer :: dimdec(reg)
  integer :: dimdec(kbox)
  integer :: dimdec(tb)
  integer :: dimdec(sb)
  real*8 :: kappar(dimv(kbox))
  real*8 :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real*8 :: ep(0:1), nu
  real*8 :: sconst(0:1), sem(0:1), sen(0:1), sep(0:1)
  real*8 :: temp(dimv(tb))
  real*8 :: state(dimv(sb), 0:SVSIZE-1)
  real*8 :: kfloor
  real*8 :: kf, teff, sct
  integer :: i, j, k
  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           teff = max(temp(i,j,k), tiny)
           teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
           kf = const(0) * &
                (state(i,j,k,DEN) ** em(0)) * &
                (teff ** (-en(0))) * &
                (nu ** (ep(0)))
           sct = sconst(0) * &
                (state(i,j,k,DEN) ** sem(0)) * &
                (teff ** (-sen(0))) * &
                (nu ** (sep(0)))
           kappar(i,j,k) = max(kf+sct, kfloor)
        enddo
     enddo
  enddo
end subroutine rosse1s

subroutine nfloor(dest, &
                  dims(dbox), &
                  dims(reg), &
                  nflr, flr, nvar)
  implicit none
  integer :: dimdec(dbox)
  integer :: dimdec(reg)
  integer :: nvar, nflr
  real*8 :: dest(dimv(dbox), 0:nvar-1)
  real*8 :: flr
  integer :: i, j, k, n
  nflr = 0
  do n = 0, nvar-1
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
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
                       dims(abox), &
                       dims(reg), &
                       kappa, &
                       dims(kbox), &
                       r, s, &
                       dt, c)
  implicit none
  integer :: dimdec(abox)
  integer :: dimdec(reg)
  integer :: dimdec(kbox)

  real*8 :: a(dimv(abox))
  real*8 :: kappa(dimv(kbox))
  real*8 :: r(regl0:regh0)
  real*8 :: s(regl1:regh1)
  real*8 :: dt, c

  integer :: i, j, k

  do k = regl2, regh2
     do j = regl1, regh1
        do i = regl0, regh0
           a(i,j,k) = c*kappa(i,j,k) + 1.d0/dt
        enddo
     enddo
  enddo

end subroutine lacoefmgfld

! *********************************
! ** END MGFLD routines          **
! *********************************

subroutine rfface(fine, &
                  dims(fbox), &
                  crse, &
                  dims(cbox), &
                  idim, irat)
  implicit none
  integer :: dimdec(fbox)
  integer :: dimdec(cbox)
  real*8 :: fine(dimv(fbox))
  real*8 :: crse(dimv(cbox))
  integer :: idim, irat(0:2)
  integer :: i, j, k
  integer :: rfac
  if (idim == 0) then
     rfac = irat(1) * irat(2)
     do k = fboxl2, fboxh2
        do j = fboxl1, fboxh1
           fine(fboxl0,j,k) = crse(cboxl0, j/irat(1), k/irat(2)) / rfac
        enddo
     enddo
  else if (idim == 1) then
     rfac = irat(0) * irat(2)
     do k = fboxl2, fboxh2
        do i = fboxl0, fboxh0
           fine(i,fboxl1,k) = crse(i/irat(0), cboxl1, k/irat(2)) / rfac
        enddo
     enddo
  else
     rfac = irat(0) * irat(1)
     do j = fboxl1, fboxh1
        do i = fboxl0, fboxh0
           fine(i,j,fboxl2) = crse(i/irat(0), j/irat(1), cboxl2) / rfac
        enddo
     enddo
  endif
end subroutine rfface


subroutine bextrp( &
                  f, fboxl0, fboxl1, fboxl2, fboxh0, fboxh1, fboxh2, &
                  regl0, regl1, regl2, regh0, regh1, regh2)

  implicit none
  integer :: fboxl0, fboxl1, fboxl2, fboxh0, fboxh1, fboxh2
  integer ::  regl0,  regl1,  regl2,  regh0,  regh1,  regh2
  real*8 :: f(fboxl0:fboxh0,fboxl1:fboxh1,fboxl2:fboxh2)
  integer :: i, j, k

  !     i direction first:
  do k = regl2, regh2
     do j = regl1, regh1
        i = regl0
        f(i-1,j,k) = 2.d0 * f(i,j,k) - f(i+1,j,k)
        i = regh0
        f(i+1,j,k) = 2.d0 * f(i,j,k) - f(i-1,j,k)
     enddo
  enddo

  !     j direction second, including edges:
  do k = regl2, regh2
     do i = regl0 - 1, regh0 + 1
        j = regl1
        f(i,j-1,k) = 2.d0 * f(i,j,k) - f(i,j+1,k)
        j = regh1
        f(i,j+1,k) = 2.d0 * f(i,j,k) - f(i,j-1,k)
     enddo
  enddo

  !     k direction third, including corners:
  do j = regl1 - 1, regh1 + 1
     do i = regl0 - 1, regh0 + 1
        k = regl2
        f(i,j,k-1) = 2.d0 * f(i,j,k) - f(i,j,k+1)
        k = regh2
        f(i,j,k+1) = 2.d0 * f(i,j,k) - f(i,j,k-1)
     enddo
  enddo

  !   corner results are the same whichever direction we extrapolate first
end subroutine bextrp


subroutine lbcoefna(bcoef, &
                    bcgrp, bboxl0, bboxl1, bboxl2, bboxh0, bboxh1, bboxh2, &
                    regl0, regl1, regl2, regh0, regh1, regh2, &
                    spec, sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2, &
                    idim)

  implicit none
  integer :: idim
  integer ::  regl0,  regl1,  regl2,  regh0,  regh1,  regh2
  integer :: bboxl0, bboxl1, bboxl2, bboxh0, bboxh1, bboxh2
  integer :: sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2
  real*8 :: bcoef(bboxl0:bboxh0,bboxl1:bboxh1,bboxl2:bboxh2)
  real*8 :: bcgrp(bboxl0:bboxh0,bboxl1:bboxh1,bboxl2:bboxh2)
  real*8 :: spec(sboxl0:sboxh0,sboxl1:sboxh1,sboxl2:sboxh2)
  integer :: i, j, k
  if (idim == 0) then
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              bcoef(i,j,k) = bcoef(i,j,k) &
                   + 0.5d0 * (spec(i-1,j,k) + spec(i,j,k)) * bcgrp(i,j,k)
           enddo
        enddo
     enddo
  else if (idim == 1) then
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              bcoef(i,j,k) = bcoef(i,j,k) &
                   + 0.5d0 * (spec(i,j-1,k) + spec(i,j,k)) * bcgrp(i,j,k)
           enddo
        enddo
     enddo
  else
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              bcoef(i,j,k) = bcoef(i,j,k) &
                   + 0.5d0 * (spec(i,j,k-1) + spec(i,j,k)) * bcgrp(i,j,k)
           enddo
        enddo
     enddo
  endif

end subroutine lbcoefna


subroutine ljupna( &
                  jnew, jboxl0, jboxl1, jboxl2, jboxh0, jboxh1, jboxh2, &
                  regl0, regl1, regl2, regh0, regh1, regh2, &
                  spec, sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2, &
                  accel, aboxl0, aboxl1, aboxl2, aboxh0, aboxh1, aboxh2, &
                  nTotal)

  implicit none
  integer :: nTotal
  integer ::  regl0,  regl1,  regl2,  regh0,  regh1,  regh2
  integer :: jboxl0, jboxl1, jboxl2, jboxh0, jboxh1, jboxh2
  integer :: sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2
  integer :: aboxl0, aboxl1, aboxl2, aboxh0, aboxh1, aboxh2
  real*8 :: jnew(jboxl0:jboxh0,jboxl1:jboxh1,jboxl2:jboxh2,0:nTotal-1)
  real*8 :: spec(sboxl0:sboxh0,sboxl1:sboxh1,sboxl2:sboxh2,0:nTotal-1)
  real*8 :: accel(aboxl0:aboxh0,aboxl1:aboxh1,aboxl2:aboxh2)

  integer :: i, j, k, n
  do n = 0, nTotal - 1
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              jnew(i,j,k,n) = jnew(i,j,k,n) + spec(i,j,k,n) * accel(i,j,k)
           enddo
        enddo
     enddo
  enddo

end subroutine ljupna
