#include "LO_BCTYPES.H"
#include "ArrayLim.H"

module rad_module

  use bl_types

  use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS, UFX, NVAR

  use rad_util_module, only : FLDlambda

  implicit none

  double precision, parameter :: tiny = 1.d-50
  double precision, parameter :: BIGKR = 1.e25_dp_t

contains

subroutine multrs(d, &
                  DIMS(dbox), &
                  DIMS(reg), &
                  r, s) bind(C, name="multrs")

  integer :: DIMDEC(dbox)
  integer :: DIMDEC(reg)
  real*8 :: d(DIMV(dbox))
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(1)
  integer :: i
  do i = reg_l1, reg_h1
     d(i) = d(i) * r(i)
  enddo
end subroutine multrs

subroutine sphc(r, s, &
                DIMS(reg), dx) bind(C, name="sphc")
  implicit none
  integer :: DIMDEC(reg)
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(1)
  real*8 :: dx(1)
  real*8 :: h1, d1
  integer :: i
  h1 = 0.5d0 * dx(1)
  d1 = 1.d0 / (3.d0 * dx(1))
  do i = reg_l1, reg_h1
     r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
  enddo
end subroutine sphc

subroutine sphe(r, s, n, &
                DIMS(reg), dx) bind(C, name="sphe")
  implicit none
  integer :: DIMDEC(reg)
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(1)
  integer :: n
  real*8 :: dx(1)
  integer :: i
  do i = reg_l1, reg_h1
     r(i) = r(i)**2
  enddo
end subroutine sphe

subroutine lacoef(a, &
                  DIMS(abox), &
                  DIMS(reg), &
                  fkp, eta, etainv, r, s, c, dt, theta) bind(C, name="lacoef")
  implicit none
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  real*8 :: a(DIMV(abox))
  real*8 :: fkp(DIMV(abox))
  real*8 :: eta(DIMV(abox))
  real*8 :: etainv(DIMV(abox))
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(1)
  real*8 :: c, dt, theta
  integer :: i
  real*8 :: dtm
  dtm = 1.d0 / dt
  do i = reg_l1, reg_h1
     a(i) = r(i) * &
          (fkp(i) * etainv(i) * c + dtm) / &
          (1.d0 - (1.d0 - theta) * eta(i))
  enddo
end subroutine lacoef

subroutine bclim(b, &
                 lambda, DIMS(bbox), &
                 DIMS(reg), &
                 n, kappar, DIMS(kbox), &
                 r, s, c, dx) bind(C, name="bclim")
  implicit none
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real*8 :: b(DIMV(bbox))
  real*8 :: lambda(DIMV(bbox))
  real*8 :: kappar(DIMV(kbox))
  real*8 :: r(reg_l1:reg_h1+1)
  real*8 :: s(1)
  real*8 :: c, dx(1)
  real*8 :: kavg
  integer :: i
  real*8 :: kap
  do i = reg_l1, reg_h1 + 1
     kap = kavg(kappar(i-1), kappar(i), dx(1),-1)
     b(i) = r(i) * c * lambda(i) / kap
  enddo
end subroutine bclim

subroutine flxlim(lambda, &
                  DIMS(rbox), &
                  DIMS(reg), limiter) bind(C, name="flxlim")
  implicit none
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: limiter
  real*8 :: lambda(DIMV(rbox))
  integer :: i
  do i = reg_l1, reg_h1
     lambda(i) = FLDlambda(lambda(i),limiter)
  enddo
end subroutine flxlim

subroutine eddfac(efact, &
                  DIMS(rbox), &
                  DIMS(reg), limiter, n)
  implicit none
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: n, limiter
  real*8 :: efact(DIMV(rbox))
  integer :: i
  real*8 :: r, lambda
  do i = reg_l1, reg_h1 + 1
     r = efact(i)
     lambda = FLDlambda(r,limiter)
     efact(i) = lambda + (lambda * r)**2
  enddo
end subroutine eddfac

subroutine scgrd1(r, &
                  DIMS(rbox), &
                  DIMS(reg), &
                  n, kappar, DIMS(kbox), &
                  er, dx) bind(C, name="scgrd1")
  implicit none
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real*8 :: r(DIMV(rbox))
  real*8 :: kappar(DIMV(kbox))
  real*8 :: er(DIMV(kbox))
  real*8 :: dx(1)
  real*8 :: kavg
  integer :: i
  real*8 :: kap
  ! gradient assembly:
  do i = reg_l1, reg_h1 + 1
     r(i) = abs(er(i) - er(i-1)) / dx(1)
  enddo
  i = reg_l1
  if (er(i-1) == -1.d0) then
     r(i) = abs(er(i+1) - er(i)) / dx(1)
  endif
  i = reg_h1 + 1
  if (er(i) == -1.d0) then
     r(i) = abs(er(i-1) - er(i-2)) / dx(1)
  endif
  ! construct scaled gradient:
  do i = reg_l1, reg_h1 + 1
     kap = kavg(kappar(i-1), kappar(i), dx(1), -1)
     r(i) = r(i) / &
          (kap * max(er(i-1), er(i), tiny))
  enddo
end subroutine scgrd1

subroutine scgrd2(r, &
                  DIMS(rbox), &
                  DIMS(reg), &
                  n, kappar, DIMS(kbox), &
                  er, dx) bind(C, name="scgrd2")
  implicit none
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real*8 :: r(DIMV(rbox))
  real*8 :: kappar(DIMV(kbox))
  real*8 :: er(DIMV(kbox))
  real*8 :: dx(1)
  real*8 :: kavg
  integer :: i
  real*8 :: kap
  ! gradient assembly:
  do i = reg_l1, reg_h1 + 1
     r(i) = abs(er(i) - er(i-1)) / dx(1)
  enddo
  i = reg_l1
  if (er(i-1) == -1.d0) then
     r(i) = abs(er(i+1) - er(i)) / dx(1)
  endif
  i = reg_h1 + 1
  if (er(i) == -1.d0) then
     r(i) = abs(er(i-1) - er(i-2)) / dx(1)
  endif
  ! construct scaled gradient:
  do i = reg_l1, reg_h1 + 1
     kap = kavg(kappar(i-1), kappar(i), dx(1), -1)
     r(i) = r(i) / &
          (kap * max(er(i-1), er(i), tiny))
  enddo
end subroutine scgrd2

subroutine scgrd3(r, &
                  DIMS(rbox), &
                  DIMS(reg), &
                  n, kappar, DIMS(kbox), &
                  er, dx) bind(C, name="scgrd3")
  implicit none
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real*8 :: r(DIMV(rbox))
  real*8 :: kappar(DIMV(kbox))
  real*8 :: er(DIMV(kbox))
  real*8 :: dx(1)
  real*8 :: kavg
  integer :: i
  real*8 :: kap
  ! gradient assembly:
  do i = reg_l1, reg_h1 + 1
     r(i) = abs(er(i) - er(i-1)) / dx(1)
  enddo
  i = reg_l1
  if (er(i-1) == -1.d0) then
     r(i) = abs(er(i+1) - er(i)) / dx(1)
  endif
  i = reg_h1 + 1
  if (er(i) == -1.d0) then
     r(i) = abs(er(i-1) - er(i-2)) / dx(1)
  endif
  ! construct scaled gradient:
  do i = reg_l1, reg_h1 + 1
     kap = kavg(kappar(i-1), kappar(i), dx(1), -1)
     r(i) = r(i) / &
          (kap * max(er(i-1), er(i), tiny))
  enddo
end subroutine scgrd3

subroutine lrhs(rhs, &
                DIMS(rbox), &
                DIMS(reg), &
                temp, fkp, eta, etainv, frhoem, frhoes, dfo, &
                ero, DIMS(ebox), edot, &
                r, s, dt, sigma, c, theta) bind(C, name="lrhs")
  implicit none
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
  real*8 :: s(1)
  real*8 :: dt, sigma, c, theta
  integer :: i
  real*8 :: dtm, ek, bs, es, ekt
  dtm = 1.d0 / dt
  do i = reg_l1, reg_h1
     ek = fkp(i) * eta(i)
     bs = etainv(i) * &
          &         4.d0 * sigma * fkp(i) * temp(i)**4
     es = eta(i) * (frhoem(i) - frhoes(i))
     ekt = (1.d0 - theta) * eta(i)
     rhs(i) = (rhs(i) + r(i) * &
          (bs + dtm * (ero(i) + es) + &
          ek * c * edot(i) - &
          ekt * dfo(i))) / &
          (1.d0 - ekt)
  enddo
end subroutine lrhs

subroutine anatw2(test, &
                  DIMS(reg), &
                  temp, p, xf, Tc, dx, xlo, lo) bind(C, name="anatw2")
  implicit none
  integer :: DIMDEC(reg)
  real*8 :: test(DIMV(reg), 0:1)
  real*8 :: temp(DIMV(reg))
  real*8 :: p, xf, Tc, dx(1), xlo(1)
  integer :: lo(1)
  integer :: i
  real*8 :: x, r2
  do i = reg_l1, reg_h1
     x  = xlo(1) + dx(1) * ((i-lo(1)) + 0.5d0)
     r2 = x*x
     test(i,0) = Tc * max((1.d0-r2/xf**2), 0.d0)**(1.d0/p)
     test(i,1) = temp(i) - test(i,0)
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
  integer :: i
  do i = reg_l1, reg_h1
     !         kin = 0.5d0 * (state(i,XMOM) ** 2) /
     !     @                 state(i,DEN)
     !         frhoe(i) = state(i,EDEN) - kin
     frhoe(i) = state(i,UEINT)
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
  real*8 :: state(DIMV(sb), NVAR)
  real*8 :: alpha, teff, ex, frhoal
  integer :: i
  if (en(0) >= 1.d0) then
     print *, "Bad exponent for cv calculation"
     stop
  endif
  ex = 1.d0 / (1.d0 - en(0))
  do i = reg_l1, reg_h1
     if (em(0) == 0.d0) then
        alpha = const(0)
     else
        alpha = const(0) * state(i,URHO) ** em(0)
     endif
     frhoal = state(i, URHO) * alpha + tiny
     if (en(0) == 0.d0) then
        temp(i) = temp(i) / frhoal
     else
        teff = max(temp(i), tiny)
        temp(i) = ((1.d0 - en(0)) * teff / frhoal) ** ex
     endif
  enddo
end subroutine gtemp

! temp contains temp on input:

subroutine gcv(DIMS(reg), &
               cv, DIMS(cbox), &
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
  real*8 :: state(DIMV(sbox),  NVAR)
  real*8 :: alpha, teff, frhoal
  integer :: i
  do i = reg_l1, reg_h1
     if (em(0) == 0.d0) then
        alpha = const(0)
     else
        alpha = const(0) * state(i, URHO) ** em(0)
     endif
     frhoal = state(i, URHO) * alpha + tiny
     if (en(0) == 0.d0) then
        cv(i) = alpha
     else
        teff = max(temp(i), tiny)
        teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
        cv(i) = alpha * teff ** (-en(0))
     endif
  enddo
end subroutine gcv

! exch contains temp on input:

subroutine cexch(DIMS(reg), &
                 exch, DIMS(xbox), &
                 er  , DIMS(ebox), &
                 fkp , DIMS(kbox), &
                 sigma, c) bind(C, name="cexch")
  implicit none
  integer :: DIMDEC(reg)
  integer :: DIMDEC(xbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(kbox)
  real*8 :: exch(DIMV(xbox))
  real*8 :: er  (DIMV(ebox))
  real*8 :: fkp (DIMV(kbox))
  real*8 :: sigma, c
  integer :: i
  do i = reg_l1, reg_h1
     exch(i) = fkp(i) * &
          (4.d0 * sigma * exch(i)**4 &
          - c * er(i))
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
  implicit none
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
  integer :: i
  fac1 = 16.d0 * sigma * dtime
  if (lagpla == 0) then
     fac0 = 0.25d0 * fac1 / dtemp
     fac2 = dtime * c / dtemp
  endif
  do i = reg_l1, reg_h1
     if (lagpla /= 0) then
        ! assume eta and fkp are the same
        d = fac1 * fkp(i) * temp(i) ** 3
     else
        d = fac0 * (eta(i) * (temp(i) + dtemp) ** 4 - &
             fkp(i) * (temp(i)        ) ** 4) - &
             fac2 * (eta(i) - fkp(i)) * er(i)
        ! alternate form, sometimes worse, sometimes better:
        !            d = fac1 * fkp(i) * temp(i) ** 3 +
        !     @          fac0 * (eta(i) - fkp(i)) * temp(i) ** 4 -
        !     @          fac2 * (eta(i) - fkp(i)) * er(i)
        ! analytic derivatives for specific test problem:
        !            d = (1.2d+6 * sigma * temp(i) ** 2 +
        !     @           1.d+5 * c * er(i) * (temp(i) + tiny) ** (-2)) * dtime
        ! another alternate form (much worse):
        !            d = fac1 * fkp(i) * (temp(i) + dtemp) ** 3 +
        !     @          fac0 * (eta(i) - fkp(i)) * (temp(i) + dtemp) ** 4 -
        !     @          fac2 * (eta(i) - fkp(i)) * er(i)
     endif
     frc = frho(i) * cv(i) + tiny
     eta(i) = d / (d + frc)
     etainv(i) = underr * frc / (d + frc)
     eta(i) = 1.d0 - etainv(i)
     !         eta(i) = 1.d0 - underr * (1.d0 - eta(i))
  enddo
end subroutine ceta2

subroutine ceup(DIMS(reg), relres, absres, &
                frhoes, DIMS(grd), &
                frhoem, eta, etainv, dfo, dfn, exch, &
                dt, theta) bind(C, name="ceup")
  implicit none
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
  integer :: i
  do i = reg_l1, reg_h1
     chg = 0.d0
     tot = 0.d0
     tmp = eta(i) * frhoes(i) + &
          etainv(i) * &
          (frhoem(i) - &
          dt * ((1.d0 - theta) * &
          (dfo(i) - dfn(i)) + &
          exch(i)))
     chg = abs(tmp - frhoes(i))
     tot = abs(frhoes(i))
     frhoes(i) = tmp
     absres = max(absres, chg)
     relres = max(relres, chg / (tot + tiny))
  enddo
end subroutine ceup

subroutine ceupdterm(DIMS(reg), relres, absres, &
                     frhoes, DIMS(grd), &
                     frhoem, eta, etainv, dfo, dfn, exch, dterm, &
                     dt, theta) bind(C, name="ceupdterm")
  implicit none
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
  integer :: i
  do i = reg_l1, reg_h1
     chg = 0.d0
     tot = 0.d0
     tmp = eta(i) * frhoes(i) + &
          etainv(i) * &
          (frhoem(i) - &
          dt * ((1.d0 - theta) * &
          (dfo(i) - dfn(i)) + &
          exch(i))) &
          + dt * dterm(i)
     chg = abs(tmp - frhoes(i))
     tot = abs(frhoes(i))
     frhoes(i) = tmp
     absres = max(absres, chg)
     relres = max(relres, chg / (tot + tiny))
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
  real*8 :: state(DIMV(sb),  NVAR)
  real*8 :: sigma, c, dt, theta, relres, absres
  real*8 :: tmp, chg, tot, exch, b, db, dbdt, frhocv
  integer :: i
  do i = reg_l1, reg_h1
     chg = 0.d0
     tot = 0.d0
     frhocv = state(i, URHO) * cv(i)
     dbdt = 16.d0 * sigma * temp(i)**3
     b = 4.d0 * sigma * temp(i)**4
     exch = fkp(i) * (b - c * er(i))
     tmp = eta(i) * frhoes(i) + etainv(i) * &
          (frhoem(i) - &
          dt * ((1.d0 - theta) * &
          (dfo(i) - dfn(i)) + &
          exch))
#if 1
     if (frhocv > tiny .AND. tmp > frhoes(i)) then
        db = (tmp - frhoes(i)) * dbdt / frhocv
        if (b + db <= 0.d0) then
           print *, i, b, db, b+db
        endif
        tmp = ((b + db) / (4.d0 * sigma))**0.25d0
        tmp = frhoes(i) + frhocv * (tmp - temp(i))
     endif
#endif
     chg = abs(tmp - frhoes(i))
     tot = abs(frhoes(i))
     frhoes(i) = tmp
     absres = max(absres, chg)
     relres = max(relres, chg / (tot + tiny))
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
  integer :: i
  do i = reg_l1, reg_h1
     !         kin = 0.5d0 * (state(i,XMOM) ** 2) /
     !     @                 state(i,DEN)
     kin = state(i, UEDEN) - state(i, UEINT)
     state(i, UEINT) = frhoe(i)
     state(i, UEDEN) = frhoe(i) + kin
  enddo
end subroutine cetot

subroutine fkpn(DIMS(reg), &
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
  real*8 :: state(DIMV(sb),  NVAR)
  real*8 :: teff
  integer :: i
  do i = reg_l1, reg_h1
     teff = max(temp(i), tiny)
     teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
     fkp(i) = const(0) * &
          (state(i, URHO) ** em(0)) * &
          (teff ** (-en(0))) * &
          (nu ** (ep(0)))
  enddo
end subroutine fkpn

subroutine rosse1( DIMS(reg), &
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
  real*8 :: state(DIMV(sb),  NVAR)
  real*8 :: kfloor
  real*8 :: kf, teff
  integer :: i
  do i = reg_l1, reg_h1
     teff = max(temp(i), tiny)
     teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
     kf = const(0) * &
          (state(i, URHO) ** em(0)) * &
          (teff ** (-en(0))) * &
          (nu ** (ep(0)))
     kappar(i) = max(kf, kfloor)
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
  real*8 :: state(DIMV(sb),  NVAR)
  real*8 :: kfloor
  real*8 :: kf, teff, sct
  integer :: i
  do i = reg_l1, reg_h1
     teff = max(temp(i), tiny)
     teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
     kf = const(0) * &
          (state(i, URHO) ** em(0)) * &
          (teff ** (-en(0))) * &
          (nu ** (ep(0)))
     sct = sconst(0) * &
          (state(i, URHO) ** sem(0)) * &
          (teff ** (-sen(0))) * &
          (nu ** (sep(0)))
     kappar(i) = max(kf+sct, kfloor)
  enddo
end subroutine rosse1s

subroutine nfloor(dest, &
                  DIMS(dbox), &
                  DIMS(reg), &
                  nflr, flr, nvar) bind(C, name="nfloor")
  implicit none
  integer :: DIMDEC(dbox)
  integer :: DIMDEC(reg)
  integer :: nvar, nflr
  real*8 :: dest(DIMV(dbox), 0:nvar-1)
  real*8 :: flr
  integer :: i, n
  nflr = 0
  do n = 0, nvar-1
     do i = reg_l1, reg_h1
        if (dest(i,n) < flr) then
           dest(i,n) = flr
           nflr = nflr + 1
        endif
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
  implicit none
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  real*8 :: a(DIMV(abox))
  real*8 :: kappa(DIMV(kbox))
  real*8 :: r(reg_l1:reg_h1)
  real*8 :: s(1)
  real*8 :: dt, c
  integer :: i
  do i = reg_l1, reg_h1
     a(i) = c*kappa(i) + 1.d0/dt
     a(i) = r(i) * a(i)
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
  implicit none
  integer :: DIMDEC(fbox)
  integer :: DIMDEC(cbox)
  real*8 :: fine(DIMV(fbox))
  real*8 :: crse(DIMV(cbox))
  integer :: idim, irat(0:0)
  fine(fbox_l1) = crse(cbox_l1)
end subroutine rfface

subroutine bextrp(f, fbox_l1, fbox_h1, reg_l1, reg_h1) bind(C, name="bextrp")
  implicit none
  integer :: fbox_l1, fbox_h1
  integer ::  reg_l1,  reg_h1
  real*8 :: f(fbox_l1:fbox_h1)
  integer :: i

  !     i direction first:
  i = reg_l1
  f(i-1) = 2.d0 * f(i) - f(i+1)
  i = reg_h1
  f(i+1) = 2.d0 * f(i) - f(i-1)
end subroutine bextrp


subroutine lbcoefna(bcoef, &
                    bcgrp, bbox_l1, bbox_h1, &
                    reg_l1, reg_h1, &
                    spec, sbox_l1, sbox_h1, &
                    idim) bind(C, name="lbcoefna")
  implicit none
  integer :: idim
  integer ::  reg_l1,  reg_h1
  integer :: bbox_l1, bbox_h1
  integer :: sbox_l1, sbox_h1
  real*8 :: bcoef(bbox_l1:bbox_h1)
  real*8 :: bcgrp(bbox_l1:bbox_h1)
  real*8 :: spec(sbox_l1:sbox_h1)
  integer :: i
  if (idim == 0) then
     do i = reg_l1, reg_h1
        bcoef(i) = bcoef(i) &
             + 0.5d0 * (spec(i-1) + spec(i)) * bcgrp(i)
     enddo
  endif
end subroutine lbcoefna


subroutine ljupna(jnew, jbox_l1, jbox_h1, &
                  reg_l1, reg_h1, &
                  spec, sbox_l1, sbox_h1, &
                  accel, abox_l1, abox_h1, &
                  nTotal) bind(C, name="ljupna")
  implicit none
  integer :: nTotal
  integer ::  reg_l1,  reg_h1
  integer :: jbox_l1, jbox_h1
  integer :: sbox_l1, sbox_h1
  integer :: abox_l1, abox_h1
  real*8 :: jnew(jbox_l1:jbox_h1, 0:nTotal-1)
  real*8 :: spec(sbox_l1:sbox_h1, 0:nTotal-1)
  real*8 :: accel(abox_l1:abox_h1)
  integer :: i, n
  do n = 0, nTotal - 1
     do i = reg_l1, reg_h1
        jnew(i,n) = jnew(i,n) + spec(i,n) * accel(i)
     enddo
  enddo
end subroutine ljupna

end module rad_module
