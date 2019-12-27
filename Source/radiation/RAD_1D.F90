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

subroutine multrs(d, &
                  DIMS(dbox), &
                  DIMS(reg), &
                  r, s) bind(C, name="multrs")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(dbox)
  integer :: DIMDEC(reg)
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(1)
  integer :: i
  do i = reg_l1, reg_h1
     d(i) = d(i) * r(i)
  enddo
end subroutine multrs

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

subroutine lacoef(a, &
                  DIMS(abox), &
                  DIMS(reg), &
                  fkp, eta, etainv, r, s, c, dt, theta) bind(C, name="lacoef")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  real(rt)         :: a(DIMV(abox))
  real(rt)         :: fkp(DIMV(abox))
  real(rt)         :: eta(DIMV(abox))
  real(rt)         :: etainv(DIMV(abox))
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(1)
  real(rt)         :: c, dt, theta
  integer :: i
  real(rt)         :: dtm
  dtm = 1.e0_rt / dt
  do i = reg_l1, reg_h1
     a(i) = r(i) * &
          (fkp(i) * etainv(i) * c + dtm) / &
          (1.e0_rt - (1.e0_rt - theta) * eta(i))
  enddo
end subroutine lacoef

subroutine bclim(b, &
                 lambda, DIMS(bbox), &
                 DIMS(reg), &
                 n, kappar, DIMS(kbox), &
                 r, s, c, dx) bind(C, name="bclim")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: lambda(DIMV(bbox))
  real(rt)         :: kappar(DIMV(kbox))
  real(rt)         :: r(reg_l1:reg_h1+1)
  real(rt)         :: s(1)
  real(rt)         :: c, dx(1)
  real(rt)         :: kavg
  integer :: i
  real(rt)         :: kap
  do i = reg_l1, reg_h1 + 1
     kap = kavg(kappar(i-1), kappar(i), dx(1),-1)
     b(i) = r(i) * c * lambda(i) / kap
  enddo
end subroutine bclim

subroutine flxlim(lambda, &
                  DIMS(rbox), &
                  DIMS(reg), limiter) bind(C, name="flxlim")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: limiter
  real(rt)         :: lambda(DIMV(rbox))
  integer :: i
  do i = reg_l1, reg_h1
     lambda(i) = FLDlambda(lambda(i),limiter)
  enddo
end subroutine flxlim

subroutine eddfac(efact, &
                  DIMS(rbox), &
                  DIMS(reg), limiter, n)
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: n, limiter
  real(rt)         :: efact(DIMV(rbox))
  integer :: i
  real(rt)         :: r, lambda
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
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real(rt)         :: r(DIMV(rbox))
  real(rt)         :: kappar(DIMV(kbox))
  real(rt)         :: er(DIMV(kbox))
  real(rt)         :: dx(1)
  real(rt)         :: kavg
  integer :: i
  real(rt)         :: kap
  ! gradient assembly:
  do i = reg_l1, reg_h1 + 1
     r(i) = abs(er(i) - er(i-1)) / dx(1)
  enddo
  i = reg_l1
  if (er(i-1) == -1.e0_rt) then
     r(i) = abs(er(i+1) - er(i)) / dx(1)
  endif
  i = reg_h1 + 1
  if (er(i) == -1.e0_rt) then
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
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real(rt)         :: r(DIMV(rbox))
  real(rt)         :: kappar(DIMV(kbox))
  real(rt)         :: er(DIMV(kbox))
  real(rt)         :: dx(1)
  real(rt)         :: kavg
  integer :: i
  real(rt)         :: kap
  ! gradient assembly:
  do i = reg_l1, reg_h1 + 1
     r(i) = abs(er(i) - er(i-1)) / dx(1)
  enddo
  i = reg_l1
  if (er(i-1) == -1.e0_rt) then
     r(i) = abs(er(i+1) - er(i)) / dx(1)
  endif
  i = reg_h1 + 1
  if (er(i) == -1.e0_rt) then
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
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(rbox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  integer :: n
  real(rt)         :: r(DIMV(rbox))
  real(rt)         :: kappar(DIMV(kbox))
  real(rt)         :: er(DIMV(kbox))
  real(rt)         :: dx(1)
  real(rt)         :: kavg
  integer :: i
  real(rt)         :: kap
  ! gradient assembly:
  do i = reg_l1, reg_h1 + 1
     r(i) = abs(er(i) - er(i-1)) / dx(1)
  enddo
  i = reg_l1
  if (er(i-1) == -1.e0_rt) then
     r(i) = abs(er(i+1) - er(i)) / dx(1)
  endif
  i = reg_h1 + 1
  if (er(i) == -1.e0_rt) then
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
  use amrex_fort_module, only : rt => amrex_real
  implicit none
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
  real(rt)         :: s(1)
  real(rt)         :: dt, sigma, c, theta
  integer :: i
  real(rt)         :: dtm, ek, bs, es, ekt
  dtm = 1.e0_rt / dt
  do i = reg_l1, reg_h1
     ek = fkp(i) * eta(i)
     bs = etainv(i) * &
          &         4.e0_rt * sigma * fkp(i) * temp(i)**4
     es = eta(i) * (frhoem(i) - frhoes(i))
     ekt = (1.e0_rt - theta) * eta(i)
     rhs(i) = (rhs(i) + r(i) * &
          (bs + dtm * (ero(i) + es) + &
          ek * c * edot(i) - &
          ekt * dfo(i))) / &
          (1.e0_rt - ekt)
  enddo
end subroutine lrhs

subroutine anatw2(test, &
                  DIMS(reg), &
                  temp, p, xf, Tc, dx, xlo, lo) bind(C, name="anatw2")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(reg)
  real(rt)         :: test(DIMV(reg), 0:1)
  real(rt)         :: temp(DIMV(reg))
  real(rt)         :: p, xf, Tc, dx(1), xlo(1)
  integer :: lo(1)
  integer :: i
  real(rt)         :: x, r2
  do i = reg_l1, reg_h1
     x  = xlo(1) + dx(1) * ((i-lo(1)) + 0.5e0_rt)
     r2 = x*x
     test(i,0) = Tc * max((1.e0_rt-r2/xf**2), 0.e0_rt)**(1.e0_rt/p)
     test(i,1) = temp(i) - test(i,0)
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
  integer :: i
  if (en(0) >= 1.e0_rt) then
     print *, "Bad exponent for cv calculation"
     stop
  endif
  ex = 1.e0_rt / (1.e0_rt - en(0))
  do i = reg_l1, reg_h1
     if (em(0) == 0.e0_rt) then
        alpha = const(0)
     else
        alpha = const(0) * state(i,URHO) ** em(0)
     endif
     frhoal = state(i, URHO) * alpha + tiny
     if (en(0) == 0.e0_rt) then
        temp(i) = temp(i) / frhoal
     else
        teff = max(temp(i), tiny)
        temp(i) = ((1.e0_rt - en(0)) * teff / frhoal) ** ex
     endif
  enddo
end subroutine gtemp

! temp contains temp on input:

subroutine gcv(DIMS(reg), &
               cv, DIMS(cbox), &
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
  real(rt)         :: state(DIMV(sbox),  NVAR)
  real(rt)         :: alpha, teff, frhoal
  integer :: i
  do i = reg_l1, reg_h1
     if (em(0) == 0.e0_rt) then
        alpha = const(0)
     else
        alpha = const(0) * state(i, URHO) ** em(0)
     endif
     frhoal = state(i, URHO) * alpha + tiny
     if (en(0) == 0.e0_rt) then
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
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(reg)
  integer :: DIMDEC(xbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(kbox)
  real(rt)         :: exch(DIMV(xbox))
  real(rt)         :: er  (DIMV(ebox))
  real(rt)         :: fkp (DIMV(kbox))
  real(rt)         :: sigma, c
  integer :: i
  do i = reg_l1, reg_h1
     exch(i) = fkp(i) * &
          (4.e0_rt * sigma * exch(i)**4 &
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
  use amrex_fort_module, only : rt => amrex_real
  implicit none
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
  integer :: i
  fac1 = 16.e0_rt * sigma * dtime
  if (lagpla == 0) then
     fac0 = 0.25e0_rt * fac1 / dtemp
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
        !            d = (1.2e+6_rt * sigma * temp(i) ** 2 +
        !     @           1.e+5_rt * c * er(i) * (temp(i) + tiny) ** (-2)) * dtime
        ! another alternate form (much worse):
        !            d = fac1 * fkp(i) * (temp(i) + dtemp) ** 3 +
        !     @          fac0 * (eta(i) - fkp(i)) * (temp(i) + dtemp) ** 4 -
        !     @          fac2 * (eta(i) - fkp(i)) * er(i)
     endif
     frc = frho(i) * cv(i) + tiny
     eta(i) = d / (d + frc)
     etainv(i) = underr * frc / (d + frc)
     eta(i) = 1.e0_rt - etainv(i)
     !         eta(i) = 1.e0_rt - underr * (1.e0_rt - eta(i))
  enddo
end subroutine ceta2

subroutine ceup(DIMS(reg), relres, absres, &
                frhoes, DIMS(grd), &
                frhoem, eta, etainv, dfo, dfn, exch, &
                dt, theta) bind(C, name="ceup")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
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
  integer :: i
  do i = reg_l1, reg_h1
     chg = 0.e0_rt
     tot = 0.e0_rt
     tmp = eta(i) * frhoes(i) + &
          etainv(i) * &
          (frhoem(i) - &
          dt * ((1.e0_rt - theta) * &
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
  use amrex_fort_module, only : rt => amrex_real
  implicit none
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
  integer :: i
  do i = reg_l1, reg_h1
     chg = 0.e0_rt
     tot = 0.e0_rt
     tmp = eta(i) * frhoes(i) + &
          etainv(i) * &
          (frhoem(i) - &
          dt * ((1.e0_rt - theta) * &
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
  real(rt)         :: state(DIMV(sb),  NVAR)
  real(rt)         :: sigma, c, dt, theta, relres, absres
  real(rt)         :: tmp, chg, tot, exch, b, db, dbdt, frhocv
  integer :: i
  do i = reg_l1, reg_h1
     chg = 0.e0_rt
     tot = 0.e0_rt
     frhocv = state(i, URHO) * cv(i)
     dbdt = 16.e0_rt * sigma * temp(i)**3
     b = 4.e0_rt * sigma * temp(i)**4
     exch = fkp(i) * (b - c * er(i))
     tmp = eta(i) * frhoes(i) + etainv(i) * &
          (frhoem(i) - &
          dt * ((1.e0_rt - theta) * &
          (dfo(i) - dfn(i)) + &
          exch))
#if 1
     if (frhocv > tiny .AND. tmp > frhoes(i)) then
        db = (tmp - frhoes(i)) * dbdt / frhocv
        if (b + db <= 0.e0_rt) then
           print *, i, b, db, b+db
        endif
        tmp = ((b + db) / (4.e0_rt * sigma))**0.25e0_rt
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

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(sb)
  integer :: DIMDEC(fb)
  real(rt)         :: state(DIMV(sb), NVAR)
  real(rt)         :: frhoe(DIMV(fb))
  real(rt)         :: kin
  integer :: i
  do i = reg_l1, reg_h1
     !         kin = 0.5e0_rt * (state(i,XMOM) ** 2) /
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

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(fb)
  integer :: DIMDEC(tb)
  integer :: DIMDEC(sb)
  real(rt)         :: fkp(DIMV(fb))
  real(rt)         :: const(0:1), em(0:1), en(0:1), tf(0:1)
  real(rt)         :: ep(0:1), nu
  real(rt)         :: temp(DIMV(tb))
  real(rt)         :: state(DIMV(sb),  NVAR)
  real(rt)         :: teff
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

subroutine nfloor(dest, &
                  DIMS(dbox), &
                  DIMS(reg), &
                  nflr, flr, nvar) bind(C, name="nfloor")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(dbox)
  integer :: DIMDEC(reg)
  integer :: nvar, nflr
  real(rt)         :: dest(DIMV(dbox), 0:nvar-1)
  real(rt)         :: flr
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
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(kbox)
  real(rt)         :: a(DIMV(abox))
  real(rt)         :: kappa(DIMV(kbox))
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: s(1)
  real(rt)         :: dt, c
  integer :: i
  do i = reg_l1, reg_h1
     a(i) = c*kappa(i) + 1.e0_rt/dt
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


subroutine lbcoefna(bcoef, &
                    bcgrp, bbox_l1, bbox_h1, &
                    reg_l1, reg_h1, &
                    spec, sbox_l1, sbox_h1, &
                    idim) bind(C, name="lbcoefna")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: idim
  integer ::  reg_l1,  reg_h1
  integer :: bbox_l1, bbox_h1
  integer :: sbox_l1, sbox_h1
  real(rt)         :: bcoef(bbox_l1:bbox_h1)
  real(rt)         :: bcgrp(bbox_l1:bbox_h1)
  real(rt)         :: spec(sbox_l1:sbox_h1)
  integer :: i
  if (idim == 0) then
     do i = reg_l1, reg_h1
        bcoef(i) = bcoef(i) &
             + 0.5e0_rt * (spec(i-1) + spec(i)) * bcgrp(i)
     enddo
  endif
end subroutine lbcoefna


subroutine ljupna(jnew, jbox_l1, jbox_h1, &
                  reg_l1, reg_h1, &
                  spec, sbox_l1, sbox_h1, &
                  accel, abox_l1, abox_h1, &
                  nTotal) bind(C, name="ljupna")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: nTotal
  integer ::  reg_l1,  reg_h1
  integer :: jbox_l1, jbox_h1
  integer :: sbox_l1, sbox_h1
  integer :: abox_l1, abox_h1
  real(rt)         :: jnew(jbox_l1:jbox_h1, 0:nTotal-1)
  real(rt)         :: spec(sbox_l1:sbox_h1, 0:nTotal-1)
  real(rt)         :: accel(abox_l1:abox_h1)
  integer :: i, n
  do n = 0, nTotal - 1
     do i = reg_l1, reg_h1
        jnew(i,n) = jnew(i,n) + spec(i,n) * accel(i)
     enddo
  enddo
end subroutine ljupna

end module rad_module
