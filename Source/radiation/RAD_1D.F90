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

end module rad_module
