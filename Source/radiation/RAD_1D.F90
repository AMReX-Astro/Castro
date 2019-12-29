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
