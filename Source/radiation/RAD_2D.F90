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

subroutine ceupdterm( DIMS(reg), relres, absres, &
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
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        chg = 0.e0_rt
        tot = 0.e0_rt
        tmp = eta(i,j) * frhoes(i,j) + &
             etainv(i,j) * &
             (frhoem(i,j) - &
             dt * ((1.e0_rt - theta) * &
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
  integer :: i, j
  do j = reg_l2, reg_h2
     do i = reg_l1, reg_h1
        chg = 0.e0_rt
        tot = 0.e0_rt
        frhocv = state(i,j, URHO) * cv(i,j)
        dbdt = 16.e0_rt * sigma * temp(i,j)**3
        b = 4.e0_rt * sigma * temp(i,j)**4
        exch = fkp(i,j) * (b - c * er(i,j))
        tmp = eta(i,j) * frhoes(i,j) + etainv(i,j) * &
             (frhoem(i,j) - &
             dt * ((1.e0_rt - theta) * &
             (dfo(i,j) - dfn(i,j)) + &
             exch))
#if 1
        if (frhocv > tiny .AND. tmp > frhoes(i,j)) then
           db = (tmp - frhoes(i,j)) * dbdt / frhocv
           if (b + db <= 0.e0_rt) then
              print *, i, j, b, db, b+db
           endif
           tmp = ((b + db) / (4.e0_rt * sigma))**0.25e0_rt
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

  use amrex_fort_module, only : rt => amrex_real
  integer :: fbox_l1, fbox_l2, fbox_h1, fbox_h2
  integer ::  reg_l1,  reg_l2,  reg_h1,  reg_h2
  real(rt)         :: f(fbox_l1:fbox_h1,fbox_l2:fbox_h2)
  integer :: i, j

  !     i direction first:
  do j = reg_l2, reg_h2
     i = reg_l1
     f(i-1,j) = 2.e0_rt * f(i,j) - f(i+1,j)
     i = reg_h1
     f(i+1,j) = 2.e0_rt * f(i,j) - f(i-1,j)
  enddo

  !     j direction second, including corners:
  do i = reg_l1 - 1, reg_h1 + 1
     j = reg_l2
     f(i,j-1) = 2.e0_rt * f(i,j) - f(i,j+1)
     j = reg_h2
     f(i,j+1) = 2.e0_rt * f(i,j) - f(i,j-1)
  enddo

  !  corner results are the same whichever direction we extrapolate first
end subroutine bextrp

end module rad_module
