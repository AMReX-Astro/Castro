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

end module rad_module

