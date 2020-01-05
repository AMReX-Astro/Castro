#include "AMReX_LO_BCTYPES.H"
#include "AMReX_ArrayLim.H"


module habec_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

subroutine hbflx(flux, &
                 DIMS(fbox), &
                 er, DIMS(ebox), &
                 DIMS(reg), &
                 cdir, bct, bho, bcl, &
                 bcval, DIMS(bcv), &
                 mask, DIMS(msk), &
                 b, DIMS(bbox), &
                 beta, dx, inhom) bind(C, name="hbflx")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(fbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: cdir, bct, bho, inhom
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: flux(DIMV(fbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: h, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i, j, k
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  if (bct == LO_DIRICHLET) then
     if (bho >= 1) then
        h2 = 0.5e0_rt * h
        th2 = 3.e0_rt * h2
        bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2))
        bfm = (beta / h) * (th2 - bcl) / (bcl + h2)
        bfm2 = (beta / h) * (bcl - h2) / (bcl + th2)
     else
        bfv = beta / (0.5e0_rt * h + bcl)
        bfm = bfv
     endif
  else
     print *, "hbflx: unsupported boundary type"
     stop
  endif
  if (inhom == 0) then
     bfv = 0.e0_rt
  endif
  if (cdir == 0) then
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              flux(i,j,k) = b(i,j,k) * &
                   (bfv * bcval(i-1,j,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k) = flux(i,j,k) - &
                      b(i,j,k) * bfm2 * er(i+1,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              flux(i+1,j,k) = -b(i+1,j,k) * &
                   (bfv * bcval(i+1,j,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i+1,j,k) = flux(i+1,j,k) + &
                      b(i+1,j,k) * bfm2 * er(i-1,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              flux(i,j,k) = b(i,j,k) * &
                   (bfv * bcval(i,j-1,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k) = flux(i,j,k) - &
                      b(i,j,k) * bfm2 * er(i,j+1,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              flux(i,j+1,k) = -b(i,j+1,k) * &
                   (bfv * bcval(i,j+1,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j+1,k) = flux(i,j+1,k) + &
                      b(i,j+1,k) * bfm2 * er(i,j-1,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              flux(i,j,k) = b(i,j,k) * &
                   (bfv * bcval(i,j,k-1) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k) = flux(i,j,k) - &
                      b(i,j,k) * bfm2 * er(i,j,k+1)
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              flux(i,j,k+1) = -b(i,j,k+1) * &
                   (bfv * bcval(i,j,k+1) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k+1) = flux(i,j,k+1) + &
                      b(i,j,k+1) * bfm2 * er(i,j,k-1)
              endif
           endif
        enddo
     enddo
  else
     print *, "hbflx: impossible face orientation"
  endif
end subroutine hbflx

subroutine hbflx3(flux, &
                  DIMS(fbox), &
                  er, DIMS(ebox), &
                  DIMS(reg), &
                  cdir, bctype, tf, bho, bcl, &
                  bcval, DIMS(bcv), &
                  mask, DIMS(msk), &
                  b, DIMS(bbox), &
                  beta, dx, c, r, inhom, &
                  spa, DIMS(spabox)) bind(C, name="hbflx3")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(fbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(spabox)
  integer :: cdir, bctype, tf(DIMV(bcv)), bho, inhom
  real(rt)         :: bcl, beta, dx(3), c
  real(rt)         :: flux(DIMV(fbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: spa(DIMV(spabox))
  real(rt)         :: r(1)
  real(rt)         :: h, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i, j, k, bct
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  if (cdir == 0) then
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i-1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = beta / (0.5e0_rt * h + bcl) * b(i,j,k)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  0.75e0_rt * beta * c
                    bfm2 = -0.25e0_rt * beta * c
                 else
                    bfm = 0.5e0_rt * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  3.0e0_rt * spa(i,j,k) * beta * c
                    bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0e0_rt * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.e0_rt
              endif
              flux(i,j,k) = (bfv * bcval(i-1,j,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k) = flux(i,j,k) - bfm2 * er(i+1,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i+1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i+1,j,k)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i+1,j,k)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i+1,j,k)
                 else
                    bfv = beta / (0.5e0_rt * h + bcl) * b(i+1,j,k)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  0.75e0_rt * beta * c
                    bfm2 = -0.25e0_rt * beta * c
                 else
                    bfm = 0.5e0_rt * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  3.0e0_rt * spa(i,j,k) * beta * c
                    bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0e0_rt * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.e0_rt
              endif
              flux(i+1,j,k) = -(bfv * bcval(i+1,j,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i+1,j,k) = flux(i+1,j,k) + bfm2 * er(i-1,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j-1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = beta / (0.5e0_rt * h + bcl) * b(i,j,k)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  0.75e0_rt * beta * c
                    bfm2 = -0.25e0_rt * beta * c
                 else
                    bfm = 0.5e0_rt * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  3.0e0_rt * spa(i,j,k) * beta * c
                    bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0e0_rt * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.e0_rt
              endif
              flux(i,j,k) = (bfv * bcval(i,j-1,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k) = flux(i,j,k) - bfm2 * er(i,j+1,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j+1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j+1,k)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j+1,k)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j+1,k)
                 else
                    bfv = beta / (0.5e0_rt * h + bcl) * b(i,j+1,k)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  0.75e0_rt * beta * c
                    bfm2 = -0.25e0_rt * beta * c
                 else
                    bfm = 0.5e0_rt * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  3.0e0_rt * spa(i,j,k) * beta * c
                    bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0e0_rt * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.e0_rt
              endif
              flux(i,j+1,k) = -(bfv * bcval(i,j+1,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j+1,k) = flux(i,j+1,k) + bfm2 * er(i,j-1,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k-1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = beta / (0.5e0_rt * h + bcl) * b(i,j,k)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  0.75e0_rt * beta * c
                    bfm2 = -0.25e0_rt * beta * c
                 else
                    bfm = 0.5e0_rt * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  3.0e0_rt * spa(i,j,k) * beta * c
                    bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0e0_rt * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.e0_rt
              endif
              flux(i,j,k) = (bfv * bcval(i,j,k-1) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k) = flux(i,j,k) - bfm2 * er(i,j,k+1)
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k+1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k+1)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k+1)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k+1)
                 else
                    bfv = beta / (0.5e0_rt * h + bcl) * b(i,j,k+1)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  0.75e0_rt * beta * c
                    bfm2 = -0.25e0_rt * beta * c
                 else
                    bfm = 0.5e0_rt * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta
                 if (bho >= 1) then
                    bfm  =  3.0e0_rt * spa(i,j,k) * beta * c
                    bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0e0_rt * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.e0_rt
              endif
              flux(i,j,k+1) = -(bfv * bcval(i,j,k+1) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k+1) = flux(i,j,k+1) + bfm2 * er(i,j,k-1)
              endif
           endif
        enddo
     enddo
  else
     print *, "hbflx3: impossible face orientation"
  endif
end subroutine hbflx3

subroutine hdterm(dterm, &
                  DIMS(dtbox), &
                  er, DIMS(ebox), &
                  DIMS(reg), &
                  cdir, bct, bcl, &
                  bcval, DIMS(bcv), &
                  mask, DIMS(msk), &
                  d, DIMS(dbox), &
                  dx) bind(C, name="hdterm")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(dtbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(dbox)
  integer :: cdir, bct
  real(rt)         :: bcl, dx(3)
  real(rt)         :: dterm(DIMV(dtbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: h, bfm, bfv
  integer :: i, j, k
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  if (bct == LO_DIRICHLET) then
     if (cdir == 0) then
        i = reg_l1
        do k = reg_l3, reg_h3
           do j = reg_l2, reg_h2
              if (mask(i-1,j,k) > 0) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i-1,j,k)) &
                      / (0.5e0_rt*h + bcl)
              endif
           enddo
        enddo
     else if (cdir == 3) then
        i = reg_h1
        do k = reg_l3, reg_h3
           do j = reg_l2, reg_h2
              if (mask(i+1,j,k) > 0) then
                 dterm(i+1,j,k) = d(i+1,j,k) * &
                      (bcval(i+1,j,k) - er(i,j,k)) &
                      / (0.5e0_rt*h + bcl)
              endif
           enddo
        enddo
     else if (cdir == 1) then
        j = reg_l2
        do k = reg_l3, reg_h3
           do i = reg_l1, reg_h1
              if (mask(i,j-1,k) > 0) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i,j-1,k)) &
                      / (0.5e0_rt*h + bcl)
              endif
           enddo
        enddo
     else if (cdir == 4) then
        j = reg_h2
        do k = reg_l3, reg_h3
           do i = reg_l1, reg_h1
              if (mask(i,j+1,k) > 0) then
                 dterm(i,j+1,k) = d(i,j+1,k) * &
                      (bcval(i,j+1,k) - er(i,j,k)) &
                      / (0.5e0_rt*h + bcl)
              endif
           enddo
        enddo
     else if (cdir == 2) then
        k = reg_l3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              if (mask(i,j,k-1) > 0) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i,j,k-1)) &
                      / (0.5e0_rt*h + bcl)
              endif
           enddo
        enddo
     else if (cdir == 5) then
        k = reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              if (mask(i,j,k+1) > 0) then
                 dterm(i,j,k+1) = d(i,j,k+1) * &
                      (bcval(i,j,k+1) - er(i,j,k)) &
                      / (0.5e0_rt*h + bcl)
              endif
           enddo
        enddo
     else
        print *, "hdterm: impossible face orientation"
     endif
  else
     print *, "hdterm: unsupported boundary type"
     stop
  endif
end subroutine hdterm

subroutine hdterm3(dterm, &
                   DIMS(dtbox), &
                   er, DIMS(ebox), &
                   DIMS(reg), &
                   cdir, bctype, tf, bcl, &
                   bcval, DIMS(bcv), &
                   mask, DIMS(msk), &
                   d, DIMS(dbox), &
                   dx) bind(C, name="hdterm3")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(dtbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(dbox)
  integer :: cdir, bctype, tf(DIMV(bcv))
  real(rt)         :: bcl, dx(3)
  real(rt)         :: dterm(DIMV(dtbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: h, bfm, bfv
  integer :: i, j, k, bct
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  if (cdir == 0) then
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i-1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i-1,j,k)) &
                      / (0.5e0_rt*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i-1,j,k) == 0.e0_rt) then
                 dterm(i,j,k) = 0.e0_rt
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i+1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i+1,j,k) = d(i+1,j,k) * &
                      (bcval(i+1,j,k) - er(i,j,k)) &
                      / (0.5e0_rt*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i+1,j,k) == 0.e0_rt) then
                 dterm(i+1,j,k) = 0.e0_rt
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j-1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i,j-1,k)) &
                      / (0.5e0_rt*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i,j-1,k) == 0.e0_rt) then
                 dterm(i,j,k) = 0.e0_rt
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j+1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i,j+1,k) = d(i,j+1,k) * &
                      (bcval(i,j+1,k) - er(i,j,k)) &
                      / (0.5e0_rt*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i,j+1,k) == 0.e0_rt) then
                 dterm(i,j+1,k) = 0.e0_rt
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k-1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i,j,k-1)) &
                      / (0.5e0_rt*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i,j,k-1) == 0.e0_rt) then
                 dterm(i,j,k) = 0.e0_rt
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k+1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i,j,k+1) = d(i,j,k+1) * &
                      (bcval(i,j,k+1) - er(i,j,k)) &
                      / (0.5e0_rt*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i,j,k+1) == 0.e0_rt) then
                 dterm(i,j,k+1) = 0.e0_rt
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else
     print *, "hdterm3: impossible face orientation"
  endif
end subroutine hdterm3

end module habec_module
