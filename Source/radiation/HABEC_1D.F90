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

subroutine hbvec(vec, &
                 DIMS(reg), &
                 cdir, bct, bho, bcl, &
                 bcval, DIMS(bcv), &
                 mask, DIMS(msk), &
                 b, DIMS(bbox), &
                 beta, dx) bind(C, name="hbvec")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: cdir, bct, bho
  real(rt)         :: bcl, beta, dx(1)
  real(rt)         :: vec(DIMV(reg))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: h, bfv
  real(rt)         :: h2, th2
  integer :: i
  h = dx(1)
  if (bct == LO_DIRICHLET) then
     if (bho >= 1) then
        h2 = 0.5e0_rt * h
        th2 = 3.e0_rt * h2
        bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
     else
        bfv = (beta / h) / (0.5e0_rt * h + bcl)
     endif
  else if (bct == LO_NEUMANN) then
     bfv = beta / h
  else
     print *, "hbvec: unsupported boundary type"
     stop
  endif
  if (cdir == 0) then
     ! Left face of grid
     i = reg_l1
     if (mask(i-1) > 0) then
        vec(i) = vec(i) + bfv * b(i) * bcval(i-1)
     endif
  else if (cdir == 1) then
     ! Right face of grid
     i = reg_h1
     if (mask(i+1) > 0) then
        vec(i) = vec(i) + bfv * b(i+1) * bcval(i+1)
     endif
  else
     print *, "hbvec: impossible face orientation"
  endif
end subroutine hbvec

subroutine hbvec3(vec, &
                  DIMS(reg), &
                  cdir, bctype, tf, bho, bcl, &
                  bcval, DIMS(bcv), &
                  mask, DIMS(msk), &
                  b, DIMS(bbox), &
                  beta, dx, r) bind(C, name="hbvec3")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: cdir, bctype, tf(DIMV(bcv)), bho
  real(rt)         :: bcl, beta, dx(1)
  real(rt)         :: vec(DIMV(reg))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: r(1)
  real(rt)         :: h, bfv, r0
  real(rt)         :: h2, th2
  integer :: i, bct
  h = dx(1)
  ! r is passed as an array but actually has only one element, which is
  ! appropriate for the face we are doing here.
  r0 = r(1)
  if (cdir == 0) then
     ! Left face of grid
     i = reg_l1
     if (mask(i-1) > 0) then
        if (bctype == -1) then
           bct = tf(i-1)
        else
           bct = bctype
        endif
        if (bct == LO_DIRICHLET) then
           if (bho >= 1) then
              h2 = 0.5e0_rt * h
              th2 = 3.e0_rt * h2
              bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
           else
              bfv = (beta / h) / (0.5e0_rt * h + bcl)
           endif
           bfv = bfv * b(i)
        else if (bct == LO_NEUMANN) then
           bfv = beta * r0 / h
        else if (bct == LO_MARSHAK .OR. &
             bct == LO_SANCHEZ_POMRANING) then
           bfv = 2.e0_rt * beta * r0 / h
        else
           print *, "hbvec3: unsupported boundary type"
           stop
        endif
        vec(i) = vec(i) + bfv * bcval(i-1)
     endif
  else if (cdir == 1) then
     ! Right face of grid
     i = reg_h1
     if (mask(i+1) > 0) then
        if (bctype == -1) then
           bct = tf(i+1)
        else
           bct = bctype
        endif
        if (bct == LO_DIRICHLET) then
           if (bho >= 1) then
              h2 = 0.5e0_rt * h
              th2 = 3.e0_rt * h2
              bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
           else
              bfv = (beta / h) / (0.5e0_rt * h + bcl)
           endif
           bfv = bfv * b(i+1)
        else if (bct == LO_NEUMANN) then
           bfv = beta * r0 / h
        else if (bct == LO_MARSHAK .OR. &
             bct == LO_SANCHEZ_POMRANING) then
           bfv = 2.e0_rt * beta * r0 / h
        else
           print *, "hbvec3: unsupported boundary type"
           stop
        endif
        vec(i) = vec(i) + bfv * bcval(i+1)
     endif
  else
     print *, "hbvec3: impossible face orientation"
  endif
end subroutine hbvec3

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
  real(rt)         :: bcl, beta, dx(1)
  real(rt)         :: flux(DIMV(fbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: h, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i
  h = dx(1)
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
     ! Left face of grid
     i = reg_l1
     if (mask(i-1) > 0) then
        flux(i) = b(i) * (bfv * bcval(i-1) - bfm * er(i))
        if (bho >= 1) then
           flux(i) = flux(i) - b(i) * bfm2 * er(i+1)
        endif
     endif
  else if (cdir == 1) then
     ! Right face of grid
     i = reg_h1
     if (mask(i+1) > 0) then
        flux(i+1) = -b(i+1) * (bfv * bcval(i+1) - bfm * er(i))
        if (bho >= 1) then
           flux(i+1) = flux(i+1) + b(i+1) * bfm2 * er(i-1)
        endif
     endif
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
  real(rt)         :: bcl, beta, dx(1), c
  real(rt)         :: flux(DIMV(fbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: spa(DIMV(spabox))
  real(rt)         :: r(1)
  real(rt)         :: h, bfm, bfv, r0
  real(rt)         :: bfm2, h2, th2
  integer :: i, bct
  h = dx(1)
  ! r is passed as an array but actually has only one element, which is
  ! appropriate for the face we are doing here.
  r0 = r(1)
  if (cdir == 0) then
     ! Left face of grid
     i = reg_l1
     if (mask(i-1) > 0) then
        if (bctype == -1) then
           bct = tf(i-1)
        else
           bct = bctype
        endif
        if (bct == LO_DIRICHLET) then
           if (bho >= 1) then
              h2 = 0.5e0_rt * h
              th2 = 3.e0_rt * h2
              bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i)
              bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i)
              bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i)
           else
              bfv = beta / (0.5e0_rt * h + bcl) * b(i)
              bfm = bfv
           endif
        else if (bct == LO_NEUMANN) then
           bfv  = beta * r0
           bfm  = 0.e0_rt
           bfm2 = 0.e0_rt
        else if (bct == LO_MARSHAK) then
           bfv = 2.e0_rt * beta * r0
           if (bho >= 1) then
              bfm  =  0.375e0_rt * c * bfv
              bfm2 = -0.125e0_rt * c * bfv
           else
              bfm = 0.25e0_rt * c * bfv
           endif
        else if (bct == LO_SANCHEZ_POMRANING) then
           bfv = 2.e0_rt * beta * r0
           if (bho >= 1) then
              bfm  =  1.5e0_rt * spa(i) * c * bfv
              bfm2 = -0.5e0_rt * spa(i) * c * bfv
           else
              bfm = spa(i) * c * bfv
           endif
        else
           print *, "hbflx3: unsupported boundary type"
           stop
        endif
        if (inhom == 0) then
           bfv = 0.e0_rt
        endif
        flux(i) = (bfv * bcval(i-1) - bfm * er(i))
        if (bho >= 1) then
           flux(i) = flux(i) - bfm2 * er(i+1)
        endif
     endif
  else if (cdir == 1) then
     ! Right face of grid
     i = reg_h1
     if (mask(i+1) > 0) then
        if (bctype == -1) then
           bct = tf(i+1)
        else
           bct = bctype
        endif
        if (bct == LO_DIRICHLET) then
           if (bho >= 1) then
              h2 = 0.5e0_rt * h
              th2 = 3.e0_rt * h2
              bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i+1)
              bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i+1)
              bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i+1)
           else
              bfv = beta / (0.5e0_rt * h + bcl) * b(i+1)
              bfm = bfv
           endif
        else if (bct == LO_NEUMANN) then
           bfv  = beta * r0
           bfm  = 0.e0_rt
           bfm2 = 0.e0_rt
        else if (bct == LO_MARSHAK) then
           bfv = 2.e0_rt * beta * r0
           if (bho >= 1) then
              bfm  =  0.375e0_rt * c * bfv
              bfm2 = -0.125e0_rt * c * bfv
           else
              bfm = 0.25e0_rt * c * bfv
           endif
        else if (bct == LO_SANCHEZ_POMRANING) then
           bfv = 2.e0_rt * beta * r0
           if (bho >= 1) then
              bfm  =  1.5e0_rt * spa(i) * c * bfv
              bfm2 = -0.5e0_rt * spa(i) * c * bfv
           else
              bfm = spa(i) * c * bfv
           endif
        else
           print *, "hbflx3: unsupported boundary type"
           stop
        endif
        if (inhom == 0) then
           bfv = 0.e0_rt
        endif
        flux(i+1) = -(bfv * bcval(i+1) - bfm * er(i))
        if (bho >= 1) then
           flux(i+1) = flux(i+1) + bfm2 * er(i-1)
        endif
     endif
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
  real(rt)         :: bcl, dx(1)
  real(rt)         :: dterm(DIMV(dtbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: h
  integer :: i
  h = dx(1)
  if (bct == LO_DIRICHLET) then
     if (cdir == 0) then
        !     Left face of grid
        i = reg_l1
        if (mask(i-1) > 0) then
           dterm(i) = d(i)*(er(i) - bcval(i-1))/(0.5e0_rt*h+bcl)
        endif
     else if (cdir == 1) then
        !     Right face of grid
        i = reg_h1
        if (mask(i+1) > 0) then
           dterm(i+1) = d(i+1)*(bcval(i+1)-er(i))/(0.5e0_rt*h+bcl)
        endif
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
  real(rt)         :: bcl, dx(1)
  real(rt)         :: dterm(DIMV(dtbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: h
  integer :: i, bct
  h = dx(1)
  if (cdir == 0) then
     ! Left face of grid
     i = reg_l1
     if (mask(i-1) > 0) then
        if (bctype == -1) then
           bct = tf(i-1)
        else
           bct = bctype
        endif
        if (bct == LO_DIRICHLET) then
           dterm(i) = d(i)*(er(i) - bcval(i-1))/(0.5e0_rt*h+bcl)
        else if (bct == LO_NEUMANN .AND. bcval(i-1) == 0.e0_rt) then
           dterm(i) = 0.e0_rt
        else
           print *, "hdterm3: unsupported boundary type"
           stop
        endif
     endif
  else if (cdir == 1) then
     ! Right face of grid
     i = reg_h1
     if (mask(i+1) > 0) then
        if (bctype == -1) then
           bct = tf(i+1)
        else
           bct = bctype
        endif
        if (bct == LO_DIRICHLET) then
           dterm(i+1) = d(i+1)*(bcval(i+1)-er(i))/(0.5e0_rt*h+bcl)
        else if (bct == LO_NEUMANN .AND. bcval(i+1) == 0.e0_rt) then
           dterm(i+1) = 0.e0_rt
        else
           print *, "hbterm3: unsupported boundary type"
           stop
        endif
     endif
  else
     print *, "hdterm3: impossible face orientation"
  endif
end subroutine hdterm3

subroutine hmac(mat, a, &
                DIMS(abox), &
                DIMS(reg), &
                alpha) bind(C, name="hmac")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  real(rt)         :: a(DIMV(abox))
  real(rt)         :: mat(0:2, DIMV(reg))
  real(rt)         :: alpha
  integer :: i
  if (alpha == 0.e0_rt) then
     do i = reg_l1, reg_h1
        mat(0,i) = 0.e0_rt
     enddo
  else
     do i = reg_l1, reg_h1
        mat(0,i) = alpha * a(i)
     enddo
  endif
end subroutine hmac

subroutine hmbc(mat, b, &
                DIMS(bbox), &
                DIMS(reg), &
                beta, dx, n) bind(C, name="hmbc")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(reg)
  integer :: n
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: mat(0:2, DIMV(reg))
  real(rt)         :: beta, dx(1)
  real(rt)         :: fac
  integer :: i
  fac = beta / (dx(1)**2)
  do i = reg_l1, reg_h1
     mat(0,i) = mat(0,i) + fac * (b(i) + b(i+1))
     mat(1,i) = - fac * b(i)
     mat(2,i) = - fac * b(i+1)
  enddo
end subroutine hmbc

subroutine hma2c(mat, a2, &
                 DIMS(bbox), &
                 DIMS(reg), &
                 alpha2, n) bind(C, name="hma2c")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(reg)
  integer :: n
  real(rt)         :: a2(DIMV(bbox))
  real(rt)         :: mat(0:2, DIMV(reg))
  real(rt)         :: alpha2
  real(rt)         :: fac
  integer :: i
  fac = 0.25e0_rt * alpha2
  do i = reg_l1, reg_h1
     mat(0,i) = mat(0,i) + fac * (a2(i) + a2(i+1))
     mat(1,i) = mat(1,i) + fac * a2(i)
     mat(2,i) = mat(2,i) + fac * a2(i+1)
  enddo
end subroutine hma2c

subroutine hmcc(mat, c, &
                DIMS(cbox), &
                DIMS(reg), &
                gamma, dx, n) bind(C, name="hmcc")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(cbox)
  integer :: DIMDEC(reg)
  integer :: n
  ! c has two components for upwinding, independent of dimension
  real(rt)         :: c(DIMV(cbox), 0:1)
  real(rt)         :: mat(0:2, DIMV(reg))
  real(rt)         :: gamma, dx(1)
  real(rt)         :: fac
  integer :: i
  fac = gamma / dx(1)
  do i = reg_l1, reg_h1
     mat(0,i) = mat(0,i) - fac * (c(i,1) - c(i+1,0))
     mat(1,i) = mat(1,i) - fac * c(i,0)
     mat(2,i) = mat(2,i) + fac * c(i+1,1)
  enddo
end subroutine hmcc

subroutine add_ccoef_flux(dir, &
                          den, DIMS(den), &
                          c, DIMS(cbox), &
                          gamma, &
                          dx, &
                          flux, DIMS(fbox))
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: DIMDEC(den)
  integer :: DIMDEC(cbox)
  integer :: DIMDEC(fbox)

  real(rt)         :: den(DIMV(den))
  real(rt)         :: c(DIMV(cbox), 0:1)
  real(rt)         :: flux(DIMV(fbox))

  integer :: dir
  real(rt)         :: gamma, dx(1)
  integer :: i

  if (dir == 0) then

     !...     x-direction

     do i = fbox_l1, fbox_h1
        flux(i) = flux(i) + gamma * (c(i,0) * den(i-1) + &
             c(i,1) * den(i))
     enddo

  endif
end subroutine add_ccoef_flux

subroutine hmd1c(mat, d1, &
                 DIMS(abox), &
                 DIMS(reg), &
                 delta1, dx, n) bind(C, name="hmd1c")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  integer :: n
  real(rt)         :: d1(DIMV(abox))
  real(rt)         :: mat(0:2, DIMV(reg))
  real(rt)         :: delta1, dx(1)
  real(rt)         :: fac
  integer :: i
  fac = 0.5e0_rt * delta1 / dx(1)
  do i = reg_l1, reg_h1
     mat(1,i) = mat(1,i) - fac * d1(i)
     mat(2,i) = mat(2,i) + fac * d1(i)
  enddo
end subroutine hmd1c

subroutine hmd2c(mat, d2, &
                 DIMS(bbox), &
                 DIMS(reg), &
                 delta2, dx, n) bind(C, name="hmd2c")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(reg)
  integer :: n
  real(rt)         :: d2(DIMV(bbox))
  real(rt)         :: mat(0:2, DIMV(reg))
  real(rt)         :: delta2, dx(1)
  real(rt)         :: fac
  integer :: i
  fac = 0.5e0_rt * delta2 / dx(1)
  do i = reg_l1, reg_h1
     mat(0,i) = mat(0,i) + fac * (d2(i) - d2(i+1))
     mat(1,i) = mat(1,i) - fac * d2(i)
     mat(2,i) = mat(2,i) + fac * d2(i+1)
  enddo
end subroutine hmd2c

subroutine hmmat(mat, &
                 DIMS(reg), &
                 cdir, bct, bho, bcl, &
                 mask, DIMS(msk), &
                 b, DIMS(bbox), &
                 beta, dx) bind(C, name="hmmat")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: DIMDEC(reg)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: cdir, bct, bho
  real(rt)         :: bcl, beta, dx(1)
  real(rt)         :: mat(0:2, DIMV(reg))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: h, fac, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i
  h = dx(1)
  fac = beta / (h**2)
  if (bct == LO_DIRICHLET) then
     if (bho >= 1) then
        h2 = 0.5e0_rt * h
        th2 = 3.e0_rt * h2
        bfm = fac * (th2 - bcl) / (bcl + h2) - fac
        bfm2 = fac * (bcl - h2) / (bcl + th2)
     else
        bfv = (beta / h) / (0.5e0_rt * h + bcl)
        bfm = bfv - fac
     endif
  else if (bct == LO_NEUMANN) then
     bfm = -fac
     bfm2 = 0.e0_rt
  else
     print *, "hmmat: unsupported boundary type"
     stop
  endif
  if (cdir == 0) then
     ! Left face of grid
     i = reg_l1
     if (mask(i-1) > 0) then
        mat(0,i) = mat(0,i) + bfm * b(i)
        mat(1,i) = 0.e0_rt
        if (bho >= 1) then
           mat(2,i) = mat(2,i) + bfm2 * b(i)
        endif
     endif
  else if (cdir == 1) then
     ! Right face of grid
     i = reg_h1
     if (mask(i+1) > 0) then
        mat(0,i) = mat(0,i) + bfm * b(i+1)
        mat(2,i) = 0.e0_rt
        if (bho >= 1) then
           mat(1,i) = mat(1,i) + bfm2 * b(i+1)
        endif
     endif
  else
     print *, "hmmat: impossible face orientation"
  endif
end subroutine hmmat

subroutine hmmat3(mat, &
                  DIMS(reg), &
                  cdir, bctype, tf, bho, bcl, &
                  DIMS(bcv), &
                  mask, DIMS(msk), &
                  b, DIMS(bbox), &
                  beta, dx, c, r, &
                  spa, DIMS(spabox)) bind(C, name="hmmat3")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(spabox)
  integer :: cdir, bctype, tf(DIMV(bcv)), bho
  real(rt)         :: bcl, beta, dx(1), c
  real(rt)         :: mat(0:2, DIMV(reg))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: spa(DIMV(spabox))
  real(rt)         :: r(1)
  real(rt)         :: h, fac, bfm, bfv, r0
  real(rt)         :: bfm2, h2, th2
  integer :: i, bct
  h = dx(1)
  ! r is passed as an array but actually has only one element, which is
  ! appropriate for the face we are doing here.
  r0 = r(1)
  ! The -fac * b(i) term applied to the matrix diagonal is the contribution
  ! from the interior stencil which must be removed at the boundary.
  fac = beta / (h**2)
  if (cdir == 0) then
     ! Left face of grid
     i = reg_l1
     if (mask(i-1) > 0) then
        if (bctype == -1) then
           bct = tf(i-1)
        else
           bct = bctype
        endif
        if (bct == LO_DIRICHLET) then
           if (bho >= 1) then
              h2 = 0.5e0_rt * h
              th2 = 3.e0_rt * h2
              bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i)
              bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i)
           else
              bfv = (beta / h) / (0.5e0_rt * h + bcl)
              bfm = bfv * b(i)
           endif
        else if (bct == LO_NEUMANN) then
           bfm  = 0.e0_rt
           bfm2 = 0.e0_rt
        else if (bct == LO_MARSHAK) then
           bfv = 2.e0_rt * beta * r0 / h
           if (bho >= 1) then
              bfm  =  0.375e0_rt * c * bfv
              bfm2 = -0.125e0_rt * c * bfv
           else
              bfm = 0.25e0_rt * c * bfv
           endif
        else if (bct == LO_SANCHEZ_POMRANING) then
           bfv = 2.e0_rt * beta * r0 / h
           if (bho >= 1) then
              bfm  =  1.5e0_rt * spa(i) * c * bfv
              bfm2 = -0.5e0_rt * spa(i) * c * bfv
           else
              bfm = spa(i) * c * bfv
           endif
        else
           print *, "hmmat3: unsupported boundary type"
           stop
        endif
        mat(0,i) = mat(0,i) + bfm - fac * b(i)
        mat(1,i) = 0.e0_rt
        if (bho >= 1) then
           mat(2,i) = mat(2,i) + bfm2
        endif
     endif
  else if (cdir == 1) then
     ! Right face of grid
     i = reg_h1
     if (mask(i+1) > 0) then
        if (bctype == -1) then
           bct = tf(i+1)
        else
           bct = bctype
        endif
        if (bct == LO_DIRICHLET) then
           if (bho >= 1) then
              h2 = 0.5e0_rt * h
              th2 = 3.e0_rt * h2
              bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i+1)
              bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i+1)
           else
              bfv = (beta / h) / (0.5e0_rt * h + bcl)
              bfm = bfv * b(i+1)
           endif
        else if (bct == LO_NEUMANN) then
           bfm  = 0.e0_rt
           bfm2 = 0.e0_rt
        else if (bct == LO_MARSHAK) then
           bfv = 2.e0_rt * beta * r0 / h
           if (bho >= 1) then
              bfm  =  0.375e0_rt * c * bfv
              bfm2 = -0.125e0_rt * c * bfv
           else
              bfm = 0.25e0_rt * c * bfv
           endif
        else if (bct == LO_SANCHEZ_POMRANING) then
           bfv = 2.e0_rt * beta * r0 / h
           if (bho >= 1) then
              bfm  =  1.5e0_rt * spa(i) * c * bfv
              bfm2 = -0.5e0_rt * spa(i) * c * bfv
           else
              bfm = spa(i) * c * bfv
           endif
        else
           print *, "hmmat3: unsupported boundary type"
           stop
        endif
        mat(0,i) = mat(0,i) + bfm - fac * b(i+1)
        mat(2,i) = 0.e0_rt
        if (bho >= 1) then
           mat(1,i) = mat(1,i) + bfm2
        endif
     endif
  else
     print *, "hmmat3: impossible face orientation"
  endif
end subroutine hmmat3

end module habec_module
