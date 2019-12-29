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

subroutine hbmat(mat, &
                 DIMS(reg), &
                 cdir, bct, bcl, &
                 mask, DIMS(msk), &
                 b, DIMS(bbox), &
                 beta, dx) bind(C, name="hbmat")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: cdir, bct
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: mat(0:3, DIMV(reg))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: h, fac, bfm, bfv
  integer :: i, j, k
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  fac = beta / (h**2)
  if (bct == LO_DIRICHLET) then
     bfv = fac * h / (0.5e0_rt * h + bcl)
     bfm = bfv - fac
  else if (bct == LO_NEUMANN) then
     bfv = beta / h
     bfm = -fac
  else
     print *, "hbmat: unsupported boundary type"
     stop
  endif
  if (cdir == 0) then
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i,j,k)
              mat(0,i,j,k) = 0.e0_rt
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i+1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i,j,k)
              mat(1,i,j,k) = 0.e0_rt
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i,j+1,k)
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i,j,k)
              mat(2,i,j,k) = 0.e0_rt
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i,j,k+1)
           endif
        enddo
     enddo
  else
     print *, "hbmat: impossible face orientation"
  endif
end subroutine hbmat

subroutine hbmat3(mat, &
                  DIMS(reg), &
                  cdir, bctype, tf, bcl, &
                  DIMS(bcv), &
                  mask, DIMS(msk), &
                  b, DIMS(bbox), &
                  beta, dx, c, r, &
                  spa, DIMS(spabox)) bind(C, name="hbmat3")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(spabox)
  integer :: cdir, bctype, tf(DIMV(bcv))
  real(rt)         :: bcl, beta, dx(3), c
  real(rt)         :: mat(0:3, DIMV(reg))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: spa(DIMV(spabox))
  real(rt)         :: r(1)
  real(rt)         :: h, fac, bfm, bfv
  integer :: i, j, k, bct
  ! The -fac * b(i,j,k) term applied to the matrix diagonal is the contribution
  ! from the interior stencil which must be removed at the boundary.
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  fac = beta / (h**2)
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
                 bfv = fac * h / (0.5e0_rt * h + bcl)
                 bfm = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = 0.25e0_rt * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = spa(i,j,k) * bfv
              else
                 print *, "hbmat3: unsupported boundary type"
                 stop
              endif
              mat(3,i,j,k) = mat(3,i,j,k) + bfm - fac * b(i,j,k)
              mat(0,i,j,k) = 0.e0_rt
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
                 bfv = fac * h / (0.5e0_rt * h + bcl)
                 bfm = bfv * b(i+1,j,k)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = 0.25e0_rt * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = spa(i,j,k) * bfv
              else
                 print *, "hbmat3: unsupported boundary type"
                 stop
              endif
              mat(3,i,j,k) = mat(3,i,j,k) + bfm - fac * b(i+1,j,k)
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
                 bfv = fac * h / (0.5e0_rt * h + bcl)
                 bfm = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = 0.25e0_rt * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = spa(i,j,k) * bfv
              else
                 print *, "hbmat3: unsupported boundary type"
                 stop
              endif
              mat(3,i,j,k) = mat(3,i,j,k) + bfm - fac * b(i,j,k)
              mat(1,i,j,k) = 0.e0_rt
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
                 bfv = fac * h / (0.5e0_rt * h + bcl)
                 bfm = bfv * b(i,j+1,k)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = 0.25e0_rt * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = spa(i,j,k) * bfv
              else
                 print *, "hbmat3: unsupported boundary type"
                 stop
              endif
              mat(3,i,j,k) = mat(3,i,j,k) + bfm - fac * b(i,j+1,k)
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
                 bfv = fac * h / (0.5e0_rt * h + bcl)
                 bfm = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = 0.25e0_rt * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = spa(i,j,k) * bfv
              else
                 print *, "hbmat3: unsupported boundary type"
                 stop
              endif
              mat(3,i,j,k) = mat(3,i,j,k) + bfm - fac * b(i,j,k)
              mat(2,i,j,k) = 0.e0_rt
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
                 bfv = fac * h / (0.5e0_rt * h + bcl)
                 bfm = bfv * b(i,j,k+1)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = 0.25e0_rt * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 bfm = spa(i,j,k) * bfv
              else
                 print *, "hbmat3: unsupported boundary type"
                 stop
              endif
              mat(3,i,j,k) = mat(3,i,j,k) + bfm - fac * b(i,j,k+1)
           endif
        enddo
     enddo
  else
     print *, "hbmat3: impossible face orientation"
  endif
end subroutine hbmat3

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
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: vec(DIMV(reg))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: h, bfv
  real(rt)         :: h2, th2
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
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k) * bcval(i-1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i+1,j,k) * bcval(i+1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k) * bcval(i,j-1,k)
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j+1,k) * bcval(i,j+1,k)
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k) * bcval(i,j,k-1)
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k+1) * bcval(i,j,k+1)
           endif
        enddo
     enddo
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
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: vec(DIMV(reg))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: r(1)
  real(rt)         :: h, bfv
  real(rt)         :: h2, th2
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
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i-1,j,k)
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
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i+1,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i+1,j,k)
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
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j-1,k)
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
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i,j+1,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j+1,k)
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
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j,k-1)
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
                    bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k+1)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j,k+1)
           endif
        enddo
     enddo
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

subroutine hmac(mat, a, &
                DIMS(abox), &
                DIMS(reg), &
                alpha) bind(C, name="hmac")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  real(rt)         :: a(DIMV(abox))
  real(rt)         :: mat(0:6, DIMV(reg))
  real(rt)         :: alpha
  integer :: i, j, k
  if (alpha == 0.e0_rt) then
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = 0.e0_rt
           enddo
        enddo
     enddo
  else
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = alpha * a(i,j,k)
           enddo
        enddo
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
  real(rt)         :: mat(0:6, DIMV(reg))
  real(rt)         :: beta, dx(3)
  real(rt)         :: fac
  integer :: i, j, k
  if (n == 0) then
     fac = beta / (dx(1)**2)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i+1,j,k))
              mat(1,i,j,k) = - fac * b(i,j,k)
              mat(2,i,j,k) = - fac * b(i+1,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     fac = beta / (dx(2)**2)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i,j+1,k))
              mat(3,i,j,k) = - fac * b(i,j,k)
              mat(4,i,j,k) = - fac * b(i,j+1,k)
           enddo
        enddo
     enddo
  else
     fac = beta / (dx(3)**2)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i,j,k+1))
              mat(5,i,j,k) = - fac * b(i,j,k)
              mat(6,i,j,k) = - fac * b(i,j,k+1)
           enddo
        enddo
     enddo
  endif
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
  real(rt)         :: mat(0:6, DIMV(reg))
  real(rt)         :: alpha2
  real(rt)         :: fac
  integer :: i, j, k
  fac = 0.25e0_rt * alpha2
  if (n == 0) then
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (a2(i,j,k) + a2(i+1,j,k))
              mat(1,i,j,k) = mat(1,i,j,k) + fac * a2(i,j,k)
              mat(2,i,j,k) = mat(2,i,j,k) + fac * a2(i+1,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (a2(i,j,k) + a2(i,j+1,k))
              mat(3,i,j,k) = mat(3,i,j,k) + fac * a2(i,j,k)
              mat(4,i,j,k) = mat(4,i,j,k) + fac * a2(i,j+1,k)
           enddo
        enddo
     enddo
  else
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (a2(i,j,k) + a2(i,j,k+1))
              mat(5,i,j,k) = mat(5,i,j,k) + fac * a2(i,j,k)
              mat(6,i,j,k) = mat(6,i,j,k) + fac * a2(i,j,k+1)
           enddo
        enddo
     enddo
  endif
end subroutine hma2c

subroutine hmcc(mat, c, &
                DIMS(bbox), &
                DIMS(reg), &
                gamma, dx, n) bind(C, name="hmcc")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(reg)
  integer :: n
  real(rt)         :: c(DIMV(bbox))
  real(rt)         :: mat(0:6, DIMV(reg))
  real(rt)         :: gamma, dx(3)
  real(rt)         :: fac
  integer :: i, j, k
  if (n == 0) then
     fac = 0.5e0_rt * gamma / dx(1)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) - fac * (c(i,j,k) - c(i+1,j,k))
              mat(1,i,j,k) = mat(1,i,j,k) - fac * c(i,j,k)
              mat(2,i,j,k) = mat(2,i,j,k) + fac * c(i+1,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     fac = 0.5e0_rt * gamma / dx(2)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) - fac * (c(i,j,k) - c(i,j+1,k))
              mat(3,i,j,k) = mat(3,i,j,k) - fac * c(i,j,k)
              mat(4,i,j,k) = mat(4,i,j,k) + fac * c(i,j+1,k)
           enddo
        enddo
     enddo
  else
     fac = 0.5e0_rt * gamma / dx(3)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) - fac * (c(i,j,k) - c(i,j,k+1))
              mat(5,i,j,k) = mat(5,i,j,k) - fac * c(i,j,k)
              mat(6,i,j,k) = mat(6,i,j,k) + fac * c(i,j,k+1)
           enddo
        enddo
     enddo
  endif
end subroutine hmcc

subroutine hmd1c(mat, d1, &
                 DIMS(abox), &
                 DIMS(reg), &
                 delta1, dx, n) bind(C, name="hmd1c")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(abox)
  integer :: DIMDEC(reg)
  integer :: n
  real(rt)         :: d1(DIMV(abox))
  real(rt)         :: mat(0:6, DIMV(reg))
  real(rt)         :: delta1, dx(3)
  real(rt)         :: fac
  integer :: i, j, k
  if (n == 0) then
     fac = 0.5e0_rt * delta1 / dx(1)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(1,i,j,k) = mat(1,i,j,k) - fac * d1(i,j,k)
              mat(2,i,j,k) = mat(2,i,j,k) + fac * d1(i,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     fac = 0.5e0_rt * delta1 / dx(2)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(3,i,j,k) = mat(3,i,j,k) - fac * d1(i,j,k)
              mat(4,i,j,k) = mat(4,i,j,k) + fac * d1(i,j,k)
           enddo
        enddo
     enddo
  else
     fac = 0.5e0_rt * delta1 / dx(3)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(5,i,j,k) = mat(5,i,j,k) - fac * d1(i,j,k)
              mat(6,i,j,k) = mat(6,i,j,k) + fac * d1(i,j,k)
           enddo
        enddo
     enddo
  endif
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
  real(rt)         :: mat(0:6, DIMV(reg))
  real(rt)         :: delta2, dx(3)
  real(rt)         :: fac
  integer :: i, j, k
  if (n == 0) then
     fac = 0.5e0_rt * delta2 / dx(1)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (d2(i,j,k) - d2(i+1,j,k))
              mat(1,i,j,k) = mat(1,i,j,k) - fac * d2(i,j,k)
              mat(2,i,j,k) = mat(2,i,j,k) + fac * d2(i+1,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     fac = 0.5e0_rt * delta2 / dx(2)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (d2(i,j,k) - d2(i,j+1,k))
              mat(3,i,j,k) = mat(3,i,j,k) - fac * d2(i,j,k)
              mat(4,i,j,k) = mat(4,i,j,k) + fac * d2(i,j+1,k)
           enddo
        enddo
     enddo
  else
     fac = 0.5e0_rt * delta2 / dx(3)
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (d2(i,j,k) - d2(i,j,k+1))
              mat(5,i,j,k) = mat(5,i,j,k) - fac * d2(i,j,k)
              mat(6,i,j,k) = mat(6,i,j,k) + fac * d2(i,j,k+1)
           enddo
        enddo
     enddo
  endif
end subroutine hmd2c

subroutine hmmat(mat, &
                 DIMS(reg), &
                 cdir, bct, bho, bcl, &
                 mask, DIMS(msk), &
                 b, DIMS(bbox), &
                 beta, dx) bind(C, name="hmmat")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(reg)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: cdir, bct, bho
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: mat(0:6, DIMV(reg))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: h, fac, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i, j, k
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
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
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(1,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(2,i,j,k) = mat(2,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i+1,j,k)
              mat(2,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(1,i,j,k) = mat(1,i,j,k) + bfm2 * b(i+1,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(3,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(4,i,j,k) = mat(4,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j+1,k)
              mat(4,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(3,i,j,k) = mat(3,i,j,k) + bfm2 * b(i,j+1,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(5,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(6,i,j,k) = mat(6,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k+1)
              mat(6,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(5,i,j,k) = mat(5,i,j,k) + bfm2 * b(i,j,k+1)
              endif
           endif
        enddo
     enddo
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
  real(rt)         :: bcl, beta, dx(3), c
  real(rt)         :: mat(0:6, DIMV(reg))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: spa(DIMV(spabox))
  real(rt)         :: r(1)
  real(rt)         :: h, fac, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i, j, k, bct
  ! The -fac * b(i,j,k) term applied to the matrix diagonal is the contribution
  ! from the interior stencil which must be removed at the boundary.
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  fac = beta / (h**2)
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
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375e0_rt * bfv
                    bfm2 = -0.125e0_rt * bfv
                 else
                    bfm = 0.25e0_rt * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                    bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(1,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(2,i,j,k) = mat(2,i,j,k) + bfm2
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
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i+1,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i+1,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i+1,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375e0_rt * bfv
                    bfm2 = -0.125e0_rt * bfv
                 else
                    bfm = 0.25e0_rt * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                    bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i+1,j,k)
              mat(2,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(1,i,j,k) = mat(1,i,j,k) + bfm2
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
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375e0_rt * bfv
                    bfm2 = -0.125e0_rt * bfv
                 else
                    bfm = 0.25e0_rt * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                    bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(3,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(4,i,j,k) = mat(4,i,j,k) + bfm2
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
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j+1,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j+1,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j+1,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375e0_rt * bfv
                    bfm2 = -0.125e0_rt * bfv
                 else
                    bfm = 0.25e0_rt * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                    bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j+1,k)
              mat(4,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(3,i,j,k) = mat(3,i,j,k) + bfm2
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
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375e0_rt * bfv
                    bfm2 = -0.125e0_rt * bfv
                 else
                    bfm = 0.25e0_rt * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                    bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(5,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(6,i,j,k) = mat(6,i,j,k) + bfm2
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
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k+1)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k+1)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k+1)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else if (bct == LO_MARSHAK) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375e0_rt * bfv
                    bfm2 = -0.125e0_rt * bfv
                 else
                    bfm = 0.25e0_rt * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.e0_rt * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                    bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k+1)
              mat(6,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(5,i,j,k) = mat(5,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else
     print *, "hmmat3: impossible face orientation"
  endif
end subroutine hmmat3

end module habec_module
