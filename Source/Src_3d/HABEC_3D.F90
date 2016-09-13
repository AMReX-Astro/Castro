#include "LO_BCTYPES.H"

#define dims(a) a/**/l0, a/**/l1, a/**/l2, a/**/h0, a/**/h1, a/**/h2
#define dimdec(a) dims(a)
#define dimv(a) a/**/l0:a/**/h0,a/**/l1:a/**/h1,a/**/l2:a/**/h2

#define tiny 1.d-50

subroutine hacoef(mat, a, &
                  dims(abox), &
                  dims(reg), &
                  alpha)
  implicit none
  integer :: dimdec(abox)
  integer :: dimdec(reg)
  real*8 :: a(dimv(abox))
  real*8 :: mat(0:3, dimv(reg))
  real*8 :: alpha
  integer :: i, j, k
  if (alpha == 0.d0) then
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(3,i,j,k) = 0.d0
           enddo
        enddo
     enddo
  else
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(3,i,j,k) = alpha * a(i,j,k)
           enddo
        enddo
     enddo
  endif
end subroutine hacoef

subroutine hbcoef(mat, b, &
                  dims(bbox), &
                  dims(reg), &
                  beta, dx, n)
  implicit none
  integer :: dimdec(bbox)
  integer :: dimdec(reg)
  integer :: n
  real*8 :: b(dimv(bbox))
  real*8 :: mat(0:3, dimv(reg))
  real*8 :: beta, dx(3)
  real*8 :: fac
  integer :: i, j, k
  if (n == 0) then
     fac = beta / (dx(1)**2)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = - fac * b(i,j,k)
              mat(3,i,j,k) = mat(3,i,j,k) + &
                   fac * (b(i,j,k) + b(i+1,j,k))
           enddo
        enddo
     enddo
  elseif (n == 1) then
     fac = beta / (dx(2)**2)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(1,i,j,k) = - fac * b(i,j,k)
              mat(3,i,j,k) = mat(3,i,j,k) + &
                   fac * (b(i,j,k) + b(i,j+1,k))
           enddo
        enddo
     enddo
  else
     fac = beta / (dx(3)**2)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(2,i,j,k) = - fac * b(i,j,k)
              mat(3,i,j,k) = mat(3,i,j,k) + &
                   fac * (b(i,j,k) + b(i,j,k+1))
           enddo
        enddo
     enddo
  endif
end subroutine hbcoef

subroutine hbmat(mat, &
                 dims(reg), &
                 cdir, bct, bcl, &
                 mask, dims(msk), &
                 b, dims(bbox), &
                 beta, dx)
  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(msk)
  integer :: dimdec(bbox)
  integer :: cdir, bct
  real*8 :: bcl, beta, dx(3)
  real*8 :: mat(0:3, dimv(reg))
  integer :: mask(dimv(msk))
  real*8 :: b(dimv(bbox))
  real*8 :: h, fac, bfm, bfv
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
     bfv = fac * h / (0.5d0 * h + bcl)
     bfm = bfv - fac
  else if (bct == LO_NEUMANN) then
     bfv = beta / h
     bfm = -fac
  else
     print *, "hbmat: unsupported boundary type"
     stop
  endif
  if (cdir == 0) then
     i = regl0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i-1,j,k) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i,j,k)
              mat(0,i,j,k) = 0.d0
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = regh0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i+1,j,k) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i+1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = regl1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j-1,k) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i,j,k)
              mat(1,i,j,k) = 0.d0
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = regh1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j+1,k) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i,j+1,k)
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = regl2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k-1) > 0) then
              mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i,j,k)
              mat(2,i,j,k) = 0.d0
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = regh2
     do j = regl1, regh1
        do i = regl0, regh0
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
                  dims(reg), &
                  cdir, bctype, tf, bcl, &
                  dims(bcv), &
                  mask, dims(msk), &
                  b, dims(bbox), &
                  beta, dx, c, r, &
                  spa, dims(spabox))
  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(bcv)
  integer :: dimdec(msk)
  integer :: dimdec(bbox)
  integer :: dimdec(spabox)
  integer :: cdir, bctype, tf(dimv(bcv))
  real*8 :: bcl, beta, dx(3), c
  real*8 :: mat(0:3, dimv(reg))
  integer :: mask(dimv(msk))
  real*8 :: b(dimv(bbox))
  real*8 :: spa(dimv(spabox))
  real*8 :: r(1)
  real*8 :: h, fac, bfm, bfv
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
     i = regl0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i-1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i-1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 bfv = fac * h / (0.5d0 * h + bcl)
                 bfm = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 bfm = 0.25d0 * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
                 bfm = spa(i,j,k) * bfv
              else
                 print *, "hbmat3: unsupported boundary type"
                 stop
              endif
              mat(3,i,j,k) = mat(3,i,j,k) + bfm - fac * b(i,j,k)
              mat(0,i,j,k) = 0.d0
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = regh0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i+1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i+1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 bfv = fac * h / (0.5d0 * h + bcl)
                 bfm = bfv * b(i+1,j,k)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 bfm = 0.25d0 * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
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
     j = regl1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j-1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j-1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 bfv = fac * h / (0.5d0 * h + bcl)
                 bfm = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 bfm = 0.25d0 * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
                 bfm = spa(i,j,k) * bfv
              else
                 print *, "hbmat3: unsupported boundary type"
                 stop
              endif
              mat(3,i,j,k) = mat(3,i,j,k) + bfm - fac * b(i,j,k)
              mat(1,i,j,k) = 0.d0
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = regh1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j+1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j+1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 bfv = fac * h / (0.5d0 * h + bcl)
                 bfm = bfv * b(i,j+1,k)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 bfm = 0.25d0 * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
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
     k = regl2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k-1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k-1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 bfv = fac * h / (0.5d0 * h + bcl)
                 bfm = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 bfm = 0.25d0 * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
                 bfm = spa(i,j,k) * bfv
              else
                 print *, "hbmat3: unsupported boundary type"
                 stop
              endif
              mat(3,i,j,k) = mat(3,i,j,k) + bfm - fac * b(i,j,k)
              mat(2,i,j,k) = 0.d0
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = regh2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k+1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k+1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 bfv = fac * h / (0.5d0 * h + bcl)
                 bfm = bfv * b(i,j,k+1)
              else if (bct == LO_NEUMANN) then
                 bfm = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 bfm = 0.25d0 * bfv
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
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
                 dims(reg), &
                 cdir, bct, bho, bcl, &
                 bcval, dims(bcv), &
                 mask, dims(msk), &
                 b, dims(bbox), &
                 beta, dx)
  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(bcv)
  integer :: dimdec(msk)
  integer :: dimdec(bbox)
  integer :: cdir, bct, bho
  real*8 :: bcl, beta, dx(3)
  real*8 :: vec(dimv(reg))
  real*8 :: bcval(dimv(bcv))
  integer :: mask(dimv(msk))
  real*8 :: b(dimv(bbox))
  real*8 :: h, bfv
  real*8 :: h2, th2
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
        h2 = 0.5d0 * h
        th2 = 3.d0 * h2
        bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
     else
        bfv = (beta / h) / (0.5d0 * h + bcl)
     endif
  else if (bct == LO_NEUMANN) then
     bfv = beta / h
  else
     print *, "hbvec: unsupported boundary type"
     stop
  endif
  if (cdir == 0) then
     i = regl0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i-1,j,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k) * bcval(i-1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = regh0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i+1,j,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i+1,j,k) * bcval(i+1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = regl1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j-1,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k) * bcval(i,j-1,k)
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = regh1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j+1,k) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j+1,k) * bcval(i,j+1,k)
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = regl2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k-1) > 0) then
              vec(i,j,k) = vec(i,j,k) + &
                   bfv * b(i,j,k) * bcval(i,j,k-1)
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = regh2
     do j = regl1, regh1
        do i = regl0, regh0
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
                  dims(reg), &
                  cdir, bctype, tf, bho, bcl, &
                  bcval, dims(bcv), &
                  mask, dims(msk), &
                  b, dims(bbox), &
                  beta, dx, r)
  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(bcv)
  integer :: dimdec(msk)
  integer :: dimdec(bbox)
  integer :: cdir, bctype, tf(dimv(bcv)), bho
  real*8 :: bcl, beta, dx(3)
  real*8 :: vec(dimv(reg))
  real*8 :: bcval(dimv(bcv))
  integer :: mask(dimv(msk))
  real*8 :: b(dimv(bbox))
  real*8 :: r(1)
  real*8 :: h, bfv
  real*8 :: h2, th2
  integer :: i, j, k, bct
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  if (cdir == 0) then
     i = regl0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i-1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i-1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i-1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = regh0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i+1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i+1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                 endif
                 bfv = bfv * b(i+1,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i+1,j,k)
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = regl1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j-1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j-1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j-1,k)
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = regh1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j+1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j+1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                 endif
                 bfv = bfv * b(i,j+1,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j+1,k)
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = regl2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k-1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k-1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta / h
              else
                 print *, "hbvec3: unsupported boundary type"
                 stop
              endif
              vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j,k-1)
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = regh2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k+1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k+1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                 endif
                 bfv = bfv * b(i,j,k+1)
              else if (bct == LO_NEUMANN) then
                 bfv = beta / h
              else if (bct == LO_MARSHAK .OR. &
                   bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta / h
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
                 dims(fbox), &
                 er, dims(ebox), &
                 dims(reg), &
                 cdir, bct, bho, bcl, &
                 bcval, dims(bcv), &
                 mask, dims(msk), &
                 b, dims(bbox), &
                 beta, dx, inhom)
  implicit none
  integer :: dimdec(fbox)
  integer :: dimdec(ebox)
  integer :: dimdec(reg)
  integer :: dimdec(bcv)
  integer :: dimdec(msk)
  integer :: dimdec(bbox)
  integer :: cdir, bct, bho, inhom
  real*8 :: bcl, beta, dx(3)
  real*8 :: flux(dimv(fbox))
  real*8 :: er(dimv(ebox))
  real*8 :: bcval(dimv(bcv))
  integer :: mask(dimv(msk))
  real*8 :: b(dimv(bbox))
  real*8 :: h, bfm, bfv
  real*8 :: bfm2, h2, th2
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
        h2 = 0.5d0 * h
        th2 = 3.d0 * h2
        bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2))
        bfm = (beta / h) * (th2 - bcl) / (bcl + h2)
        bfm2 = (beta / h) * (bcl - h2) / (bcl + th2)
     else
        bfv = beta / (0.5d0 * h + bcl)
        bfm = bfv
     endif
  else
     print *, "hbflx: unsupported boundary type"
     stop
  endif
  if (inhom == 0) then
     bfv = 0.d0
  endif
  if (cdir == 0) then
     i = regl0
     do k = regl2, regh2
        do j = regl1, regh1
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
     i = regh0
     do k = regl2, regh2
        do j = regl1, regh1
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
     j = regl1
     do k = regl2, regh2
        do i = regl0, regh0
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
     j = regh1
     do k = regl2, regh2
        do i = regl0, regh0
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
     k = regl2
     do j = regl1, regh1
        do i = regl0, regh0
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
     k = regh2
     do j = regl1, regh1
        do i = regl0, regh0
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
                  dims(fbox), &
                  er, dims(ebox), &
                  dims(reg), &
                  cdir, bctype, tf, bho, bcl, &
                  bcval, dims(bcv), &
                  mask, dims(msk), &
                  b, dims(bbox), &
                  beta, dx, c, r, inhom, &
                  spa, dims(spabox))
  implicit none
  integer :: dimdec(fbox)
  integer :: dimdec(ebox)
  integer :: dimdec(reg)
  integer :: dimdec(bcv)
  integer :: dimdec(msk)
  integer :: dimdec(bbox)
  integer :: dimdec(spabox)
  integer :: cdir, bctype, tf(dimv(bcv)), bho, inhom
  real*8 :: bcl, beta, dx(3), c
  real*8 :: flux(dimv(fbox))
  real*8 :: er(dimv(ebox))
  real*8 :: bcval(dimv(bcv))
  integer :: mask(dimv(msk))
  real*8 :: b(dimv(bbox))
  real*8 :: spa(dimv(spabox))
  real*8 :: r(1)
  real*8 :: h, bfm, bfv
  real*8 :: bfm2, h2, th2
  integer :: i, j, k, bct
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  if (cdir == 0) then
     i = regl0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i-1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i-1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = beta / (0.5d0 * h + bcl) * b(i,j,k)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  0.75d0 * beta * c
                    bfm2 = -0.25d0 * beta * c
                 else
                    bfm = 0.5d0 * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  3.0d0 * spa(i,j,k) * beta * c
                    bfm2 = -1.0d0 * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0d0 * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.d0
              endif
              flux(i,j,k) = (bfv * bcval(i-1,j,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k) = flux(i,j,k) - bfm2 * er(i+1,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = regh0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i+1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i+1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i+1,j,k)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i+1,j,k)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i+1,j,k)
                 else
                    bfv = beta / (0.5d0 * h + bcl) * b(i+1,j,k)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  0.75d0 * beta * c
                    bfm2 = -0.25d0 * beta * c
                 else
                    bfm = 0.5d0 * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  3.0d0 * spa(i,j,k) * beta * c
                    bfm2 = -1.0d0 * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0d0 * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.d0
              endif
              flux(i+1,j,k) = -(bfv * bcval(i+1,j,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i+1,j,k) = flux(i+1,j,k) + bfm2 * er(i-1,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = regl1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j-1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j-1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = beta / (0.5d0 * h + bcl) * b(i,j,k)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  0.75d0 * beta * c
                    bfm2 = -0.25d0 * beta * c
                 else
                    bfm = 0.5d0 * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  3.0d0 * spa(i,j,k) * beta * c
                    bfm2 = -1.0d0 * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0d0 * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.d0
              endif
              flux(i,j,k) = (bfv * bcval(i,j-1,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k) = flux(i,j,k) - bfm2 * er(i,j+1,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = regh1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j+1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j+1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j+1,k)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j+1,k)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j+1,k)
                 else
                    bfv = beta / (0.5d0 * h + bcl) * b(i,j+1,k)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  0.75d0 * beta * c
                    bfm2 = -0.25d0 * beta * c
                 else
                    bfm = 0.5d0 * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  3.0d0 * spa(i,j,k) * beta * c
                    bfm2 = -1.0d0 * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0d0 * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.d0
              endif
              flux(i,j+1,k) = -(bfv * bcval(i,j+1,k) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j+1,k) = flux(i,j+1,k) + bfm2 * er(i,j-1,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = regl2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k-1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k-1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = beta / (0.5d0 * h + bcl) * b(i,j,k)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  0.75d0 * beta * c
                    bfm2 = -0.25d0 * beta * c
                 else
                    bfm = 0.5d0 * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  3.0d0 * spa(i,j,k) * beta * c
                    bfm2 = -1.0d0 * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0d0 * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.d0
              endif
              flux(i,j,k) = (bfv * bcval(i,j,k-1) - bfm * er(i,j,k))
              if (bho >= 1) then
                 flux(i,j,k) = flux(i,j,k) - bfm2 * er(i,j,k+1)
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = regh2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k+1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k+1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k+1)
                    bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k+1)
                    bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k+1)
                 else
                    bfv = beta / (0.5d0 * h + bcl) * b(i,j,k+1)
                    bfm = bfv
                 endif
              else if (bct == LO_NEUMANN) then
                 bfv  = beta
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  0.75d0 * beta * c
                    bfm2 = -0.25d0 * beta * c
                 else
                    bfm = 0.5d0 * beta * c
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * beta
                 if (bho >= 1) then
                    bfm  =  3.0d0 * spa(i,j,k) * beta * c
                    bfm2 = -1.0d0 * spa(i,j,k) * beta * c
                 else
                    bfm = 2.0d0 * spa(i,j,k) * beta * c
                 endif
              else
                 print *, "hbflx3: unsupported boundary type"
                 stop
              endif
              if (inhom == 0) then
                 bfv = 0.d0
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
                  dims(dtbox), &
                  er, dims(ebox), &
                  dims(reg), &
                  cdir, bct, bcl, &
                  bcval, dims(bcv), &
                  mask, dims(msk), &
                  d, dims(dbox), &
                  dx)
  implicit none
  integer :: dimdec(dtbox)
  integer :: dimdec(ebox)
  integer :: dimdec(reg)
  integer :: dimdec(bcv)
  integer :: dimdec(msk)
  integer :: dimdec(dbox)
  integer :: cdir, bct
  real*8 :: bcl, dx(3)
  real*8 :: dterm(dimv(dtbox))
  real*8 :: er(dimv(ebox))
  real*8 :: bcval(dimv(bcv))
  integer :: mask(dimv(msk))
  real*8 :: d(dimv(dbox))
  real*8 :: h, bfm, bfv
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
        i = regl0
        do k = regl2, regh2
           do j = regl1, regh1
              if (mask(i-1,j,k) > 0) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i-1,j,k)) &
                      / (0.5d0*h + bcl)
              endif
           enddo
        enddo
     else if (cdir == 3) then
        i = regh0
        do k = regl2, regh2
           do j = regl1, regh1
              if (mask(i+1,j,k) > 0) then
                 dterm(i+1,j,k) = d(i+1,j,k) * &
                      (bcval(i+1,j,k) - er(i,j,k)) &
                      / (0.5d0*h + bcl)
              endif
           enddo
        enddo
     else if (cdir == 1) then
        j = regl1
        do k = regl2, regh2
           do i = regl0, regh0
              if (mask(i,j-1,k) > 0) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i,j-1,k)) &
                      / (0.5d0*h + bcl)
              endif
           enddo
        enddo
     else if (cdir == 4) then
        j = regh1
        do k = regl2, regh2
           do i = regl0, regh0
              if (mask(i,j+1,k) > 0) then
                 dterm(i,j+1,k) = d(i,j+1,k) * &
                      (bcval(i,j+1,k) - er(i,j,k)) &
                      / (0.5d0*h + bcl)
              endif
           enddo
        enddo
     else if (cdir == 2) then
        k = regl2
        do j = regl1, regh1
           do i = regl0, regh0
              if (mask(i,j,k-1) > 0) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i,j,k-1)) &
                      / (0.5d0*h + bcl)
              endif
           enddo
        enddo
     else if (cdir == 5) then
        k = regh2
        do j = regl1, regh1
           do i = regl0, regh0
              if (mask(i,j,k+1) > 0) then
                 dterm(i,j,k+1) = d(i,j,k+1) * &
                      (bcval(i,j,k+1) - er(i,j,k)) &
                      / (0.5d0*h + bcl)
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
                   dims(dtbox), &
                   er, dims(ebox), &
                   dims(reg), &
                   cdir, bctype, tf, bcl, &
                   bcval, dims(bcv), &
                   mask, dims(msk), &
                   d, dims(dbox), &
                   dx)
  implicit none
  integer :: dimdec(dtbox)
  integer :: dimdec(ebox)
  integer :: dimdec(reg)
  integer :: dimdec(bcv)
  integer :: dimdec(msk)
  integer :: dimdec(dbox)
  integer :: cdir, bctype, tf(dimv(bcv))
  real*8 :: bcl, dx(3)
  real*8 :: dterm(dimv(dtbox))
  real*8 :: er(dimv(ebox))
  real*8 :: bcval(dimv(bcv))
  integer :: mask(dimv(msk))
  real*8 :: d(dimv(dbox))
  real*8 :: h, bfm, bfv
  integer :: i, j, k, bct
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  if (cdir == 0) then
     i = regl0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i-1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i-1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i-1,j,k)) &
                      / (0.5d0*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i-1,j,k) == 0.d0) then
                 dterm(i,j,k) = 0.d0
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = regh0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i+1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i+1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i+1,j,k) = d(i+1,j,k) * &
                      (bcval(i+1,j,k) - er(i,j,k)) &
                      / (0.5d0*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i+1,j,k) == 0.d0) then
                 dterm(i+1,j,k) = 0.d0
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = regl1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j-1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j-1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i,j-1,k)) &
                      / (0.5d0*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i,j-1,k) == 0.d0) then
                 dterm(i,j,k) = 0.d0
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = regh1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j+1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j+1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i,j+1,k) = d(i,j+1,k) * &
                      (bcval(i,j+1,k) - er(i,j,k)) &
                      / (0.5d0*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i,j+1,k) == 0.d0) then
                 dterm(i,j+1,k) = 0.d0
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = regl2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k-1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k-1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i,j,k) = d(i,j,k) * &
                      (er(i,j,k) - bcval(i,j,k-1)) &
                      / (0.5d0*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i,j,k-1) == 0.d0) then
                 dterm(i,j,k) = 0.d0
              else
                 print *, "hdterm3: unsupported boundary type"
                 stop
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = regh2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k+1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k+1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 dterm(i,j,k+1) = d(i,j,k+1) * &
                      (bcval(i,j,k+1) - er(i,j,k)) &
                      / (0.5d0*h + bcl)
              else if (bct == LO_NEUMANN &
                   .AND. bcval(i,j,k+1) == 0.d0) then
                 dterm(i,j,k+1) = 0.d0
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
                dims(abox), &
                dims(reg), &
                alpha)
  implicit none
  integer :: dimdec(abox)
  integer :: dimdec(reg)
  real*8 :: a(dimv(abox))
  real*8 :: mat(0:6, dimv(reg))
  real*8 :: alpha
  integer :: i, j, k
  if (alpha == 0.d0) then
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = 0.d0
           enddo
        enddo
     enddo
  else
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = alpha * a(i,j,k)
           enddo
        enddo
     enddo
  endif
end subroutine hmac

subroutine hmbc(mat, b, &
                dims(bbox), &
                dims(reg), &
                beta, dx, n)
  implicit none
  integer :: dimdec(bbox)
  integer :: dimdec(reg)
  integer :: n
  real*8 :: b(dimv(bbox))
  real*8 :: mat(0:6, dimv(reg))
  real*8 :: beta, dx(3)
  real*8 :: fac
  integer :: i, j, k
  if (n == 0) then
     fac = beta / (dx(1)**2)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i+1,j,k))
              mat(1,i,j,k) = - fac * b(i,j,k)
              mat(2,i,j,k) = - fac * b(i+1,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     fac = beta / (dx(2)**2)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i,j+1,k))
              mat(3,i,j,k) = - fac * b(i,j,k)
              mat(4,i,j,k) = - fac * b(i,j+1,k)
           enddo
        enddo
     enddo
  else
     fac = beta / (dx(3)**2)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i,j,k+1))
              mat(5,i,j,k) = - fac * b(i,j,k)
              mat(6,i,j,k) = - fac * b(i,j,k+1)
           enddo
        enddo
     enddo
  endif
end subroutine hmbc

subroutine hma2c(mat, a2, &
                 dims(bbox), &
                 dims(reg), &
                 alpha2, n)
  implicit none
  integer :: dimdec(bbox)
  integer :: dimdec(reg)
  integer :: n
  real*8 :: a2(dimv(bbox))
  real*8 :: mat(0:6, dimv(reg))
  real*8 :: alpha2
  real*8 :: fac
  integer :: i, j, k
  fac = 0.25d0 * alpha2
  if (n == 0) then
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (a2(i,j,k) + a2(i+1,j,k))
              mat(1,i,j,k) = mat(1,i,j,k) + fac * a2(i,j,k)
              mat(2,i,j,k) = mat(2,i,j,k) + fac * a2(i+1,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (a2(i,j,k) + a2(i,j+1,k))
              mat(3,i,j,k) = mat(3,i,j,k) + fac * a2(i,j,k)
              mat(4,i,j,k) = mat(4,i,j,k) + fac * a2(i,j+1,k)
           enddo
        enddo
     enddo
  else
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (a2(i,j,k) + a2(i,j,k+1))
              mat(5,i,j,k) = mat(5,i,j,k) + fac * a2(i,j,k)
              mat(6,i,j,k) = mat(6,i,j,k) + fac * a2(i,j,k+1)
           enddo
        enddo
     enddo
  endif
end subroutine hma2c

subroutine hmcc(mat, c, &
                dims(bbox), &
                dims(reg), &
                gamma, dx, n)
  implicit none
  integer :: dimdec(bbox)
  integer :: dimdec(reg)
  integer :: n
  real*8 :: c(dimv(bbox))
  real*8 :: mat(0:6, dimv(reg))
  real*8 :: gamma, dx(3)
  real*8 :: fac
  integer :: i, j, k
  if (n == 0) then
     fac = 0.5d0 * gamma / dx(1)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) - fac * (c(i,j,k) - c(i+1,j,k))
              mat(1,i,j,k) = mat(1,i,j,k) - fac * c(i,j,k)
              mat(2,i,j,k) = mat(2,i,j,k) + fac * c(i+1,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     fac = 0.5d0 * gamma / dx(2)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) - fac * (c(i,j,k) - c(i,j+1,k))
              mat(3,i,j,k) = mat(3,i,j,k) - fac * c(i,j,k)
              mat(4,i,j,k) = mat(4,i,j,k) + fac * c(i,j+1,k)
           enddo
        enddo
     enddo
  else
     fac = 0.5d0 * gamma / dx(3)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) - fac * (c(i,j,k) - c(i,j,k+1))
              mat(5,i,j,k) = mat(5,i,j,k) - fac * c(i,j,k)
              mat(6,i,j,k) = mat(6,i,j,k) + fac * c(i,j,k+1)
           enddo
        enddo
     enddo
  endif
end subroutine hmcc

subroutine hmd1c(mat, d1, &
                 dims(abox), &
                 dims(reg), &
                 delta1, dx, n)
  implicit none
  integer :: dimdec(abox)
  integer :: dimdec(reg)
  integer :: n
  real*8 :: d1(dimv(abox))
  real*8 :: mat(0:6, dimv(reg))
  real*8 :: delta1, dx(3)
  real*8 :: fac
  integer :: i, j, k
  if (n == 0) then
     fac = 0.5d0 * delta1 / dx(1)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(1,i,j,k) = mat(1,i,j,k) - fac * d1(i,j,k)
              mat(2,i,j,k) = mat(2,i,j,k) + fac * d1(i,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     fac = 0.5d0 * delta1 / dx(2)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(3,i,j,k) = mat(3,i,j,k) - fac * d1(i,j,k)
              mat(4,i,j,k) = mat(4,i,j,k) + fac * d1(i,j,k)
           enddo
        enddo
     enddo
  else
     fac = 0.5d0 * delta1 / dx(3)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(5,i,j,k) = mat(5,i,j,k) - fac * d1(i,j,k)
              mat(6,i,j,k) = mat(6,i,j,k) + fac * d1(i,j,k)
           enddo
        enddo
     enddo
  endif
end subroutine hmd1c

subroutine hmd2c(mat, d2, &
                 dims(bbox), &
                 dims(reg), &
                 delta2, dx, n)
  implicit none
  integer :: dimdec(bbox)
  integer :: dimdec(reg)
  integer :: n
  real*8 :: d2(dimv(bbox))
  real*8 :: mat(0:6, dimv(reg))
  real*8 :: delta2, dx(3)
  real*8 :: fac
  integer :: i, j, k
  if (n == 0) then
     fac = 0.5d0 * delta2 / dx(1)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (d2(i,j,k) - d2(i+1,j,k))
              mat(1,i,j,k) = mat(1,i,j,k) - fac * d2(i,j,k)
              mat(2,i,j,k) = mat(2,i,j,k) + fac * d2(i+1,j,k)
           enddo
        enddo
     enddo
  elseif (n == 1) then
     fac = 0.5d0 * delta2 / dx(2)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (d2(i,j,k) - d2(i,j+1,k))
              mat(3,i,j,k) = mat(3,i,j,k) - fac * d2(i,j,k)
              mat(4,i,j,k) = mat(4,i,j,k) + fac * d2(i,j+1,k)
           enddo
        enddo
     enddo
  else
     fac = 0.5d0 * delta2 / dx(3)
     do k = regl2, regh2
        do j = regl1, regh1
           do i = regl0, regh0
              mat(0,i,j,k) = mat(0,i,j,k) + fac * (d2(i,j,k) - d2(i,j,k+1))
              mat(5,i,j,k) = mat(5,i,j,k) - fac * d2(i,j,k)
              mat(6,i,j,k) = mat(6,i,j,k) + fac * d2(i,j,k+1)
           enddo
        enddo
     enddo
  endif
end subroutine hmd2c

subroutine hmmat(mat, &
                 dims(reg), &
                 cdir, bct, bho, bcl, &
                 mask, dims(msk), &
                 b, dims(bbox), &
                 beta, dx)
  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(msk)
  integer :: dimdec(bbox)
  integer :: cdir, bct, bho
  real*8 :: bcl, beta, dx(3)
  real*8 :: mat(0:6, dimv(reg))
  integer :: mask(dimv(msk))
  real*8 :: b(dimv(bbox))
  real*8 :: h, fac, bfm, bfv
  real*8 :: bfm2, h2, th2
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
        h2 = 0.5d0 * h
        th2 = 3.d0 * h2
        bfm = fac * (th2 - bcl) / (bcl + h2) - fac
        bfm2 = fac * (bcl - h2) / (bcl + th2)
     else
        bfv = (beta / h) / (0.5d0 * h + bcl)
        bfm = bfv - fac
     endif
  else if (bct == LO_NEUMANN) then
     bfm = -fac
     bfm2 = 0.d0
  else
     print *, "hmmat: unsupported boundary type"
     stop
  endif
  if (cdir == 0) then
     i = regl0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i-1,j,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(1,i,j,k) = 0.d0
              if (bho >= 1) then
                 mat(2,i,j,k) = mat(2,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = regh0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i+1,j,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i+1,j,k)
              mat(2,i,j,k) = 0.d0
              if (bho >= 1) then
                 mat(1,i,j,k) = mat(1,i,j,k) + bfm2 * b(i+1,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = regl1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j-1,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(3,i,j,k) = 0.d0
              if (bho >= 1) then
                 mat(4,i,j,k) = mat(4,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = regh1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j+1,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j+1,k)
              mat(4,i,j,k) = 0.d0
              if (bho >= 1) then
                 mat(3,i,j,k) = mat(3,i,j,k) + bfm2 * b(i,j+1,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = regl2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k-1) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(5,i,j,k) = 0.d0
              if (bho >= 1) then
                 mat(6,i,j,k) = mat(6,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = regh2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k+1) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k+1)
              mat(6,i,j,k) = 0.d0
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
                  dims(reg), &
                  cdir, bctype, tf, bho, bcl, &
                  dims(bcv), &
                  mask, dims(msk), &
                  b, dims(bbox), &
                  beta, dx, c, r, &
                  spa, dims(spabox))
  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(bcv)
  integer :: dimdec(msk)
  integer :: dimdec(bbox)
  integer :: dimdec(spabox)
  integer :: cdir, bctype, tf(dimv(bcv)), bho
  real*8 :: bcl, beta, dx(3), c
  real*8 :: mat(0:6, dimv(reg))
  integer :: mask(dimv(msk))
  real*8 :: b(dimv(bbox))
  real*8 :: spa(dimv(spabox))
  real*8 :: r(1)
  real*8 :: h, fac, bfm, bfv
  real*8 :: bfm2, h2, th2
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
     i = regl0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i-1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i-1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375d0 * bfv
                    bfm2 = -0.125d0 * bfv
                 else
                    bfm = 0.25d0 * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5d0 * spa(i,j,k) * bfv
                    bfm2 = -0.5d0 * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(1,i,j,k) = 0.d0
              if (bho >= 1) then
                 mat(2,i,j,k) = mat(2,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = regh0
     do k = regl2, regh2
        do j = regl1, regh1
           if (mask(i+1,j,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i+1,j,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i+1,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i+1,j,k)
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                    bfm = bfv * b(i+1,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375d0 * bfv
                    bfm2 = -0.125d0 * bfv
                 else
                    bfm = 0.25d0 * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5d0 * spa(i,j,k) * bfv
                    bfm2 = -0.5d0 * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i+1,j,k)
              mat(2,i,j,k) = 0.d0
              if (bho >= 1) then
                 mat(1,i,j,k) = mat(1,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = regl1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j-1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j-1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375d0 * bfv
                    bfm2 = -0.125d0 * bfv
                 else
                    bfm = 0.25d0 * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5d0 * spa(i,j,k) * bfv
                    bfm2 = -0.5d0 * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(3,i,j,k) = 0.d0
              if (bho >= 1) then
                 mat(4,i,j,k) = mat(4,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = regh1
     do k = regl2, regh2
        do i = regl0, regh0
           if (mask(i,j+1,k) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j+1,k)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j+1,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j+1,k)
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                    bfm = bfv * b(i,j+1,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375d0 * bfv
                    bfm2 = -0.125d0 * bfv
                 else
                    bfm = 0.25d0 * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5d0 * spa(i,j,k) * bfv
                    bfm2 = -0.5d0 * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j+1,k)
              mat(4,i,j,k) = 0.d0
              if (bho >= 1) then
                 mat(3,i,j,k) = mat(3,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = regl2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k-1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k-1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375d0 * bfv
                    bfm2 = -0.125d0 * bfv
                 else
                    bfm = 0.25d0 * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5d0 * spa(i,j,k) * bfv
                    bfm2 = -0.5d0 * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(5,i,j,k) = 0.d0
              if (bho >= 1) then
                 mat(6,i,j,k) = mat(6,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = regh2
     do j = regl1, regh1
        do i = regl0, regh0
           if (mask(i,j,k+1) > 0) then
              if (bctype == -1) then
                 bct = tf(i,j,k+1)
              else
                 bct = bctype
              endif
              if (bct == LO_DIRICHLET) then
                 if (bho >= 1) then
                    h2 = 0.5d0 * h
                    th2 = 3.d0 * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k+1)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k+1)
                 else
                    bfv = (beta / h) / (0.5d0 * h + bcl)
                    bfm = bfv * b(i,j,k+1)
                 endif
              else if (bct == LO_NEUMANN) then
                 bfm  = 0.d0
                 bfm2 = 0.d0
              else if (bct == LO_MARSHAK) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  0.375d0 * bfv
                    bfm2 = -0.125d0 * bfv
                 else
                    bfm = 0.25d0 * bfv
                 endif
              else if (bct == LO_SANCHEZ_POMRANING) then
                 bfv = 2.d0 * c * beta / h
                 if (bho >= 1) then
                    bfm  =  1.5d0 * spa(i,j,k) * bfv
                    bfm2 = -0.5d0 * spa(i,j,k) * bfv
                 else
                    bfm = spa(i,j,k) * bfv
                 endif
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k+1)
              mat(6,i,j,k) = 0.d0
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

subroutine setabecflux( &
                       dims(reg), dir, &
                       density, dims(density), &
                       dcoef, dims(dcoef), &
                       beta, &
                       dx, &
                       flux, dims(flux))

  implicit none
  integer :: dimdec(reg)
  integer :: dimdec(density)
  integer :: dimdec(dcoef)
  integer :: dimdec(flux)

  real*8 :: density(dimv(density))
  real*8 :: dcoef(dimv(dcoef))
  real*8 :: flux(dimv(flux))

  integer :: dir,i,j,k
  real*8 :: beta, dx(BL_SPACEDIM), fac


  if( dir == 0 ) then

     !...     x-direction

     fac = - beta / dx(1)

     do k = regl2,regh2
        do j = regl1,regh1
           do i = regl0,regh0
              flux(i,j,k) = dcoef(i,j,k) * &
                   (density(i,j,k) - density(i-1,j,k)) * fac
           end do
        end do
     end do

  else if( dir == 1 ) then

     !...     y-direction

     fac = - beta / dx(2)

     do k = regl2,regh2
        do j = regl1,regh1
           do i = regl0,regh0
              flux(i,j,k) = dcoef(i,j,k) * &
                   (density(i,j,k) - density(i,j-1,k)) * fac
           end do
        end do
     end do

  else if( dir == 2 ) then

     !...     z-direction

     fac = - beta / dx(3)

     do k = regl2,regh2
        do j = regl1,regh1
           do i = regl0,regh0
              flux(i,j,k) = dcoef(i,j,k) * &
                   (density(i,j,k) - density(i,j,k-1)) * fac
           end do
        end do
     end do
  endif

  return
end subroutine setabecflux
