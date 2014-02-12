
#include "LO_BCTYPES.H"

#define dims(a) a l0, a l1, a h0, a h1
#define dimdec(a) dims(a)
#define dimv(a) a l0:a h0,a l1:a h1

#define tiny 1.d-50

      subroutine hacoef(mat, a,
     @  dims(abox),
     @  dims(reg),
     @  alpha)
      implicit none
      integer dimdec(abox)
      integer dimdec(reg)
      real*8 a(dimv(abox))
      real*8 mat(0:2, dimv(reg))
      real*8 alpha
      integer i, j
      if (alpha .eq. 0.d0) then
         do j = regl1, regh1
            do i = regl0, regh0
               mat(2,i,j) = 0.d0
            enddo
         enddo
      else
         do j = regl1, regh1
            do i = regl0, regh0
               mat(2,i,j) = alpha * a(i,j)
            enddo
         enddo
      endif
      end

      subroutine hbcoef(mat, b,
     @  dims(bbox),
     @  dims(reg),
     @  beta, dx, n)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer n
      real*8 b(dimv(bbox))
      real*8 mat(0:2, dimv(reg))
      real*8 beta, dx(2)
      real*8 fac
      integer i, j
      if (n .eq. 0) then
         fac = beta / (dx(1)**2)
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = - fac * b(i,j)
               mat(2,i,j) = mat(2,i,j) + fac * (b(i,j) + b(i+1,j))
            enddo
         enddo
      else
         fac = beta / (dx(2)**2)
         do j = regl1, regh1
            do i = regl0, regh0
               mat(1,i,j) = - fac * b(i,j)
               mat(2,i,j) = mat(2,i,j) + fac * (b(i,j) + b(i,j+1))
            enddo
         enddo
      endif
      end

      subroutine hbmat(mat,
     @  dims(reg),
     @  cdir, bct, bcl,
     @  mask, dims(msk),
     @  b, dims(bbox),
     @  beta, dx)
      implicit none
      integer dimdec(reg)
      integer dimdec(msk)
      integer dimdec(bbox)
      integer cdir, bct
      real*8 bcl, beta, dx(2)
      real*8 mat(0:2, dimv(reg))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 h, fac, bfm, bfv
      integer i, j
      if (cdir .eq. 0 .or. cdir .eq. 2) then
         h = dx(1)
      else
         h = dx(2)
      endif
      fac = beta / (h**2)
      if (bct .eq. LO_DIRICHLET) then
         bfv = fac * h / (0.5d0 * h + bcl)
         bfm = bfv - fac
      else
         print *, "hbmat: unsupported boundary type"
         stop
      endif
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         do j = regl1, regh1
            if (mask(i-1,j) .gt. 0) then
               mat(2,i,j) = mat(2,i,j) + bfm * b(i,j)
               mat(0,i,j) = 0.d0
            endif
         enddo
      else if (cdir .eq. 2) then
c Right face of grid
         i = regh0
         do j = regl1, regh1
            if (mask(i+1,j) .gt. 0) then
               mat(2,i,j) = mat(2,i,j) + bfm * b(i+1,j)
            endif
         enddo
      else if (cdir .eq. 1) then
c Bottom face of grid
         j = regl1
         do i = regl0, regh0
            if (mask(i,j-1) .gt. 0) then
               mat(2,i,j) = mat(2,i,j) + bfm * b(i,j)
               mat(1,i,j) = 0.d0
            endif
         enddo
      else if (cdir .eq. 3) then
c Top face of grid
         j = regh1
         do i = regl0, regh0
            if (mask(i,j+1) .gt. 0) then
               mat(2,i,j) = mat(2,i,j) + bfm * b(i,j+1)
            endif
         enddo
      else
         print *, "hbmat: impossible face orientation"
      endif
      end

      subroutine hbmat3(mat,
     @  dims(reg),
     @  cdir, bctype, tf, bcl,
     @  dims(bcv),
     @  mask, dims(msk),
     @  b, dims(bbox),
     @  beta, dx, c, r,
     @  spa, dims(spabox))
      implicit none
      integer dimdec(reg)
      integer dimdec(bcv)
      integer dimdec(msk)
      integer dimdec(bbox)
      integer dimdec(spabox)
      integer cdir, bctype, tf(dimv(bcv))
      real*8 bcl, beta, dx(2), c
      real*8 mat(0:2, dimv(reg))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 spa(dimv(spabox))
      real*8 r(regl0:regh0)
      real*8 h, fac, bfm, bfv, r0
      integer i, j, bct
      if (cdir .eq. 0 .or. cdir .eq. 2) then
         h = dx(1)
c For the left and right faces, r is constant and the array actually has
c only one element.  This following "r(regl0)" looks wrong, but what it
c really means is that the array dimensions are meaningless for this case.
         r0 = r(regl0)
      else
         h = dx(2)
      endif
c The -fac * b(i,j) term applied to the matrix diagonal is the contribution
c from the interior stencil which must be removed at the boundary.
      fac = beta / (h**2)
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         do j = regl1, regh1
            if (mask(i-1,j) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i-1,j)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  bfv = fac * h / (0.5d0 * h + bcl)
                  bfm = bfv * b(i,j)
               else if (bct .eq. LO_NEUMANN) then
                  bfm = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r0 / h
                  bfm = 0.25d0 * c * bfv
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r0 / h
                  bfm = spa(i,j) * c * bfv
               else
                  print *, "hbmat3: unsupported boundary type"
                  stop
               endif
               mat(2,i,j) = mat(2,i,j) + bfm - fac * b(i,j)
               mat(0,i,j) = 0.d0
            endif
         enddo
      else if (cdir .eq. 2) then
c Right face of grid
         i = regh0
         do j = regl1, regh1
            if (mask(i+1,j) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i+1,j)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  bfv = fac * h / (0.5d0 * h + bcl)
                  bfm = bfv * b(i+1,j)
               else if (bct .eq. LO_NEUMANN) then
                  bfm = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r0 / h
                  bfm = 0.25d0 * c * bfv
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r0 / h
                  bfm = spa(i,j) * c * bfv
               else
                  print *, "hbmat3: unsupported boundary type"
                  stop
               endif
               mat(2,i,j) = mat(2,i,j) + bfm - fac * b(i+1,j)
            endif
         enddo
      else if (cdir .eq. 1) then
c Bottom face of grid
         j = regl1
         do i = regl0, regh0
            if (mask(i,j-1) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i,j-1)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  bfv = fac * h / (0.5d0 * h + bcl)
                  bfm = bfv * b(i,j)
               else if (bct .eq. LO_NEUMANN) then
                  bfm = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r(i) / h
                  bfm = 0.25d0 * c * bfv
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r(i) / h
                  bfm = spa(i,j) * c * bfv
               else
                  print *, "hbmat3: unsupported boundary type"
                  stop
               endif
               mat(2,i,j) = mat(2,i,j) + bfm - fac * b(i,j)
               mat(1,i,j) = 0.d0
            endif
         enddo
      else if (cdir .eq. 3) then
c Top face of grid
         j = regh1
         do i = regl0, regh0
            if (mask(i,j+1) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i,j+1)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  bfv = fac * h / (0.5d0 * h + bcl)
                  bfm = bfv * b(i,j+1)
               else if (bct .eq. LO_NEUMANN) then
                  bfm = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r(i) / h
                  bfm = 0.25d0 * c * bfv
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r(i) / h
                  bfm = spa(i,j) * c * bfv
               else
                  print *, "hbmat3: unsupported boundary type"
                  stop
               endif
               mat(2,i,j) = mat(2,i,j) + bfm - fac * b(i,j+1)
            endif
         enddo
      else
         print *, "hbmat3: impossible face orientation"
      endif
      end

      subroutine hbvec(vec,
     @  dims(reg),
     @  cdir, bct, bho, bcl,
     @  bcval, dims(bcv),
     @  mask, dims(msk),
     @  b, dims(bbox),
     @  beta, dx)
      implicit none
      integer dimdec(reg)
      integer dimdec(bcv)
      integer dimdec(msk)
      integer dimdec(bbox)
      integer cdir, bct, bho
      real*8 bcl, beta, dx(2)
      real*8 vec(dimv(reg))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 h, bfv
      real*8 h2, th2
      integer i, j
      if (cdir .eq. 0 .or. cdir .eq. 2) then
         h = dx(1)
      else
         h = dx(2)
      endif
      if (bct .eq. LO_DIRICHLET) then
         if (bho .ge. 1) then
            h2 = 0.5d0 * h
            th2 = 3.d0 * h2
            bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
         else
            bfv = (beta / h) / (0.5d0 * h + bcl)
         endif
      else
         print *, "hbvec: unsupported boundary type"
         stop
      endif
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         do j = regl1, regh1
            if (mask(i-1,j) .gt. 0) then
               vec(i,j) = vec(i,j) + bfv * b(i,j) * bcval(i-1,j)
            endif
         enddo
      else if (cdir .eq. 2) then
c Right face of grid
         i = regh0
         do j = regl1, regh1
            if (mask(i+1,j) .gt. 0) then
               vec(i,j) = vec(i,j) + bfv * b(i+1,j) * bcval(i+1,j)
            endif
         enddo
      else if (cdir .eq. 1) then
c Bottom face of grid
         j = regl1
         do i = regl0, regh0
            if (mask(i,j-1) .gt. 0) then
               vec(i,j) = vec(i,j) + bfv * b(i,j) * bcval(i,j-1)
            endif
         enddo
      else if (cdir .eq. 3) then
c Top face of grid
         j = regh1
         do i = regl0, regh0
            if (mask(i,j+1) .gt. 0) then
               vec(i,j) = vec(i,j) + bfv * b(i,j+1) * bcval(i,j+1)
            endif
         enddo
      else
         print *, "hbvec: impossible face orientation"
      endif
      end

      subroutine hbvec3(vec,
     @  dims(reg),
     @  cdir, bctype, tf, bho, bcl,
     @  bcval, dims(bcv),
     @  mask, dims(msk),
     @  b, dims(bbox),
     @  beta, dx, r)
      implicit none
      integer dimdec(reg)
      integer dimdec(bcv)
      integer dimdec(msk)
      integer dimdec(bbox)
      integer cdir, bctype, tf(dimv(bcv)), bho
      real*8 bcl, beta, dx(2)
      real*8 vec(dimv(reg))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 r(regl0:regh0)
      real*8 h, bfv, r0
      real*8 h2, th2
      integer i, j, bct
      if (cdir .eq. 0 .or. cdir .eq. 2) then
         h = dx(1)
c For the left and right faces, r is constant and the array actually has
c only one element.  This following "r(regl0)" looks wrong, but what it
c really means is that the array dimensions are meaningless for this case.
         r0 = r(regl0)
      else
         h = dx(2)
      endif
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         do j = regl1, regh1
            if (mask(i-1,j) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i-1,j)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
                  else
                     bfv = (beta / h) / (0.5d0 * h + bcl)
                  endif
                  bfv = bfv * b(i,j)
               else if (bct .eq. LO_NEUMANN) then
                  bfv = beta / h
               else if (bct .eq. LO_MARSHAK .or.
     @                  bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r0 / h
               else
                  print *, "hbvec3: unsupported boundary type"
                  stop
               endif
               vec(i,j) = vec(i,j) + bfv * bcval(i-1,j)
            endif
         enddo
      else if (cdir .eq. 2) then
c Right face of grid
         i = regh0
         do j = regl1, regh1
            if (mask(i+1,j) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i+1,j)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
                  else
                     bfv = (beta / h) / (0.5d0 * h + bcl)
                  endif
                  bfv = bfv * b(i+1,j)
               else if (bct .eq. LO_NEUMANN) then
                  bfv = beta / h
               else if (bct .eq. LO_MARSHAK .or.
     @                  bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r0 / h
               else
                  print *, "hbvec3: unsupported boundary type"
                  stop
               endif
               vec(i,j) = vec(i,j) + bfv * bcval(i+1,j)
            endif
         enddo
      else if (cdir .eq. 1) then
c Bottom face of grid
         j = regl1
         do i = regl0, regh0
            if (mask(i,j-1) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i,j-1)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
                  else
                     bfv = (beta / h) / (0.5d0 * h + bcl)
                  endif
                  bfv = bfv * b(i,j)
               else if (bct .eq. LO_NEUMANN) then
                  bfv = beta / h
               else if (bct .eq. LO_MARSHAK .or.
     @                  bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r(i) / h
               else
                  print *, "hbvec3: unsupported boundary type"
                  stop
               endif
               vec(i,j) = vec(i,j) + bfv * bcval(i,j-1)
            endif
         enddo
      else if (cdir .eq. 3) then
c Top face of grid
         j = regh1
         do i = regl0, regh0
            if (mask(i,j+1) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i,j+1)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
                  else
                     bfv = (beta / h) / (0.5d0 * h + bcl)
                  endif
                  bfv = bfv * b(i,j+1)
               else if (bct .eq. LO_NEUMANN) then
                  bfv = beta / h
               else if (bct .eq. LO_MARSHAK .or.
     @                  bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r(i) / h
               else
                  print *, "hbvec3: unsupported boundary type"
                  stop
               endif
               vec(i,j) = vec(i,j) + bfv * bcval(i,j+1)
            endif
         enddo
      else
         print *, "hbvec3: impossible face orientation"
      endif
      end

      subroutine hbflx(flux,
     @  dims(fbox),
     @  er, dims(ebox),
     @  dims(reg),
     @  cdir, bct, bho, bcl,
     @  bcval, dims(bcv),
     @  mask, dims(msk),
     @  b, dims(bbox),
     @  beta, dx, inhom)
      implicit none
      integer dimdec(fbox)
      integer dimdec(ebox)
      integer dimdec(reg)
      integer dimdec(bcv)
      integer dimdec(msk)
      integer dimdec(bbox)
      integer cdir, bct, bho, inhom
      real*8 bcl, beta, dx(2)
      real*8 flux(dimv(fbox))
      real*8 er(dimv(ebox))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 h, bfm, bfv
      real*8 bfm2, h2, th2
      integer i, j
      if (cdir .eq. 0 .or. cdir .eq. 2) then
         h = dx(1)
      else
         h = dx(2)
      endif
      if (bct .eq. LO_DIRICHLET) then
         if (bho .ge. 1) then
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
      if (inhom .eq. 0) then
         bfv = 0.d0
      endif
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         do j = regl1, regh1
            if (mask(i-1,j) .gt. 0) then
               flux(i,j) = b(i,j) * (bfv * bcval(i-1,j) - bfm * er(i,j))
               if (bho .ge. 1) then
                  flux(i,j) = flux(i,j) - b(i,j) * bfm2 * er(i+1,j)
               endif
            endif
         enddo
      else if (cdir .eq. 2) then
c Right face of grid
         i = regh0
         do j = regl1, regh1
            if (mask(i+1,j) .gt. 0) then
               flux(i+1,j) = -b(i+1,j) * (bfv * bcval(i+1,j) - bfm * er(i,j))
               if (bho .ge. 1) then
                  flux(i+1,j) = flux(i+1,j) + b(i+1,j) * bfm2 * er(i-1,j)
               endif
            endif
         enddo
      else if (cdir .eq. 1) then
c Bottom face of grid
         j = regl1
         do i = regl0, regh0
            if (mask(i,j-1) .gt. 0) then
               flux(i,j) = b(i,j) * (bfv * bcval(i,j-1) - bfm * er(i,j))
               if (bho .ge. 1) then
                  flux(i,j) = flux(i,j) - b(i,j) * bfm2 * er(i,j+1)
               endif
            endif
         enddo
      else if (cdir .eq. 3) then
c Top face of grid
         j = regh1
         do i = regl0, regh0
            if (mask(i,j+1) .gt. 0) then
               flux(i,j+1) = -b(i,j+1) * (bfv * bcval(i,j+1) - bfm * er(i,j))
               if (bho .ge. 1) then
                  flux(i,j+1) = flux(i,j+1) + b(i,j+1) * bfm2 * er(i,j-1)
               endif
            endif
         enddo
      else
         print *, "hbflx: impossible face orientation"
      endif
      end

      subroutine hbflx3(flux,
     @  dims(fbox),
     @  er, dims(ebox),
     @  dims(reg),
     @  cdir, bctype, tf, bho, bcl,
     @  bcval, dims(bcv),
     @  mask, dims(msk),
     @  b, dims(bbox),
     @  beta, dx, c, r, inhom,
     @  spa, dims(spabox))
      implicit none
      integer dimdec(fbox)
      integer dimdec(ebox)
      integer dimdec(reg)
      integer dimdec(bcv)
      integer dimdec(msk)
      integer dimdec(bbox)
      integer dimdec(spabox)
      integer cdir, bctype, tf(dimv(bcv)), bho, inhom
      real*8 bcl, beta, dx(2), c
      real*8 flux(dimv(fbox))
      real*8 er(dimv(ebox))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 spa(dimv(spabox))
      real*8 r(regl0:regh0)
      real*8 h, bfm, bfv, r0
      real*8 bfm2, h2, th2
      integer i, j, bct
      if (cdir .eq. 0 .or. cdir .eq. 2) then
         h = dx(1)
c For the left and right faces, r is constant and the array actually has
c only one element.  This following "r(regl0)" looks wrong, but what it
c really means is that the array dimensions are meaningless for this case.
         r0 = r(regl0)
      else
         h = dx(2)
      endif
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         do j = regl1, regh1
            if (mask(i-1,j) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i-1,j)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j)
                     bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j)
                     bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j)
                  else
                     bfv = beta / (0.5d0 * h + bcl) * b(i,j)
                     bfm = bfv
                  endif
               else if (bct .eq. LO_NEUMANN) then
                  bfv  = beta
                  bfm  = 0.d0
                  bfm2 = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r0
                  if (bho .ge. 1) then
                     bfm  =  0.375d0 * c * bfv
                     bfm2 = -0.125d0 * c * bfv
                  else
                     bfm = 0.25d0 * c * bfv
                  endif
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r0
                  if (bho .ge. 1) then
                     bfm  =  1.5d0 * spa(i,j) * c * bfv
                     bfm2 = -0.5d0 * spa(i,j) * c * bfv
                  else
                     bfm = spa(i,j) * c * bfv
                  endif
               else
                  print *, "hbflx3: unsupported boundary type"
                  stop
               endif
               if (inhom .eq. 0) then
                  bfv = 0.d0
               endif
               flux(i,j) = (bfv * bcval(i-1,j) - bfm * er(i,j))
               if (bho .ge. 1) then
                  flux(i,j) = flux(i,j) - bfm2 * er(i+1,j)
               endif
            endif
         enddo
      else if (cdir .eq. 2) then
c Right face of grid
         i = regh0
         do j = regl1, regh1
            if (mask(i+1,j) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i+1,j)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i+1,j)
                     bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i+1,j)
                     bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i+1,j)
                  else
                     bfv = beta / (0.5d0 * h + bcl) * b(i+1,j)
                     bfm = bfv
                  endif
               else if (bct .eq. LO_NEUMANN) then
                  bfv  = beta
                  bfm  = 0.d0
                  bfm2 = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r0
                  if (bho .ge. 1) then
                     bfm  =  0.375d0 * c * bfv
                     bfm2 = -0.125d0 * c * bfv
                  else
                     bfm = 0.25d0 * c * bfv
                  endif
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r0
                  if (bho .ge. 1) then
                     bfm  =  1.5d0 * spa(i,j) * c * bfv
                     bfm2 = -0.5d0 * spa(i,j) * c * bfv
                  else
                     bfm = spa(i,j) * c * bfv
                  endif
               else
                  print *, "hbflx3: unsupported boundary type"
                  stop
               endif
               if (inhom .eq. 0) then
                  bfv = 0.d0
               endif
               flux(i+1,j) = -(bfv * bcval(i+1,j) - bfm * er(i,j))
               if (bho .ge. 1) then
                  flux(i+1,j) = flux(i+1,j) + bfm2 * er(i-1,j)
               endif
            endif
         enddo
      else if (cdir .eq. 1) then
c Bottom face of grid
         j = regl1
         do i = regl0, regh0
            if (mask(i,j-1) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i,j-1)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j)
                     bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j)
                     bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j)
                  else
                     bfv = beta / (0.5d0 * h + bcl) * b(i,j)
                     bfm = bfv
                  endif
               else if (bct .eq. LO_NEUMANN) then
                  bfv  = beta
                  bfm  = 0.d0
                  bfm2 = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r(i)
                  if (bho .ge. 1) then
                     bfm  =  0.375d0 * c * bfv
                     bfm2 = -0.125d0 * c * bfv
                  else
                     bfm = 0.25d0 * c * bfv
                  endif
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r(i)
                  if (bho .ge. 1) then
                     bfm  =  1.5d0 * spa(i,j) * c * bfv
                     bfm2 = -0.5d0 * spa(i,j) * c * bfv
                  else
                     bfm = spa(i,j) * c * bfv
                  endif
               else
                  print *, "hbflx3: unsupported boundary type"
                  stop
               endif
               if (inhom .eq. 0) then
                  bfv = 0.d0
               endif
               flux(i,j) = (bfv * bcval(i,j-1) - bfm * er(i,j))
               if (bho .ge. 1) then
                  flux(i,j) = flux(i,j) - bfm2 * er(i,j+1)
               endif
            endif
         enddo
      else if (cdir .eq. 3) then
c Top face of grid
         j = regh1
         do i = regl0, regh0
            if (mask(i,j+1) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i,j+1)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j+1)
                     bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j+1)
                     bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j+1)
                  else
                     bfv = beta / (0.5d0 * h + bcl) * b(i,j+1)
                     bfm = bfv
                  endif
               else if (bct .eq. LO_NEUMANN) then
                  bfv  = beta
                  bfm  = 0.d0
                  bfm2 = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r(i)
                  if (bho .ge. 1) then
                     bfm  =  0.375d0 * c * bfv
                     bfm2 = -0.125d0 * c * bfv
                  else
                     bfm = 0.25d0 * c * bfv
                  endif
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r(i)
                  if (bho .ge. 1) then
                     bfm  =  1.5d0 * spa(i,j) * c * bfv
                     bfm2 = -0.5d0 * spa(i,j) * c * bfv
                  else
                     bfm = spa(i,j) * c * bfv
                  endif
               else
                  print *, "hbflx3: unsupported boundary type"
                  stop
               endif
               if (inhom .eq. 0) then
                  bfv = 0.d0
               endif
               flux(i,j+1) = -(bfv * bcval(i,j+1) - bfm * er(i,j))
               if (bho .ge. 1) then
                  flux(i,j+1) = flux(i,j+1) + bfm2 * er(i,j-1)
               endif
            endif
         enddo
      else
         print *, "hbflx3: impossible face orientation"
      endif
      end

      subroutine hdterm(dterm,
     @  dims(dtbox),
     @  er, dims(ebox),
     @  dims(reg),
     @  cdir, bct, bcl,
     @  bcval, dims(bcv),
     @  mask, dims(msk),
     @  d, dims(dbox),
     @  dx)
      implicit none
      integer dimdec(dtbox)
      integer dimdec(ebox)
      integer dimdec(reg)
      integer dimdec(bcv)
      integer dimdec(msk)
      integer dimdec(dbox)
      integer cdir, bct
      real*8 bcl, dx(2)
      real*8 dterm(dimv(dtbox))
      real*8 er(dimv(ebox))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 d(dimv(dbox))
      real*8 h 
      integer i, j
      if (cdir .eq. 0 .or. cdir .eq. 2) then
         h = dx(1)
      else
         h = dx(2)
      endif
      if (bct .eq. LO_DIRICHLET) then
         if (cdir .eq. 0) then
c     Left face of grid
            i = regl0
            do j = regl1, regh1
               if (mask(i-1,j) .gt. 0) then
                  dterm(i,j) = d(i,j)*(er(i,j)-bcval(i-1,j))/(0.5d0*h+bcl)
               endif
            enddo
         else if (cdir .eq. 2) then
c     Right face of grid
            i = regh0
            do j = regl1, regh1
               if (mask(i+1,j) .gt. 0) then
                  dterm(i+1,j) = d(i+1,j)*(bcval(i+1,j)-er(i,j))/(0.5d0*h+bcl)
               endif
            enddo
         else if (cdir .eq. 1) then
c     Bottom face of grid
            j = regl1
            do i = regl0, regh0
               if (mask(i,j-1) .gt. 0) then
                  dterm(i,j) = d(i,j)*(er(i,j)-bcval(i,j-1))/(0.5d0*h+bcl)
               endif
            enddo
         else if (cdir .eq. 3) then
c     Top face of grid
            j = regh1
            do i = regl0, regh0
               if (mask(i,j+1) .gt. 0) then
                  dterm(i,j+1) = d(i,j+1)*(bcval(i,j+1)-er(i,j))/(0.5d0*h+bcl)
               endif
            enddo
         else
            print *, "hdterm: impossible face orientation"
         endif
      else
         print *, "hdterm: unsupported boundary type"
         stop
      endif
      end

      subroutine hdterm3(dterm,
     @  dims(dtbox),
     @  er, dims(ebox),
     @  dims(reg),
     @  cdir, bctype, tf, bcl,
     @  bcval, dims(bcv),
     @  mask, dims(msk),
     @  d, dims(dbox),
     @  dx)
      implicit none
      integer dimdec(dtbox)
      integer dimdec(ebox)
      integer dimdec(reg)
      integer dimdec(bcv)
      integer dimdec(msk)
      integer dimdec(dbox)
      integer cdir, bctype, tf(dimv(bcv))
      real*8 bcl, dx(2)
      real*8 dterm(dimv(dtbox))
      real*8 er(dimv(ebox))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 d(dimv(dbox))
      real*8 h
      integer i, j, bct
      if (cdir .eq. 0 .or. cdir .eq. 2) then
         h = dx(1)
      else
         h = dx(2)
      endif
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         do j = regl1, regh1
            if (mask(i-1,j) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i-1,j)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  dterm(i,j) = d(i,j)*(er(i,j)-bcval(i-1,j))/(0.5d0*h+bcl)
               else if (bct.eq.LO_NEUMANN .and. bcval(i-1,j).eq.0.d0) then
                  dterm(i,j) = 0.d0
               else
                  print *, "hdterm3: unsupported boundary type"
                  stop
               endif
            endif
         enddo
      else if (cdir .eq. 2) then
c Right face of grid
         i = regh0
         do j = regl1, regh1
            if (mask(i+1,j) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i+1,j)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  dterm(i+1,j) = d(i+1,j)*(bcval(i+1,j)-er(i,j))/(0.5d0*h+bcl)
               else if (bct.eq.LO_NEUMANN .and. bcval(i+1,j).eq.0.d0) then
                  dterm(i+1,j) = 0.d0
               else
                  print *, "hdterm3: unsupported boundary type"
                  stop
               endif
            endif
         enddo
      else if (cdir .eq. 1) then
c Bottom face of grid
         j = regl1
         do i = regl0, regh0
            if (mask(i,j-1) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i,j-1)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  dterm(i,j) = d(i,j)*(er(i,j)-bcval(i,j-1))/(0.5d0*h+bcl)
               else if (bct.eq.LO_NEUMANN .and. bcval(i,j-1).eq.0.d0) then
                  dterm(i,j) = 0.d0
               else
                  print *, "hdterm3: unsupported boundary type"
                  stop
               endif
            endif
         enddo
      else if (cdir .eq. 3) then
c Top face of grid
         j = regh1
         do i = regl0, regh0
            if (mask(i,j+1) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i,j+1)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  dterm(i,j+1) = d(i,j+1)*(bcval(i,j+1)-er(i,j))/(0.5d0*h+bcl)
               else if (bct.eq.LO_NEUMANN .and. bcval(i,j+1).eq.0.d0) then
                  dterm(i,j+1) = 0.d0
               else
                  print *, "hdterm3: unsupported boundary type"
                  stop
               endif
            endif
         enddo
      else
         print *, "hdterm3: impossible face orientation"
      endif
      end

      subroutine hmac(mat, a,
     @  dims(abox),
     @  dims(reg),
     @  alpha)
      implicit none
      integer dimdec(abox)
      integer dimdec(reg)
      real*8 a(dimv(abox))
      real*8 mat(0:4, dimv(reg))
      real*8 alpha
      integer i, j
      if (alpha .eq. 0.d0) then
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = 0.d0
            enddo
         enddo
      else
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = alpha * a(i,j)
            enddo
         enddo
      endif
      end

      subroutine hmbc(mat, b,
     @  dims(bbox),
     @  dims(reg),
     @  beta, dx, n)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer n
      real*8 b(dimv(bbox))
      real*8 mat(0:4, dimv(reg))
      real*8 beta, dx(2)
      real*8 fac
      integer i, j
      if (n .eq. 0) then
         fac = beta / (dx(1)**2)
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = mat(0,i,j) + fac * (b(i,j) + b(i+1,j))
               mat(1,i,j) = - fac * b(i,j)
               mat(2,i,j) = - fac * b(i+1,j)
            enddo
         enddo
      else
         fac = beta / (dx(2)**2)
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = mat(0,i,j) + fac * (b(i,j) + b(i,j+1))
               mat(3,i,j) = - fac * b(i,j)
               mat(4,i,j) = - fac * b(i,j+1)
            enddo
         enddo
      endif
      end

      subroutine hma2c(mat, a2,
     @  dims(bbox),
     @  dims(reg),
     @  alpha2, n)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer n
      real*8 a2(dimv(bbox))
      real*8 mat(0:4, dimv(reg))
      real*8 alpha2
      real*8 fac
      integer i, j
      fac = 0.25d0 * alpha2
      if (n .eq. 0) then
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = mat(0,i,j) + fac * (a2(i,j) + a2(i+1,j))
               mat(1,i,j) = mat(1,i,j) + fac * a2(i,j)
               mat(2,i,j) = mat(2,i,j) + fac * a2(i+1,j)
            enddo
         enddo
      else
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = mat(0,i,j) + fac * (a2(i,j) + a2(i,j+1))
               mat(3,i,j) = mat(3,i,j) + fac * a2(i,j)
               mat(4,i,j) = mat(4,i,j) + fac * a2(i,j+1)
            enddo
         enddo
      endif
      end

      subroutine hmcc(mat, c,
     @  dims(bbox),
     @  dims(reg),
     @  gamma, dx, n)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer n
      real*8 c(dimv(bbox))
      real*8 mat(0:4, dimv(reg))
      real*8 gamma, dx(2)
      real*8 fac
      integer i, j
      if (n .eq. 0) then
         fac = 0.5d0 * gamma / dx(1)
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = mat(0,i,j) - fac * (c(i,j) - c(i+1,j))
               mat(1,i,j) = mat(1,i,j) - fac * c(i,j)
               mat(2,i,j) = mat(2,i,j) + fac * c(i+1,j)
            enddo
         enddo
      else
         fac = 0.5d0 * gamma / dx(2)
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = mat(0,i,j) - fac * (c(i,j) - c(i,j+1))
               mat(3,i,j) = mat(3,i,j) - fac * c(i,j)
               mat(4,i,j) = mat(4,i,j) + fac * c(i,j+1)
            enddo
         enddo
      endif
      end

      subroutine hmd1c(mat, d1,
     @  dims(abox),
     @  dims(reg),
     @  delta1, dx, n)
      implicit none
      integer dimdec(abox)
      integer dimdec(reg)
      integer n
      real*8 d1(dimv(abox))
      real*8 mat(0:4, dimv(reg))
      real*8 delta1, dx(2)
      real*8 fac
      integer i, j
      if (n .eq. 0) then
         fac = 0.5d0 * delta1 / dx(1)
         do j = regl1, regh1
            do i = regl0, regh0
               mat(1,i,j) = mat(1,i,j) - fac * d1(i,j)
               mat(2,i,j) = mat(2,i,j) + fac * d1(i,j)
            enddo
         enddo
      else
         fac = 0.5d0 * delta1 / dx(2)
         do j = regl1, regh1
            do i = regl0, regh0
               mat(3,i,j) = mat(3,i,j) - fac * d1(i,j)
               mat(4,i,j) = mat(4,i,j) + fac * d1(i,j)
            enddo
         enddo
      endif
      end

      subroutine hmd2c(mat, d2,
     @  dims(bbox),
     @  dims(reg),
     @  delta2, dx, n)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer n
      real*8 d2(dimv(bbox))
      real*8 mat(0:4, dimv(reg))
      real*8 delta2, dx(2)
      real*8 fac
      integer i, j
      if (n .eq. 0) then
         fac = 0.5d0 * delta2 / dx(1)
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = mat(0,i,j) + fac * (d2(i,j) - d2(i+1,j))
               mat(1,i,j) = mat(1,i,j) - fac * d2(i,j)
               mat(2,i,j) = mat(2,i,j) + fac * d2(i+1,j)
            enddo
         enddo
      else
         fac = 0.5d0 * delta2 / dx(2)
         do j = regl1, regh1
            do i = regl0, regh0
               mat(0,i,j) = mat(0,i,j) + fac * (d2(i,j) - d2(i,j+1))
               mat(3,i,j) = mat(3,i,j) - fac * d2(i,j)
               mat(4,i,j) = mat(4,i,j) + fac * d2(i,j+1)
            enddo
         enddo
      endif
      end

      subroutine hmmat(mat,
     @  dims(reg),
     @  cdir, bct, bho, bcl,
     @  mask, dims(msk),
     @  b, dims(bbox),
     @  beta, dx)
      implicit none
      integer dimdec(reg)
      integer dimdec(msk)
      integer dimdec(bbox)
      integer cdir, bct, bho
      real*8 bcl, beta, dx(2)
      real*8 mat(0:4, dimv(reg))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 h, fac, bfm, bfv
      real*8 bfm2, h2, th2
      integer i, j
      if (cdir .eq. 0 .or. cdir .eq. 2) then
         h = dx(1)
      else
         h = dx(2)
      endif
      fac = beta / (h**2)
      if (bct .eq. LO_DIRICHLET) then
         if (bho .ge. 1) then
            h2 = 0.5d0 * h
            th2 = 3.d0 * h2
            bfm = fac * (th2 - bcl) / (bcl + h2) - fac
            bfm2 = fac * (bcl - h2) / (bcl + th2)
         else
            bfv = (beta / h) / (0.5d0 * h + bcl)
            bfm = bfv - fac
         endif
      else
         print *, "hmmat: unsupported boundary type"
         stop
      endif
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         do j = regl1, regh1
            if (mask(i-1,j) .gt. 0) then
               mat(0,i,j) = mat(0,i,j) + bfm * b(i,j)
               mat(1,i,j) = 0.d0
               if (bho .ge. 1) then
                  mat(2,i,j) = mat(2,i,j) + bfm2 * b(i,j)
               endif
            endif
         enddo
      else if (cdir .eq. 2) then
c Right face of grid
         i = regh0
         do j = regl1, regh1
            if (mask(i+1,j) .gt. 0) then
               mat(0,i,j) = mat(0,i,j) + bfm * b(i+1,j)
               mat(2,i,j) = 0.d0
               if (bho .ge. 1) then
                  mat(1,i,j) = mat(1,i,j) + bfm2 * b(i+1,j)
               endif
            endif
         enddo
      else if (cdir .eq. 1) then
c Bottom face of grid
         j = regl1
         do i = regl0, regh0
            if (mask(i,j-1) .gt. 0) then
               mat(0,i,j) = mat(0,i,j) + bfm * b(i,j)
               mat(3,i,j) = 0.d0
               if (bho .ge. 1) then
                  mat(4,i,j) = mat(4,i,j) + bfm2 * b(i,j)
               endif
            endif
         enddo
      else if (cdir .eq. 3) then
c Top face of grid
         j = regh1
         do i = regl0, regh0
            if (mask(i,j+1) .gt. 0) then
               mat(0,i,j) = mat(0,i,j) + bfm * b(i,j+1)
               mat(4,i,j) = 0.d0
               if (bho .ge. 1) then
                  mat(3,i,j) = mat(3,i,j) + bfm2 * b(i,j+1)
               endif
            endif
         enddo
      else
         print *, "hmmat: impossible face orientation"
      endif
      end

      subroutine hmmat3(mat,
     @  dims(reg),
     @  cdir, bctype, tf, bho, bcl,
     @  dims(bcv),
     @  mask, dims(msk),
     @  b, dims(bbox),
     @  beta, dx, c, r,
     @  spa, dims(spabox))
      implicit none
      integer dimdec(reg)
      integer dimdec(bcv)
      integer dimdec(msk)
      integer dimdec(bbox)
      integer dimdec(spabox)
      integer cdir, bctype, tf(dimv(bcv)), bho
      real*8 bcl, beta, dx(2), c
      real*8 mat(0:4, dimv(reg))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 spa(dimv(spabox))
      real*8 r(regl0:regh0)
      real*8 h, fac, bfm, bfv, r0
      real*8 bfm2, h2, th2
      integer i, j, bct
      if (cdir .eq. 0 .or. cdir .eq. 2) then
         h = dx(1)
c For the left and right faces, r is constant and the array actually has
c only one element.  This following "r(regl0)" looks wrong, but what it
c really means is that the array dimensions are meaningless for this case.
         r0 = r(regl0)
      else
         h = dx(2)
      endif
c The -fac * b(i,j) term applied to the matrix diagonal is the contribution
c from the interior stencil which must be removed at the boundary.
      fac = beta / (h**2)
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         do j = regl1, regh1
            if (mask(i-1,j) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i-1,j)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j)
                     bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j)
                  else
                     bfv = (beta / h) / (0.5d0 * h + bcl)
                     bfm = bfv * b(i,j)
                  endif
               else if (bct .eq. LO_NEUMANN) then
                  bfm  = 0.d0
                  bfm2 = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r0 / h
                  if (bho .ge. 1) then
                     bfm  =  0.375d0 * c * bfv
                     bfm2 = -0.125d0 * c * bfv
                  else
                     bfm = 0.25d0 * c * bfv
                  endif
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r0 / h
                  if (bho .ge. 1) then
                     bfm  =  1.5d0 * spa(i,j) * c * bfv
                     bfm2 = -0.5d0 * spa(i,j) * c * bfv
                  else
                     bfm = spa(i,j) * c * bfv
                  endif
               else
                  print *, "hmmat3: unsupported boundary type"
                  stop
               endif
               mat(0,i,j) = mat(0,i,j) + bfm - fac * b(i,j)
               mat(1,i,j) = 0.d0
               if (bho .ge. 1) then
                  mat(2,i,j) = mat(2,i,j) + bfm2
               endif
            endif
         enddo
      else if (cdir .eq. 2) then
c Right face of grid
         i = regh0
         do j = regl1, regh1
            if (mask(i+1,j) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i+1,j)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i+1,j)
                     bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i+1,j)
                  else
                     bfv = (beta / h) / (0.5d0 * h + bcl)
                     bfm = bfv * b(i+1,j)
                  endif
               else if (bct .eq. LO_NEUMANN) then
                  bfm  = 0.d0
                  bfm2 = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r0 / h
                  if (bho .ge. 1) then
                     bfm  =  0.375d0 * c * bfv
                     bfm2 = -0.125d0 * c * bfv
                  else
                     bfm = 0.25d0 * c * bfv
                  endif
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r0 / h
                  if (bho .ge. 1) then
                     bfm  =  1.5d0 * spa(i,j) * c * bfv
                     bfm2 = -0.5d0 * spa(i,j) * c * bfv
                  else
                     bfm = spa(i,j) * c * bfv
                  endif
               else
                  print *, "hmmat3: unsupported boundary type"
                  stop
               endif
               mat(0,i,j) = mat(0,i,j) + bfm - fac * b(i+1,j)
               mat(2,i,j) = 0.d0
               if (bho .ge. 1) then
                  mat(1,i,j) = mat(1,i,j) + bfm2
               endif
            endif
         enddo
      else if (cdir .eq. 1) then
c Bottom face of grid
         j = regl1
         do i = regl0, regh0
            if (mask(i,j-1) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i,j-1)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j)
                     bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j)
                  else
                     bfv = (beta / h) / (0.5d0 * h + bcl)
                     bfm = bfv * b(i,j)
                  endif
               else if (bct .eq. LO_NEUMANN) then
                  bfm  = 0.d0
                  bfm2 = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r(i) / h
                  if (bho .ge. 1) then
                     bfm  =  0.375d0 * c * bfv
                     bfm2 = -0.125d0 * c * bfv
                  else
                     bfm = 0.25d0 * c * bfv
                  endif
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r(i) / h
                  if (bho .ge. 1) then
                     bfm  =  1.5d0 * spa(i,j) * c * bfv
                     bfm2 = -0.5d0 * spa(i,j) * c * bfv
                  else
                     bfm = spa(i,j) * c * bfv
                  endif
               else
                  print *, "hmmat3: unsupported boundary type"
                  stop
               endif
               mat(0,i,j) = mat(0,i,j) + bfm - fac * b(i,j)
               mat(3,i,j) = 0.d0
               if (bho .ge. 1) then
                  mat(4,i,j) = mat(4,i,j) + bfm2
               endif
            endif
         enddo
      else if (cdir .eq. 3) then
c Top face of grid
         j = regh1
         do i = regl0, regh0
            if (mask(i,j+1) .gt. 0) then
               if (bctype .eq. -1) then
                  bct = tf(i,j+1)
               else
                  bct = bctype
               endif
               if (bct .eq. LO_DIRICHLET) then
                  if (bho .ge. 1) then
                     h2 = 0.5d0 * h
                     th2 = 3.d0 * h2
                     bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j+1)
                     bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j+1)
                  else
                     bfv = (beta / h) / (0.5d0 * h + bcl)
                     bfm = bfv * b(i,j+1)
                  endif
               else if (bct .eq. LO_NEUMANN) then
                  bfm  = 0.d0
                  bfm2 = 0.d0
               else if (bct .eq. LO_MARSHAK) then
                  bfv = 2.d0 * beta * r(i) / h
                  if (bho .ge. 1) then
                     bfm  =  0.375d0 * c * bfv
                     bfm2 = -0.125d0 * c * bfv
                  else
                     bfm = 0.25d0 * c * bfv
                  endif
               else if (bct .eq. LO_SANCHEZ_POMRANING) then
                  bfv = 2.d0 * beta * r(i) / h
                  if (bho .ge. 1) then
                     bfm  =  1.5d0 * spa(i,j) * c * bfv
                     bfm2 = -0.5d0 * spa(i,j) * c * bfv
                  else
                     bfm = spa(i,j) * c * bfv
                  endif
               else
                  print *, "hmmat3: unsupported boundary type"
                  stop
               endif
               mat(0,i,j) = mat(0,i,j) + bfm - fac * b(i,j+1)
               mat(4,i,j) = 0.d0
               if (bho .ge. 1) then
                  mat(3,i,j) = mat(3,i,j) + bfm2
               endif
            endif
         enddo
      else
         print *, "hmmat3: impossible face orientation"
      endif
      end
