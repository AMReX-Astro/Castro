
#include "LO_BCTYPES.H"

#define dims(a) a l0, a h0
#define dimdec(a) dims(a)
#define dimv(a) a l0:a h0

#define tiny 1.d-50

      subroutine hacoef(mat, a,
     @  dims(abox),
     @  dims(reg),
     @  alpha)
      implicit none
      integer dimdec(abox)
      integer dimdec(reg)
      real*8 a(dimv(abox))
      real*8 mat(0:1, dimv(reg))
      real*8 alpha
      integer i
      if (alpha .eq. 0.d0) then
         do i = regl0, regh0
            mat(1,i) = 0.d0
         enddo
      else
         do i = regl0, regh0
            mat(1,i) = alpha * a(i)
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
      real*8 mat(0:1, dimv(reg))
      real*8 beta, dx(1)
      real*8 fac
      integer i
      fac = beta / (dx(1)**2)
      do i = regl0, regh0
         mat(0,i) = - fac * b(i)
         mat(1,i) = mat(1,i) + fac * (b(i) + b(i+1))
      enddo
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
      real*8 bcl, beta, dx(1)
      real*8 mat(0:1, dimv(reg))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 h, fac, bfm, bfv
      integer i
c Boundary conditions generalized and condensed for maximum obscurity.
c Some of the b(i) (= c) factors in the Neumann and Marshak boundary
c conditions are to compensate for the 1/c in bcval.  The -fac term in
c each matrix boundary condition is the contribution from the interior
c stencil which should not be included at the boundary.
      h = dx(1)
      fac = beta / (h**2)
      if (bct .eq. LO_DIRICHLET) then
         bfv = fac * h / (0.5d0 * h + bcl)
         bfm = bfv - fac
      else if (bct .eq. LO_NEUMANN) then
         bfv = beta / h
         bfm = -fac
      else if (bct .eq. LO_MARSHAK) then
         bfv = 2.d0 * beta / h
         bfm = 0.25d0 * bfv - fac
      else
         print *, "hbmat: unsupported boundary type"
         stop
      endif
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         if (mask(i-1) .gt. 0) then
            mat(1,i) = mat(1,i) + bfm * b(i)
            mat(0,i) = 0.d0
         endif
      else if (cdir .eq. 1) then
c Right face of grid
         i = regh0
         if (mask(i+1) .gt. 0) then
            mat(1,i) = mat(1,i) + bfm * b(i+1)
         endif
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
      real*8 bcl, beta, dx(1), c
      real*8 mat(0:1, dimv(reg))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 spa(dimv(spabox))
      real*8 r(1)
      real*8 h, fac, bfm, bfv, r0
      integer i, bct
      h = dx(1)
c r is passed as an array but actually has only one element, which is
c appropriate for the face we are doing here.
      r0 = r(1)
c The -fac * b(i) term applied to the matrix diagonal is the contribution
c from the interior stencil which must be removed at the boundary.
      fac = beta / (h**2)
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         if (mask(i-1) .gt. 0) then
            if (bctype .eq. -1) then
               bct = tf(i-1)
            else
               bct = bctype
            endif
            if (bct .eq. LO_DIRICHLET) then
               bfv = fac * h / (0.5d0 * h + bcl)
               bfm = bfv * b(i)
            else if (bct .eq. LO_NEUMANN) then
               bfm = 0.d0
            else if (bct .eq. LO_MARSHAK) then
               bfv = 2.d0 * beta * r0 / h
               bfm = 0.25d0 * c * bfv
            else if (bct .eq. LO_SANCHEZ_POMRANING) then
               bfv = 2.d0 * beta * r0 / h
               bfm = spa(i) * c * bfv
            else
               print *, "hbmat3: unsupported boundary type"
               stop
            endif
            mat(1,i) = mat(1,i) + bfm - fac * b(i)
            mat(0,i) = 0.d0
         endif
      else if (cdir .eq. 1) then
c Right face of grid
         i = regh0
         if (mask(i+1) .gt. 0) then
            if (bctype .eq. -1) then
               bct = tf(i+1)
            else
               bct = bctype
            endif
            if (bct .eq. LO_DIRICHLET) then
               bfv = fac * h / (0.5d0 * h + bcl)
               bfm = bfv * b(i+1)
            else if (bct .eq. LO_NEUMANN) then
               bfm = 0.d0
            else if (bct .eq. LO_MARSHAK) then
               bfv = 2.d0 * beta * r0 / h
               bfm = 0.25d0 * c * bfv
            else if (bct .eq. LO_SANCHEZ_POMRANING) then
               bfv = 2.d0 * beta * r0 / h
               bfm = spa(i) * c * bfv
            else
               print *, "hbmat3: unsupported boundary type"
               stop
            endif
            mat(1,i) = mat(1,i) + bfm - fac * b(i+1)
         endif
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
      real*8 bcl, beta, dx(1)
      real*8 vec(dimv(reg))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 h, bfv
      real*8 h2, th2
      integer i
c Boundary conditions generalized and condensed for maximum obscurity.
c Some of the b(i) (= c) factors in the Neumann and Marshak boundary
c conditions are to compensate for the 1/c in bcval.
      h = dx(1)
      if (bct .eq. LO_DIRICHLET) then
         if (bho .ge. 1) then
            h2 = 0.5d0 * h
            th2 = 3.d0 * h2
            bfv = 2.d0 * beta / ((bcl + h2) * (bcl + th2))
         else
            bfv = (beta / h) / (0.5d0 * h + bcl)
         endif
      else if (bct .eq. LO_NEUMANN) then
         bfv = beta / h
      else if (bct .eq. LO_MARSHAK) then
         bfv = 2.d0 * beta / h
      else
         print *, "hbvec: unsupported boundary type"
         stop
      endif
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         if (mask(i-1) .gt. 0) then
            vec(i) = vec(i) + bfv * b(i) * bcval(i-1)
         endif
      else if (cdir .eq. 1) then
c Right face of grid
         i = regh0
         if (mask(i+1) .gt. 0) then
            vec(i) = vec(i) + bfv * b(i+1) * bcval(i+1)
         endif
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
      real*8 bcl, beta, dx(1)
      real*8 vec(dimv(reg))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 r(1)
      real*8 h, bfv, r0
      real*8 h2, th2
      integer i, bct
      h = dx(1)
c r is passed as an array but actually has only one element, which is
c appropriate for the face we are doing here.
      r0 = r(1)
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         if (mask(i-1) .gt. 0) then
            if (bctype .eq. -1) then
               bct = tf(i-1)
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
               bfv = bfv * b(i)
            else if (bct .eq. LO_NEUMANN) then
               bfv = beta * r0 / h
            else if (bct .eq. LO_MARSHAK .or.
     @               bct .eq. LO_SANCHEZ_POMRANING) then
               bfv = 2.d0 * beta * r0 / h
            else
               print *, "hbvec3: unsupported boundary type"
               stop
            endif
            vec(i) = vec(i) + bfv * bcval(i-1)
         endif
      else if (cdir .eq. 1) then
c Right face of grid
         i = regh0
         if (mask(i+1) .gt. 0) then
            if (bctype .eq. -1) then
               bct = tf(i+1)
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
               bfv = bfv * b(i+1)
            else if (bct .eq. LO_NEUMANN) then
               bfv = beta * r0 / h
            else if (bct .eq. LO_MARSHAK .or.
     @               bct .eq. LO_SANCHEZ_POMRANING) then
               bfv = 2.d0 * beta * r0 / h
            else
               print *, "hbvec3: unsupported boundary type"
               stop
            endif
            vec(i) = vec(i) + bfv * bcval(i+1)
         endif
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
      real*8 bcl, beta, dx(1)
      real*8 flux(dimv(fbox))
      real*8 er(dimv(ebox))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 h, bfm, bfv
      real*8 bfm2, h2, th2
      integer i
c Boundary conditions generalized and condensed for maximum obscurity.
c At present this routine is only called for internal coarse-fine
c Dirichlet conditions, not for physical boundaries.  The support for
c Neumann and Marshak boundaries is therefore never used, and the
c implementation here is not correct for the way coefficients are
c currently handled at the boundaries.  (The b coefficient should not
c be multiplied by bcval at Neumann and Marshak boundaries.  The
c implemetation here is left over from a version where b coefficients
c were modified these boundaries, and cannot be trivially fixed.
c See hbflx3 for the current correct implementation.)
      h = dx(1)
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
      else if (bct .eq. LO_NEUMANN) then
         print *, "hbflx: unsupported boundary type"
         stop
c Not correct.  See above.
         bfv = beta
         bfm = 0.d0
         bfm2 = 0.d0
      else if (bct .eq. LO_MARSHAK) then
         print *, "hbflx: unsupported boundary type"
         stop
c Not correct.  See above.
         bfv = 2.d0 * beta
         if (bho .ge. 1) then
            bfm  =  0.75d0 * beta
            bfm2 = -0.25d0 * beta
         else
            bfm = 0.5d0 * beta
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
         if (mask(i-1) .gt. 0) then
            flux(i) = b(i) * (bfv * bcval(i-1) - bfm * er(i))
            if (bho .ge. 1) then
               flux(i) = flux(i) - b(i) * bfm2 * er(i+1)
            endif
         endif
      else if (cdir .eq. 1) then
c Right face of grid
         i = regh0
         if (mask(i+1) .gt. 0) then
            flux(i+1) = -b(i+1) * (bfv * bcval(i+1) - bfm * er(i))
            if (bho .ge. 1) then
               flux(i+1) = flux(i+1) + b(i+1) * bfm2 * er(i-1)
            endif
         endif
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
      real*8 bcl, beta, dx(1), c
      real*8 flux(dimv(fbox))
      real*8 er(dimv(ebox))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 spa(dimv(spabox))
      real*8 r(1)
      real*8 h, bfm, bfv, r0
      real*8 bfm2, h2, th2
      integer i, bct
      h = dx(1)
c r is passed as an array but actually has only one element, which is
c appropriate for the face we are doing here.
      r0 = r(1)
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         if (mask(i-1) .gt. 0) then
            if (bctype .eq. -1) then
               bct = tf(i-1)
            else
               bct = bctype
            endif
            if (bct .eq. LO_DIRICHLET) then
               if (bho .ge. 1) then
                  h2 = 0.5d0 * h
                  th2 = 3.d0 * h2
                  bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i)
                  bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i)
                  bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i)
               else
                  bfv = beta / (0.5d0 * h + bcl) * b(i)
                  bfm = bfv
               endif
            else if (bct .eq. LO_NEUMANN) then
               bfv  = beta * r0
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
                  bfm  =  1.5d0 * spa(i) * c * bfv
                  bfm2 = -0.5d0 * spa(i) * c * bfv
               else
                  bfm = spa(i) * c * bfv
               endif
            else
               print *, "hbflx3: unsupported boundary type"
               stop
            endif
            if (inhom .eq. 0) then
               bfv = 0.d0
            endif
            flux(i) = (bfv * bcval(i-1) - bfm * er(i))
            if (bho .ge. 1) then
               flux(i) = flux(i) - bfm2 * er(i+1)
            endif
         endif
      else if (cdir .eq. 1) then
c Right face of grid
         i = regh0
         if (mask(i+1) .gt. 0) then
            if (bctype .eq. -1) then
               bct = tf(i+1)
            else
               bct = bctype
            endif
            if (bct .eq. LO_DIRICHLET) then
               if (bho .ge. 1) then
                  h2 = 0.5d0 * h
                  th2 = 3.d0 * h2
                  bfv = 2.d0 * beta * h / ((bcl + h2) * (bcl + th2)) * b(i+1)
                  bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i+1)
                  bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i+1)
               else
                  bfv = beta / (0.5d0 * h + bcl) * b(i+1)
                  bfm = bfv
               endif
            else if (bct .eq. LO_NEUMANN) then
               bfv  = beta * r0
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
                  bfm  =  1.5d0 * spa(i) * c * bfv
                  bfm2 = -0.5d0 * spa(i) * c * bfv
               else
                  bfm = spa(i) * c * bfv
               endif
            else
               print *, "hbflx3: unsupported boundary type"
               stop
            endif
            if (inhom .eq. 0) then
               bfv = 0.d0
            endif
            flux(i+1) = -(bfv * bcval(i+1) - bfm * er(i))
            if (bho .ge. 1) then
               flux(i+1) = flux(i+1) + bfm2 * er(i-1)
            endif
         endif
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
      real*8 bcl, dx(1)
      real*8 dterm(dimv(dtbox))
      real*8 er(dimv(ebox))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 d(dimv(dbox))
      real*8 h 
      integer i
      h = dx(1)
      if (bct .eq. LO_DIRICHLET) then
         if (cdir .eq. 0) then
c     Left face of grid
            i = regl0
            if (mask(i-1) .gt. 0) then
               dterm(i) = d(i)*(er(i) - bcval(i-1))/(0.5d0*h+bcl)
            endif
         else if (cdir .eq. 1) then
c     Right face of grid
            i = regh0
            if (mask(i+1) .gt. 0) then
               dterm(i+1) = d(i+1)*(bcval(i+1)-er(i))/(0.5d0*h+bcl)
            endif
         else
            print *, "hdterm: impossible face orientation"
         endif
      else
         print *, "hdterm: unsupported boundary type"
         stop
      endif
      end

      subroutine hdterm3(dterm,
     @     dims(dtbox),
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
      real*8 bcl, dx(1) 
      real*8 dterm(dimv(dtbox))
      real*8 er(dimv(ebox))
      real*8 bcval(dimv(bcv))
      integer mask(dimv(msk))
      real*8 d(dimv(dbox))
      real*8 h 
      integer i, bct
      h = dx(1)
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         if (mask(i-1) .gt. 0) then
            if (bctype .eq. -1) then
               bct = tf(i-1)
            else
               bct = bctype
            endif
            if (bct .eq. LO_DIRICHLET) then
               dterm(i) = d(i)*(er(i) - bcval(i-1))/(0.5d0*h+bcl)
            else if (bct.eq.LO_NEUMANN .and. bcval(i-1).eq.0.d0) then
               dterm(i) = 0.d0
            else
               print *, "hdterm3: unsupported boundary type"
               stop
            endif
         endif
      else if (cdir .eq. 1) then
c Right face of grid
         i = regh0
         if (mask(i+1) .gt. 0) then
            if (bctype .eq. -1) then
               bct = tf(i+1)
            else
               bct = bctype
            endif
            if (bct .eq. LO_DIRICHLET) then
               dterm(i+1) = d(i+1)*(bcval(i+1)-er(i))/(0.5d0*h+bcl)
            else if (bct.eq.LO_NEUMANN .and. bcval(i+1).eq.0.d0) then
               dterm(i+1) = 0.d0
            else
               print *, "hbterm3: unsupported boundary type"
               stop
            endif
         endif
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
      real*8 mat(0:2, dimv(reg))
      real*8 alpha
      integer i
      if (alpha .eq. 0.d0) then
         do i = regl0, regh0
            mat(0,i) = 0.d0
         enddo
      else
         do i = regl0, regh0
            mat(0,i) = alpha * a(i)
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
      real*8 mat(0:2, dimv(reg))
      real*8 beta, dx(1)
      real*8 fac
      integer i
      fac = beta / (dx(1)**2)
      do i = regl0, regh0
         mat(0,i) = mat(0,i) + fac * (b(i) + b(i+1))
         mat(1,i) = - fac * b(i)
         mat(2,i) = - fac * b(i+1)
      enddo
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
      real*8 mat(0:2, dimv(reg))
      real*8 alpha2
      real*8 fac
      integer i
      fac = 0.25d0 * alpha2
      do i = regl0, regh0
         mat(0,i) = mat(0,i) + fac * (a2(i) + a2(i+1))
         mat(1,i) = mat(1,i) + fac * a2(i)
         mat(2,i) = mat(2,i) + fac * a2(i+1)
      enddo
      end

      subroutine hmcc(mat, c,
     @  dims(cbox),
     @  dims(reg),
     @  gamma, dx, n)
      implicit none
      integer dimdec(cbox)
      integer dimdec(reg)
      integer n
c c has two components for upwinding, independent of dimension
      real*8 c(dimv(cbox), 0:1)
      real*8 mat(0:2, dimv(reg))
      real*8 gamma, dx(1)
      real*8 fac
      integer i
      fac = gamma / dx(1)
      do i = regl0, regh0
         mat(0,i) = mat(0,i) - fac * (c(i,1) - c(i+1,0))
         mat(1,i) = mat(1,i) - fac * c(i,0)
         mat(2,i) = mat(2,i) + fac * c(i+1,1)
      enddo
      end

      subroutine add_ccoef_flux(dir,
     @  den, dims(den),
     @  c, dims(cbox),
     @  gamma,
     @  dx,
     @  flux, dims(fbox))
      implicit none

      integer dimdec(den)
      integer dimdec(cbox)
      integer dimdec(fbox)

      real*8 den(dimv(den))
      real*8 c(dimv(cbox), 0:1)
      real*8 flux(dimv(fbox))

      integer dir
      real*8 gamma, dx(1)
      integer i

      if (dir .eq. 0) then

c...     x-direction

         do i = fboxl0, fboxh0
            flux(i) = flux(i) + gamma * (c(i,0) * den(i-1) +
     @                                   c(i,1) * den(i))
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
      real*8 mat(0:2, dimv(reg))
      real*8 delta1, dx(1)
      real*8 fac
      integer i
      fac = 0.5d0 * delta1 / dx(1)
      do i = regl0, regh0
         mat(1,i) = mat(1,i) - fac * d1(i)
         mat(2,i) = mat(2,i) + fac * d1(i)
      enddo
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
      real*8 mat(0:2, dimv(reg))
      real*8 delta2, dx(1)
      real*8 fac
      integer i
      fac = 0.5d0 * delta2 / dx(1)
      do i = regl0, regh0
         mat(0,i) = mat(0,i) + fac * (d2(i) - d2(i+1))
         mat(1,i) = mat(1,i) - fac * d2(i)
         mat(2,i) = mat(2,i) + fac * d2(i+1)
      enddo
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
      real*8 bcl, beta, dx(1)
      real*8 mat(0:2, dimv(reg))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 h, fac, bfm, bfv
      real*8 bfm2, h2, th2
      integer i
c Boundary conditions generalized and condensed for maximum obscurity.
c Some of the b(i) (= c) factors in the Neumann and Marshak boundary
c conditions are to compensate for the 1/c in bcval.  The -fac term in
c each matrix boundary condition is the contribution from the interior
c stencil which should not be included at the boundary.
      h = dx(1)
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
      else if (bct .eq. LO_NEUMANN) then
         bfm = -fac
         bfm2 = 0.d0
      else if (bct .eq. LO_MARSHAK) then
         bfv = 2.d0 * beta / h
         if (bho .ge. 1) then
            bfm  =  0.375d0 * bfv - fac
            bfm2 = -0.125d0 * bfv
         else
            bfm = 0.25d0 * bfv - fac
         endif
      else
         print *, "hmmat: unsupported boundary type"
         stop
      endif
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         if (mask(i-1) .gt. 0) then
            mat(0,i) = mat(0,i) + bfm * b(i)
            mat(1,i) = 0.d0
            if (bho .ge. 1) then
               mat(2,i) = mat(2,i) + bfm2 * b(i)
            endif
         endif
      else if (cdir .eq. 1) then
c Right face of grid
         i = regh0
         if (mask(i+1) .gt. 0) then
            mat(0,i) = mat(0,i) + bfm * b(i+1)
            mat(2,i) = 0.d0
            if (bho .ge. 1) then
               mat(1,i) = mat(1,i) + bfm2 * b(i+1)
            endif
         endif
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
      real*8 bcl, beta, dx(1), c
      real*8 mat(0:2, dimv(reg))
      integer mask(dimv(msk))
      real*8 b(dimv(bbox))
      real*8 spa(dimv(spabox))
      real*8 r(1)
      real*8 h, fac, bfm, bfv, r0
      real*8 bfm2, h2, th2
      integer i, bct
      h = dx(1)
c r is passed as an array but actually has only one element, which is
c appropriate for the face we are doing here.
      r0 = r(1)
c The -fac * b(i) term applied to the matrix diagonal is the contribution
c from the interior stencil which must be removed at the boundary.
      fac = beta / (h**2)
      if (cdir .eq. 0) then
c Left face of grid
         i = regl0
         if (mask(i-1) .gt. 0) then
            if (bctype .eq. -1) then
               bct = tf(i-1)
            else
               bct = bctype
            endif
            if (bct .eq. LO_DIRICHLET) then
               if (bho .ge. 1) then
                  h2 = 0.5d0 * h
                  th2 = 3.d0 * h2
                  bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i)
                  bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i)
               else
                  bfv = (beta / h) / (0.5d0 * h + bcl)
                  bfm = bfv * b(i)
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
                  bfm  =  1.5d0 * spa(i) * c * bfv
                  bfm2 = -0.5d0 * spa(i) * c * bfv
               else
                  bfm = spa(i) * c * bfv
               endif
            else
               print *, "hmmat3: unsupported boundary type"
               stop
            endif
            mat(0,i) = mat(0,i) + bfm - fac * b(i)
            mat(1,i) = 0.d0
            if (bho .ge. 1) then
               mat(2,i) = mat(2,i) + bfm2
            endif
         endif
      else if (cdir .eq. 1) then
c Right face of grid
         i = regh0
         if (mask(i+1) .gt. 0) then
            if (bctype .eq. -1) then
               bct = tf(i+1)
            else
               bct = bctype
            endif
            if (bct .eq. LO_DIRICHLET) then
               if (bho .ge. 1) then
                  h2 = 0.5d0 * h
                  th2 = 3.d0 * h2
                  bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i+1)
                  bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i+1)
               else
                  bfv = (beta / h) / (0.5d0 * h + bcl)
                  bfm = bfv * b(i+1)
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
                  bfm  =  1.5d0 * spa(i) * c * bfv
                  bfm2 = -0.5d0 * spa(i) * c * bfv
               else
                  bfm = spa(i) * c * bfv
               endif
            else
               print *, "hmmat3: unsupported boundary type"
               stop
            endif
            mat(0,i) = mat(0,i) + bfm - fac * b(i+1)
            mat(2,i) = 0.d0
            if (bho .ge. 1) then
               mat(1,i) = mat(1,i) + bfm2
            endif
         endif
      else
         print *, "hmmat3: impossible face orientation"
      endif
      end
