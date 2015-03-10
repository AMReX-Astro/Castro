
#include "LO_BCTYPES.H"

#define dims(a) a l0, a l1, a l2, a h0, a h1, a h2
#define dimdec(a) dims(a)
#define dimv(a) a l0:a h0,a l1:a h1,a l2:a h2

#define tiny 1.d-50
c#define big  1.d+50

#define BIGKR 1.d25
c#define BIGKP 1.d0

#define REAL_T real*8

#include "FLD_limiter.F"

      subroutine setlayout(Den_in, Xmom_in, Eden_in, Eint_in, Temp_in, 
     @     FirstSpec, FirstAux, NUM_STATE)
      implicit none
      INCLUDE 'StateLayout.H'
      integer Den_in, Xmom_in, Eden_in, Eint_in, Temp_in
      integer FirstSpec, FirstAux, NUM_STATE
      DEN = Den_in
      XMOM = Xmom_in
      YMOM = Xmom_in+1
      ZMOM = Xmom_in+2
      EDEN = Eden_in
      EINT = Eint_in
      ITEMP = Temp_in
      IFS = FirstSpec
      IFX = FirstAux
      SVSIZE = NUM_STATE
      end

c The following routines implement metric terms in 2D and are included
c in the 3D source only to enable the code to link.

      subroutine multrs(d,
     @  dims(dbox),
     @  dims(reg),
     @  r, s)
      implicit none
      integer dimdec(dbox)
      integer dimdec(reg)
      real*8 d(dimv(dbox))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               d(i,j,k) = d(i,j,k) * r(i) * s(j)
            enddo
         enddo
      enddo
      end

      subroutine divr(d,
     @  dims(dbox),
     @  dims(reg),
     @  r)
      implicit none
      integer dimdec(dbox)
      integer dimdec(reg)
      real*8 d(dimv(dbox))
      real*8 r(dboxl0:dboxh0)
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               d(i,j,k) = d(i,j,k) / r(i)
            enddo
         enddo
      enddo
      end

      subroutine divrs(d,
     @  dims(dbox),
     @  dims(reg),
     @  r, s)
      implicit none
      integer dimdec(dbox)
      integer dimdec(reg)
      real*8 d(dimv(dbox))
      real*8 r(dboxl0:dboxh0)
      real*8 s(dboxl1:dboxh1)
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               d(i,j,k) = d(i,j,k) / (r(i) * s(j))
            enddo
         enddo
      enddo
      end

      subroutine sphc(r, s,
     @  dims(reg), dx)
      implicit none
      integer dimdec(reg)
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 dx(2)
      real*8 h1, h2, d1, d2
      integer i, j
      h1 = 0.5d0 * dx(1)
      h2 = 0.5d0 * dx(2)
      d1 = 1.d0 / (3.d0 * dx(1))
      d2 = 1.d0 / dx(2)
      do i = regl0, regh0
         r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
      enddo
      do j = regl1, regh1
         s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
      enddo
      end

      subroutine sphe(r, s, n,
     @  dims(reg), dx)
      implicit none
      integer dimdec(reg)
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      integer n
      real*8 dx(2)
      real*8 h1, h2, d1, d2
      integer i, j
      if (n .eq. 0) then
         do i = regl0, regh0
            r(i) = r(i)**2
         enddo
         h2 = 0.5d0 * dx(2)
         d2 = 1.d0 / dx(2)
         do j = regl1, regh1
            s(j) = d2 * (cos(s(j) - h2) - cos(s(j) + h2))
         enddo
      else
         h1 = 0.5d0 * dx(1)
         d1 = 1.d0 / (3.d0 * dx(1))
         do i = regl0, regh0
            r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
         enddo
         do j = regl1, regh1
            s(j) = sin(s(j))
         enddo
      endif
      end

      subroutine macoef(a,
     @  dims(abox),
     @  dims(reg),
     @  fkp, r, s, c)
      implicit none
      integer dimdec(abox)
      integer dimdec(reg)
      real*8 a(dimv(abox))
      real*8 fkp(dimv(reg))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 c
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               a(i,j,k) = r(i) * s(j) * fkp(i,j,k) * c
            enddo
         enddo
      enddo
      end

      subroutine mrhs(rhs,
     @  dims(rbox),
     @  dims(reg),
     @  temp, fkp, ert, r, s,
     @  sigma)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      real*8 rhs(dimv(rbox))
      real*8 temp(dimv(reg))
      real*8 fkp(dimv(reg))
      real*8 ert(dimv(reg))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 sigma
      integer i, j, k
      real*8 bs
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               bs = 4.d0 * sigma * fkp(i,j,k) * temp(i,j,k)**4
               rhs(i,j,k) = r(i) * s(j) * (bs - ert(i,j,k))
            enddo
         enddo
      enddo
      end

      subroutine sacoef(a,
     @  dims(abox),
     @  dims(reg),
     @  fkp, etainv, r, s, c, dt)
      implicit none
      integer dimdec(abox)
      integer dimdec(reg)
      real*8 a(dimv(abox))
      real*8 fkp(dimv(reg))
      real*8 etainv(dimv(reg))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 c, dt
      integer i, j, k
      real*8 dtm, eki
      dtm = 1.d0 / dt
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               eki = fkp(i,j,k) * etainv(i,j,k)
               a(i,j,k) = r(i) * s(j) * (eki * c + dtm)
            enddo
         enddo
      enddo
      end

      subroutine lacoef(a,
     @  dims(abox),
     @  dims(reg),
     @  fkp, eta, etainv, r, s, c, dt, theta)
      implicit none
      integer dimdec(abox)
      integer dimdec(reg)
      real*8 a(dimv(abox))
      real*8 fkp(dimv(abox))
      real*8 eta(dimv(abox))
      real*8 etainv(dimv(abox))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 c, dt, theta
      integer i, j, k
      real*8 dtm
      dtm = 1.d0 / dt
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               a(i,j,k) = r(i) * s(j) *
     @                    (fkp(i,j,k) * etainv(i,j,k) * c + dtm) /
     @                    (1.d0 - (1.d0 - theta) * eta(i,j,k))
            enddo
         enddo
      enddo
      end

c arithmetic average, geometrically correct(?) but underestimates surface flux

#define KAVG0(a,b) (0.5d0 * (a + b + tiny))

c#define KAVG(a,b,d) KAVG0(a,b)

c harmonic average, overestimates surface flux

#define KAVG1(a,b) ((2.d0 * a * b) / (a + b + tiny) + tiny)

c#define KAVG(a,b,d) KAVG1(a,b)

c chooses arithmetic where optically thin, harmonic where optically thick,
c surface flux approximation at a thick/thin boundary

#define KAVG2(a,b,d) min(KAVG0(a,b), max(KAVG1(a,b), 4.d0 / (3.d0 * d)))

c#define KAVG(a,b,d) KAVG2(a,b,d)

#define KAVG(a,b,d) kavg(a,b,d,-1)

      real*8 function kavg(a, b, d, iopt)
      implicit none
      real*8 a, b, d
      integer iopt, opt
      save opt
      if (iopt .ge. 0) then
         opt = iopt
         if (opt .gt. 2) then
            print *, "Fortran KAVG: invalid averaging option"
         endif
      endif
      if (opt .eq. 0) then
         kavg = KAVG0(a,b)
      else if (opt .eq. 1) then
         kavg = KAVG1(a,b)
      else
         kavg = KAVG2(a,b,d)
      endif
      end

      subroutine bclim(b,
     @  lambda, dims(bbox),
     @  dims(reg),
     @  n, kappar, dims(kbox),
     @  r, s, c, dx)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer n
      real*8 b(dimv(bbox))
      real*8 lambda(dimv(bbox))
      real*8 kappar(dimv(kbox))
      real*8 r(regl0:regh0+1)
      real*8 s(regl1:regh1+1)
      real*8 c, dx(3)
      real*8 kavg
      integer i, j, k
      real*8 kap
      if (n .eq. 0) then
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0 + 1
                  kap = KAVG(kappar(i-1,j,k), kappar(i,j,k), dx(1))
                  b(i,j,k) = r(i) * s(j) * c * lambda(i,j,k) / kap
               enddo
            enddo
         enddo
      else if (n .eq. 1) then
         do k = regl2, regh2
            do j = regl1, regh1 + 1
               do i = regl0, regh0
                  kap = KAVG(kappar(i,j-1,k), kappar(i,j,k), dx(2))
                  b(i,j,k) = r(i) * s(j) * c * lambda(i,j,k) / kap
               enddo
            enddo
         enddo
      else 
         do k = regl2, regh2 + 1
            do j = regl1, regh1 
               do i = regl0, regh0
                  kap = KAVG(kappar(i,j,k-1), kappar(i,j,k), dx(3))
                  b(i,j,k) = r(i) * s(j) * c * lambda(i,j,k) / kap
               enddo
            enddo
         enddo
      endif
      end

      subroutine flxlim(lambda,
     @  dims(rbox),
     @  dims(reg), limiter)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer limiter
      real*8 lambda(dimv(rbox)), FLDlambda
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               lambda(i,j,k) = FLDlambda(lambda(i,j,k),limiter)
            enddo
         enddo
      enddo
      end

      subroutine eddfac(efact,
     @  dims(rbox),
     @  dims(reg), limiter, n)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer n, limiter
      real*8 efact(dimv(rbox))
      integer i, j, k
      real*8 r, lambda, FLDlambda
      integer dir(0:2)
      dir(0) = 0
      dir(1) = 0
      dir(2) = 0
      dir(n) = 1
      do k = regl2, regh2 + dir(2)
         do j = regl1, regh1 + dir(1)
            do i = regl0, regh0 + dir(0)
               r = efact(i,j,k)
               lambda = FLDlambda(r,limiter)
               efact(i,j,k) = lambda + (lambda * r)**2
            enddo
         enddo
      enddo
      end

      subroutine scgrd1(r,
     @  dims(rbox),
     @  dims(reg),
     @  n, kappar, dims(kbox),
     @  er, dx)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer n
      real*8 r(dimv(rbox))
      real*8 kappar(dimv(kbox))
      real*8 er(dimv(kbox))
      real*8 dx(3)
      real*8 kavg
      integer i, j, k
      real*8 kap
      if (n .eq. 0) then
         do k = regl2, regh2
            do j = regl1, regh1
c     x derivatives, gradient assembly:
               do i = regl0, regh0 + 1
                  r(i,j,k) = abs(er(i,j,k) - er(i-1,j,k)) / dx(1)
               enddo
               i = regl0
               if (er(i-1,j,k) .eq. -1.d0) then
                  r(i,j,k) = abs(er(i+1,j,k) - er(i,j,k)) / dx(1)
               endif
               i = regh0 + 1
               if (er(i,j,k) .eq. -1.d0) then
                  r(i,j,k) = abs(er(i-1,j,k) - er(i-2,j,k)) / dx(1)
               endif
c     construct coefficients:
               do i = regl0, regh0 + 1
                  kap = KAVG(kappar(i-1,j,k), kappar(i,j,k), dx(1))
                  r(i,j,k) = r(i,j,k) /
     @                 (kap * max(er(i-1,j,k), er(i,j,k), tiny))
               enddo
            enddo
         enddo
      else if (n .eq. 1) then
         do k = regl2, regh2 
c     y derivatives, gradient assembly:
            do j = regl1, regh1 + 1
               do i = regl0, regh0
                  r(i,j,k) = abs(er(i,j,k) - er(i,j-1,k)) / dx(2)
               enddo
            enddo
            do i = regl0, regh0
               j = regl1
               if (er(i,j-1,k) .eq. -1.d0) then
                  r(i,j,k) = abs(er(i,j+1,k) - er(i,j,k)) / dx(2)
               endif
               j = regh1 + 1
               if (er(i,j,k) .eq. -1.d0) then
                  r(i,j,k) = abs(er(i,j-1,k) - er(i,j-2,k)) / dx(2)
               endif
            enddo
c construct coefficients:
            do j = regl1, regh1 + 1
               do i = regl0, regh0
                  kap = KAVG(kappar(i,j-1,k), kappar(i,j,k), dx(2))
                  r(i,j,k) = r(i,j,k) /
     @                 (kap * max(er(i,j-1,k), er(i,j,k), tiny))
               enddo
            enddo
         enddo
      else 
c     z derivatives, gradient assembly:
         do k = regl2, regh2 +1
            do j = regl1, regh1 
               do i = regl0, regh0
                  r(i,j,k) = abs(er(i,j,k-1) - er(i,j,k)) / dx(3)
               enddo
            enddo
         enddo
         do j = regl1, regh1 
            do i = regl0, regh0
               k = regl2
               if (er(i,j,k-1) .eq. -1.d0) then
                  r(i,j,k) = abs(er(i,j,k+1) - er(i,j,k)) / dx(3)
               endif
               k = regh2 + 1
               if (er(i,j,k) .eq. -1.d0) then
                  r(i,j,k) = abs(er(i,j,k-2) - er(i,j,k-1)) / dx(3)
               endif
            enddo
         enddo
c construct coefficients:
         do k = regl2, regh2 + 1
            do j = regl1, regh1 
               do i = regl0, regh0
                  kap = KAVG(kappar(i,j,k-1), kappar(i,j,k), dx(3))
                  r(i,j,k) = r(i,j,k) /
     @                 (kap * max(er(i,j,k-1), er(i,j,k), tiny))
               enddo
            enddo
         enddo
      endif
      end

      subroutine scgrd2(r,
     @  dims(rbox),
     @  dims(reg),
     @  n, kappar, dims(kbox), er,
     @  dims(dbox), da, db, dx)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer dimdec(dbox)
      integer n
      real*8 r(dimv(rbox))
      real*8 kappar(dimv(kbox))
      real*8 er(dimv(kbox))
      real*8 da(dimv(dbox))
      real*8 db(dimv(dbox))
      real*8 dx(3)
      real*8 kavg
      integer i, j, k
      real*8 kap
      if (n .eq. 0) then
c     y & z derivatives:
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0 - 1, regh0 + 1
                  da(i,j,k) = er(i,j+1,k  ) - er(i,j-1,k  )
                  db(i,j,k) = er(i,j  ,k+1) - er(i,j  ,k-1)                  
               enddo
            enddo
         enddo                  

c     check y derivatives
         do k = regl2, regh2
            do i = regl0 - 1, regh0 + 1
               j = regl1
               if (er(i,j-1,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
               endif
               j = regh1
               if (er(i,j+1,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
               endif
            enddo
            do j = regl1, regh1
               i = regl0 - 1
c     (check at j-1 and j+1.  if value at j is bad it will not be used at all.)
               if (er(i,j-1,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
               else if (er(i,j+1,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
               endif
               i = regh0 + 1
               if (er(i,j-1,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
               else if (er(i,j+1,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
               endif
            enddo
         enddo                  

c     check z derivatives
         do j = regl1, regh1
            do i = regl0-1, regh0 + 1
               k = regl2
               if (er(i,j,k-1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
               endif
               k = regh2
               if (er(i,j,k+1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
               endif
            enddo
         enddo
         do k = regl2, regh2
            do j = regl1, regh1
               i = regl0 - 1
c     (check at k-1 and k+1.  if value at k is bad it will not be used at all.)
               if (er(i,j,k-1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
               else if (er(i,j,k+1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
               endif
               i = regh0 + 1
               if (er(i,j,k-1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
               else if (er(i,j,k+1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
               endif
            enddo
         enddo
c     x derivatives 
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0 + 1
                  r(i,j,k) = ((er(i,j,k) - er(i-1,j,k)) / dx(1)) ** 2 
     @                 + ((da(i-1,j,k) + da(i,j,k)) / 
     @                                 (4.d0 * dx(2))) ** 2
     @                 + ((db(i-1,j,k) + db(i,j,k)) / 
     @                                 (4.d0 * dx(3))) ** 2
               enddo
               i = regl0
               if (er(i-1,j,k) .eq. -1.d0) then
                  r(i,j,k) = ((er(i+1,j,k) - er(i,j,k)) / dx(1)) ** 2
     @                 + (da(i,j,k) / (2.d0 * dx(2))) ** 2
     @                 + (db(i,j,k) / (2.d0 * dx(3))) ** 2
               endif
               i = regh0 + 1
               if (er(i,j,k) .eq. -1.d0) then
                  r(i,j,k) = ((er(i-1,j,k) - er(i-2,j,k)) / dx(1)) ** 2
     @                 + (da(i-1,j,k) / (2.d0 * dx(2))) ** 2
     @                 + (db(i-1,j,k) / (2.d0 * dx(3))) ** 2
               endif
c     construct scaled gradient
               do i = regl0, regh0 + 1
                  kap = KAVG(kappar(i-1,j,k), kappar(i,j,k), dx(1))
                  r(i,j,k) = sqrt(r(i,j,k)) /
     @                 (kap * max(er(i-1,j,k), er(i,j,k), tiny))
               enddo
            enddo
         enddo
      else if (n .eq. 1) then
c     x & z derivatives:
         do k = regl2, regh2
            do j = regl1-1, regh1+1
               do i = regl0, regh0
                  da(i,j,k) = er(i+1,j,k) - er(i-1,j,k)
                  db(i,j,k) = er(i,j,k+1) - er(i,j,k-1)
               enddo
            enddo
         enddo

c     check x derivatives
         do k = regl2, regh2
            do j = regl1-1, regh1+1
               i = regl0
               if (er(i-1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
               endif
               i = regh0
               if (er(i+1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i,j,k) - er(i-1,j,k))
               endif
            enddo
            do i = regl0, regh0
               j = regl1-1
c     (check at i-1 and i+1.  if value at i is bad it will not be used at all.)
               if (er(i-1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
               else if (er(i+1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i,j,k) - er(i-1,j,k))
               endif
               j = regh1+1
               if (er(i-1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
               else if (er(i+1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i,j,k) - er(i-1,j,k))
               endif
            enddo
         enddo

c     check z derivatives
         do j = regl1-1, regh1+1
            do i = regl0, regh0
               k = regl2
               if (er(i,j,k-1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
               endif
               k = regh2
               if (er(i,j,k+1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
               endif               
            enddo
         enddo
         do k = regl2, regh2
            do i = regl0, regh0
               j = regl1-1
c     (check at k-1 and k+1.  if value at k is bad it will not be used at all.)
               if (er(i,j,k-1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
               else if (er(i,j,k+1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
               endif
               j = regh1+1
               if (er(i,j,k-1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k+1) - er(i,j,k))
               else if (er(i,j,k+1) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j,k-1))
               endif
            enddo
         enddo

c     y derivatives 
         do k = regl2, regh2
            do j = regl1, regh1 + 1
               do i = regl0, regh0
                  r(i,j,k) = ((er(i,j,k) - er(i,j-1,k)) / dx(2)) ** 2
     @                 + ((da(i,j-1,k) + da(i,j,k)) /
     @                                 (4.d0 * dx(1))) ** 2
     @                 + ((db(i,j-1,k) + db(i,j,k)) /
     @                                 (4.d0 * dx(3))) ** 2
               enddo
            enddo
            do i = regl0, regh0
               j = regl1
               if (er(i,j-1,k) .eq. -1.d0) then
                  r(i,j,k) = ((er(i,j+1,k) - er(i,j,k)) / dx(2)) ** 2
     @                 + (da(i,j,k) / (2.d0 * dx(1))) ** 2
     @                 + (db(i,j,k) / (2.d0 * dx(3))) ** 2
               endif
               j = regh1 + 1
               if (er(i,j,k) .eq. -1.d0) then
                  r(i,j,k) = ((er(i,j-1,k) - er(i,j-2,k)) / dx(2)) ** 2
     @                 + (da(i,j-1,k) / (2.d0 * dx(1))) ** 2
     @                 + (db(i,j-1,k) / (2.d0 * dx(3))) ** 2
               endif
            enddo
c     construct scaled gradient
            do j = regl1, regh1 + 1
               do i = regl0, regh0
                  kap = KAVG(kappar(i,j-1,k), kappar(i,j,k), dx(2))
                  r(i,j,k) = sqrt(r(i,j,k)) /
     @                 (kap * max(er(i,j-1,k), er(i,j,k), tiny))
               enddo
            enddo
         enddo
      else
c     x & y derivatives
         do k = regl2-1, regh2+1
            do j = regl1, regh1
               do i = regl0, regh0
                  da(i,j,k) = er(i+1,j,k) - er(i-1,j,k)
                  db(i,j,k) = er(i,j+1,k) - er(i,j-1,k)
               enddo
            enddo
         enddo                  

c     check x derivatives
         do k = regl2-1, regh2+1
            do j = regl1, regh1
               i = regl0
               if (er(i-1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
               endif
               i = regh0
               if (er(i+1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.0d0 * (er(i,j,k) - er(i-1,j,k))
               endif
            enddo
         enddo     
         do j = regl1, regh1
            do i = regl0, regh0
               k = regl2 - 1
c     (check at i-1 and i+1.  if value at i is bad it will not be used at all.)
               if (er(i-1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
               else if (er(i+1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.0d0 * (er(i,j,k) - er(i-1,j,k))
               endif
               k = regh2 + 1
               if (er(i-1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.d0 * (er(i+1,j,k) - er(i,j,k))
               else if (er(i+1,j,k) .eq. -1.d0) then
                  da(i,j,k) = 2.0d0 * (er(i,j,k) - er(i-1,j,k))
               endif
            enddo
         enddo

c     check y derivative
         do k = regl2-1, regh2+1
            do i = regl0, regh0
               j = regl1
               if (er(i,j-1,k) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
               endif
               j = regh1
               if (er(i,j+1,k) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
               endif
            enddo
         enddo
         do j = regl1, regh1
            do i = regl0, regh0
               k = regl2 - 1
c     (check at j-1 and j+1.  if value at j is bad it will not be used at all.)
               if (er(i,j-1,k) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
               else if (er(i,j+1,k) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
               endif
               k = regh2 + 1
               if (er(i,j-1,k) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j+1,k) - er(i,j,k))
               else if (er(i,j+1,k) .eq. -1.d0) then
                  db(i,j,k) = 2.d0 * (er(i,j,k) - er(i,j-1,k))
               endif
            enddo
         enddo

c z derivatives
         do k = regl2, regh2+1
            do j = regl1, regh1
               do i = regl0, regh0 
                  r(i,j,k) = ((er(i,j,k) - er(i,j,k-1)) / dx(3)) ** 2
     @                 + ((da(i,j,k) + da(i,j,k-1)) /
     @                               (4.d0 * dx(1))) ** 2
     @                 + ((db(i,j,k) + db(i,j,k-1)) /
     @                               (4.d0 * dx(2))) ** 2
               enddo
            enddo
         enddo

         do j = regl1, regh1
            do i = regl0, regh0 
               k = regl2
               if (er(i,j,k-1) .eq. -1.d0) then
                  r(i,j,k) = ((er(i,j,k+1) - er(i,j,k)) / dx(3)) ** 2
     @                 + (da(i,j,k) / (2.d0 * dx(1))) ** 2
     @                 + (db(i,j,k) / (2.d0 * dx(2))) ** 2
               endif
               k = regh2+1
               if (er(i,j,k) .eq. -1.d0) then
                  r(i,j,k) = ((er(i,j,k-1) - er(i,j,k-2)) / dx(3)) ** 2
     @                 + (da(i,j,k-1) / (2.d0 * dx(1))) ** 2
     @                 + (db(i,j,k-1) / (2.d0 * dx(2))) ** 2
               endif
            enddo
         enddo

c     gradient assembly
         do k = regl2, regh2+1
            do j = regl1, regh1
               do i = regl0, regh0 
                  kap = KAVG(kappar(i,j,k), kappar(i,j,k-1), dx(3))
                  r(i,j,k) = sqrt(r(i,j,k)) /
     @                 (kap * max(er(i,j,k-1), er(i,j,k), tiny))
               enddo
            enddo
         enddo
      endif
      end

      subroutine scgrd3(r,
     @  dims(rbox),
     @  dims(reg),
     @  n, kappar, dims(kbox), er,
     @  dims(dbox), da, db, dx)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer dimdec(dbox)
      integer n
      real*8 r(dimv(rbox))
      real*8 kappar(dimv(kbox))
      real*8 er(dimv(kbox))
      real*8 da(dimv(dbox))
      real*8 db(dimv(dbox))
      real*8 dx(3)
      real*8 kavg
      integer i
      real*8 kap
      print *, "scgrd3 not implemented in 3d"
      stop
      end

      subroutine lrhs(rhs,
     @  dims(rbox),
     @  dims(reg),
     @  temp, fkp, eta, etainv, frhoem, frhoes, dfo, 
     @  ero, dims(ebox), edot,
     @  r, s, dt, sigma, c, theta)
      implicit none
      integer dimdec(rbox)
      integer dimdec(ebox)
      integer dimdec(reg)
      real*8 rhs(dimv(rbox))
      real*8 temp(dimv(reg))
      real*8 fkp(dimv(reg))
      real*8 eta(dimv(reg))
      real*8 etainv(dimv(reg))
      real*8 frhoem(dimv(reg))
      real*8 frhoes(dimv(reg))
      real*8 dfo(dimv(reg))
      real*8 ero(dimv(ebox))
      real*8 edot(dimv(reg))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 dt, sigma, c, theta
      integer i, j, k
      real*8 dtm, ek, bs, es, ekt
      dtm = 1.d0 / dt
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               ek = fkp(i,j,k) * eta(i,j,k)
               bs = etainv(i,j,k) *
     @              4.d0 * sigma * fkp(i,j,k) * temp(i,j,k)**4
               es = eta(i,j,k) * (frhoem(i,j,k) - frhoes(i,j,k))
               ekt = (1.d0 - theta) * eta(i,j,k)
               rhs(i,j,k) = (rhs(i,j,k) + r(i) * s(j) *
     @                       (bs + dtm * (ero(i,j,k) + es) +
     @                        ek * c * edot(i,j,k) -
     @                        ekt * dfo(i,j,k))) /
     @                      (1.d0 - ekt)
            enddo
         enddo
      enddo
      end

      subroutine anatw2(test,
     @  dims(reg),
     @  temp, p, xf, Tc, dx, xlo, lo)
      implicit none
      integer dimdec(reg)
      real*8 test(dimv(reg), 0:1)
      real*8 temp(dimv(reg))
      real*8 p, xf, Tc, dx(3), xlo(3)
      integer lo(3)
      integer i, j, k
      real*8 x, y, z, r2
      do k = regl2, regh2
         z = xlo(3) + dx(3) * ((k-lo(3)) + 0.5d0)
         do j = regl1, regh1
            y = xlo(2) + dx(2) * ((j-lo(2)) + 0.5d0)
            do i = regl0, regh0
               x  = xlo(1) + dx(1) * ((i-lo(1)) + 0.5d0)
               r2 = x*x + y*y + z*z
               test(i,j,k,0) = Tc * max((1.d0-r2/xf**2), 0.d0)**(1.d0/p)
               test(i,j,k,1) = temp(i,j,k) - test(i,j,k,0)
            enddo
         enddo
      enddo
      end

      subroutine cfrhoe( dims(reg),
     @  frhoe,
     @  dims(fb),
     @  state,
     @  dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(fb)
      integer dimdec(sb)
      real*8 frhoe(dimv(fb))
      real*8 state(dimv(sb), 0:SVSIZE-1)
c      real*8 kin
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
c               kin = 0.5d0 * (state(i,j,k,XMOM)   ** 2 +
c     @                        state(i,j,k,XMOM+1) ** 2 +
c     @                        state(i,j,k,XMOM+2) ** 2) /
c     @                       state(i,j,k,DEN)
c               frhoe(i,j,k) = state(i,j,k,EDEN) - kin
               frhoe(i,j,k) = state(i,j,k,EINT)
            enddo
         enddo
      enddo
      end

c temp contains frhoe on input:

      subroutine gtemp(dims(reg),
     @  temp, dims(tb),
     @  const, em, en, 
     @  state, dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(tb)
      integer dimdec(sb)
      real*8 temp(dimv(tb))
      real*8 const(0:1), em(0:1), en(0:1)
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 alpha, teff, ex, frhoal
      integer i, j, k
      if (en(0) .ge. 1.d0) then
         print *, "Bad exponent for cv calculation"
         stop
      endif
      ex = 1.d0 / (1.d0 - en(0))
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               if (em(0) .eq. 0.d0) then
                  alpha = const(0)
               else
                  alpha = const(0) * state(i,j,k,DEN) ** em(0)
               endif
               frhoal = state(i,j,k,DEN) * alpha + tiny
               if (en(0) .eq. 0.d0) then
                  temp(i,j,k) = temp(i,j,k) / frhoal
               else
                  teff = max(temp(i,j,k), tiny)
                  temp(i,j,k) = ((1.d0 - en(0)) * teff / frhoal) ** ex
               endif
            enddo
         enddo
      enddo
      end

c temp contains temp on input:

      subroutine gcv(dims(reg),
     @  cv,dims(cbox),
     @  temp, dims(tbox),
     @  const, em, en, tf,
     @  state, dims(sbox))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(cbox)
      integer dimdec(tbox)
      integer dimdec(sbox)
      real*8 cv(dimv(cbox))
      real*8 temp(dimv(tbox))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 state(dimv(sbox), 0:SVSIZE-1)
      real*8 alpha, teff, frhoal
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               if (em(0) .eq. 0.d0) then
                  alpha = const(0)
               else
                  alpha = const(0) * state(i,j,k,DEN) ** em(0)
               endif
               frhoal = state(i,j,k,DEN) * alpha + tiny
               if (en(0) .eq. 0.d0) then
                  cv(i,j,k) = alpha
               else
                  teff = max(temp(i,j,k), tiny)
                  teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
                  cv(i,j,k) = alpha * teff ** (-en(0))
               endif
            enddo
         enddo
      enddo
      end

c exch contains temp on input:

      subroutine cexch( dims(reg),
     @  exch, dims(xbox),
     @  er  , dims(ebox),
     @  fkp , dims(kbox),
     @  sigma, c)
      implicit none
      integer dimdec(reg)
      integer dimdec(xbox)
      integer dimdec(ebox)
      integer dimdec(kbox)
      real*8 exch(dimv(xbox))
      real*8 er  (dimv(ebox))
      real*8 fkp (dimv(kbox))
      real*8 sigma, c
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               exch(i,j,k) = fkp(i,j,k) *
     @                  (4.d0 * sigma * exch(i,j,k)**4
     @                    - c * er(i,j,k))
            enddo
         enddo
      enddo
      end

      subroutine extrp1(frhoes,
     @  dims(reg),
     @  frhoem, exch, state,
     @  dims(sb), dt)
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(sb)
      real*8 frhoes(dimv(reg))
      real*8 frhoem(dimv(reg))
      real*8 exch(dimv(reg))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 dt
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               frhoes(i,j,k) = frhoem(i,j,k) - dt * exch(i,j,k)
            enddo
         enddo
      enddo
      end

c exch contains temp on input:

      subroutine extrp3(exch,
     @  dims(reg),
     @  er, fkp, cv, state,
     @  dims(sb),
     @  sigma, c, dt)
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(sb)
      real*8 exch(dimv(reg))
      real*8 er(dimv(reg))
      real*8 fkp(dimv(reg))
      real*8 cv(dimv(reg))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 sigma, c, dt
      real*8 a, b, d, t, t3
      integer i, j, k, it
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               b = state(i,j,k,DEN)
               t = max(exch(i,j,k), tiny)
               a = fkp(i,j,k) * 4.d0 * sigma
               b = b * cv(i,j,k) / dt
               d = b * t + fkp(i,j,k) * c * er(i,j,k)
               do it = 1, 10
                  t3 = t**3
                  t = t - (t * (a * t3 + b) - d) / (4.d0 * a * t3 + b)
               enddo
               exch(i,j,k) = fkp(i,j,k) *
     @            (4.d0 * sigma * t**4 - c * er(i,j,k))
            enddo
         enddo
      enddo
      end

c exch contains temp on input:

      subroutine extrp4(exch,
     @  dims(reg),
     @  er, fkp, kpen, cv, state,
     @  dims(sb),
     @  sigma, c, dt)
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(sb)
      real*8 exch(dimv(reg))
      real*8 er(dimv(reg))
      real*8 fkp(dimv(reg))
      real*8 kpen(0:1)
      real*8 cv(dimv(reg))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 sigma, c, dt
      real*8 a, b, d, e, en, t, t3, tn
      integer i, j, k, it
      en = kpen(0)
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               b = state(i,j,k,DEN)
               t = max(exch(i,j,k), tiny)
               e = fkp(i,j,k) * t**en
               a = e * 4.d0 * sigma
               e = e * c * er(i,j,k)
               b = b * cv(i,j,k) / dt
               d = b * t
               do it = 1, 20
                  tn = t**(-en)
                  t3 = t**3 * tn
                  tn = t - (t * (a * t3 + b) - e * tn - d) /
     @               ((4.d0 - en) * a * t3 + b + en * e * tn / t)
                  if (tn .le. tiny) then
                     t = 0.5d0 * t
                  else
                     t = tn
                  endif
               enddo
               exch(i,j,k) = (a * t**4 - e) * t**(-en)
            enddo
         enddo
      enddo
      end

      subroutine ceta2(dims(reg),
     @  eta, etainv, dims(etab),
     @  frho, dims(sb),
     @  temp, dims(tb),
     @  cv, dims(cb),
     @  fkp, dims(fb),
     @  er, dims(ebox),
     @  dtemp, dtime, sigma, c, underr, lagpla)
      implicit none
      integer dimdec(reg)
      integer dimdec(etab)
      integer dimdec(sb)
      integer dimdec(tb)
      integer dimdec(cb)
      integer dimdec(fb)
      integer dimdec(ebox)
      real*8 eta(dimv(etab))
      real*8 etainv(dimv(etab))
      real*8 frho(dimv(sb))
      real*8 temp(dimv(tb))
      real*8 cv(dimv(cb))
      real*8 fkp(dimv(fb))
      real*8 er(dimv(ebox))
      real*8 dtemp, dtime, sigma, c, underr
      integer lagpla
      real*8 d, frc, fac0, fac1, fac2
      integer i, j, k
      fac1 = 16.d0 * sigma * dtime
      if (lagpla .eq. 0) then
         fac0 = 0.25d0 * fac1 / dtemp
         fac2 = dtime * c / dtemp
      endif
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               if (lagpla .ne. 0) then
c assume eta and fkp are the same
                  d = fac1 * fkp(i,j,k) * temp(i,j,k) ** 3
               else
                  d = fac0 * (eta(i,j,k) * (temp(i,j,k) + dtemp) ** 4 -
     @                        fkp(i,j,k) * (temp(i,j,k)        ) ** 4) -
     @                fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
c alternate form, sometimes worse, sometimes better:
c                  d = fac1 * fkp(i,j,k) * temp(i,j,k) ** 3 +
c     @                fac0 * (eta(i,j,k) - fkp(i,j,k)) * temp(i,j,k) ** 4 -
c     @                fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
c another alternate form (much worse):
c                  d = fac1 * fkp(i,j,k) * (temp(i,j,k) + dtemp) ** 3 +
c     @                fac0 * (eta(i,j,k) - fkp(i,j,k))
c     @                     * (temp(i,j,k) + dtemp) ** 4 -
c     @                fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
               endif
               frc = frho(i,j,k) * cv(i,j,k) + tiny
               eta(i,j,k) = d / (d + frc)
               etainv(i,j,k) = underr * frc / (d + frc)
               eta(i,j,k) = 1.d0 - etainv(i,j,k)
c               eta(i,j,k) = 1.d0 - underr * (1.d0 - eta(i,j,k))
            enddo
         enddo
      enddo
      end

      subroutine ssup(ptemp1, ptemp0,
     @  dims(reg),
     @  ecorr, eta, etainv, fkp, c, dt)
      implicit none
      integer dimdec(reg)
      real*8 ptemp1(dimv(reg))
      real*8 ptemp0(dimv(reg))
      real*8 ecorr(dimv(reg))
      real*8 eta(dimv(reg))
      real*8 etainv(dimv(reg))
      real*8 fkp(dimv(reg))
      real*8 c, dt
      real*8 fac, up
      integer i, j, k
      fac = c * dt
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               up = fac * ecorr(i,j,k) * fkp(i,j,k) * etainv(i,j,k)
               ptemp1(i,j,k) = ptemp0(i,j,k) + up
c            state(i,j,k,EDEN) = state(i,j,k,EDEN) + up
            enddo
         enddo
      enddo
      end

      subroutine ceup( dims(reg), relres, absres, 
     @  frhoes, dims(grd),
     @  frhoem, eta, etainv, dfo, dfn, exch,
     @  dt, theta)
      implicit none
      integer dimdec(reg)
      integer dimdec(grd)
      real*8 frhoes(dimv(grd))
      real*8 frhoem(dimv(grd))
      real*8 eta(dimv(grd))
      real*8 etainv(dimv(grd))
      real*8 dfo(dimv(grd))
      real*8 dfn(dimv(grd))
      real*8 exch(dimv(grd))
      real*8 dt, theta, relres, absres
      real*8 tmp, chg, tot
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               chg = 0.d0
               tot = 0.d0
               tmp = eta(i,j,k) * frhoes(i,j,k) +
     @           etainv(i,j,k) *
     @           (frhoem(i,j,k) -
     @            dt * ((1.d0 - theta) *
     @                  (dfo(i,j,k) - dfn(i,j,k)) +
     @                  exch(i,j,k)))
               chg = abs(tmp - frhoes(i,j,k))
               tot = abs(frhoes(i,j,k))
               frhoes(i,j,k) = tmp
               absres = max(absres, chg)
               relres = max(relres, chg / (tot + tiny))
            enddo
         enddo
      enddo
      end

      subroutine ceupdterm( dims(reg), relres, absres, 
     @  frhoes, dims(grd),
     @  frhoem, eta, etainv, dfo, dfn, exch, dterm,
     @  dt, theta)
      implicit none
      integer dimdec(reg)
      integer dimdec(grd)
      real*8 frhoes(dimv(grd))
      real*8 frhoem(dimv(grd))
      real*8 eta(dimv(grd))
      real*8 etainv(dimv(grd))
      real*8 dfo(dimv(grd))
      real*8 dfn(dimv(grd))
      real*8 exch(dimv(grd))
      real*8 dterm(dimv(grd))
      real*8 dt, theta, relres, absres
      real*8 tmp, chg, tot
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               chg = 0.d0
               tot = 0.d0
               tmp = eta(i,j,k) * frhoes(i,j,k) +
     @           etainv(i,j,k) *
     @           (frhoem(i,j,k) -
     @            dt * ((1.d0 - theta) *
     @                  (dfo(i,j,k) - dfn(i,j,k)) +
     @                  exch(i,j,k)))
     @           + dt * dterm(i,j,k)
               chg = abs(tmp - frhoes(i,j,k))
               tot = abs(frhoes(i,j,k))
               frhoes(i,j,k) = tmp
               absres = max(absres, chg)
               relres = max(relres, chg / (tot + tiny))
            enddo
         enddo
      enddo
      end

c convergence check for multigroup
      subroutine cconv(relres, absres, frhoes, frhoess,
     @   dims(fbox),
     @   dims(reg))
      implicit none
      integer dimdec(fbox)
      integer dimdec(reg)
      real*8 frhoes(dimv(fbox))
      real*8 frhoess(dimv(fbox))
      real*8 relres, absres
      real*8 chg, tot
      integer i, j, k

      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0

               chg = abs(frhoes(i,j,k) - frhoess(i,j,k))
               tot = abs(frhoes(i,j,k))

               absres = max(absres, chg)
               relres = max(relres, chg / (tot + tiny))

            enddo
         enddo
      enddo
      end

c nonconservative form based on delta B
      subroutine nceup(dims(reg), relres, absres, 
     @  frhoes, dims(grd),
     @  frhoem, eta, etainv,
     @  er, dims(ebox),
     @  dfo, dfn, temp, fkp, cv,
     @  state, dims(sb),
     @  sigma, c, dt, theta)
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(grd)
      integer dimdec(sb)
      integer dimdec(ebox)
      real*8 frhoes(dimv(grd))
      real*8 frhoem(dimv(grd))
      real*8 eta(dimv(grd))
      real*8 etainv(dimv(grd))
      real*8 er(dimv(ebox))
      real*8 dfo(dimv(grd))
      real*8 dfn(dimv(grd))
      real*8 temp(dimv(grd))
      real*8 fkp(dimv(grd))
      real*8 cv(dimv(reg))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 sigma, c, dt, theta, relres, absres
      real*8 tmp, chg, tot, exch, b, db, dbdt, frhocv
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               chg = 0.d0
               tot = 0.d0
               frhocv = state(i,j,k,DEN) * cv(i,j,k)
               dbdt = 16.d0 * sigma * temp(i,j,k)**3
               b = 4.d0 * sigma * temp(i,j,k)**4
               exch = fkp(i,j,k) * (b - c * er(i,j,k))
               tmp = eta(i,j,k) * frhoes(i,j,k) + etainv(i,j,k) *
     @               (frhoem(i,j,k) -
     @                dt * ((1.d0 - theta) *
     @                      (dfo(i,j,k) - dfn(i,j,k)) +
     @                      exch))
#if 1
               if (frhocv .gt. tiny .and. tmp .gt. frhoes(i,j,k)) then
                  db = (tmp - frhoes(i,j,k)) * dbdt / frhocv
                  if (b + db .le. 0.d0) then
                     print *, i, j, k, b, db, b+db
                  endif
                  tmp = ((b + db) / (4.d0 * sigma))**0.25d0
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
      end

      subroutine cetot(dims(reg),
     @  state, dims(sb),
     @  frhoe, dims(fb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(sb)
      integer dimdec(fb)
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 frhoe(dimv(fb))
      real*8 kin
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
c               kin = 0.5d0 * (state(i,j,k,XMOM)   ** 2 +
c     @                        state(i,j,k,XMOM+1) ** 2 +
c     @                        state(i,j,k,XMOM+2) ** 2) /
c     @                       state(i,j,k,DEN)
               kin = state(i,j,k,EDEN) - state(i,j,k,EINT)
               state(i,j,k,EINT) = frhoe(i,j,k)
               state(i,j,k,EDEN) = frhoe(i,j,k) + kin
            enddo
         enddo
      enddo
      end

      subroutine fkpn( dims(reg),
     @  fkp, dims(fb),
     @  const, em, en,
     @  ep, nu, tf, 
     @  temp, dims(tb),
     @  state, dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(fb)
      integer dimdec(tb)
      integer dimdec(sb)
      real*8 fkp(dimv(fb))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 ep(0:1), nu
      real*8 temp(dimv(tb))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 teff
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               teff = max(temp(i,j,k), tiny)
               teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
               fkp(i,j,k) = const(0) *
     @                        (state(i,j,k,DEN) ** em(0)) *
     @                        (teff ** (-en(0))) *
     @                        (nu ** (ep(0)))
            enddo
         enddo
      enddo
      end

      subroutine rosse1( dims(reg),
     @  kappar, dims(kbox),
     @  const, em, en,
     @  ep, nu,
     @  tf, kfloor, 
     @  temp, dims(tb),
     @  state, dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(kbox)
      integer dimdec(tb)
      integer dimdec(sb)
      real*8 kappar(dimv(kbox))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 ep(0:1), nu
      real*8 temp(dimv(tb))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 kfloor
      real*8 kf, teff
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               teff = max(temp(i,j,k), tiny)
               teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
               kf = const(0) *
     @                (state(i,j,k,DEN) ** em(0)) *
     @                (teff ** (-en(0))) *
     @             (nu ** (ep(0)))
               kappar(i,j,k) = max(kf, kfloor)
            enddo
         enddo
      enddo
      end

      subroutine rosse1s( dims(reg),
     @  kappar, dims(kbox),
     @  const, em, en,
     @  ep, 
     @  sconst, sem, sen,
     @  sep, 
     @  nu,
     @  tf, kfloor, 
     @  temp, dims(tb),
     @  state, dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(kbox)
      integer dimdec(tb)
      integer dimdec(sb)
      real*8 kappar(dimv(kbox))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 ep(0:1), nu
      real*8 sconst(0:1), sem(0:1), sen(0:1), sep(0:1)
      real*8 temp(dimv(tb))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 kfloor
      real*8 kf, teff, sct
      integer i, j, k
      do k = regl2, regh2
         do j = regl1, regh1
            do i = regl0, regh0
               teff = max(temp(i,j,k), tiny)
               teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
               kf = const(0) *
     @                (state(i,j,k,DEN) ** em(0)) *
     @                (teff ** (-en(0))) *
     @                (nu ** (ep(0)))
               sct = sconst(0) *
     @                (state(i,j,k,DEN) ** sem(0)) *
     @                (teff ** (-sen(0))) *
     @                (nu ** (sep(0)))
               kappar(i,j,k) = max(kf+sct, kfloor)
            enddo
         enddo
      enddo
      end

      subroutine nfloor(dest,
     @  dims(dbox),
     @  dims(reg),
     @  nflr, flr, nvar)
      implicit none
      integer dimdec(dbox)
      integer dimdec(reg)
      integer nvar, nflr
      real*8 dest(dimv(dbox), 0:nvar-1)
      real*8 flr
      integer i, j, k, n
      nflr = 0
      do n = 0, nvar-1
         do k = regl2, regh2
            do j = regl1, regh1
               do i = regl0, regh0
                  if (dest(i,j,k,n) .lt. flr) then
                     dest(i,j,k,n) = flr
                     nflr = nflr + 1
                  endif
               enddo
            enddo
         enddo
      enddo
      end

c *********************************
c ** BEGIN MGFLD routines        **
c *********************************

      subroutine lacoefmgfld(a,
     +     dims(abox),
     +     dims(reg),
     +     kappa,
     +     dims(kbox),
     +     r, s,
     +     dt, c)
      implicit none
      integer dimdec(abox)
      integer dimdec(reg)
      integer dimdec(kbox)

      real*8 a(dimv(abox))
      real*8 kappa(dimv(kbox))
      real*8 r(regl0:regh0)
      real*8 s(regl1:regh1)
      real*8 dt, c

      integer i, j, k

      do k = regl2, regh2
      do j = regl1, regh1
      do i = regl0, regh0
         a(i,j,k) = c*kappa(i,j,k) + 1.d0/dt
      enddo
      enddo
      enddo

      end

c *********************************
c ** END MGFLD routines          **
c *********************************

      subroutine rfface(fine,
     @     dims(fbox),
     @     crse,
     @     dims(cbox),
     @     idim, irat)
      implicit none
      integer dimdec(fbox)
      integer dimdec(cbox)
      real*8 fine(dimv(fbox))
      real*8 crse(dimv(cbox))
      integer idim, irat(0:2)
      integer i, j, k
      integer rfac
      if (idim .eq. 0) then
         rfac = irat(1) * irat(2)
         do k = fboxl2, fboxh2
         do j = fboxl1, fboxh1
            fine(fboxl0,j,k) = crse(cboxl0, j/irat(1), k/irat(2)) / rfac
         enddo
         enddo
      else if (idim .eq. 1) then
         rfac = irat(0) * irat(2)
         do k = fboxl2, fboxh2
         do i = fboxl0, fboxh0
            fine(i,fboxl1,k) = crse(i/irat(0), cboxl1, k/irat(2)) / rfac
         enddo
         enddo
      else
         rfac = irat(0) * irat(1)
         do j = fboxl1, fboxh1
         do i = fboxl0, fboxh0
            fine(i,j,fboxl2) = crse(i/irat(0), j/irat(1), cboxl2) / rfac
         enddo
         enddo
      endif
      end

