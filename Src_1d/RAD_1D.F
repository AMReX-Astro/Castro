
#include "LO_BCTYPES.H"

#define dims(a) a l0, a h0
#define dimdec(a) dims(a)
#define dimv(a) a l0:a h0

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
      YMOM = -1
      ZMOM = -1
      EDEN = Eden_in
      EINT = Eint_in
      ITEMP = Temp_in
      IFS = FirstSpec
      IFX = FirstAux
      SVSIZE = NUM_STATE
      end

      subroutine multrs(d,
     @  dims(dbox),
     @  dims(reg),
     @  r, s)
      implicit none
      integer dimdec(dbox)
      integer dimdec(reg)
      real*8 d(dimv(dbox))
      real*8 r(dboxl0:dboxh0)
      real*8 s(1)
      integer i
      do i = regl0, regh0
         d(i) = d(i) * r(i)
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
      integer i
      do i = regl0, regh0
         d(i) = d(i) / r(i)
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
      real*8 s(1)
      integer i
      do i = regl0, regh0
         d(i) = d(i) / (r(i))
      enddo
      end

      subroutine sphc(r, s,
     @  dims(reg), dx)
      implicit none
      integer dimdec(reg)
      real*8 r(regl0:regh0)
      real*8 s(1)
      real*8 dx(1)
      real*8 h1, d1
      integer i
      h1 = 0.5d0 * dx(1)
      d1 = 1.d0 / (3.d0 * dx(1))
      do i = regl0, regh0
         r(i) = d1 * ((r(i) + h1)**3 - (r(i) - h1)**3)
      enddo
      end

      subroutine sphe(r, s, n,
     @  dims(reg), dx)
      implicit none
      integer dimdec(reg)
      real*8 r(regl0:regh0)
      real*8 s(1)
      integer n
      real*8 dx(1)
      integer i
      do i = regl0, regh0
         r(i) = r(i)**2
      enddo
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
      real*8 s(1)
      real*8 c
      integer i
      do i = regl0, regh0
         a(i) = r(i) * fkp(i) * c
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
      real*8 s(1)
      real*8 sigma
      integer i
      real*8 bs
      do i = regl0, regh0
         bs = 4.d0 * sigma * fkp(i) * temp(i)**4
         rhs(i) = r(i) * (bs - ert(i))
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
      real*8 s(1)
      real*8 c, dt
      integer i
      real*8 dtm, eki
      dtm = 1.d0 / dt
      do i = regl0, regh0
         eki = fkp(i) * etainv(i)
         a(i) = r(i) * (eki * c + dtm)
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
      real*8 fkp(dimv(reg))
      real*8 eta(dimv(reg))
      real*8 etainv(dimv(reg))
      real*8 r(regl0:regh0)
      real*8 s(1)
      real*8 c, dt, theta
      integer i
      real*8 dtm
      dtm = 1.d0 / dt
      do i = regl0, regh0
         a(i) = r(i) *
     @          (fkp(i) * etainv(i) * c + dtm) /
     @          (1.d0 - (1.d0 - theta) * eta(i))
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

      subroutine bclim0(b,
     @  dims(bbox),
     @  dims(reg),
     @  n, kappar, dims(kbox),
     @  q, r, s, c, dx)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer n
      real*8 b(dimv(bbox))
      real*8 kappar(dimv(kbox))
      real*8 q(regl0:regh0+1)
      real*8 r(regl0:regh0+1)
      real*8 s(1)
      real*8 c, dx(1)
      real*8 kavg
      integer i
      do i = regl0, regh0 + 1
         b(i) = r(i) * (c / 3.d0) / KAVG(kappar(i-1), kappar(i), dx(1))
      enddo
      end

      subroutine bclim1(b,
     @  dims(bbox),
     @  dims(reg),
     @  n, kappar, dims(kbox),
     @  er, q, r, s, c, dx, limiter)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer n, limiter
      real*8 b(dimv(bbox))
      real*8 kappar(dimv(kbox))
      real*8 er(dimv(kbox))
      real*8 q(regl0:regh0+1)
      real*8 r(regl0:regh0+1)
      real*8 s(1)
      real*8 c, dx(1)
      real*8 kavg
      integer i
      real*8 kap, par, lambda, FLDlambda
c gradient assembly:
      do i = regl0, regh0 + 1
         b(i) = abs(er(i) - er(i-1)) / dx(1)
      enddo
      i = regl0
      if (er(i-1) .eq. -1.d0) then
         b(i) = abs(er(i+1) - er(i)) / dx(1)
      endif
      i = regh0 + 1
      if (er(i) .eq. -1.d0) then
         b(i) = abs(er(i-1) - er(i-2)) / dx(1)
      endif
c construct coefficients:
      do i = regl0, regh0 + 1
         kap = KAVG(kappar(i-1), kappar(i), dx(1))
         par = b(i) /
     @          (kap * max(er(i-1), er(i), tiny))
         lambda = FLDlambda(par, limiter) 
         b(i) = r(i) * c * lambda / kap
      enddo
      end

      subroutine bclim2(b,
     @  dims(bbox),
     @  dims(reg),
     @  n, kappar, dims(kbox),
     @  er, q, r, s, c, dx, limiter)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer n, limiter
      real*8 b(dimv(bbox))
      real*8 kappar(dimv(kbox))
      real*8 er(dimv(kbox))
      real*8 q(regl0:regh0+1)
      real*8 r(regl0:regh0+1)
      real*8 s(1)
      real*8 c, dx(1)
      real*8 kavg
      integer i
      real*8 kap, par, lambda, FLDlambda
c gradient assembly:
      do i = regl0, regh0 + 1
         b(i) = abs(er(i) - er(i-1)) / dx(1)
      enddo
      i = regl0
      if (er(i-1) .eq. -1.d0) then
         b(i) = abs(er(i+1) - er(i)) / dx(1)
      endif
      i = regh0 + 1
      if (er(i) .eq. -1.d0) then
         b(i) = abs(er(i-1) - er(i-2)) / dx(1)
      endif
c construct coefficients:
      do i = regl0, regh0 + 1
         kap = KAVG(kappar(i-1), kappar(i), dx(1))
         par = b(i) /
c     @          (kap * max(er(i-1), er(i), tiny))
     @          (kap * max(abs(er(i-1)), abs(er(i)), tiny))
         lambda = FLDlambda(par, limiter)
         b(i) = r(i) * c * lambda / kap
      enddo
      end

      subroutine bclim3(b,
     @  dims(bbox),
     @  dims(reg),
     @  n, kappar, dims(kbox),
     @  er, q, r, s, c, dx, limiter)
      implicit none
      integer dimdec(bbox)
      integer dimdec(reg)
      integer dimdec(kbox)
      integer n, limiter
      real*8 b(dimv(bbox))
      real*8 kappar(dimv(kbox))
      real*8 er(dimv(kbox))
      real*8 q(regl0:regh0+1)
      real*8 r(regl0:regh0+1)
      real*8 s(1)
      real*8 c, dx(1)
      real*8 kavg
      integer i
      real*8 kap, par, lambda, FLDlambda
c gradient assembly:
      do i = regl0, regh0 + 1
         b(i) = abs(er(i) - er(i-1)) / dx(1)
      enddo
      i = regl0
      if (er(i-1) .eq. -1.d0) then
         b(i) = abs(er(i+1) - er(i)) / dx(1)
      endif
      i = regh0 + 1
      if (er(i) .eq. -1.d0) then
         b(i) = abs(er(i-1) - er(i-2)) / dx(1)
      endif
c construct coefficients:
      do i = regl0, regh0 + 1
         kap = KAVG(kappar(i-1), kappar(i), dx(1))
         par = b(i) /
     @          (kap * max(er(i-1), er(i), tiny))
         lambda = FLDlambda(par,limiter)
         b(i) = r(i) * c * lambda / kap
      enddo
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
      real*8 s(1)
      real*8 c, dx(1)
      real*8 kavg
      integer i
      real*8 kap
      do i = regl0, regh0 + 1
         kap = KAVG(kappar(i-1), kappar(i), dx(1))
         b(i) = r(i) * c * lambda(i) / kap
      enddo
      end

      subroutine flxlim(lambda,
     @  dims(rbox),
     @  dims(reg), limiter, n)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer n, limiter
      real*8 lambda(dimv(rbox)), FLDlambda
      integer i
      real*8 r
      do i = regl0, regh0 + 1
         r = lambda(i)
         lambda(i) = FLDlambda(r,limiter)
      enddo
      end

      subroutine eddfac(efact,
     @  dims(rbox),
     @  dims(reg), limiter, n)
      implicit none
      integer dimdec(rbox)
      integer dimdec(reg)
      integer n, limiter
      real*8 efact(dimv(rbox)), FLDlambda
      integer i
      real*8 r, lambda
      do i = regl0, regh0 + 1
         r = efact(i)
         lambda = FLDlambda(r,limiter)
         efact(i) = lambda + (lambda * r)**2
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
      real*8 dx(1)
      real*8 kavg
      integer i
      real*8 kap
c gradient assembly:
      do i = regl0, regh0 + 1
         r(i) = abs(er(i) - er(i-1)) / dx(1)
      enddo
      i = regl0
      if (er(i-1) .eq. -1.d0) then
         r(i) = abs(er(i+1) - er(i)) / dx(1)
      endif
      i = regh0 + 1
      if (er(i) .eq. -1.d0) then
         r(i) = abs(er(i-1) - er(i-2)) / dx(1)
      endif
c construct scaled gradient:
      do i = regl0, regh0 + 1
         kap = KAVG(kappar(i-1), kappar(i), dx(1))
         r(i) = r(i) /
     @          (kap * max(er(i-1), er(i), tiny))
      enddo
      end

      subroutine scgrd2(r,
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
      real*8 dx(1)
      real*8 kavg
      integer i
      real*8 kap
c gradient assembly:
      do i = regl0, regh0 + 1
         r(i) = abs(er(i) - er(i-1)) / dx(1)
      enddo
      i = regl0
      if (er(i-1) .eq. -1.d0) then
         r(i) = abs(er(i+1) - er(i)) / dx(1)
      endif
      i = regh0 + 1
      if (er(i) .eq. -1.d0) then
         r(i) = abs(er(i-1) - er(i-2)) / dx(1)
      endif
c construct scaled gradient:
      do i = regl0, regh0 + 1
         kap = KAVG(kappar(i-1), kappar(i), dx(1))
         r(i) = r(i) /
     @          (kap * max(er(i-1), er(i), tiny))
      enddo
      end

      subroutine scgrd3(r,
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
      real*8 dx(1)
      real*8 kavg
      integer i
      real*8 kap
c gradient assembly:
      do i = regl0, regh0 + 1
         r(i) = abs(er(i) - er(i-1)) / dx(1)
      enddo
      i = regl0
      if (er(i-1) .eq. -1.d0) then
         r(i) = abs(er(i+1) - er(i)) / dx(1)
      endif
      i = regh0 + 1
      if (er(i) .eq. -1.d0) then
         r(i) = abs(er(i-1) - er(i-2)) / dx(1)
      endif
c construct scaled gradient:
      do i = regl0, regh0 + 1
         kap = KAVG(kappar(i-1), kappar(i), dx(1))
         r(i) = r(i) /
     @          (kap * max(er(i-1), er(i), tiny))
      enddo
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
      real*8 s(1)
      real*8 dt, sigma, c, theta
      integer i
      real*8 dtm, ek, bs, es, ekt
      dtm = 1.d0 / dt
      do i = regl0, regh0
         ek = fkp(i) * eta(i)
         bs = etainv(i) *
     @        4.d0 * sigma * fkp(i) * temp(i)**4
         es = eta(i) * (frhoem(i) - frhoes(i))
         ekt = (1.d0 - theta) * eta(i)
         rhs(i) = (rhs(i) + r(i) *
     @               (bs + dtm * (ero(i) + es) +
     @                ek * c * edot(i) -
     @                ekt * dfo(i))) /
     @              (1.d0 - ekt)
      enddo
      end

      subroutine anatw2(test,
     @  dims(reg),
     @  temp, p, xf, Tc, dx, xlo, lo)
      implicit none
      integer dimdec(reg)
      real*8 test(dimv(reg), 0:1)
      real*8 temp(dimv(reg))
      real*8 p, xf, Tc, dx(1), xlo(1)
      integer lo(1)
      integer i
      real*8 x, r2
      do i = regl0, regh0
         x  = xlo(1) + dx(1) * ((i-lo(1)) + 0.5d0)
         r2 = x*x 
         test(i,0) = Tc * max((1.d0-r2/xf**2), 0.d0)**(1.d0/p)
         test(i,1) = temp(i) - test(i,0)
      enddo
      end

      subroutine cfrhoe(frhoe,
     @  dims(reg),
     @  state,
     @  dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(sb)
      real*8 frhoe(dimv(reg))
      real*8 state(dimv(sb), 0:SVSIZE-1)
c      real*8 kin
      integer i
      do i = regl0, regh0
c         kin = 0.5d0 * (state(i,XMOM) ** 2) /
c     @                 state(i,DEN)
c         frhoe(i) = state(i,EDEN) - kin
         frhoe(i) = state(i,EINT)
      enddo
      end

c temp contains frhoe on input:

      subroutine gtemp(temp,
     @  dims(reg),
     @  const, em, en, state,
     @  dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(sb)
      real*8 temp(dimv(reg))
      real*8 const(0:1), em(0:1), en(0:1)
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 alpha, teff, ex, frhoal
      integer i
      if (en(0) .ge. 1.d0) then
         print *, "Bad exponent for cv calculation"
         stop
      endif
      ex = 1.d0 / (1.d0 - en(0))
      do i = regl0, regh0
         if (em(0) .eq. 0.d0) then
            alpha = const(0)
         else
            alpha = const(0) * state(i,DEN) ** em(0)
         endif
         frhoal = state(i,DEN) * alpha + tiny
         if (en(0) .eq. 0.d0) then
            temp(i) = temp(i) / frhoal
         else
            teff = max(temp(i), tiny)
            temp(i) = ((1.d0 - en(0)) * teff / frhoal) ** ex
         endif
      enddo
      end

c temp contains temp on input:

      subroutine gcv(dims(reg),
     @  cv, dims(cbox),
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
      integer i
      do i = regl0, regh0
         if (em(0) .eq. 0.d0) then
            alpha = const(0)
         else
            alpha = const(0) * state(i,DEN) ** em(0)
         endif
         frhoal = state(i,DEN) * alpha + tiny
         if (en(0) .eq. 0.d0) then
            cv(i) = alpha
         else
            teff = max(temp(i), tiny)
            teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
            cv(i) = alpha * teff ** (-en(0))
         endif
      enddo
      end

c exch contains temp on input:

      subroutine cexch(exch,
     @  dims(reg),
     @  er, dims(ebox),
     @  fkp, sigma, c)
      implicit none
      integer dimdec(reg)
      integer dimdec(ebox)
      real*8 exch(dimv(reg))
      real*8 er(dimv(ebox))
      real*8 fkp(dimv(reg))
      real*8 sigma, c
      integer i
      do i = regl0, regh0
         exch(i) = fkp(i) *
     @             (4.d0 * sigma * exch(i)**4
     @               - c * er(i))
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
      integer i
      do i = regl0, regh0
         frhoes(i) = frhoem(i) - dt * exch(i)
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
      integer i, it
      do i = regl0, regh0
         b = state(i,DEN)
         t = max(exch(i), tiny)
         a = fkp(i) * 4.d0 * sigma
         b = b * cv(i) / dt
         d = b * t + fkp(i) * c * er(i)
         do it = 1, 10
            t3 = t**3
            t = t - (t * (a * t3 + b) - d) / (4.d0 * a * t3 + b)
         enddo
         exch(i) = fkp(i) *
     @      (4.d0 * sigma * t**4 - c * er(i))
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
      integer i, it
      en = kpen(0)
      do i = regl0, regh0
         b = state(i,DEN)
         t = max(exch(i), tiny)
         e = fkp(i) * t**en
         a = e * 4.d0 * sigma
         e = e * c * er(i)
         b = b * cv(i) / dt
         d = b * t
         do it = 1, 20
            tn = t**(-en)
            t3 = t**3 * tn
            tn = t - (t * (a * t3 + b) - e * tn - d) /
     @         ((4.d0 - en) * a * t3 + b + en * e * tn / t)
            if (tn .le. tiny) then
               t = 0.5d0 * t
            else
               t = tn
            endif
         enddo
         exch(i) = (a * t**4 - e) * t**(-en)
      enddo
      end

      subroutine ceta2(eta, etainv,
     @  dims(reg),
     @  frho,
     @  dims(sb),
     @  temp, cv, fkp,
     @  er, dims(ebox),
     @  dtemp, dtime, sigma, c, underr, lagpla)
      implicit none
      integer dimdec(reg)
      integer dimdec(sb)
      integer dimdec(ebox)
      real*8 eta(dimv(reg))
      real*8 etainv(dimv(reg))
      real*8 frho(dimv(sb))
      real*8 temp(dimv(reg))
      real*8 cv(dimv(reg))
      real*8 fkp(dimv(reg))
      real*8 er(dimv(ebox))
      real*8 dtemp, dtime, sigma, c, underr
      integer lagpla
      real*8 d, frc, fac0, fac1, fac2
      integer i
      fac1 = 16.d0 * sigma * dtime
      if (lagpla .eq. 0) then
         fac0 = 0.25d0 * fac1 / dtemp
         fac2 = dtime * c / dtemp
      endif
      do i = regl0, regh0
         if (lagpla .ne. 0) then
c assume eta and fkp are the same
            d = fac1 * fkp(i) * temp(i) ** 3
         else
            d = fac0 * (eta(i) * (temp(i) + dtemp) ** 4 -
     @                  fkp(i) * (temp(i)        ) ** 4) -
     @          fac2 * (eta(i) - fkp(i)) * er(i)
c alternate form, sometimes worse, sometimes better:
c            d = fac1 * fkp(i) * temp(i) ** 3 +
c     @          fac0 * (eta(i) - fkp(i)) * temp(i) ** 4 -
c     @          fac2 * (eta(i) - fkp(i)) * er(i)
c analytic derivatives for specific test problem:
c            d = (1.2d+6 * sigma * temp(i) ** 2 +
c     @           1.d+5 * c * er(i) * (temp(i) + tiny) ** (-2)) * dtime
c another alternate form (much worse):
c            d = fac1 * fkp(i) * (temp(i) + dtemp) ** 3 +
c     @          fac0 * (eta(i) - fkp(i)) * (temp(i) + dtemp) ** 4 -
c     @          fac2 * (eta(i) - fkp(i)) * er(i)
         endif
         frc = frho(i) * cv(i) + tiny
         eta(i) = d / (d + frc)
         etainv(i) = underr * frc / (d + frc)
         eta(i) = 1.d0 - etainv(i)
c         eta(i) = 1.d0 - underr * (1.d0 - eta(i))
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
      integer i
      fac = c * dt
      do i = regl0, regh0
         up = fac * ecorr(i) * fkp(i) * etainv(i)
         ptemp1(i) = ptemp0(i) + up
c            state(i,EDEN) = state(i,EDEN) + up
      enddo
      end

      subroutine ceup(relres, absres, frhoes,
     @  dims(reg),
     @  frhoem, eta, etainv, dfo, dfn, exch,
     @  dt, theta)
      implicit none
      integer dimdec(reg)
      real*8 frhoes(dimv(reg))
      real*8 frhoem(dimv(reg))
      real*8 eta(dimv(reg))
      real*8 etainv(dimv(reg))
      real*8 dfo(dimv(reg))
      real*8 dfn(dimv(reg))
      real*8 exch(dimv(reg))
      real*8 dt, theta, relres, absres
      real*8 tmp, chg, tot
      integer i
      do i = regl0, regh0
         chg = 0.d0
         tot = 0.d0
         tmp = eta(i) * frhoes(i) +
     @     etainv(i) *
     @     (frhoem(i) -
     @      dt * ((1.d0 - theta) *
     @            (dfo(i) - dfn(i)) +
     @            exch(i)))
         chg = abs(tmp - frhoes(i))
         tot = abs(frhoes(i))
         frhoes(i) = tmp
         absres = max(absres, chg)
         relres = max(relres, chg / (tot + tiny))
      enddo
      end

      subroutine ceupdterm(relres, absres, frhoes,
     @  dims(reg),
     @  frhoem, eta, etainv, dfo, dfn, exch, dterm,
     @  dt, theta)
      implicit none
      integer dimdec(reg)
      real*8 frhoes(dimv(reg))
      real*8 frhoem(dimv(reg))
      real*8 eta(dimv(reg))
      real*8 etainv(dimv(reg))
      real*8 dfo(dimv(reg))
      real*8 dfn(dimv(reg))
      real*8 exch(dimv(reg))
      real*8 dterm(dimv(reg))
      real*8 dt, theta, relres, absres
      real*8 tmp, chg, tot
      integer i
      do i = regl0, regh0
         chg = 0.d0
         tot = 0.d0
         tmp = eta(i) * frhoes(i) +
     @     etainv(i) *
     @     (frhoem(i) -
     @      dt * ((1.d0 - theta) *
     @            (dfo(i) - dfn(i)) +
     @            exch(i)))
     @     + dt * dterm(i)
         chg = abs(tmp - frhoes(i))
         tot = abs(frhoes(i))
         frhoes(i) = tmp
         absres = max(absres, chg)
         relres = max(relres, chg / (tot + tiny))
      enddo
      end

c#define PRINTMAX

c convergence check for multigroup
      subroutine cconv(relres, absres, frhoes, frhoess,
     +   dims(fbox),
     +   dims(reg))
      implicit none
      integer dimdec(fbox)
      integer dimdec(reg)
      real*8 frhoes(dimv(fbox))
      real*8 frhoess(dimv(fbox))
      real*8 relres, absres
      real*8 chg, tot
      integer i

#ifdef PRINTMAX
      integer isave
      isave = regl0 - 1
#endif

      do i = regl0, regh0

         chg = abs(frhoes(i) - frhoess(i))
         tot = abs(frhoes(i))

#ifdef PRINTMAX
         if (chg / (tot + tiny) .gt. relres) then
            isave = i
         endif
#endif

         absres = max(absres, chg)
         relres = max(relres, chg / (tot + tiny))

      enddo

#ifdef PRINTMAX
      if (isave .ge. regl0) then
         i = isave
         print *, "MAX IS ", i, abs(frhoes(i) - frhoess(i)), abs(frhoes(i))
      endif
#endif

      end

c nonconservative form based on delta B
      subroutine nceup(relres, absres, frhoes,
     @  dims(reg),
     @  frhoem, eta, etainv,
     @  er, dims(ebox),
     @  dfo, dfn, temp, fkp, cv,
     @  state, dims(sb),
     @  sigma, c, dt, theta)
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(sb)
      integer dimdec(ebox)
      real*8 frhoes(dimv(reg))
      real*8 frhoem(dimv(reg))
      real*8 eta(dimv(reg))
      real*8 etainv(dimv(reg))
      real*8 er(dimv(ebox))
      real*8 dfo(dimv(reg))
      real*8 dfn(dimv(reg))
      real*8 temp(dimv(reg))
      real*8 fkp(dimv(reg))
      real*8 cv(dimv(reg))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 sigma, c, dt, theta, relres, absres
      real*8 tmp, chg, tot, exch, b, db, dbdt, frhocv
      integer i
      do i = regl0, regh0
         chg = 0.d0
         tot = 0.d0
         frhocv = state(i,DEN) * cv(i)
         dbdt = 16.d0 * sigma * temp(i)**3
         b = 4.d0 * sigma * temp(i)**4
         exch = fkp(i) * (b - c * er(i))
         tmp = eta(i) * frhoes(i) + etainv(i) *
     @         (frhoem(i) -
     @          dt * ((1.d0 - theta) *
     @                (dfo(i) - dfn(i)) +
     @                exch))
#if 1
         if (frhocv .gt. tiny .and. tmp .gt. frhoes(i)) then
            db = (tmp - frhoes(i)) * dbdt / frhocv
            if (b + db .le. 0.d0) then
               print *, i, b, db, b+db
            endif
            tmp = ((b + db) / (4.d0 * sigma))**0.25d0
            tmp = frhoes(i) + frhocv * (tmp - temp(i))
         endif
#endif
         chg = abs(tmp - frhoes(i))
         tot = abs(frhoes(i))
         frhoes(i) = tmp
         absres = max(absres, chg)
         relres = max(relres, chg / (tot + tiny))
      enddo
      end

      subroutine cetot(state,
     @  dims(sb),
     @  dims(reg),
     @  frhoe)
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(sb)
      integer dimdec(reg)
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 frhoe(dimv(reg))
      real*8 kin
      integer i
      do i = regl0, regh0
c         kin = 0.5d0 * (state(i,XMOM) ** 2) /
c     @                 state(i,DEN)
         kin = state(i,EDEN) - state(i,EINT)
         state(i,EINT) = frhoe(i)
         state(i,EDEN) = frhoe(i) + kin
      enddo
      end

      subroutine fkpn(fkp,
     @  dims(reg),
     @  const, em, en,
     @  ep, nu,
     @  tf, temp, state,
     @  dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(reg)
      integer dimdec(sb)
      real*8 fkp(dimv(reg))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 ep(0:1), nu
      real*8 temp(dimv(reg))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 teff
      integer i
      do i = regl0, regh0
         teff = max(temp(i), tiny)
         teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
         fkp(i) = const(0) *
     @              (state(i,DEN) ** em(0)) *
     @              (teff ** (-en(0))) *
     @              (nu ** (ep(0)))
      enddo
      end

      subroutine rosse1(kappar,
     @  dims(kbox),
     @  dims(reg),
     @  const, em, en,
     @  ep, nu,
     @  tf, kfloor, temp, state,
     @  dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(kbox)
      integer dimdec(reg)
      integer dimdec(sb)
      real*8 kappar(dimv(kbox))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 ep(0:1), nu
      real*8 temp(dimv(reg))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 kfloor
      real*8 kf, teff
      integer i
      do i = regl0, regh0
         teff = max(temp(i), tiny)
         teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
         kf = const(0) *
     @          (state(i,DEN) ** em(0)) *
     @          (teff ** (-en(0))) *
     @          (nu ** (ep(0)))
         kappar(i) = max(kf, kfloor)
      enddo
      end

      subroutine rosse1s(kappar,
     @  dims(kbox),
     @  dims(reg),
     @  const, em, en,
     @  ep, 
     @  sconst, sem, sen,
     @  sep, 
     @  nu,
     @  tf, kfloor, temp, state,
     @  dims(sb))
      implicit none
      INCLUDE 'StateLayout.H'
      integer dimdec(kbox)
      integer dimdec(reg)
      integer dimdec(sb)
      real*8 kappar(dimv(kbox))
      real*8 const(0:1), em(0:1), en(0:1), tf(0:1)
      real*8 ep(0:1), nu
      real*8 sconst(0:1), sem(0:1), sen(0:1), sep(0:1)
      real*8 temp(dimv(reg))
      real*8 state(dimv(sb), 0:SVSIZE-1)
      real*8 kfloor
      real*8 kf, teff, sct
      integer i
      do i = regl0, regh0
         teff = max(temp(i), tiny)
         teff = teff + tf(0) * exp(-teff / (tf(0) + tiny))
         kf = const(0) *
     @          (state(i,DEN) ** em(0)) *
     @          (teff ** (-en(0))) *
     @          (nu ** (ep(0)))
         sct = sconst(0) *
     @          (state(i,DEN) ** sem(0)) *
     @          (teff ** (-sen(0))) *
     @          (nu ** (sep(0)))
         kappar(i) = max(kf+sct, kfloor)
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
      integer i, n
      nflr = 0
      do n = 0, nvar-1
         do i = regl0, regh0
            if (dest(i,n) .lt. flr) then
               dest(i,n) = flr
               nflr = nflr + 1
            endif
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
      real*8 s(1)
      real*8 dt, c
      integer i
      do i = regl0, regh0
         a(i) = c*kappa(i) + 1.d0/dt
         a(i) = r(i) * a(i)
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
      integer idim, irat(0:0)
      fine(fboxl0) = crse(cboxl0)
      end

