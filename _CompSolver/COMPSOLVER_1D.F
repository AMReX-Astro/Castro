
#include <BC_TYPES.H>
#include <COMPSOLVER_F.H>

#define dims(a) a loi1, a hii1
#define dimdec(a) dims(a)
#define dimv(a) a loi1:a hii1

       subroutine FORT_EDGE_AVG_DOWN(
     &   dir,
     &   fine, dims(fine),
     &   nref,
     &   coarse, dims(coarse))
c
       integer dir, nref(1)
       integer dimdec(fine)
       integer dimdec(coarse)
       real*8 fine(dimv(fine))
       real*8 coarse(dimv(coarse))
c
       integer i
       integer imin,imax
       integer nrefx
c
       nrefx = nref(1)

       imin = coarseloi1
       imax = coarsehii1

       if( dir .eq. 0 ) then
c...
c... X-Direction
c...
c
          do i=imin,imax
             coarse(i) = fine(nrefx*i)
          end do
c
       endif
       return
       end

      subroutine FORT_AVG_DOWN (crse, dims(crse),nvar,
     &                          fine, dims(fine),lo,hi,lrat)
      integer  lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer  nvar, lrat(BL_SPACEDIM)
      integer dimdec(crse)
      integer dimdec(fine)
      real*8 crse(dimv(crse),nvar)
      real*8 fine(dimv(fine),nvar)
      real*8 fac
      
      integer  i, n, ic, ioff
      integer  lratx

      lratx = lrat(1)

      fac = 1.0d0 / lratx

      do 100 n = 1, nvar

c     ::::: set coarse grid to zero on overlap
         do ic = lo(1), hi(1)
            crse(ic,n) = 0.0d0
         enddo

c     ::::: sum fine data
         do ioff = 0, lratx-1
            do ic = lo(1), hi(1)
               i = ic*lratx + ioff
               crse(ic,n) = crse(ic,n) + fine(i,n)
            enddo
         enddo
         
c     ::::: divide out by volume weight
         do ic = lo(1), hi(1)
            crse(ic,n) = crse(ic,n) * fac
         enddo
	 
 100  continue

      return
      end

      subroutine FORT_SET_ABEC_FLUX(
     &     dir,
     &     density, dims(density),
     &     dcoef, dims(dcoef),
     &     beta,
     &     dx,
     &     flux, dims(flux))

      integer dimdec(density)
      integer dimdec(dcoef)
      integer dimdec(flux)

      real*8 density(dimv(density))
      real*8 dcoef(dimv(dcoef))
      real*8 flux(dimv(flux))

      integer dir,i
      real*8 beta, dx(BL_SPACEDIM), fac

      if( dir .eq. 0 ) then

c...     x-direction

         fac = - beta / dx(1)

         do i = fluxloi1,fluxhii1
            flux(i) = dcoef(i) * (density(i) - density(i-1)) * fac
         end do

      endif

      return
      end

      subroutine FORT_EDGE_HARM_AVG(
     &   dir,
     &   fine, dims(fine),
     &   nref,
     &   coarse, dims(coarse))
c
       integer dir, nref(1)
       integer dimdec(fine)
       integer dimdec(coarse)
       real*8 fine(dimv(fine))
       real*8 coarse(dimv(coarse))
c
       integer i
       integer imin,imax
       integer nrefx
       real*8 center, left, right, denom
c
       nrefx = nref(1)

       imin = coarseloi1
       imax = coarsehii1

       if( dir .eq. 0 ) then
c...
c... X-Direction
c...
c
          do i=imin,imax

             left   = fine(nrefx*i-1)
             center = fine(nrefx*i)
             right  = fine(nrefx*i+1)

             denom = center*(left+right) + 2.0d0*left*right
             if( denom .ne. 0.0d0 ) then
                coarse(i) = 4.0d0 * left * center * right / denom
             else
                coarse(i) = 0.0d0
             endif

          end do

       endif
c
       return
       end

      subroutine FORT_CONS_INTERP (
     $     f, dims(f),
     $	   c, dims(c),
     $     lo, hi, nref,
     $     dlo, dhi
     $     )
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer dlo(BL_SPACEDIM)
      integer dhi(BL_SPACEDIM)

      integer dimdec(f)
      integer dimdec(c)
      real*8 f(dimv(f))
      real*8 c(dimv(c))

      integer nref(BL_SPACEDIM)
c
      integer ic, if, ir
      real*8 xoff, sx

      ir = nref(1)
c
      if (lo(1) .eq. dlo(1)) then
         c(lo(1)-1) = 2.d0 * c(lo(1)) - c(lo(1)+1)
      endif
      if (hi(1) .eq. dhi(1)) then
         c(hi(1)+1) = 2.d0 * c(hi(1)) - c(hi(1)-1)
      endif

      do if = lo(1) * ir, (hi(1) + 1) * ir - 1
         ic = if / ir
         xoff = (mod(if,ir) + 0.5d0) / ir - 0.5d0
         sx = 0.5d0 * (c(ic+1) - c(ic-1))
         f(if) = c(ic) + xoff * sx
      enddo
c
      return
      end

      subroutine FORT_PC_INTERP (
     $	   f, dims(f),
     $	   c, dims(c),
     $     lo, hi, nref
     $     )

      integer dimdec(f)
      integer dimdec(c)
      real*8 f(dimv(f))
      real*8 c(dimv(c))

      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer nref(BL_SPACEDIM)
c
      integer ic, if, nrefx

      nrefx = nref(1)
c
      do ic = lo(1), hi(1)
         do if = nrefx*ic, nrefx*(ic+1)-1
            f(if) = c(ic) + f(if)
         enddo
      enddo
c
      return
      end

      subroutine FORT_PL_INTERP (
     $	   f, dims(f),
     $	   c, dims(c),
     $     lo, hi, ratio
     $     )
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer ratio(BL_SPACEDIM)

      integer dimdec(f)
      integer dimdec(c)
      real*8 f(dimv(f))
      real*8 c(dimv(c))

      real*8 e,n,w,s,co16,co4
c
      integer i, j

      print *, "FORT_PL_INTERP not implemented in 1D"
      stop
c
#if 0
      if( ratio(1) .eq. 2 ) then

         if( ratio(2) .eq. 2 ) then

            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  co16 = c(i,j) * 0.0625d0

                  e = 1.0d0
                  w = 1.0d0
                  n = 1.0d0
                  s = 1.0d0

                  if( i .eq. lo(1) ) w = 0.0d0
                  if( i .eq. hi(1) ) e = 0.0d0
                  if( j .eq. lo(2) ) s = 0.0d0
                  if( j .eq. hi(2) ) n = 0.0d0

                  f(2*i-1,2*j+2) = f(2*i-1,2*j+2) + co16*n*w
                  f(2*i  ,2*j+2) = f(2*i  ,2*j+2) + co16*n*(2.0d0+w)
                  f(2*i+1,2*j+2) = f(2*i+1,2*j+2) + co16*n*(2.0d0+e)
                  f(2*i+2,2*j+2) = f(2*i+2,2*j+2) + co16*n*e

                  f(2*i-1,2*j+1) = f(2*i-1,2*j+1) + co16*(2.0d0+n)*w
                  f(2*i  ,2*j+1) = f(2*i  ,2*j+1) + co16*(2.0d0+n)*(2.0d0+w)
                  f(2*i+1,2*j+1) = f(2*i+1,2*j+1) + co16*(2.0d0+n)*(2.0d0+e)
                  f(2*i+2,2*j+1) = f(2*i+2,2*j+1) + co16*(2.0d0+n)*e

                  f(2*i-1,2*j  ) = f(2*i-1,2*j  ) + co16*(2.0d0+s)*w
                  f(2*i  ,2*j  ) = f(2*i  ,2*j  ) + co16*(2.0d0+s)*(2.0d0+w)
                  f(2*i+1,2*j  ) = f(2*i+1,2*j  ) + co16*(2.0d0+s)*(2.0d0+e)
                  f(2*i+2,2*j  ) = f(2*i+2,2*j  ) + co16*(2.0d0+s)*e

                  f(2*i-1,2*j-1) = f(2*i-1,2*j-1) + co16*s*w
                  f(2*i  ,2*j-1) = f(2*i  ,2*j-1) + co16*s*(2.0d0+w)
                  f(2*i+1,2*j-1) = f(2*i+1,2*j-1) + co16*s*(2.0d0+e)
                  f(2*i+2,2*j-1) = f(2*i+2,2*j-1) + co16*s*e

               enddo
            enddo

         else if ( ratio(2) .eq. 1 ) then

            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  co4 = c(i,j) * 0.25d0

                  e = 1.0d0
                  w = 1.0d0
                  if( i .eq. lo(1) ) w = 0.0d0
                  if( i .eq. hi(1) ) e = 0.0d0

                  f(2*i-1,j) = f(2*i-1,j) + co4*w
                  f(2*i  ,j) = f(2*i  ,j) + co4*(2.0d0+w)
                  f(2*i+1,j) = f(2*i+1,j) + co4*(2.0d0+e)
                  f(2*i+2,j) = f(2*i+2,j) + co4*e

               enddo
            enddo

         endif

      else if ( ratio(1) .eq. 1 ) then

         if( ratio(2) .eq. 2 ) then

            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  co4 = c(i,j) * 0.25d0

                  n = 1.0d0
                  s = 1.0d0
                  if( j .eq. lo(2) ) s = 0.0d0
                  if( j .eq. hi(2) ) n = 0.0d0

                  f(i,2*j+2) = f(i,2*j+2) + co4*n
                  f(i,2*j+1) = f(i,2*j+1) + co4*(2.0d0+n)
                  f(i,2*j  ) = f(i,2*j  ) + co4*(2.0d0+s)
                  f(i,2*j-1) = f(i,2*j-1) + co4*s

               enddo
            enddo

         else if( ratio(2) .eq. 1 ) then

            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  f(i,j) = f(i,j) + c(i,j)
               enddo
            enddo

         endif

      endif
#endif
c
      end

      subroutine FORT_PC_RESTRICT (
     $	   f, dims(f),
     $	   c, dims(c),
     $     lo, hi, nref
     $     )

      integer dimdec(f)
      integer dimdec(c)
      real*8 f(dimv(f))
      real*8 c(dimv(c))

      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer nref(BL_SPACEDIM)
c
      integer ic, if, nrefx
      real*8 fac

      nrefx = nref(1)
      fac = 1.0d0 / nrefx
c
      do ic = lo(1), hi(1)
         c(ic) = 0.0d0
         do if = nrefx*ic, nrefx*(ic+1)-1
            c(ic) = c(ic) + f(if)
         enddo
         c(ic) = c(ic) * fac
      enddo
c
      return
      end

      subroutine FORT_PL_RESTRICT (
     $	   f, dims(f),
     $	   c, dims(c),
     $     lo, hi, ratio
     $     )

      integer dimdec(f)
      integer dimdec(c)
      real*8 f(dimv(f))
      real*8 c(dimv(c))

      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer ratio(BL_SPACEDIM)
      real*8 e,n,w,s
c
      integer i, j
c
      print *, "FORT_PL_RESTRICT not implemented in 1D"
      stop
#if 0

      if( ratio(1) .eq. 2 ) then

         if( ratio(2) .eq. 2 ) then

            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  e = 1.0d0
                  w = 1.0d0
                  n = 1.0d0
                  s = 1.0d0
                  if( i .eq. lo(1) ) w = 0.0d0
                  if( i .eq. hi(1) ) e = 0.0d0
                  if( j .eq. lo(2) ) s = 0.0d0
                  if( j .eq. hi(2) ) n = 0.0d0

       c(i,j) = (
     &  f(2*i-1,2*j+2)*n*w                 + f(2*i  ,2*j+2)*n*(2.0d0+w)
     &+ f(2*i+1,2*j+2)*n*(2.0d0+e)         + f(2*i+2,2*j+2)*n*e
     &+ f(2*i-1,2*j+1)*(2.0d0+n)*w         + f(2*i  ,2*j+1)*(2.0d0+n)*(2.0d0+w)
     &+ f(2*i+1,2*j+1)*(2.0d0+n)*(2.0d0+e) + f(2*i+2,2*j+1)*(2.0d0+n)*e
     &+ f(2*i-1,2*j  )*(2.0d0+s)*w         + f(2*i  ,2*j  )*(2.0d0+s)*(2.0d0+w)
     &+ f(2*i+1,2*j  )*(2.0d0+s)*(2.0d0+e) + f(2*i+2,2*j  )*(2.0d0+s)*e
     &+ f(2*i-1,2*j-1)*s*w                 + f(2*i  ,2*j-1)*s*(2.0d0+w)
     &+ f(2*i+1,2*j-1)*s*(2.0d0+e)         + f(2*i+2,2*j-1)*s*e 
     &          ) * 0.015625d0

               enddo
            enddo

         else if ( ratio(2) .eq. 1 ) then

            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  e = 1.0d0
                  w = 1.0d0
                  if( i .eq. lo(1) ) w = 0.0d0
                  if( i .eq. hi(1) ) e = 0.0d0

                  c(i,j) = (
     &                       f(2*i-1,j)*w         + f(2*i  ,j)*(2.0d0+w)
     &                     + f(2*i+1,j)*(2.0d0+e) + f(2*i+2,j)*e
     &                     ) * 0.125d0

               enddo
            enddo

         endif

      else if ( ratio(1) .eq. 1 ) then

         if( ratio(2) .eq. 2 ) then

            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  n = 1.0d0
                  s = 1.0d0
                  if( j .eq. lo(2) ) s = 0.0d0
                  if( j .eq. hi(2) ) n = 0.0d0

                  c(i,j) = (
     &                       f(i,2*j+2)*n         + f(i,2*j+1)*(2.0d0+n)
     &                     + f(i,2*j  )*(2.0d0+s) + f(i,2*j-1)*s
     &                     ) * 0.125d0

               enddo
            enddo

         else if( ratio(2) .eq. 1 ) then

            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j) = f(i,j)
               enddo
            enddo

         endif

      endif
#endif
c
      end
