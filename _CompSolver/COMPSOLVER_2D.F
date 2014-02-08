
#include <BC_TYPES.H>
#include <COMPSOLVER_F.H>

#define dims(a) a loi1, a loi2, a hii1, a hii2
#define dimdec(a) dims(a)
#define dimv(a) a loi1:a hii1,a loi2:a hii2

       subroutine FORT_EDGE_AVG_DOWN(
     &   dir,
     &   fine, dims(fine),
     &   nref,
     &   coarse, dims(coarse))
c
       integer dir, nref(2)
       integer dimdec(fine)
       integer dimdec(coarse)
       real*8 fine(dimv(fine))
       real*8 coarse(dimv(coarse))
c
       integer i,j,k
       integer imin,imax,jmin,jmax,kmin,kmax
       integer nrefx, nrefy
       real*8 sum
c
       nrefx = nref(1)
       nrefy = nref(2)

       imin = coarseloi1
       imax = coarsehii1
       jmin = coarseloi2
       jmax = coarsehii2

       if( dir .eq. 0 ) then
c...
c... X-Direction
c...
c
          do j=jmin,jmax
             do i=imin,imax
                kmin = j*nrefy
                kmax = (j+1)*nrefy-1
                sum = 0.0
                do k=kmin,kmax
                   sum = sum + fine(nrefx*i,k)
                end do
                coarse(i,j) = sum / dble(nrefy)
             end do
          end do

       else if( dir .eq. 1 ) then
c...
c... Y-Direction
c...
c
          do j=jmin,jmax
             do i=imin,imax
                kmin = i*nrefx
                kmax = (i+1)*nrefx-1
                sum = 0.0
                do k=kmin,kmax
                   sum = sum + fine(k,nrefy*j)
                end do
                coarse(i,j) = sum / dble(nrefx)
             end do
          end do
       endif
c
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
      
      integer  i, j, n, ic, jc, ioff, joff
      integer  lratx, lraty

      lratx = lrat(1)
      lraty = lrat(2)

      fac = 1.0d0 / (lratx * lraty)

      do 100 n = 1, nvar

c     ::::: set coarse grid to zero on overlap
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,n) = 0.0d0
            enddo
         enddo

c     ::::: sum fine data
         do joff = 0, lraty-1
            do jc = lo(2), hi(2)
               j = jc*lraty + joff
               do ioff = 0, lratx-1
                  do ic = lo(1), hi(1)
                     i = ic*lratx + ioff
                     crse(ic,jc,n) = crse(ic,jc,n) + fine(i,j,n)
                  enddo
               enddo
            enddo
         enddo
         
c     ::::: divide out by volume weight
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,n) = crse(ic,jc,n) * fac
            enddo
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

      integer dir,i,j
      real*8 beta, dx(BL_SPACEDIM), fac


      if( dir .eq. 0 ) then

c...     x-direction

         fac = - beta / dx(1)

         do j = fluxloi2,fluxhii2
            do i = fluxloi1,fluxhii1
               flux(i,j) = dcoef(i,j) * (density(i,j) - density(i-1,j)) * fac
            end do
         end do

      else if( dir .eq. 1 ) then

c...     y-direction

         fac = - beta / dx(2)

         do j = fluxloi2,fluxhii2
            do i = fluxloi1,fluxhii1
               flux(i,j) = dcoef(i,j) * (density(i,j) - density(i,j-1)) * fac
            end do
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
       integer dir, nref(2)
       integer dimdec(fine)
       integer dimdec(coarse)
       real*8 fine(dimv(fine))
       real*8 coarse(dimv(coarse))
c
       integer i,j,k
       integer imin,imax,jmin,jmax,kmin,kmax
       integer nrefx, nrefy
       real*8 center, left, right, top, bottom, denom
c
       nrefx = nref(1)
       nrefy = nref(2)

       imin = coarseloi1
       imax = coarsehii1
       jmin = coarseloi2
       jmax = coarsehii2

       if( dir .eq. 0 ) then
c...
c... X-Direction
c...
c
          do j=jmin,jmax
             do i=imin,imax
                kmin = j*nrefy
                kmax = (j+1)*nrefy-1

                left = 0.0
                do k=kmin,kmax
                   left = left + fine(nrefx*i-1,k)
                end do
                left = left / dble(nrefy)

                center = 0.0d0
                do k=kmin,kmax
                   center = center + fine(nrefx*i,k)
                end do
                center = center / dble(nrefy)

                right = 0.0d0
                do k=kmin,kmax
                   right = right + fine(nrefx*i+1,k)
                end do
                right = right / dble(nrefy)

                denom = center*(left+right) + 2.0d0*left*right
                if( denom .ne. 0.0d0 ) then
                   coarse(i,j) = 4.0d0 * left * center * right / denom
                else
                   coarse(i,j) = 0.0d0
                endif

             end do
          end do

       else if( dir .eq. 1 ) then
c...
c... Y-Direction
c...
c
          do j=jmin,jmax
             do i=imin,imax
                kmin = i*nrefx
                kmax = (i+1)*nrefx-1

                top = 0.0d0
                do k=kmin,kmax
                   top = top + fine(k,nrefy*j+1)
                end do
                top = top / dble(nrefx)

                center = 0.0d0
                do k=kmin,kmax
                   center = center + fine(k,nrefy*j)
                end do
                center = center / dble(nrefx)

                bottom = 0.0d0
                do k=kmin,kmax
                   bottom = bottom + fine(k,nrefy*j-1)
                end do
                bottom = bottom / dble(nrefx)

                denom = center*(top+bottom) + 2.0d0*top*bottom
                if( denom .ne. 0.0d0 ) then
                   coarse(i,j) = 4.0d0 * top * center * bottom / denom
                else
                   coarse(i,j) = 0.0d0
                endif

             end do
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
      integer ic, jc, if, jf, ir, jr
      real*8 xoff, yoff, sx, sy

      ir = nref(1)
      jr = nref(2)
c
      if (lo(1) .eq. dlo(1)) then
         do jc = lo(2), hi(2)
            c(lo(1)-1,jc) = 2.d0 * c(lo(1),jc) - c(lo(1)+1,jc)
         enddo
      endif
      if (hi(1) .eq. dhi(1)) then
         do jc = lo(2), hi(2)
            c(hi(1)+1,jc) = 2.d0 * c(hi(1),jc) - c(hi(1)-1,jc)
         enddo
      endif
      if (lo(2) .eq. dlo(2)) then
         do ic = lo(1), hi(1)
            c(ic,lo(2)-1) = 2.d0 * c(ic,lo(2)) - c(ic,lo(2)+1)
         enddo
      endif
      if (hi(2) .eq. dhi(2)) then
         do ic = lo(1), hi(1)
            c(ic,hi(2)+1) = 2.d0 * c(ic,hi(2)) - c(ic,hi(2)-1)
         enddo
      endif

      do jf = lo(2) * jr, (hi(2) + 1) * jr - 1
         jc = jf / jr
         yoff = (mod(jf,jr) + 0.5d0) / jr - 0.5d0
         do if = lo(1) * ir, (hi(1) + 1) * ir - 1
            ic = if / ir
            xoff = (mod(if,ir) + 0.5d0) / ir - 0.5d0
            sy = 0.5d0 * (c(ic,jc+1) - c(ic,jc-1))
            sx = 0.5d0 * (c(ic+1,jc) - c(ic-1,jc))
            f(if,jf) = c(ic,jc) + xoff * sx + yoff * sy
         enddo
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
      integer ic, jc, if, jf, nrefx, nrefy

      nrefx = nref(1)
      nrefy = nref(2)
c
      do jc = lo(2), hi(2)
         do ic = lo(1), hi(1)
            do jf = nrefy*jc, nrefy*(jc+1)-1
               do if = nrefx*ic, nrefx*(ic+1)-1
                  f(if,jf) = c(ic,jc) + f(if,jf)
               enddo
            enddo
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


c
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
      integer ic, jc, if, jf, nrefx, nrefy
      real*8 fac

      nrefx = nref(1)
      nrefy = nref(2)
      fac = 1.0d0 / (nrefx * nrefy)
c
      do jc = lo(2), hi(2)
         do ic = lo(1), hi(1)
            c(ic,jc) = 0.0d0
            do jf = nrefy*jc, nrefy*(jc+1)-1
               do if = nrefx*ic, nrefx*(ic+1)-1
                  c(ic,jc) = c(ic,jc) + f(if,jf)
               enddo
            enddo
            c(ic,jc) = c(ic,jc) * fac
         enddo
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
c
      end
