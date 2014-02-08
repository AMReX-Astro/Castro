
#include <BC_TYPES.H>
#include <COMPSOLVER_F.H>

#define dims(a) a loi1, a loi2, a loi3, a hii1, a hii2, a hii3
#define dimdec(a) dims(a)
#define dimv(a) a loi1:a hii1,a loi2:a hii2,a loi3:a hii3

       subroutine FORT_EDGE_AVG_DOWN(
     &   dir,
     &   fine, dims(fine),
     &   nref,
     &   coarse, dims(coarse))
c
       integer dir, nref(3)
       integer dimdec(fine)
       integer dimdec(coarse)
       real*8 fine(dimv(fine))
       real*8 coarse(dimv(coarse))
c
       integer i,j,k,m,n
       integer imin,imax,jmin,jmax,kmin,kmax,mmin,mmax,nmin,nmax
       integer nrefx, nrefy, nrefz
       real*8 sum
c
       nrefx = nref(1)
       nrefy = nref(2)
       nrefz = nref(3)

       imin = coarseloi1
       imax = coarsehii1
       jmin = coarseloi2
       jmax = coarsehii2
       kmin = coarseloi3
       kmax = coarsehii3

       if( dir .eq. 0 ) then
c...
c... X-Direction
c...
c
          do k=kmin,kmax
             do j=jmin,jmax
                do i=imin,imax
                   mmin = j*nrefy
                   mmax = (j+1)*nrefy-1
                   nmin = k*nrefz
                   nmax = (k+1)*nrefz-1
                   sum = 0.0
                   do m=mmin,mmax
                      do n=nmin,nmax
                         sum = sum + fine(nrefx*i,m,n)
                      end do
                   end do
                   coarse(i,j,k) = sum / dble(nrefy*nrefz)
                end do
             end do
          end do

       else if( dir .eq. 1 ) then
c...
c... Y-Direction
c...
c
          do k=kmin,kmax
             do j=jmin,jmax
                do i=imin,imax
                   mmin = i*nrefx
                   mmax = (i+1)*nrefx-1
                   nmin = k*nrefz
                   nmax = (k+1)*nrefz-1
                   sum = 0.0
                   do m=mmin,mmax
                      do n=nmin,nmax
                         sum = sum + fine(m,nrefy*j,n)
                      end do
                   end do
                   coarse(i,j,k) = sum / dble(nrefx*nrefz)
                end do
             end do
          end do

       else if( dir .eq. 2 ) then
c...
c... Z-Direction
c...
c
          do k=kmin,kmax
             do j=jmin,jmax
                do i=imin,imax
                   mmin = i*nrefx
                   mmax = (i+1)*nrefx-1
                   nmin = j*nrefy
                   nmax = (j+1)*nrefy-1
                   sum = 0.0
                   do m=mmin,mmax
                      do n=nmin,nmax
                         sum = sum + fine(m,n,nrefz*k)
                      end do
                   end do
                   coarse(i,j,k) = sum / dble(nrefx*nrefy)
                end do
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
      
      integer  i, j, k, n, ic, jc, kc, ioff, joff, koff
      integer  lratx, lraty, lratz

      lratx = lrat(1)
      lraty = lrat(2)
      lratz = lrat(3)

      fac = 1.0d0 / (lratx * lraty * lratz)

      do 100 n = 1, nvar

c     ::::: set coarse grid to zero on overlap
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = 0.0d0
               enddo
            enddo
         enddo

c     ::::: sum fine data
         do koff = 0, lratz-1
            do kc = lo(3), hi(3)
               k = kc*lratz + koff
               do joff = 0, lraty-1
                  do jc = lo(2), hi(2)
                     j = jc*lraty + joff
                     do ioff = 0, lratx-1
                        do ic = lo(1), hi(1)
                           i = ic*lratx + ioff
                           crse(ic,jc,kc,n) = crse(ic,jc,kc,n) +
     &                                        fine(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

c     ::::: divide out by volume weight
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = crse(ic,jc,kc,n) * fac
               enddo
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

      integer dir,i,j,k
      real*8 beta, dx(BL_SPACEDIM), fac


      if( dir .eq. 0 ) then

c...     x-direction

         fac = - beta / dx(1)

         do k = fluxloi3,fluxhii3
            do j = fluxloi2,fluxhii2
               do i = fluxloi1,fluxhii1
                  flux(i,j,k) = dcoef(i,j,k) *
     &               (density(i,j,k) - density(i-1,j,k)) * fac
               end do
            end do
         end do

      else if( dir .eq. 1 ) then

c...     y-direction

         fac = - beta / dx(2)

         do k = fluxloi3,fluxhii3
            do j = fluxloi2,fluxhii2
               do i = fluxloi1,fluxhii1
                  flux(i,j,k) = dcoef(i,j,k) *
     &               (density(i,j,k) - density(i,j-1,k)) * fac
               end do
            end do
         end do

      else if( dir .eq. 2 ) then

c...     z-direction

         fac = - beta / dx(3)

         do k = fluxloi3,fluxhii3
            do j = fluxloi2,fluxhii2
               do i = fluxloi1,fluxhii1
                  flux(i,j,k) = dcoef(i,j,k) *
     &               (density(i,j,k) - density(i,j,k-1)) * fac
               end do
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
       integer dir, nref(3)
       integer dimdec(fine)
       integer dimdec(coarse)
       real*8 fine(dimv(fine))
       real*8 coarse(dimv(coarse))
c
       integer i,j,k,m,n
       integer imin,imax,jmin,jmax,kmin,kmax,mmin,mmax,nmin,nmax
       integer nrefx, nrefy, nrefz
       real*8 center, left, right, denom
c
       nrefx = nref(1)
       nrefy = nref(2)
       nrefz = nref(3)

       imin = coarseloi1
       imax = coarsehii1
       jmin = coarseloi2
       jmax = coarsehii2
       kmin = coarseloi3
       kmax = coarsehii3

       if( dir .eq. 0 ) then
c...
c... X-Direction
c...
c
          do k=kmin,kmax
             do j=jmin,jmax
                do i=imin,imax
                   mmin = j*nrefy
                   mmax = (j+1)*nrefy-1
                   nmin = k*nrefz
                   nmax = (k+1)*nrefz-1

                   left = 0.0d0
                   center = 0.0d0
                   right = 0.0d0
                   do m=mmin,mmax
                      do n=nmin,nmax
                         left   = left   + fine(nrefx*i-1,m,n)
                         center = center + fine(nrefx*i,  m,n)
                         right  = right  + fine(nrefx*i+1,m,n)
                      end do
                   end do
                   left   = left   / dble(nrefy*nrefz)
                   center = center / dble(nrefy*nrefz)
                   right  = right  / dble(nrefy*nrefz)

                   denom = center*(left+right) + 2.0d0*left*right
                   if( denom .ne. 0.0d0 ) then
                      coarse(i,j,k) = 4.0d0 * left*center*right / denom
                   else
                      coarse(i,j,k) = 0.0d0
                   endif

                end do
             end do
          end do

       else if( dir .eq. 1 ) then
c...
c... Y-Direction
c...
c
          do k=kmin,kmax
             do j=jmin,jmax
                do i=imin,imax
                   mmin = i*nrefx
                   mmax = (i+1)*nrefx-1
                   nmin = k*nrefz
                   nmax = (k+1)*nrefz-1

                   left = 0.0d0
                   center = 0.0d0
                   right = 0.0d0
                   do m=mmin,mmax
                      do n=nmin,nmax
                         left   = left   + fine(m,nrefy*j-1,n)
                         center = center + fine(m,nrefy*j,  n)
                         right  = right  + fine(m,nrefy*j+1,n)
                      end do
                   end do
                   left   = left   / dble(nrefx*nrefz)
                   center = center / dble(nrefx*nrefz)
                   right  = right  / dble(nrefx*nrefz)

                   denom = center*(left+right) + 2.0d0*left*right
                   if( denom .ne. 0.0d0 ) then
                      coarse(i,j,k) = 4.0d0 * left*center*right / denom
                   else
                      coarse(i,j,k) = 0.0d0
                   endif

                end do
             end do
          end do
       else if( dir .eq. 2 ) then
c...
c... Z-Direction
c...
c
          do k=kmin,kmax
             do j=jmin,jmax
                do i=imin,imax
                   mmin = i*nrefx
                   mmax = (i+1)*nrefx-1
                   nmin = j*nrefy
                   nmax = (j+1)*nrefy-1

                   left = 0.0d0
                   center = 0.0d0
                   right = 0.0d0
                   do m=mmin,mmax
                      do n=nmin,nmax
                         left   = left   + fine(m,n,nrefz*k-1)
                         center = center + fine(m,n,nrefz*k)
                         right  = right  + fine(m,n,nrefz*k+1)
                      end do
                   end do
                   left   = left   / dble(nrefx*nrefy)
                   center = center / dble(nrefx*nrefy)
                   right  = right  / dble(nrefx*nrefy)

                   denom = center*(left+right) + 2.0d0*left*right
                   if( denom .ne. 0.0d0 ) then
                      coarse(i,j,k) = 4.0d0 * left*center*right / denom
                   else
                      coarse(i,j,k) = 0.0d0
                   endif

                end do
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
      integer ic, jc, kc, if, jf, kf, ir, jr, kr
      real*8 xoff, yoff, zoff, sx, sy, sz

      ir = nref(1)
      jr = nref(2)
      kr = nref(3)
c
      if (lo(1) .eq. dlo(1)) then
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               c(lo(1)-1,jc,kc) =
     &            2.d0 * c(lo(1),jc,kc) - c(lo(1)+1,jc,kc)
            enddo
         enddo
      endif
      if (hi(1) .eq. dhi(1)) then
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               c(hi(1)+1,jc,kc) =
     &            2.d0 * c(hi(1),jc,kc) - c(hi(1)-1,jc,kc)
            enddo
         enddo
      endif
      if (lo(2) .eq. dlo(2)) then
         do kc = lo(3), hi(3)
            do ic = lo(1), hi(1)
               c(ic,lo(2)-1,kc) =
     &            2.d0 * c(ic,lo(2),kc) - c(ic,lo(2)+1,kc)
            enddo
         enddo
      endif
      if (hi(2) .eq. dhi(2)) then
         do kc = lo(3), hi(3)
            do ic = lo(1), hi(1)
               c(ic,hi(2)+1,kc) =
     &            2.d0 * c(ic,hi(2),kc) - c(ic,hi(2)-1,kc)
            enddo
         enddo
      endif
      if (lo(3) .eq. dlo(3)) then
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               c(ic,jc,lo(3)-1) =
     &            2.d0 * c(ic,jc,lo(3)) - c(ic,jc,lo(3)+1)
            enddo
         enddo
      endif
      if (hi(3) .eq. dhi(3)) then
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               c(ic,jc,hi(3)+1) =
     &            2.d0 * c(ic,jc,hi(3)) - c(ic,jc,hi(3)-1)
            enddo
         enddo
      endif

      do kf = lo(3) * kr, (hi(3) + 1) * kr - 1
         kc = kf / kr
         zoff = (mod(kf,kr) + 0.5d0) / kr - 0.5d0
         do jf = lo(2) * jr, (hi(2) + 1) * jr - 1
            jc = jf / jr
            yoff = (mod(jf,jr) + 0.5d0) / jr - 0.5d0
            do if = lo(1) * ir, (hi(1) + 1) * ir - 1
               ic = if / ir
               xoff = (mod(if,ir) + 0.5d0) / ir - 0.5d0
               sx = 0.5d0 * (c(ic+1,jc,kc) - c(ic-1,jc,kc))
               sy = 0.5d0 * (c(ic,jc+1,kc) - c(ic,jc-1,kc))
               sz = 0.5d0 * (c(ic,jc,kc+1) - c(ic,jc,kc-1))
               f(if,jf,kf) = c(ic,jc,kc) +
     &                       xoff * sx + yoff * sy + zoff * sz
            enddo
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
      integer ic, jc, kc, if, jf, kf, nrefx, nrefy, nrefz

      nrefx = nref(1)
      nrefy = nref(2)
      nrefz = nref(3)
c
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               do kf = nrefz*kc, nrefz*(kc+1)-1
                  do jf = nrefy*jc, nrefy*(jc+1)-1
                     do if = nrefx*ic, nrefx*(ic+1)-1
                        f(if,jf,kf) = c(ic,jc,kc) + f(if,jf,kf)
                     enddo
                  enddo
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

      integer dimdec(f)
      integer dimdec(c)
      real*8 f(dimv(f))
      real*8 c(dimv(c))

      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer ratio(BL_SPACEDIM)
      real*8 e,n,w,s,co16,co4
c
      integer i, j

      print *, "FORT_PL_INTERP not implemented in 3D"
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
      integer ic, jc, kc, if, jf, kf, nrefx, nrefy, nrefz
      real*8 fac

      nrefx = nref(1)
      nrefy = nref(2)
      nrefz = nref(3)
      fac = 1.0d0 / (nrefx * nrefy * nrefz)
c
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               c(ic,jc,kc) = 0.0d0
               do kf = nrefz*kc, nrefz*(kc+1)-1
                  do jf = nrefy*jc, nrefy*(jc+1)-1
                     do if = nrefx*ic, nrefx*(ic+1)-1
                        c(ic,jc,kc) = c(ic,jc,kc) + f(if,jf,kf)
                     enddo
                  enddo
               enddo
               c(ic,jc,kc) = c(ic,jc,kc) * fac
            enddo
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
      print *, "FORT_PL_RESTRICT not implemented in 3D"
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

            do j = lo(2), h(2)
               do i = lo(1), hi(1)
                  c(i,j) = f(i,j)
               enddo
            enddo

         endif

      endif
#endif
c
      end
