
#include "EDGEVAR_F.H"

#define dims(a) a loi1, a hii1
#define dimdec(a) dims(a)
#define dimv(a) a loi1:a hii1

       subroutine FORT_SET_GRADIENT(
     &   scalar, dims(scalar),
     &   dx,
     &   ncomp,
     &   gradx, dims(gradx) )
c
       real*8 dx(BL_SPACEDIM)
       real*8 dxinv
       integer i,n
       integer ncomp
       integer dimdec(scalar)
       real*8 scalar(dimv(scalar),ncomp)
       integer dimdec(gradx)
       real*8 gradx(dimv(gradx),ncomp)
c
       dxinv = 1.0d0 / dx(1)
c              
       do n=1,ncomp
          do i=gradxloi1,gradxhii1
             gradx(i,n) = dxinv * ( scalar(i,n)-scalar(i-1,n) )
          end do
       end do
c
       return
       end


       subroutine FORT_EDGE_AVG_DOWN(
     &   fine_x, dims(fine_x),
     &   nref,
     &   coarse_x, dims(coarse_x) )
c
       integer nref(1)
       integer dimdec(fine_x)
       integer dimdec(coarse_x)
       real*8 fine_x(dimv(fine_x))
       real*8 coarse_x(dimv(coarse_x))
c
       integer i
       integer imin,imax
       integer nrefx
c
       nrefx = nref(1)
c...
c... X-Direction
c...
       imin = coarse_xloi1
       imax = coarse_xhii1
c
       do i=imin,imax
          coarse_x(i) = fine_x(nrefx*i)
       end do
c
       return
       end


       subroutine FORT_AVG_DIV(
     &   gradx, dims(gradx),
     &   dx,
     &   ncomp,
     &   div, dims(div) ) 
c
       real*8 dx(BL_SPACEDIM)
       real*8 dxinv
       integer i,n
       integer ncomp
       integer dimdec(div)
       real*8 div(dimv(div), ncomp)
       integer dimdec(gradx)
       real*8 gradx(dimv(gradx), ncomp)
c
       dxinv = 1.0d0 / dx(1)
c
       do n=1,ncomp
          do i=divloi1,divhii1
             div(i,n) = dxinv * ( gradx(i+1,n) - gradx(i,n) )
          end do
       end do
c
       return
       end

       subroutine FORT_EDGE_INTERP(
     &   coarse_x, dims(coarse_x),
     &   nref,
     &   fine_x, dims(fine_x) )
c
       integer nref(1)
       integer dimdec(fine_x)
       integer dimdec(coarse_x)
       real*8 fine_x(dimv(fine_x))
       real*8 coarse_x(dimv(coarse_x))
c
       integer i
       integer imin,imax
       integer nrefx
c
       nrefx = nref(1)
c...
c... X-Direction
c...
       imin = coarse_xloi1
       imax = coarse_xhii1
c
       do i=imin,imax
          fine_x(nrefx*i) = coarse_x(i)
       end do
c
       return
       end


