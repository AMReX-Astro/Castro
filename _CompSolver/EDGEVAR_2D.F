
#include "EDGEVAR_F.H"

#define dims(a) a loi1, a loi2, a hii1, a hii2
#define dimdec(a) dims(a)
#define dimv(a) a loi1:a hii1,a loi2:a hii2

       subroutine FORT_SET_GRADIENT(
     &   scalar, dims(scalar),
     &   dx,
     &   ncomp,
     &   gradx, dims(gradx),
     &   grady, dims(grady) ) 
c
       real*8 dx(BL_SPACEDIM)
       real*8 dxinv, dyinv
       integer i,j,n
       integer ncomp
       integer dimdec(scalar)
       real*8 scalar(dimv(scalar),ncomp)
       integer dimdec(gradx)
       real*8 gradx(dimv(gradx),ncomp)
       integer dimdec(grady)
       real*8 grady(dimv(grady),ncomp)
c
       dxinv = 1.0d0 / dx(1)
       dyinv = 1.0d0 / dx(2)
c              
       do n=1,ncomp
          do j=gradxloi2,gradxhii2
             do i=gradxloi1,gradxhii1
                gradx(i,j,n) = dxinv * ( scalar(i,j,n)-scalar(i-1,j,n) )
             end do
          end do
       end do
c
       do n=1,ncomp
          do j=gradyloi2,gradyhii2
             do i=gradyloi1,gradyhii1
                grady(i,j,n) = dyinv * (scalar(i,j,n)-scalar(i,j-1,n) )
             end do
          end do
       end do
c
       return
       end


       subroutine FORT_EDGE_AVG_DOWN(
     &   fine_x, dims(fine_x),
     &   fine_y, dims(fine_y),
     &   nref,
     &   coarse_x, dims(coarse_x),
     &   coarse_y, dims(coarse_y) )
c
       integer nref(2)
       integer dimdec(fine_x)
       integer dimdec(fine_y)
       integer dimdec(coarse_x)
       integer dimdec(coarse_y)
       real*8 fine_x(dimv(fine_x))
       real*8 fine_y(dimv(fine_y))
       real*8 coarse_x(dimv(coarse_x))
       real*8 coarse_y(dimv(coarse_y))
c
       integer i,j,k
       integer imin,imax,jmin,jmax,kmin,kmax
       integer nrefx, nrefy
       real*8 sum
c
       nrefx = nref(1)
       nrefy = nref(2)
c...
c... X-Direction
c...
       imin = coarse_xloi1
       imax = coarse_xhii1
       jmin = coarse_xloi2
       jmax = coarse_xhii2
c
       do j=jmin,jmax
          do i=imin,imax
             kmin = j*nrefy
             kmax = (j+1)*nrefy-1
             sum = 0.0
             do k=kmin,kmax
                sum = sum + fine_x(nrefx*i,k)
             end do
             coarse_x(i,j) = sum / dble(nrefy)
          end do
       end do
c
c...
c... Y-Direction
c...
       imin = coarse_yloi1
       imax = coarse_yhii1
       jmin = coarse_yloi2
       jmax = coarse_yhii2
c
       do j=jmin,jmax
          do i=imin,imax
             kmin = i*nrefx
             kmax = (i+1)*nrefx-1
             sum = 0.0
             do k=kmin,kmax
                sum = sum + fine_y(k,nrefy*j)
             end do
             coarse_y(i,j) = sum / dble(nrefx)
          end do
       end do
c
       return
       end


       subroutine FORT_AVG_DIV(
     &   gradx, dims(gradx),
     &   grady, dims(grady),
     &   dx,
     &   ncomp,
     &   div, dims(div) ) 
c
       real*8 dx(BL_SPACEDIM)
       real*8 dxinv, dyinv
       integer i,j,n
       integer ncomp
       integer dimdec(div)
       real*8 div(dimv(div), ncomp)
       integer dimdec(gradx)
       real*8 gradx(dimv(gradx), ncomp)
       integer dimdec(grady)
       real*8 grady(dimv(grady), ncomp)
c
       dxinv = 1.0d0 / dx(1)
       dyinv = 1.0d0 / dx(2)
c
       do n=1,ncomp
          do j=divloi2,divhii2
             do i=divloi1,divhii1
                div(i,j,n) =              dxinv * ( gradx(i+1,j,n) - gradx(i,j,n) )
                div(i,j,n) = div(i,j,n) + dyinv * ( grady(i,j+1,n) - grady(i,j,n) )
             end do
          end do
       end do
c
       return
       end

       subroutine FORT_EDGE_INTERP(
     &   coarse_x, dims(coarse_x),
     &   coarse_y, dims(coarse_y),
     &   nref,
     &   fine_x, dims(fine_x),
     &   fine_y, dims(fine_y) )
c
       integer nref(2)
       integer dimdec(fine_x)
       integer dimdec(fine_y)
       integer dimdec(coarse_x)
       integer dimdec(coarse_y)
       real*8 fine_x(dimv(fine_x))
       real*8 fine_y(dimv(fine_y))
       real*8 coarse_x(dimv(coarse_x))
       real*8 coarse_y(dimv(coarse_y))
c
       integer i,j,k
       integer imin,imax,jmin,jmax,kmin,kmax
       integer nrefx, nrefy
c
       nrefx = nref(1)
       nrefy = nref(2)
c...
c... X-Direction
c...
       imin = coarse_xloi1
       imax = coarse_xhii1
       jmin = coarse_xloi2
       jmax = coarse_xhii2
c
       do j=jmin,jmax
          do i=imin,imax
             kmin = j*nrefy
             kmax = (j+1)*nrefy-1
             do k=kmin,kmax
                fine_x(nrefx*i,k) = coarse_x(i,j)
             enddo
          end do
       end do
c
c...
c... Y-Direction
c...
       imin = coarse_yloi1
       imax = coarse_yhii1
       jmin = coarse_yloi2
       jmax = coarse_yhii2
c
       do j=jmin,jmax
          do i=imin,imax
             kmin = i*nrefx
             kmax = (i+1)*nrefx-1
             do k=kmin,kmax
                fine_y(k,nrefy*j) = coarse_y(i,j)
             end do
          end do
       end do
c
       return
       end


