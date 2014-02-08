
#include "EDGEVAR_F.H"

#define dims(a) a loi1, a loi2, a loi3, a hii1, a hii2, a hii3
#define dimdec(a) dims(a)
#define dimv(a) a loi1:a hii1,a loi2:a hii2,a loi3:a hii3

       subroutine FORT_SET_GRADIENT(
     &   scalar, dims(scalar),
     &   dx,
     &   ncomp,
     &   gradx, dims(gradx),
     &   grady, dims(grady),
     &   gradz, dims(gradz) ) 
c
       real*8 dx(BL_SPACEDIM)
       real*8 dxinv, dyinv, dzinv
       integer i,j,k,n
       integer ncomp
       integer dimdec(scalar)
       real*8 scalar(dimv(scalar),ncomp)
       integer dimdec(gradx)
       real*8 gradx(dimv(gradx),ncomp)
       integer dimdec(grady)
       real*8 grady(dimv(grady),ncomp)
       integer dimdec(gradz)
       real*8 gradz(dimv(gradz),ncomp)
c
       dxinv = 1.0d0 / dx(1)
       dyinv = 1.0d0 / dx(2)
       dzinv = 1.0d0 / dx(3)
c              
       do n=1,ncomp
          do k=gradxloi3,gradxhii3
             do j=gradxloi2,gradxhii2
                do i=gradxloi1,gradxhii1
                   gradx(i,j,k,n) =
     @                dxinv * ( scalar(i,j,k,n)-scalar(i-1,j,k,n) )
                end do
             end do
          end do
       end do
c
       do n=1,ncomp
          do k=gradxloi3,gradxhii3
             do j=gradyloi2,gradyhii2
                do i=gradyloi1,gradyhii1
                   grady(i,j,k,n) =
     @                dyinv * (scalar(i,j,k,n)-scalar(i,j-1,k,n) )
                end do
             end do
          end do
       end do
c
       do n=1,ncomp
          do k=gradxloi3,gradxhii3
             do j=gradyloi2,gradyhii2
                do i=gradyloi1,gradyhii1
                   gradz(i,j,k,n) =
     @                dzinv * (scalar(i,j,k,n)-scalar(i,j,k-1,n) )
                end do
             end do
          end do
       end do
c
       return
       end


       subroutine FORT_EDGE_AVG_DOWN(
     &   fine_x, dims(fine_x),
     &   fine_y, dims(fine_y),
     &   fine_z, dims(fine_z),
     &   nref,
     &   coarse_x, dims(coarse_x),
     &   coarse_y, dims(coarse_y),
     &   coarse_z, dims(coarse_z) )
c
       integer nref(3)
       integer dimdec(fine_x)
       integer dimdec(fine_y)
       integer dimdec(fine_z)
       integer dimdec(coarse_x)
       integer dimdec(coarse_y)
       integer dimdec(coarse_z)
       real*8 fine_x(dimv(fine_x))
       real*8 fine_y(dimv(fine_y))
       real*8 fine_z(dimv(fine_z))
       real*8 coarse_x(dimv(coarse_x))
       real*8 coarse_y(dimv(coarse_y))
       real*8 coarse_z(dimv(coarse_z))
c
       integer i,j,k,m,n
       integer imin,imax,jmin,jmax,kmin,kmax,mmin,mmax,nmin,nmax
       integer nrefx, nrefy, nrefz
       real*8 sum
c
       nrefx = nref(1)
       nrefy = nref(2)
       nrefz = nref(3)
c...
c... X-Direction
c...
       imin = coarse_xloi1
       imax = coarse_xhii1
       jmin = coarse_xloi2
       jmax = coarse_xhii2
       kmin = coarse_xloi3
       kmax = coarse_xhii3
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
                      sum = sum + fine_x(nrefx*i,m,n)
                   end do
                end do
                coarse_x(i,j,k) = sum / dble(nrefy*nrefz)
             end do
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
       kmin = coarse_yloi3
       kmax = coarse_yhii3
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
                      sum = sum + fine_y(m,nrefy*j,n)
                   end do
                end do
                coarse_y(i,j,k) = sum / dble(nrefx*nrefz)
             end do
          end do
       end do
c
c...
c... Z-Direction
c...
       imin = coarse_zloi1
       imax = coarse_zhii1
       jmin = coarse_zloi2
       jmax = coarse_zhii2
       kmin = coarse_zloi3
       kmax = coarse_zhii3
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
                      sum = sum + fine_z(m,n,nrefz*k)
                   end do
                end do
                coarse_z(i,j,k) = sum / dble(nrefx*nrefy)
             end do
          end do
       end do
c
       return
       end


       subroutine FORT_AVG_DIV(
     &   gradx, dims(gradx),
     &   grady, dims(grady),
     &   gradz, dims(gradz),
     &   dx,
     &   ncomp,
     &   div, dims(div) ) 
c
       real*8 dx(BL_SPACEDIM)
       real*8 dxinv, dyinv, dzinv
       integer i,j,k,n
       integer ncomp
       integer dimdec(div)
       real*8 div(dimv(div),ncomp)
       integer dimdec(gradx)
       real*8 gradx(dimv(gradx),ncomp)
       integer dimdec(grady)
       real*8 grady(dimv(grady),ncomp)
       integer dimdec(gradz)
       real*8 gradz(dimv(gradz),ncomp)
c
       dxinv = 1.0d0 / dx(1)
       dyinv = 1.0d0 / dx(2)
       dzinv = 1.0d0 / dx(3)
c
       do n=1,ncomp
          do k=divloi3,divhii3
             do j=divloi2,divhii2
                do i=divloi1,divhii1
                   div(i,j,k,n) =
     @                dxinv * ( gradx(i+1,j,k,n) - gradx(i,j,k,n) ) +
     @                dyinv * ( grady(i,j+1,k,n) - grady(i,j,k,n) ) +
     @                dzinv * ( gradz(i,j,k+1,n) - gradz(i,j,k,n) )
                end do
             end do
          end do
       end do
c
       return
       end

       subroutine FORT_EDGE_INTERP(
     &   coarse_x, dims(coarse_x),
     &   coarse_y, dims(coarse_y),
     &   coarse_z, dims(coarse_z),
     &   nref,
     &   fine_x, dims(fine_x),
     &   fine_y, dims(fine_y),
     &   fine_z, dims(fine_z) )
c
       integer nref(3)
       integer dimdec(fine_x)
       integer dimdec(fine_y)
       integer dimdec(fine_z)
       integer dimdec(coarse_x)
       integer dimdec(coarse_y)
       integer dimdec(coarse_z)
       real*8 fine_x(dimv(fine_x))
       real*8 fine_y(dimv(fine_y))
       real*8 fine_z(dimv(fine_z))
       real*8 coarse_x(dimv(coarse_x))
       real*8 coarse_y(dimv(coarse_y))
       real*8 coarse_z(dimv(coarse_z))
c
       integer i,j,k,m,n
       integer imin,imax,jmin,jmax,kmin,kmax,mmin,mmax,nmin,nmax
       integer nrefx, nrefy, nrefz
c
       nrefx = nref(1)
       nrefy = nref(2)
       nrefz = nref(3)
c...
c... X-Direction
c...
       imin = coarse_xloi1
       imax = coarse_xhii1
       jmin = coarse_xloi2
       jmax = coarse_xhii2
       kmin = coarse_xloi3
       kmax = coarse_xhii3
c
       do k=kmin,kmax
          do j=jmin,jmax
             do i=imin,imax
                mmin = j*nrefy
                mmax = (j+1)*nrefy-1
                nmin = k*nrefz
                nmax = (k+1)*nrefz-1
                do m=mmin,mmax
                   do n=nmin,nmax
                      fine_x(nrefx*i,m,n) = coarse_x(i,j,k)
                   end do
                end do
             end do
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
       kmin = coarse_yloi3
       kmax = coarse_yhii3
c
       do k=kmin,kmax
          do j=jmin,jmax
             do i=imin,imax
                mmin = i*nrefx
                mmax = (i+1)*nrefx-1
                nmin = k*nrefz
                nmax = (k+1)*nrefz-1
                do m=mmin,mmax
                   do n=nmin,nmax
                      fine_y(m,nrefy*j,n) = coarse_y(i,j,k)
                   end do
                end do
             end do
          end do
       end do
c
c...
c... Z-Direction
c...
       imin = coarse_zloi1
       imax = coarse_zhii1
       jmin = coarse_zloi2
       jmax = coarse_zhii2
       kmin = coarse_zloi3
       kmax = coarse_zhii3
c
       do k=kmin,kmax
          do j=jmin,jmax
             do i=imin,imax
                mmin = i*nrefx
                mmax = (i+1)*nrefx-1
                nmin = j*nrefy
                nmax = (j+1)*nrefy-1
                do m=mmin,mmax
                   do n=nmin,nmax
                      fine_z(m,n,nrefz*k) = coarse_z(i,j,k)
                   end do
                end do
             end do
          end do
       end do
c
       return
       end
