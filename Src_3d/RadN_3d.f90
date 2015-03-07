
      subroutine bextrp(f, fboxl0, fboxl1, fboxl2, fboxh0, fboxh1, fboxh2, &
           regl0, regl1, regl2, regh0, regh1, regh2)

        implicit none
        integer fboxl0, fboxl1, fboxl2, fboxh0, fboxh1, fboxh2
        integer  regl0,  regl1,  regl2,  regh0,  regh1,  regh2
        real*8 f(fboxl0:fboxh0,fboxl1:fboxh1,fboxl2:fboxh2)
        integer i, j, k

! i direction first:
        do k = regl2, regh2
        do j = regl1, regh1
           i = regl0
           f(i-1,j,k) = 2.d0 * f(i,j,k) - f(i+1,j,k)
           i = regh0
           f(i+1,j,k) = 2.d0 * f(i,j,k) - f(i-1,j,k)
        enddo
        enddo

! j direction second, including edges:
        do k = regl2, regh2
        do i = regl0 - 1, regh0 + 1
           j = regl1
           f(i,j-1,k) = 2.d0 * f(i,j,k) - f(i,j+1,k)
           j = regh1
           f(i,j+1,k) = 2.d0 * f(i,j,k) - f(i,j-1,k)
        enddo
        enddo

! k direction third, including corners:
        do j = regl1 - 1, regh1 + 1
        do i = regl0 - 1, regh0 + 1
           k = regl2
           f(i,j,k-1) = 2.d0 * f(i,j,k) - f(i,j,k+1)
           k = regh2
           f(i,j,k+1) = 2.d0 * f(i,j,k) - f(i,j,k-1)
        enddo
        enddo

! corner results are the same whichever direction we extrapolate first
      end subroutine bextrp


      subroutine lbcoefna(bcoef, &
           bcgrp, bboxl0, bboxl1, bboxl2, bboxh0, bboxh1, bboxh2, &
           regl0, regl1, regl2, regh0, regh1, regh2, &
           spec, sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2, &
           idim)

        implicit none
        integer idim
        integer  regl0,  regl1,  regl2,  regh0,  regh1,  regh2
        integer bboxl0, bboxl1, bboxl2, bboxh0, bboxh1, bboxh2
        integer sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2
        real*8 bcoef(bboxl0:bboxh0,bboxl1:bboxh1,bboxl2:bboxh2)
        real*8 bcgrp(bboxl0:bboxh0,bboxl1:bboxh1,bboxl2:bboxh2)
        real*8 spec(sboxl0:sboxh0,sboxl1:sboxh1,sboxl2:sboxh2)
        integer i, j, k
        if (idim .eq. 0) then
           do k = regl2, regh2
           do j = regl1, regh1
           do i = regl0, regh0
              bcoef(i,j,k) = bcoef(i,j,k) &
                 + 0.5d0 * (spec(i-1,j,k) + spec(i,j,k)) * bcgrp(i,j,k)
           enddo
           enddo
           enddo
        else if (idim .eq. 1) then
           do k = regl2, regh2
           do j = regl1, regh1
           do i = regl0, regh0
              bcoef(i,j,k) = bcoef(i,j,k) &
                 + 0.5d0 * (spec(i,j-1,k) + spec(i,j,k)) * bcgrp(i,j,k)
           enddo
           enddo
           enddo
        else
           do k = regl2, regh2
           do j = regl1, regh1
           do i = regl0, regh0
              bcoef(i,j,k) = bcoef(i,j,k) &
                 + 0.5d0 * (spec(i,j,k-1) + spec(i,j,k)) * bcgrp(i,j,k)
           enddo
           enddo
           enddo
        endif

      end subroutine lbcoefna


      subroutine ljupna(jnew, jboxl0, jboxl1, jboxl2, jboxh0, jboxh1, jboxh2, &
           regl0, regl1, regl2, regh0, regh1, regh2, &
           spec, sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2, &
           accel, nTotal)

        implicit none
        integer nTotal
        integer  regl0,  regl1,  regl2,  regh0,  regh1,  regh2
        integer jboxl0, jboxl1, jboxl2, jboxh0, jboxh1, jboxh2
        integer sboxl0, sboxl1, sboxl2, sboxh0, sboxh1, sboxh2
        real*8 jnew(jboxl0:jboxh0,jboxl1:jboxh1,jboxl2:jboxh2,0:nTotal-1)
        real*8 spec(sboxl0:sboxh0,sboxl1:sboxh1,sboxl2:sboxh2,0:nTotal-1)
        real*8 accel(regl0:regh0,regl1:regh1,regl2:regh2)

        integer i, j, k, n
        do n = 0, nTotal - 1
        do k = regl2, regh2
        do j = regl1, regh1
        do i = regl0, regh0
           jnew(i,j,k,n) = jnew(i,j,k,n) + spec(i,j,k,n) * accel(i,j,k)
        enddo
        enddo
        enddo
        enddo

      end subroutine ljupna
