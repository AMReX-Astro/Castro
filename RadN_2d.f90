
      subroutine bextrp(f, fboxl0, fboxl1, fboxh0, fboxh1, &
           regl0, regl1, regh0, regh1)

        implicit none
        integer fboxl0, fboxl1, fboxh0, fboxh1
        integer  regl0,  regl1,  regh0,  regh1
        real*8 f(fboxl0:fboxh0,fboxl1:fboxh1)
        integer i, j

! i direction first:
        do j = regl1, regh1
           i = regl0
           f(i-1,j) = 2.d0 * f(i,j) - f(i+1,j)
           i = regh0
           f(i+1,j) = 2.d0 * f(i,j) - f(i-1,j)
        enddo

! j direction second, including corners:
        do i = regl0 - 1, regh0 + 1
           j = regl1
           f(i,j-1) = 2.d0 * f(i,j) - f(i,j+1)
           j = regh1
           f(i,j+1) = 2.d0 * f(i,j) - f(i,j-1)
        enddo

! corner results are the same whichever direction we extrapolate first
      end subroutine bextrp


      subroutine lbcoefna(bcoef, &
           bcgrp, bboxl0, bboxl1, bboxh0, bboxh1, &
           regl0, regl1, regh0, regh1, &
           spec, sboxl0, sboxl1, sboxh0, sboxh1, &
           idim)

        implicit none
        integer idim
        integer  regl0,  regl1,  regh0,  regh1
        integer bboxl0, bboxl1, bboxh0, bboxh1
        integer sboxl0, sboxl1, sboxh0, sboxh1
        real*8 bcoef(bboxl0:bboxh0,bboxl1:bboxh1)
        real*8 bcgrp(bboxl0:bboxh0,bboxl1:bboxh1)
        real*8 spec(sboxl0:sboxh0,sboxl1:sboxh1)
        integer i, j
        if (idim .eq. 0) then
           do j = regl1, regh1
           do i = regl0, regh0 + 1
              bcoef(i,j) = bcoef(i,j) &
                 + 0.5d0 * (spec(i-1,j) + spec(i,j)) * bcgrp(i,j)
           enddo
           enddo
        else
           do j = regl1, regh1 + 1
           do i = regl0, regh0
              bcoef(i,j) = bcoef(i,j) &
                 + 0.5d0 * (spec(i,j-1) + spec(i,j)) * bcgrp(i,j)
           enddo
           enddo
        endif

      end subroutine lbcoefna


      subroutine ljupna(jnew, jboxl0, jboxl1, jboxh0, jboxh1, &
           regl0, regl1, regh0, regh1, &
           spec, sboxl0, sboxl1, sboxh0, sboxh1, &
           accel, nTotal)

        implicit none
        integer nTotal
        integer  regl0,  regl1,  regh0, regh1
        integer jboxl0, jboxl1, jboxh0, jboxh1
        integer sboxl0, sboxl1, sboxh0, sboxh1
        real*8 jnew(jboxl0:jboxh0,jboxl1:jboxh1,0:nTotal-1)
        real*8 spec(sboxl0:sboxh0,sboxl1:sboxh1,0:nTotal-1)
        real*8 accel(regl0:regh0,regl1:regh1)

        integer i, j, n
        do n = 0, nTotal - 1
        do j = regl1, regh1
        do i = regl0, regh0
           jnew(i,j,n) = jnew(i,j,n) + spec(i,j,n) * accel(i,j)
        enddo
        enddo
        enddo

      end subroutine ljupna
