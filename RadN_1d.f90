
      subroutine bextrp(f, fboxl0, fboxh0, &
           regl0, regh0)

        implicit none
        integer fboxl0, fboxh0
        integer  regl0,  regh0
        real*8 f(fboxl0:fboxh0)
        integer i

! i direction first:
        i = regl0
        f(i-1) = 2.d0 * f(i) - f(i+1)
        i = regh0
        f(i+1) = 2.d0 * f(i) - f(i-1)

      end subroutine bextrp


      subroutine lbcoefna(bcoef, &
           bcgrp, bboxl0, bboxh0, &
           regl0, regh0, &
           spec, sboxl0, sboxh0, &
           idim)

        implicit none
        integer idim
        integer  regl0,  regh0
        integer bboxl0, bboxh0
        integer sboxl0, sboxh0
        real*8 bcoef(bboxl0:bboxh0)
        real*8 bcgrp(bboxl0:bboxh0)
        real*8 spec(sboxl0:sboxh0)
        integer i
        if (idim .eq. 0) then
           do i = regl0, regh0 + 1
              bcoef(i) = bcoef(i) &
                 + 0.5d0 * (spec(i-1) + spec(i)) * bcgrp(i)
           enddo
        endif

      end subroutine lbcoefna


      subroutine ljupna(jnew, jboxl0, jboxh0, &
           regl0, regh0, &
           spec, sboxl0, sboxh0, &
           accel, nTotal)

        implicit none
        integer nTotal
        integer  regl0,  regh0
        integer jboxl0, jboxh0
        integer sboxl0, sboxh0
        real*8 jnew(jboxl0:jboxh0, 0:nTotal-1)
        real*8 spec(sboxl0:sboxh0, 0:nTotal-1)
        real*8 accel(regl0:regh0)

        integer i, n
        do n = 0, nTotal - 1
        do i = regl0, regh0
           jnew(i,n) = jnew(i,n) + spec(i,n) * accel(i)
        enddo
        enddo

      end subroutine ljupna

