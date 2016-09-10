      real*8 function FLDlambda(r, limiter)
      implicit none
      real*8 r
      integer limiter
      if (limiter .eq. 0) then
c     no limiter
         FLDlambda = 1.d0/3.d0
      else if (limiter < 10) then
c     approximate LP
         FLDlambda = (2.d0 + r) / (6.d0 + r * (3.d0 + r))
      else if (limiter < 20) then
c     Bruenn
         FLDlambda = 1.d0 / (3.d0 + r)
      else if (limiter < 30) then
c     Larsen's square root
         FLDlambda = 1.d0 / sqrt(9.d0 + r**2)
      else if (limiter < 40) then 
c     Minerbo
         if (r .lt. 1.5d0) then
            FLDlambda = 2.d0/(3.d0 + sqrt(9.d0+12.d0*r**2))
         else 
            FLDlambda = 1.d0/(1.d0+r+sqrt(1.d0+2.d0*r))
         end if
      else
         print *, "Unknown limiter ", limiter
         stop
      endif
      end

