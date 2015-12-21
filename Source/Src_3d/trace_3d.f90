module trace_module

  implicit none

  private

  public tracexy, tracez
  
contains

      subroutine tracexy(q,c,qd_lo,qd_hi, &
                         dqx,dqy,dq_lo,dq_hi, &
                         qxm,qxp,qym,qyp,qpd_lo,qpd_hi, &
                         ilo1,ilo2,ihi1,ihi2,dx,dt,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QESGS, QPRES, QFS, QFX, &
                                     npassive, qpass_map, small_dens, ppm_type
      use bl_constants_module
      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: dq_lo(3), dq_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer          :: ilo1, ilo2, ihi1, ihi2
      integer          :: kc, k3d

      double precision ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
      double precision ::     c(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))

      double precision ::  dqx(dq_lo(1):dq_hi(1),dq_lo(2):dq_hi(2),dq_lo(3):dq_hi(3),QVAR)
      double precision ::  dqy(dq_lo(1):dq_hi(1),dq_lo(2):dq_hi(2),dq_lo(3):dq_hi(3),QVAR)

      double precision :: qxm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
      double precision :: qxp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
      double precision :: qym(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
      double precision :: qyp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
      double precision :: dx(3), dt

      ! Local variables
      integer i, j, n, ipassive

      double precision dtdx, dtdy
      double precision cc, csq, rho, u, v, w, p, rhoe
      double precision drho, du, dv, dw, dp, drhoe

      double precision enth, alpham, alphap, alpha0r, alpha0e
      double precision alpha0u, alpha0v, alpha0w
      double precision spminus, spplus, spzero
      double precision apright, amright, azrright, azeright
      double precision azu1rght, azv1rght, azw1rght
      double precision apleft, amleft, azrleft, azeleft
      double precision azu1left, azv1left, azw1left
      double precision acmprght, acmpleft, acmpbot, acmptop
      double precision ascmprght, ascmpleft, ascmpbot, ascmptop

      dtdx = dt/dx(1)
      dtdy = dt/dx(2)

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in tracexy with ppm_type != 0'
        call bl_error("Error:: Castro_advection_3d.f90 :: tracexy")
      end if

      !!!!!!!!!!!!!!!
      ! NON-PPM CODE
      !!!!!!!!!!!!!!!
      
      ! Compute left and right traced states
      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,QRHO)
            u = q(i,j,k3d,QU)
            v = q(i,j,k3d,QV)
            w = q(i,j,k3d,QW)
            p = q(i,j,k3d,QPRES)
            rhoe = q(i,j,k3d,QREINT)
            enth = (rhoe+p)/(rho*csq)

            drho = dqx(i,j,kc,QRHO)
            du = dqx(i,j,kc,QU)
            dv = dqx(i,j,kc,QV)
            dw = dqx(i,j,kc,QW)
            dp = dqx(i,j,kc,QPRES)
            drhoe = dqx(i,j,kc,QREINT)

            alpham = HALF*((dp/(rho*cc)) - du)*(rho/cc)
            alphap = HALF*((dp/(rho*cc)) + du)*(rho/cc)
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0v = dv
            alpha0w = dw

            if (u-cc .gt. ZERO) then
               spminus = -ONE
            else
               spminus = (u-cc)*dtdx
            endif
            if (u+cc .gt. ZERO) then
               spplus = -ONE
            else
               spplus = (u+cc)*dtdx
            endif
            if (u .gt. ZERO) then
               spzero = -ONE
            else
               spzero = u*dtdx
            endif

            apright = HALF*(-ONE - spplus )*alphap
            amright = HALF*(-ONE - spminus)*alpham
            azrright = HALF*(-ONE - spzero )*alpha0r
            azeright = HALF*(-ONE - spzero )*alpha0e
            azv1rght = HALF*(-ONE - spzero )*alpha0v
            azw1rght = HALF*(-ONE - spzero )*alpha0w

            if (i .ge. ilo1) then
               qxp(i,j,kc,QRHO) = rho + apright + amright + azrright
               qxp(i,j,kc,QRHO) = max(small_dens,qxp(i,j,kc,QRHO))
               qxp(i,j,kc,QU) = u + (apright - amright)*(cc/rho)
               qxp(i,j,kc,QV) = v + azv1rght
               qxp(i,j,kc,QW) = w + azw1rght
               qxp(i,j,kc,QPRES) = p + (apright + amright)*csq
               qxp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
            end if

            if (u-cc .ge. ZERO) then
               spminus = (u-cc)*dtdx
            else
               spminus = ONE
            endif
            if (u+cc .ge. ZERO) then
               spplus = (u+cc)*dtdx
            else
               spplus = ONE
            endif
            if (u .ge. ZERO) then
               spzero = u*dtdx
            else
               spzero = ONE
            endif

            apleft = HALF*(ONE - spplus )*alphap
            amleft = HALF*(ONE - spminus)*alpham
            azrleft = HALF*(ONE - spzero )*alpha0r
            azeleft = HALF*(ONE - spzero )*alpha0e
            azv1left = HALF*(ONE - spzero )*alpha0v
            azw1left = HALF*(ONE - spzero )*alpha0w

            if (i .le. ihi1) then
               qxm(i+1,j,kc,QRHO) = rho + apleft + amleft + azrleft
               qxm(i+1,j,kc,QRHO) = max(small_dens, qxm(i+1,j,kc,QRHO))
               qxm(i+1,j,kc,QU) = u + (apleft - amleft)*(cc/rho)
               qxm(i+1,j,kc,QV) = v + azv1left
               qxm(i+1,j,kc,QW) = w + azw1left
               qxm(i+1,j,kc,QPRES) = p + (apleft + amleft)*csq
               qxm(i+1,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
            endif

         enddo
      enddo

      ! Treat K as a passively advected quantity
      if (QESGS .gt. -1) then
         n = QESGS
         do j = ilo2-1, ihi2+1
            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. ZERO) then
                  spzero = -ONE
               else
                  spzero = u*dtdx
               endif
               acmprght = HALF*(-ONE - spzero )*dqx(i,j,kc,n)
               qxp(i,j,kc,n) = q(i,j,k3d,n) + acmprght
            enddo
 
            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. ZERO) then
                  spzero = u*dtdx
               else
                  spzero = ONE
               endif
               acmpleft = HALF*(ONE - spzero )*dqx(i,j,kc,n)
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + acmpleft
            enddo
         enddo
      endif

      do ipassive = 1, npassive
         n = qpass_map(ipassive)

         do j = ilo2-1, ihi2+1

            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. ZERO) then
                  spzero = -ONE
               else
                  spzero = u*dtdx
               endif
               acmprght = HALF*(-ONE - spzero )*dqx(i,j,kc,n)
               qxp(i,j,kc,n) = q(i,j,k3d,n) + acmprght
            enddo

            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. ZERO) then
                  spzero = u*dtdx
               else
                  spzero = ONE
               endif
               acmpleft = HALF*(ONE - spzero )*dqx(i,j,kc,n)
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + acmpleft
            enddo

         enddo
      enddo

      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,QRHO)
            u = q(i,j,k3d,QU)
            v = q(i,j,k3d,QV)
            w = q(i,j,k3d,QW)
            p = q(i,j,k3d,QPRES)
            rhoe = q(i,j,k3d,QREINT)
            enth = (rhoe+p)/(rho*csq)

            drho = dqy(i,j,kc,QRHO)
            du = dqy(i,j,kc,QU)
            dv = dqy(i,j,kc,QV)
            dw = dqy(i,j,kc,QW)
            dp = dqy(i,j,kc,QPRES)
            drhoe = dqy(i,j,kc,QREINT)

            alpham = HALF*((dp/(rho*cc)) - dv)*(rho/cc)
            alphap = HALF*((dp/(rho*cc)) + dv)*(rho/cc)
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0w = dw

            if (v-cc .gt. ZERO) then
               spminus = -ONE
            else
               spminus = (v-cc)*dtdy
            endif
            if (v+cc .gt. ZERO) then
               spplus = -ONE
            else
               spplus = (v+cc)*dtdy
            endif
            if (v .gt. ZERO) then
               spzero = -ONE
            else
               spzero = v*dtdy
            endif

            apright = HALF*(-ONE - spplus )*alphap
            amright = HALF*(-ONE - spminus)*alpham
            azrright = HALF*(-ONE - spzero )*alpha0r
            azeright = HALF*(-ONE - spzero )*alpha0e
            azu1rght = HALF*(-ONE - spzero )*alpha0u
            azw1rght = HALF*(-ONE - spzero )*alpha0w

            if (j .ge. ilo2) then
               qyp(i,j,kc,QRHO) = rho + apright + amright + azrright
               qyp(i,j,kc,QRHO) = max(small_dens, qyp(i,j,kc,QRHO))
               qyp(i,j,kc,QV) = v + (apright - amright)*(cc/rho)
               qyp(i,j,kc,QU) = u + azu1rght
               qyp(i,j,kc,QW) = w + azw1rght
               qyp(i,j,kc,QPRES) = p + (apright + amright)*csq
               qyp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
            end if

            if (v-cc .ge. ZERO) then
               spminus = (v-cc)*dtdy
            else
               spminus = ONE
            endif
            if (v+cc .ge. ZERO) then
               spplus = (v+cc)*dtdy
            else
               spplus = ONE
            endif
            if (v .ge. ZERO) then
               spzero = v*dtdy
            else
               spzero = ONE
            endif

            apleft = HALF*(ONE - spplus )*alphap
            amleft = HALF*(ONE - spminus)*alpham
            azrleft = HALF*(ONE - spzero )*alpha0r
            azeleft = HALF*(ONE - spzero )*alpha0e
            azu1left = HALF*(ONE - spzero )*alpha0u
            azw1left = HALF*(ONE - spzero )*alpha0w

            if (j .le. ihi2) then
               qym(i,j+1,kc,QRHO) = rho + apleft + amleft + azrleft
               qym(i,j+1,kc,QRHO) = max(small_dens, qym(i,j+1,kc,QRHO))
               qym(i,j+1,kc,QV) = v + (apleft - amleft)*(cc/rho)
               qym(i,j+1,kc,QU) = u + azu1left
               qym(i,j+1,kc,QW) = w + azw1left
               qym(i,j+1,kc,QPRES) = p + (apleft + amleft)*csq
               qym(i,j+1,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
            endif

         enddo
      enddo

      ! Treat K as a passively advected quantity
      if (QESGS .gt. -1) then
         n = QESGS
         do i = ilo1-1, ihi1+1
 
            ! Top state
            do j = ilo2, ihi2+1
               v = q(i,j,k3d,QV)
               if (v .gt. ZERO) then
                  spzero = -ONE
               else
                  spzero = v*dtdy
               endif
               acmptop = HALF*(-ONE - spzero )*dqy(i,j,kc,n)
               qyp(i,j,kc,n) = q(i,j,k3d,n) + acmptop
            enddo
 
            ! Bottom state
            do j = ilo2-1, ihi2
               v = q(i,j,k3d,QV)
               if (v .ge. ZERO) then
                  spzero = v*dtdy
               else
                  spzero = ONE
               endif
               acmpbot = HALF*(ONE - spzero )*dqy(i,j,kc,n)
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + acmpbot
            enddo
         enddo
      endif

      do ipassive = 1, npassive
         n = qpass_map(ipassive)

         do i = ilo1-1, ihi1+1

            ! Top state
            do j = ilo2, ihi2+1
               v = q(i,j,k3d,QV)
               if (v .gt. ZERO) then
                  spzero = -ONE
               else
                  spzero = v*dtdy
               endif
               acmptop = HALF*(-ONE - spzero )*dqy(i,j,kc,n)
               qyp(i,j,kc,n) = q(i,j,k3d,n) + acmptop
            enddo

            ! Bottom state
            do j = ilo2-1, ihi2
               v = q(i,j,k3d,QV)
               if (v .ge. ZERO) then
                  spzero = v*dtdy
               else
                  spzero = ONE
               endif
               acmpbot = HALF*(ONE - spzero )*dqy(i,j,kc,n)
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + acmpbot
            enddo

         enddo
      enddo

    end subroutine tracexy

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine tracez(q,c,qd_lo,qd_hi, &
                        dqz,dq_lo,dq_hi, &
                        qzm,qzp,qpd_lo,qpd_hi, &
                        ilo1,ilo2,ihi1,ihi2,dx,dt,km,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QESGS, QPRES, QFS, QFX, &
                                     npassive, qpass_map, small_dens, ppm_type
      use bl_constants_module

      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: dq_lo(3), dq_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer          :: ilo1, ilo2, ihi1, ihi2
      integer          :: km, kc, k3d

      double precision ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
      double precision ::     c(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))

      double precision :: dqz( dq_lo(1):dq_hi(1),  dq_lo(2):dq_hi(2),  dq_lo(3):dq_hi(3), QVAR)
      double precision :: qzm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
      double precision :: qzp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
      double precision :: dx(3), dt

      ! Local variables
      integer i, j
      integer n, ipassive

      double precision dtdz
      double precision cc, csq, rho, u, v, w, p, rhoe

      double precision drho, du, dv, dw, dp, drhoe
      double precision enth, alpham, alphap, alpha0r, alpha0e
      double precision alpha0u, alpha0v
      double precision spminus, spplus, spzero
      double precision apright, amright, azrright, azeright
      double precision azu1rght, azv1rght
      double precision apleft, amleft, azrleft, azeleft
      double precision azu1left, azv1left
      double precision acmpbot, acmptop
      double precision ascmpbot, ascmptop

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in tracez with ppm_type != 0'
        call bl_error("Error:: Castro_advection_3d.f90 :: tracez")
      end if

      dtdz = dt/dx(3)
      
      !!!!!!!!!!!!!!!
      ! NON-PPM CODE
      !!!!!!!!!!!!!!!
      
      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,QRHO)
            u = q(i,j,k3d,QU)
            v = q(i,j,k3d,QV)
            w = q(i,j,k3d,QW)
            p = q(i,j,k3d,QPRES)
            rhoe = q(i,j,k3d,QREINT)
            enth = (rhoe+p)/(rho*csq)

            drho = dqz(i,j,kc,QRHO)
            du = dqz(i,j,kc,QU)
            dv = dqz(i,j,kc,QV)
            dw = dqz(i,j,kc,QW)
            dp = dqz(i,j,kc,QPRES)
            drhoe = dqz(i,j,kc,QREINT)

            alpham = HALF*((dp/(rho*cc)) - dw)*(rho/cc)
            alphap = HALF*((dp/(rho*cc)) + dw)*(rho/cc)
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0v = dv

            if (w-cc .gt. ZERO) then
               spminus = -ONE
            else
               spminus = (w-cc)*dtdz
            endif
            if (w+cc .gt. ZERO) then
               spplus = -ONE
            else
               spplus = (w+cc)*dtdz
            endif
            if (w .gt. ZERO) then
               spzero = -ONE
            else
               spzero = w*dtdz
            endif

            apright = HALF*(-ONE - spplus )*alphap
            amright = HALF*(-ONE - spminus)*alpham
            azrright = HALF*(-ONE - spzero )*alpha0r
            azeright = HALF*(-ONE - spzero )*alpha0e
            azu1rght = HALF*(-ONE - spzero )*alpha0u
            azv1rght = HALF*(-ONE - spzero )*alpha0v

            qzp(i,j,kc,QRHO) = rho + apright + amright + azrright
            qzp(i,j,kc,QRHO) = max(small_dens, qzp(i,j,kc,QRHO))
            qzp(i,j,kc,QW) = w + (apright - amright)*(cc/rho)
            qzp(i,j,kc,QU) = u + azu1rght
            qzp(i,j,kc,QV) = v + azv1rght
            qzp(i,j,kc,QPRES) = p + (apright + amright)*csq
            qzp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

            ! repeat above with km (k3d-1) to get qzm at kc
            cc = c(i,j,k3d-1)
            csq = cc**2
            rho = q(i,j,k3d-1,QRHO)
            u = q(i,j,k3d-1,QU)
            v = q(i,j,k3d-1,QV)
            w = q(i,j,k3d-1,QW)
            p = q(i,j,k3d-1,QPRES)
            rhoe = q(i,j,k3d-1,QREINT)
            enth = (rhoe+p)/(rho*csq)

            drho = dqz(i,j,km,QRHO)
            du = dqz(i,j,km,QU)
            dv = dqz(i,j,km,QV)
            dw = dqz(i,j,km,QW)
            dp = dqz(i,j,km,QPRES)
            drhoe = dqz(i,j,km,QREINT)

            alpham = HALF*((dp/(rho*cc)) - dw)*(rho/cc)
            alphap = HALF*((dp/(rho*cc)) + dw)*(rho/cc)
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0v = dv

            if (w-cc .ge. ZERO) then
               spminus = (w-cc)*dtdz
            else
               spminus = ONE
            endif
            if (w+cc .ge. ZERO) then
               spplus = (w+cc)*dtdz
            else
               spplus = ONE
            endif
            if (w .ge. ZERO) then
               spzero = w*dtdz
            else
               spzero = ONE
            endif

            apleft = HALF*(ONE - spplus )*alphap
            amleft = HALF*(ONE - spminus)*alpham
            azrleft = HALF*(ONE - spzero )*alpha0r
            azeleft = HALF*(ONE - spzero )*alpha0e
            azu1left = HALF*(ONE - spzero )*alpha0u
            azv1left = HALF*(ONE - spzero )*alpha0v

            qzm(i,j,kc,QRHO) = rho + apleft + amleft + azrleft
            qzm(i,j,kc,QRHO) = max(small_dens, qzm(i,j,kc,QRHO))
            qzm(i,j,kc,QW) = w + (apleft - amleft)*(cc/rho)
            qzm(i,j,kc,QU) = u + azu1left
            qzm(i,j,kc,QV) = v + azv1left
            qzm(i,j,kc,QPRES) = p + (apleft + amleft)*csq
            qzm(i,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

         enddo
      enddo

      ! Treat K as a passively advected quantity
      if (QESGS .gt. -1) then
         n = QESGS
         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1
 
               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. ZERO) then
                  spzero = -ONE
               else
                  spzero = w*dtdz
               endif
               acmptop = HALF*(-ONE - spzero )*dqz(i,j,kc,n)
               qzp(i,j,kc,n) = q(i,j,k3d,n) + acmptop
 
               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. ZERO) then
                  spzero = w*dtdz
               else
                  spzero = ONE
               endif
               acmpbot = HALF*(ONE - spzero )*dqz(i,j,km,n)
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + acmpbot
 
            enddo
         enddo
      endif

      do ipassive = 1, npassive
         n = qpass_map(ipassive)

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. ZERO) then
                  spzero = -ONE
               else
                  spzero = w*dtdz
               endif
               acmptop = HALF*(-ONE - spzero )*dqz(i,j,kc,n)
               qzp(i,j,kc,n) = q(i,j,k3d,n) + acmptop

               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. ZERO) then
                  spzero = w*dtdz
               else
                  spzero = ONE
               endif
               acmpbot = HALF*(ONE - spzero )*dqz(i,j,km,n)
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + acmpbot
            enddo
         enddo
      enddo

    end subroutine tracez

end module trace_module
