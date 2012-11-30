
      subroutine tracexy(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         dqx,dqy,dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3, &
                         qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QESGS, QPRES, QFA, QFS, QFX, nadv, small_dens, &
                                     ppm_type
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2
      integer kc,k3d

      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      double precision  dqx(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)
      double precision  dqy(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)

      double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision dx, dy, dt

      ! Local variables
      integer i, j, n
      integer iadv, ispec, iaux

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

      dtdx = dt/dx
      dtdy = dt/dy

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in tracexy with ppm_type != 0'
        call bl_error("Error:: Castro_advection_3d.f90 :: tracexy")
      end if

      !!!!!!!!!!!!!!!
      ! NON-PPM CODE
      !!!!!!!!!!!!!!!
      
      ! Compute left and right traced states
      !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth,drho,du,dv,dw,dp,drhoe,alpham,alphap,alpha0r,alpha0e) &
      !$OMP PRIVATE(alpha0v,alpha0w,spminus,spplus,spzero,apright,amright,azrright,azeright,azv1rght,azw1rght,apleft) &
      !$OMP PRIVATE(amleft,azrleft,azeleft,azv1left,azw1left)
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
            enth = ( (rhoe+p)/rho )/csq

            drho = dqx(i,j,kc,QRHO)
            du = dqx(i,j,kc,QU)
            dv = dqx(i,j,kc,QV)
            dw = dqx(i,j,kc,QW)
            dp = dqx(i,j,kc,QPRES)
            drhoe = dqx(i,j,kc,QREINT)

            alpham = 0.5d0*(dp/(rho*cc) - du)*rho/cc
            alphap = 0.5d0*(dp/(rho*cc) + du)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0v = dv
            alpha0w = dw

            if (u-cc .gt. 0.d0) then
               spminus = -1.d0
            else
               spminus = (u-cc)*dtdx
            endif
            if (u+cc .gt. 0.d0) then
               spplus = -1.d0
            else
               spplus = (u+cc)*dtdx
            endif
            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif

            apright = 0.5d0*(-1.d0 - spplus )*alphap
            amright = 0.5d0*(-1.d0 - spminus)*alpham
            azrright = 0.5d0*(-1.d0 - spzero )*alpha0r
            azeright = 0.5d0*(-1.d0 - spzero )*alpha0e
            azv1rght = 0.5d0*(-1.d0 - spzero )*alpha0v
            azw1rght = 0.5d0*(-1.d0 - spzero )*alpha0w

            if (i .ge. ilo1) then
               qxp(i,j,kc,QRHO) = rho + apright + amright + azrright
               qxp(i,j,kc,QRHO) = max(small_dens,qxp(i,j,kc,QRHO))
               qxp(i,j,kc,QU) = u + (apright - amright)*cc/rho
               qxp(i,j,kc,QV) = v + azv1rght
               qxp(i,j,kc,QW) = w + azw1rght
               qxp(i,j,kc,QPRES) = p + (apright + amright)*csq
               qxp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
            end if

            if (u-cc .ge. 0.d0) then
               spminus = (u-cc)*dtdx
            else
               spminus = 1.d0
            endif
            if (u+cc .ge. 0.d0) then
               spplus = (u+cc)*dtdx
            else
               spplus = 1.d0
            endif
            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif

            apleft = 0.5d0*(1.d0 - spplus )*alphap
            amleft = 0.5d0*(1.d0 - spminus)*alpham
            azrleft = 0.5d0*(1.d0 - spzero )*alpha0r
            azeleft = 0.5d0*(1.d0 - spzero )*alpha0e
            azv1left = 0.5d0*(1.d0 - spzero )*alpha0v
            azw1left = 0.5d0*(1.d0 - spzero )*alpha0w

            if (i .le. ihi1) then
               qxm(i+1,j,kc,QRHO) = rho + apleft + amleft + azrleft
               qxm(i+1,j,kc,QRHO) = max(small_dens, qxm(i+1,j,kc,QRHO))
               qxm(i+1,j,kc,QU) = u + (apleft - amleft)*cc/rho
               qxm(i+1,j,kc,QV) = v + azv1left
               qxm(i+1,j,kc,QW) = w + azw1left
               qxm(i+1,j,kc,QPRES) = p + (apleft + amleft)*csq
               qxm(i+1,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
            endif

         enddo
      enddo
      !$OMP END PARALLEL DO

      ! Treat K as a passively advected quantity
      if (QESGS .gt. -1) then
         n = QESGS
         do j = ilo2-1, ihi2+1
            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = u*dtdx
               endif
               acmprght = 0.5d0*(-1.d0 - spzero )*dqx(i,j,kc,n)
               qxp(i,j,kc,n) = q(i,j,k3d,n) + acmprght
            enddo
 
            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. 0.d0) then
                  spzero = u*dtdx
               else
                  spzero = 1.d0
               endif
               acmpleft = 0.5d0*(1.d0 - spzero )*dqx(i,j,kc,n)
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + acmpleft
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,i,j,u,spzero,acmprght,acmpleft) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = QFA + iadv - 1

         do j = ilo2-1, ihi2+1

            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = u*dtdx
               endif
               acmprght = 0.5d0*(-1.d0 - spzero )*dqx(i,j,kc,n)
               qxp(i,j,kc,n) = q(i,j,k3d,n) + acmprght
            enddo

            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. 0.d0) then
                  spzero = u*dtdx
               else
                  spzero = 1.d0
               endif
               acmpleft = 0.5d0*(1.d0 - spzero )*dqx(i,j,kc,n)
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + acmpleft
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,i,j,u,spzero,ascmprght,ascmpleft) IF(nspec.gt.1)
      do ispec = 1, nspec
         n = QFS + ispec - 1

         do j = ilo2-1, ihi2+1

            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = u*dtdx
               endif
               ascmprght = 0.5d0*(-1.d0 - spzero )*dqx(i,j,kc,n)
               qxp(i,j,kc,n) = q(i,j,k3d,n) + ascmprght
            enddo

            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. 0.d0) then
                  spzero = u*dtdx
               else
                  spzero = 1.d0
               endif
               ascmpleft = 0.5d0*(1.d0 - spzero )*dqx(i,j,kc,n)
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + ascmpleft
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,i,j,u,spzero,ascmprght,ascmpleft) IF(naux.gt.1)
      do iaux = 1, naux
         n = QFX + iaux - 1

         do j = ilo2-1, ihi2+1

            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = u*dtdx
               endif
               ascmprght = 0.5d0*(-1.d0 - spzero )*dqx(i,j,kc,n)
               qxp(i,j,kc,n) = q(i,j,k3d,n) + ascmprght
            enddo

            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. 0.d0) then
                  spzero = u*dtdx
               else
                  spzero = 1.d0
               endif
               ascmpleft = 0.5d0*(1.d0 - spzero )*dqx(i,j,kc,n)
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + ascmpleft
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth,drho,du,dv,dw,dp,drhoe,alpham,alphap,alpha0r) &
      !$OMP PRIVATE(alpha0e,alpha0v,alpha0w,spminus,spplus,spzero,apright,amright,azrright,azeright,azv1rght) &
      !$OMP PRIVATE(azw1rght,apleft,amleft,azrleft,azeleft,azv1left,azw1left) &
      !$OMP PRIVATE(alpha0u,azu1rght,azu1left)
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
            enth = ( (rhoe+p)/rho )/csq

            drho = dqy(i,j,kc,QRHO)
            du = dqy(i,j,kc,QU)
            dv = dqy(i,j,kc,QV)
            dw = dqy(i,j,kc,QW)
            dp = dqy(i,j,kc,QPRES)
            drhoe = dqy(i,j,kc,QREINT)

            alpham = 0.5d0*(dp/(rho*cc) - dv)*rho/cc
            alphap = 0.5d0*(dp/(rho*cc) + dv)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0w = dw

            if (v-cc .gt. 0.d0) then
               spminus = -1.d0
            else
               spminus = (v-cc)*dtdy
            endif
            if (v+cc .gt. 0.d0) then
               spplus = -1.d0
            else
               spplus = (v+cc)*dtdy
            endif
            if (v .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = v*dtdy
            endif

            apright = 0.5d0*(-1.d0 - spplus )*alphap
            amright = 0.5d0*(-1.d0 - spminus)*alpham
            azrright = 0.5d0*(-1.d0 - spzero )*alpha0r
            azeright = 0.5d0*(-1.d0 - spzero )*alpha0e
            azu1rght = 0.5d0*(-1.d0 - spzero )*alpha0u
            azw1rght = 0.5d0*(-1.d0 - spzero )*alpha0w

            if (j .ge. ilo2) then
               qyp(i,j,kc,QRHO) = rho + apright + amright + azrright
               qyp(i,j,kc,QRHO) = max(small_dens, qyp(i,j,kc,QRHO))
               qyp(i,j,kc,QV) = v + (apright - amright)*cc/rho
               qyp(i,j,kc,QU) = u + azu1rght
               qyp(i,j,kc,QW) = w + azw1rght
               qyp(i,j,kc,QPRES) = p + (apright + amright)*csq
               qyp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
            end if

            if (v-cc .ge. 0.d0) then
               spminus = (v-cc)*dtdy
            else
               spminus = 1.d0
            endif
            if (v+cc .ge. 0.d0) then
               spplus = (v+cc)*dtdy
            else
               spplus = 1.d0
            endif
            if (v .ge. 0.d0) then
               spzero = v*dtdy
            else
               spzero = 1.d0
            endif

            apleft = 0.5d0*(1.d0 - spplus )*alphap
            amleft = 0.5d0*(1.d0 - spminus)*alpham
            azrleft = 0.5d0*(1.d0 - spzero )*alpha0r
            azeleft = 0.5d0*(1.d0 - spzero )*alpha0e
            azu1left = 0.5d0*(1.d0 - spzero )*alpha0u
            azw1left = 0.5d0*(1.d0 - spzero )*alpha0w

            if (j .le. ihi2) then
               qym(i,j+1,kc,QRHO) = rho + apleft + amleft + azrleft
               qym(i,j+1,kc,QRHO) = max(small_dens, qym(i,j+1,kc,QRHO))
               qym(i,j+1,kc,QV) = v + (apleft - amleft)*cc/rho
               qym(i,j+1,kc,QU) = u + azu1left
               qym(i,j+1,kc,QW) = w + azw1left
               qym(i,j+1,kc,QPRES) = p + (apleft + amleft)*csq
               qym(i,j+1,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
            endif

         enddo
      enddo
      !$OMP END PARALLEL DO

      ! Treat K as a passively advected quantity
      if (QESGS .gt. -1) then
         n = QESGS
         do i = ilo1-1, ihi1+1
 
            ! Top state
            do j = ilo2, ihi2+1
               v = q(i,j,k3d,QV)
               if (v .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = v*dtdy
               endif
               acmptop = 0.5d0*(-1.d0 - spzero )*dqy(i,j,kc,n)
               qyp(i,j,kc,n) = q(i,j,k3d,n) + acmptop
            enddo
 
            ! Bottom state
            do j = ilo2-1, ihi2
               v = q(i,j,k3d,QV)
               if (v .ge. 0.d0) then
                  spzero = v*dtdy
               else
                  spzero = 1.d0
               endif
               acmpbot = 0.5d0*(1.d0 - spzero )*dqy(i,j,kc,n)
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + acmpbot
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,i,j,v,spzero,acmptop,acmpbot) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = QFA + iadv - 1

         do i = ilo1-1, ihi1+1

            ! Top state
            do j = ilo2, ihi2+1
               v = q(i,j,k3d,QV)
               if (v .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = v*dtdy
               endif
               acmptop = 0.5d0*(-1.d0 - spzero )*dqy(i,j,kc,n)
               qyp(i,j,kc,n) = q(i,j,k3d,n) + acmptop
            enddo

            ! Bottom state
            do j = ilo2-1, ihi2
               v = q(i,j,k3d,QV)
               if (v .ge. 0.d0) then
                  spzero = v*dtdy
               else
                  spzero = 1.d0
               endif
               acmpbot = 0.5d0*(1.d0 - spzero )*dqy(i,j,kc,n)
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + acmpbot
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,i,j,v,spzero,ascmptop,ascmpbot) IF(nspec.gt.1)
      do ispec = 1, nspec
         n = QFS + ispec - 1

         do i = ilo1-1, ihi1+1

            ! Top state
            do j = ilo2, ihi2+1
               v = q(i,j,k3d,QV)
               if (v .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = v*dtdy
               endif
               ascmptop = 0.5d0*(-1.d0 - spzero )*dqy(i,j,kc,n)
               qyp(i,j,kc,n) = q(i,j,k3d,n) + ascmptop
            enddo

            ! Bottom state
            do j = ilo2-1, ihi2
               v = q(i,j,k3d,QV)
               if (v .ge. 0.d0) then
                  spzero = v*dtdy
               else
                  spzero = 1.d0
               endif
               ascmpbot = 0.5d0*(1.d0 - spzero )*dqy(i,j,kc,n)
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + ascmpbot
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,i,j,v,spzero,ascmptop,ascmpbot) IF(naux.gt.1)
      do iaux = 1, naux
         n = QFX + iaux - 1

         do i = ilo1-1, ihi1+1

            ! Top state
            do j = ilo2, ihi2+1
               v = q(i,j,k3d,QV)
               if (v .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = v*dtdy
               endif
               ascmptop = 0.5d0*(-1.d0 - spzero )*dqy(i,j,kc,n)
               qyp(i,j,kc,n) = q(i,j,k3d,n) + ascmptop
            enddo

            ! Bottom state
            do j = ilo2-1, ihi2
               v = q(i,j,k3d,QV)
               if (v .ge. 0.d0) then
                  spzero = v*dtdy
               else
                  spzero = 1.d0
               endif
               ascmpbot = 0.5d0*(1.d0 - spzero )*dqy(i,j,kc,n)
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + ascmpbot
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

    end subroutine tracexy

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine tracez(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
           dqz,dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3, &
           qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
           ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QESGS, QPRES, QFA, QFS, QFX, nadv, small_dens, &
                                     ppm_type

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2
      integer km,kc,k3d

      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      double precision  dqz(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)
      double precision qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision dz, dt

      ! Local variables
      integer i, j
      integer n
      integer iadv, ispec, iaux

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

      dtdz = dt/dz
      
      !!!!!!!!!!!!!!!
      ! NON-PPM CODE
      !!!!!!!!!!!!!!!
      
      !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth,drho,du,dv,dw,dp,drhoe,alpham,alphap,alpha0r,alpha0e) &
      !$OMP PRIVATE(alpha0u,alpha0v,spminus,spplus,spzero,apright,amright,azrright,azeright,azu1rght,azv1rght,apleft) &
      !$OMP PRIVATE(amleft,azrleft,azeleft,azu1left,azv1left)
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
            enth = ( (rhoe+p)/rho )/csq

            drho = dqz(i,j,kc,QRHO)
            du = dqz(i,j,kc,QU)
            dv = dqz(i,j,kc,QV)
            dw = dqz(i,j,kc,QW)
            dp = dqz(i,j,kc,QPRES)
            drhoe = dqz(i,j,kc,QREINT)

            alpham = 0.5d0*(dp/(rho*cc) - dw)*rho/cc
            alphap = 0.5d0*(dp/(rho*cc) + dw)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0v = dv

            if (w-cc .gt. 0.d0) then
               spminus = -1.d0
            else
               spminus = (w-cc)*dtdz
            endif
            if (w+cc .gt. 0.d0) then
               spplus = -1.d0
            else
               spplus = (w+cc)*dtdz
            endif
            if (w .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = w*dtdz
            endif

            apright = 0.5d0*(-1.d0 - spplus )*alphap
            amright = 0.5d0*(-1.d0 - spminus)*alpham
            azrright = 0.5d0*(-1.d0 - spzero )*alpha0r
            azeright = 0.5d0*(-1.d0 - spzero )*alpha0e
            azu1rght = 0.5d0*(-1.d0 - spzero )*alpha0u
            azv1rght = 0.5d0*(-1.d0 - spzero )*alpha0v

            qzp(i,j,kc,QRHO) = rho + apright + amright + azrright
            qzp(i,j,kc,QRHO) = max(small_dens, qzp(i,j,kc,QRHO))
            qzp(i,j,kc,QW) = w + (apright - amright)*cc/rho
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
            enth = ( (rhoe+p)/rho )/csq

            drho = dqz(i,j,km,QRHO)
            du = dqz(i,j,km,QU)
            dv = dqz(i,j,km,QV)
            dw = dqz(i,j,km,QW)
            dp = dqz(i,j,km,QPRES)
            drhoe = dqz(i,j,km,QREINT)

            alpham = 0.5d0*(dp/(rho*cc) - dw)*rho/cc
            alphap = 0.5d0*(dp/(rho*cc) + dw)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0v = dv

            if (w-cc .ge. 0.d0) then
               spminus = (w-cc)*dtdz
            else
               spminus = 1.d0
            endif
            if (w+cc .ge. 0.d0) then
               spplus = (w+cc)*dtdz
            else
               spplus = 1.d0
            endif
            if (w .ge. 0.d0) then
               spzero = w*dtdz
            else
               spzero = 1.d0
            endif

            apleft = 0.5d0*(1.d0 - spplus )*alphap
            amleft = 0.5d0*(1.d0 - spminus)*alpham
            azrleft = 0.5d0*(1.d0 - spzero )*alpha0r
            azeleft = 0.5d0*(1.d0 - spzero )*alpha0e
            azu1left = 0.5d0*(1.d0 - spzero )*alpha0u
            azv1left = 0.5d0*(1.d0 - spzero )*alpha0v

            qzm(i,j,kc,QRHO) = rho + apleft + amleft + azrleft
            qzm(i,j,kc,QRHO) = max(small_dens, qzm(i,j,kc,QRHO))
            qzm(i,j,kc,QW) = w + (apleft - amleft)*cc/rho
            qzm(i,j,kc,QU) = u + azu1left
            qzm(i,j,kc,QV) = v + azv1left
            qzm(i,j,kc,QPRES) = p + (apleft + amleft)*csq
            qzm(i,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

         enddo
      enddo
      !$OMP END PARALLEL DO

      ! Treat K as a passively advected quantity
      if (QESGS .gt. -1) then
         n = QESGS
         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1
 
               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = w*dtdz
               endif
               acmptop = 0.5d0*(-1.d0 - spzero )*dqz(i,j,kc,n)
               qzp(i,j,kc,n) = q(i,j,k3d,n) + acmptop
 
               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. 0.d0) then
                  spzero = w*dtdz
               else
                  spzero = 1.d0
               endif
               acmpbot = 0.5d0*(1.d0 - spzero )*dqz(i,j,km,n)
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + acmpbot
 
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,i,j,w,acmptop,acmpbot,spzero) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = QFA + iadv - 1

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = w*dtdz
               endif
               acmptop = 0.5d0*(-1.d0 - spzero )*dqz(i,j,kc,n)
               qzp(i,j,kc,n) = q(i,j,k3d,n) + acmptop

               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. 0.d0) then
                  spzero = w*dtdz
               else
                  spzero = 1.d0
               endif
               acmpbot = 0.5d0*(1.d0 - spzero )*dqz(i,j,km,n)
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + acmpbot
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,i,j,w,ascmptop,ascmpbot,spzero) IF(nspec.gt.1)
      do ispec = 1, nspec
         n = QFS + ispec - 1

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = w*dtdz
               endif
               ascmptop = 0.5d0*(-1.d0 - spzero )*dqz(i,j,kc,n)
               qzp(i,j,kc,n) = q(i,j,k3d,n) + ascmptop

               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. 0.d0) then
                  spzero = w*dtdz
               else
                  spzero = 1.d0
               endif
               ascmpbot = 0.5d0*(1.d0 - spzero )*dqz(i,j,km,n)
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + ascmpbot
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,i,j,w,ascmptop,ascmpbot,spzero) IF(naux.gt.1)
      do iaux = 1, naux
         n = QFX + iaux - 1

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = w*dtdz
               endif
               ascmptop = 0.5d0*(-1.d0 - spzero )*dqz(i,j,kc,n)
               qzp(i,j,kc,n) = q(i,j,k3d,n) + ascmptop

               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. 0.d0) then
                  spzero = w*dtdz
               else
                  spzero = 1.d0
               endif
               ascmpbot = 0.5d0*(1.d0 - spzero )*dqz(i,j,km,n)
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + ascmpbot
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

    end subroutine tracez
