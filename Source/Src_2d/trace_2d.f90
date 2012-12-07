module trace_module

  implicit none

  private

  public trace

contains

      subroutine trace(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                       dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                       dq,qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       grav,gv_l1,gv_l2,gv_h1,gv_h2, &
                       ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

      use network, only : nspec, naux
      use meth_params_module, only : iorder, QVAR, QRHO, QU, QV, &
                                     QREINT, QPRES, QFA, QFS, QFX, &
                                     nadv, small_dens, ppm_type, use_pslope
      use slope_module, only : uslope, pslope

      implicit none

      integer ilo1,ilo2,ihi1,ihi2
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
      integer gv_l1,gv_l2,gv_h1,gv_h2

      double precision dx, dy, dt
      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
      double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
      double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
      double precision  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,2)

      ! Local variables
      integer i, j
      integer n, iadv, ispec, iaux

      double precision dtdx, dtdy
      double precision cc, csq, rho, u, v, p, rhoe
      double precision drho, du, dv, dp, drhoe

      double precision enth, alpham, alphap, alpha0r, alpha0e
      double precision alpha0u, alpha0v
      double precision spminus, spplus, spzero
      double precision apright, amright, azrright, azeright
      double precision azu1rght, azv1rght
      double precision apleft, amleft, azrleft, azeleft
      double precision azu1left, azv1left
      double precision acmprght, acmpleft, acmpbot, acmptop
      double precision ascmprght, ascmpleft, ascmpbot, ascmptop
      double precision sourcr,sourcp,source,courn,eta,dlogatmp

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in trace with ppm_type != 0'
        call bl_error("Error:: Castro_2d.f90 :: trace")
      end if

      dtdx = dt/dx
      dtdy = dt/dy

         ! Compute slopes
         if (iorder .eq. 1) then

            dq(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:QVAR) = 0.d0

         else
            
            call uslope(q, &
                        flatn, qd_l1, qd_l2, qd_h1, qd_h2, &
                        dq   ,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        ilo1,ilo2,ihi1,ihi2,QVAR,1)

            if (use_pslope .eq. 1) &
               call pslope(q(:,:,QPRES),q(:,:,QRHO),  &
                           flatn        , qd_l1, qd_l2, qd_h1, qd_h2, &
                           dq(:,:,QPRES),qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                           grav         , gv_l1, gv_l2, gv_h1, gv_h2, &
                           ilo1,ilo2,ihi1,ihi2,dx,dy,1)

         endif
      
         ! Compute left and right traced states
         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1
               
               cc = c(i,j)
               csq = cc**2
               rho = q(i,j,QRHO)
               u = q(i,j,QU)
               v = q(i,j,QV)
               p = q(i,j,QPRES)
               rhoe = q(i,j,QREINT)
               enth = ( (rhoe+p)/rho )/csq
               
               drho = dq(i,j,QRHO)
               du = dq(i,j,QU)
               dv = dq(i,j,QV)
               dp = dq(i,j,QPRES)
               drhoe = dq(i,j,QREINT)
               
               alpham = 0.5d0*(dp/(rho*cc) - du)*rho/cc
               alphap = 0.5d0*(dp/(rho*cc) + du)*rho/cc
               alpha0r = drho - dp/csq
               alpha0e = drhoe - dp*enth
               alpha0v = dv
               
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

               if (i .ge. ilo1) then
                  qxp(i,j,QRHO) = rho + apright + amright + azrright
                  qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))
                  qxp(i,j,QU) = u + (apright - amright)*cc/rho
                  qxp(i,j,QV) = v + azv1rght
                  qxp(i,j,QPRES) = p + (apright + amright)*csq
                  qxp(i,j,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
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

               if (i .le. ihi1) then
                  qxm(i+1,j,QRHO) = rho + apleft + amleft + azrleft
                  qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
                  qxm(i+1,j,QU) = u + (apleft - amleft)*cc/rho
                  qxm(i+1,j,QV) = v + azv1left
                  qxm(i+1,j,QPRES) = p + (apleft + amleft)*csq
                  qxm(i+1,j,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
               end if

               if(dloga(i,j).ne.0)then
                  courn = dtdx*(cc+abs(u))
                  eta = (1.d0-courn)/(cc*dt*abs(dloga(i,j)))
                  dlogatmp = min(eta,1.d0)*dloga(i,j)
                  sourcr = -0.5d0*dt*rho*dlogatmp*u
                  sourcp = sourcr*csq
                  source = sourcp*enth
                  if (i .le. ihi1) then
                     qxm(i+1,j,QRHO) = qxm(i+1,j,QRHO) + sourcr
                     qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
                     qxm(i+1,j,QPRES) = qxm(i+1,j,QPRES) + sourcp
                     qxm(i+1,j,QREINT) = qxm(i+1,j,QREINT) + source
                  end if
                  if (i .ge. ilo1) then
                     qxp(i,j,QRHO) = qxp(i,j,QRHO) + sourcr
                     qxp(i,j,QRHO) = max(qxp(i,j,QRHO),small_dens)
                     qxp(i,j,QPRES) = qxp(i,j,QPRES) + sourcp
                     qxp(i,j,QREINT) = qxp(i,j,QREINT) + source
                  end if
               endif

            enddo
         enddo

         do iadv = 1, nadv
            n = QFA + iadv - 1
            do j = ilo2-1, ihi2+1

               ! Right state
               do i = ilo1, ihi1+1
                  u = q(i,j,QU)
                  if (u .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = u*dtdx
                  endif
                  acmprght = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qxp(i,j,n) = q(i,j,n) + acmprght
               enddo

               ! Left state
               do i = ilo1-1, ihi1
                  u = q(i,j,QU)
                  if (u .ge. 0.d0) then
                     spzero = u*dtdx
                  else
                     spzero = 1.d0
                  endif
                  acmpleft = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qxm(i+1,j,n) = q(i,j,n) + acmpleft
               enddo

            enddo
         enddo

         do ispec = 1, nspec
            n = QFS + ispec - 1
            do j = ilo2-1, ihi2+1

               ! Right state
               do i = ilo1, ihi1+1
                  u = q(i,j,QU)
                  if (u .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = u*dtdx
                  endif
                  ascmprght = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qxp(i,j,n) = q(i,j,n) + ascmprght
               enddo

               ! Left state
               do i = ilo1-1, ihi1
                  u = q(i,j,QU)
                  if (u .ge. 0.d0) then
                     spzero = u*dtdx
                  else
                     spzero = 1.d0
                  endif
                  ascmpleft = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qxm(i+1,j,n) = q(i,j,n) + ascmpleft
               enddo

            enddo
         enddo

         do iaux = 1, naux
            n = QFX + iaux - 1
            do j = ilo2-1, ihi2+1

               ! Right state
               do i = ilo1, ihi1+1
                  u = q(i,j,QU)
                  if (u .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = u*dtdx
                  endif
                  ascmprght = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qxp(i,j,n) = q(i,j,n) + ascmprght
               enddo

               ! Left state
               do i = ilo1-1, ihi1
                  u = q(i,j,QU)
                  if (u .ge. 0.d0) then
                     spzero = u*dtdx
                  else
                     spzero = 1.d0
                  endif
                  ascmpleft = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qxm(i+1,j,n) = q(i,j,n) + ascmpleft
               enddo

            enddo
         enddo

         ! Compute slopes
         if (iorder .ne. 1) then

            call uslope(q,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                        dq,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        ilo1,ilo2,ihi1,ihi2,QVAR,2)

            call pslope(q(:,:,QPRES),q(:,:,QRHO),  &
                        flatn        , qd_l1, qd_l2, qd_h1, qd_h2, &
                        dq(:,:,QPRES),qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        grav         , gv_l1, gv_l2, gv_h1, gv_h2, &
                        ilo1,ilo2,ihi1,ihi2,dx,dy,2)

         endif

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               cc = c(i,j)
               csq = cc**2
               rho = q(i,j,QRHO)
               u = q(i,j,QU)
               v = q(i,j,QV)
               p = q(i,j,QPRES)
               rhoe = q(i,j,QREINT)
               enth = ( (rhoe+p)/rho )/csq

               drho = dq(i,j,QRHO)
               du = dq(i,j,QU)
               dv = dq(i,j,QV)
               dp = dq(i,j,QPRES)
               drhoe = dq(i,j,QREINT)

               alpham = 0.5d0*(dp/(rho*cc) - dv)*rho/cc
               alphap = 0.5d0*(dp/(rho*cc) + dv)*rho/cc
               alpha0r = drho - dp/csq
               alpha0e = drhoe - dp*enth
               alpha0u = du

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

               if (j .ge. ilo2) then
                  qyp(i,j,QRHO) = rho + apright + amright + azrright
                  qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))
                  qyp(i,j,QV) = v + (apright - amright)*cc/rho
                  qyp(i,j,QU) = u + azu1rght
                  qyp(i,j,QPRES) = p + (apright + amright)*csq
                  qyp(i,j,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
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

               if (j .le. ihi2) then
                  qym(i,j+1,QRHO) = rho + apleft + amleft + azrleft
                  qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))
                  qym(i,j+1,QV) = v + (apleft - amleft)*cc/rho
                  qym(i,j+1,QU) = u + azu1left
                  qym(i,j+1,QPRES) = p + (apleft + amleft)*csq
                  qym(i,j+1,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
               end if

            enddo
         enddo

         do iadv = 1, nadv
            n = QFA + iadv - 1
            do i = ilo1-1, ihi1+1

               ! Top state
               do j = ilo2, ihi2+1
                  v = q(i,j,QV)
                  if (v .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = v*dtdy
                  endif
                  acmptop = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qyp(i,j,n) = q(i,j,n) + acmptop
               enddo

               ! Bottom state
               do j = ilo2-1, ihi2
                  v = q(i,j,QV)
                  if (v .ge. 0.d0) then
                     spzero = v*dtdy
                  else
                     spzero = 1.d0
                  endif
                  acmpbot = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qym(i,j+1,n) = q(i,j,n) + acmpbot
               enddo

            enddo
         enddo

         do ispec = 1, nspec
            n = QFS + ispec - 1
            do i = ilo1-1, ihi1+1

               ! Top state
               do j = ilo2, ihi2+1
                  v = q(i,j,QV)
                  if (v .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = v*dtdy
                  endif
                  ascmptop = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qyp(i,j,n) = q(i,j,n) + ascmptop
               enddo

               ! Bottom state
               do j = ilo2-1, ihi2
                  v = q(i,j,QV)
                  if (v .ge. 0.d0) then
                     spzero = v*dtdy
                  else
                     spzero = 1.d0
                  endif
                  ascmpbot = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qym(i,j+1,n) = q(i,j,n) + ascmpbot
               enddo

            enddo
         enddo

         do iaux = 1, naux
            n = QFX + iaux - 1
            do i = ilo1-1, ihi1+1

               ! Top state
               do j = ilo2, ihi2+1
                  v = q(i,j,QV)
                  if (v .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = v*dtdy
                  endif
                  ascmptop = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qyp(i,j,n) = q(i,j,n) + ascmptop
               enddo

               ! Bottom state
               do j = ilo2-1, ihi2
                  v = q(i,j,QV)
                  if (v .ge. 0.d0) then
                     spzero = v*dtdy
                  else
                     spzero = 1.d0
                  endif
                  ascmpbot = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qym(i,j+1,n) = q(i,j,n) + ascmpbot
               enddo

            enddo
         enddo

      end subroutine trace

end module trace_module
