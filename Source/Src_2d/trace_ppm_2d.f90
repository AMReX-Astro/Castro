
      subroutine trace_ppm(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
           dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
           qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
           ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

        use network, only : nspec, naux
        use meth_params_module, only : iorder, QVAR, QRHO, QU, QV, &
             QREINT, QPRES, QFA, QFS, QFX, &
             nadv, small_dens, ppm_type

        implicit none

        integer ilo1,ilo2,ihi1,ihi2
        integer qd_l1,qd_l2,qd_h1,qd_h2
        integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
        integer qpd_l1,qpd_l2,qpd_h1,qpd_h2

        double precision dx, dy, dt
        double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
        double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
        double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
        double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
        double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
        double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
        double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
        double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)

        ! Local variables
        integer i, j
        integer n, iadv
        integer ns, ispec, iaux

        double precision dtdx, dtdy
        double precision cc, csq, rho, u, v, p, rhoe
        double precision drho, du, dv, dp, drhoe
        double precision drhop, dup, dvp, dpp, drhoep
        double precision drhom, dum, dvm, dpm, drhoem

        double precision enth, alpham, alphap, alpha0r, alpha0e
        double precision alpha0u, alpha0v
        double precision apright, amright, azrright, azeright
        double precision azu1rght, azv1rght
        double precision apleft, amleft, azrleft, azeleft
        double precision azu1left, azv1left
        double precision sourcr,sourcp,source,courn,eta,dlogatmp

        double precision, allocatable :: Ip(:,:,:,:,:)
        double precision, allocatable :: Im(:,:,:,:,:)

        dtdx = dt/dx
        dtdy = dt/dy

        if (ppm_type .eq. 0) then
           print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
           call bl_error("Error:: ppm_2d.f90 :: trace_ppm")
        end if

        allocate(Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QVAR))
        allocate(Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QVAR))

        ! Compute Ip and Im
        do n=1,QVAR
           call ppm(q(:,:,n),qd_l1,qd_l2,qd_h1,qd_h2,q(:,:,QU:),c, &
                Ip(:,:,:,:,n),Im(:,:,:,:,n),ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
        end do

        ! Trace to left and right edges using upwind PPM
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

              ! plus state on face i
              drhom  = flatn(i,j)*(rho  - Im(i,j,1,1,QRHO))
              dum    = flatn(i,j)*(u    - Im(i,j,1,1,QU))
              dvm    = flatn(i,j)*(v    - Im(i,j,1,1,QV))
              dpm    = flatn(i,j)*(p    - Im(i,j,1,1,QPRES))
              drhoem = flatn(i,j)*(rhoe - Im(i,j,1,1,QREINT))

              drho  = flatn(i,j)*(rho  - Im(i,j,1,2,QRHO))
              du    = flatn(i,j)*(u    - Im(i,j,1,2,QU))
              dv    = flatn(i,j)*(v    - Im(i,j,1,2,QV))
              dp    = flatn(i,j)*(p    - Im(i,j,1,2,QPRES))
              drhoe = flatn(i,j)*(rhoe - Im(i,j,1,2,QREINT))

              drhop  = flatn(i,j)*(rho  - Im(i,j,1,3,QRHO))
              dup    = flatn(i,j)*(u    - Im(i,j,1,3,QU))
              dvp    = flatn(i,j)*(v    - Im(i,j,1,3,QV))
              dpp    = flatn(i,j)*(p    - Im(i,j,1,3,QPRES))
              drhoep = flatn(i,j)*(rhoe - Im(i,j,1,3,QREINT))

              alpham = 0.5d0*(dpm/(rho*cc) - dum)*rho/cc
              alphap = 0.5d0*(dpp/(rho*cc) + dup)*rho/cc
              alpha0r = drho - dp/csq
              alpha0e = drhoe - dp*enth
              alpha0v = dv

              if (u-cc .gt. 0.d0) then
                 amright = 0.d0
              else if (u-cc .lt. 0.d0) then
                 amright = -alpham
              else
                 amright = -0.5d0*alpham
              endif
              if (u+cc .gt. 0.d0) then
                 apright = 0.d0
              else if (u+cc .lt. 0.d0) then
                 apright = -alphap
              else
                 apright = -0.5d0*alphap
              endif
              if (u .gt. 0.d0) then
                 azrright = 0.d0
                 azeright = 0.d0
                 azv1rght = 0.d0
              else if (u .lt. 0.d0) then
                 azrright = -alpha0r
                 azeright = -alpha0e
                 azv1rght = -alpha0v
              else
                 azrright = -0.5d0*alpha0r
                 azeright = -0.5d0*alpha0e
                 azv1rght = -0.5d0*alpha0v
              endif

              if (i .ge. ilo1) then
                 qxp(i,j,QRHO) = rho + apright + amright + azrright
                 qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))
                 qxp(i,j,QU) = u + (apright - amright)*cc/rho
                 qxp(i,j,QV) = v + azv1rght
                 qxp(i,j,QPRES) = p + (apright + amright)*csq
                 qxp(i,j,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
              end if

              ! minus state on face i+1
              drhom  = flatn(i,j)*(rho  - Ip(i,j,1,1,QRHO))
              dum    = flatn(i,j)*(u    - Ip(i,j,1,1,QU))
              dvm    = flatn(i,j)*(v    - Ip(i,j,1,1,QV))
              dpm    = flatn(i,j)*(p    - Ip(i,j,1,1,QPRES))
              drhoem = flatn(i,j)*(rhoe - Ip(i,j,1,1,QREINT))

              drho  = flatn(i,j)*(rho  - Ip(i,j,1,2,QRHO))
              du    = flatn(i,j)*(u    - Ip(i,j,1,2,QU))
              dv    = flatn(i,j)*(v    - Ip(i,j,1,2,QV))
              dp    = flatn(i,j)*(p    - Ip(i,j,1,2,QPRES))
              drhoe = flatn(i,j)*(rhoe - Ip(i,j,1,2,QREINT))

              drhop  = flatn(i,j)*(rho  - Ip(i,j,1,3,QRHO))
              dup    = flatn(i,j)*(u    - Ip(i,j,1,3,QU))
              dvp    = flatn(i,j)*(v    - Ip(i,j,1,3,QV))
              dpp    = flatn(i,j)*(p    - Ip(i,j,1,3,QPRES))
              drhoep = flatn(i,j)*(rhoe - Ip(i,j,1,3,QREINT))

              alpham = 0.5d0*(dpm/(rho*cc) - dum)*rho/cc
              alphap = 0.5d0*(dpp/(rho*cc) + dup)*rho/cc
              alpha0r = drho - dp/csq
              alpha0e = drhoe - dp*enth
              alpha0v = dv

              if (u-cc .gt. 0.d0) then
                 amleft = -alpham
              else if (u-cc .lt. 0.d0) then
                 amleft = 0.d0
              else
                 amleft = -0.5d0*alpham
              endif
              if (u+cc .gt. 0.d0) then
                 apleft = -alphap
              else if (u+cc .lt. 0.d0) then
                 apleft = 0.d0
              else
                 apleft = -0.5d0*alphap
              endif
              if (u .gt. 0.d0) then
                 azrleft = -alpha0r
                 azeleft = -alpha0e
                 azv1left = -alpha0v
              else if (u .lt. 0.d0) then
                 azrleft = 0.d0
                 azeleft = 0.d0
                 azv1left = 0.d0
              else
                 azrleft = -0.5d0*alpha0r
                 azeleft = -0.5d0*alpha0e
                 azv1left = -0.5d0*alpha0v
              endif

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

           end do
        end do

        ! Now do the passively advected quantities
        do iadv = 1, nadv
           n = QFA + iadv - 1
           do j = ilo2-1, ihi2+1

              ! plus state on face i
              do i = ilo1, ihi1+1
                 u = q(i,j,QU)
                 if (u .gt. 0.d0) then
                    qxp(i,j,n) = q(i,j,n)
                 else if (u .lt. 0.d0) then
                    qxp(i,j,n) = q(i,j,n) + flatn(i,j)*(Im(i,j,1,2,n) - q(i,j,n))
                 else
                    qxp(i,j,n) = q(i,j,n) + 0.5d0*flatn(i,j)*(Im(i,j,1,2,n) - q(i,j,n))
                 endif
              enddo

              ! minus state on face i+1
              do i = ilo1-1, ihi1
                 u = q(i,j,QU)
                 if (u .gt. 0.d0) then
                    qxm(i+1,j,n) = q(i,j,n) + flatn(i,j)*(Ip(i,j,1,2,n) - q(i,j,n))
                 else if (u .lt. 0.d0) then
                    qxm(i+1,j,n) = q(i,j,n)
                 else
                    qxm(i+1,j,n) = q(i,j,n) + 0.5d0*flatn(i,j)*(Ip(i,j,1,2,n) - q(i,j,n))
                 endif
              enddo

           enddo
        enddo

        do ispec = 1, nspec
           ns = QFS + ispec - 1
           do j = ilo2-1, ihi2+1

              ! plus state on face i
              do i = ilo1, ihi1+1
                 u = q(i,j,QU)
                 if (u .gt. 0.d0) then
                    qxp(i,j,ns) = q(i,j,ns)
                 else if (u .lt. 0.d0) then
                    qxp(i,j,ns) = q(i,j,ns) + flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
                 else
                    qxp(i,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
                 endif
              enddo

              ! minus state on face i+1
              do i = ilo1-1, ihi1
                 u = q(i,j,QU)
                 if (u .gt. 0.d0) then
                    qxm(i+1,j,ns) = q(i,j,ns) + flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
                 else if (u .lt. 0.d0) then
                    qxm(i+1,j,ns) = q(i,j,ns)
                 else
                    qxm(i+1,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
                 endif
              enddo

           enddo
        enddo

        do iaux = 1, naux
           ns = QFX + iaux - 1
           do j = ilo2-1, ihi2+1

              ! plus state on face i
              do i = ilo1, ihi1+1
                 u = q(i,j,QU)
                 if (u .gt. 0.d0) then
                    qxp(i,j,ns) = q(i,j,ns)
                 else if (u .lt. 0.d0) then
                    qxp(i,j,ns) = q(i,j,ns) + flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
                 else
                    qxp(i,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
                 endif
              enddo

              ! minus state on face i+1
              do i = ilo1-1, ihi1
                 u = q(i,j,QU)
                 if (u .gt. 0.d0) then
                    qxm(i+1,j,ns) = q(i,j,ns) + flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
                 else if (u .lt. 0.d0) then
                    qxm(i+1,j,ns) = q(i,j,ns)
                 else
                    qxm(i+1,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
                 endif
              enddo

           enddo
        enddo

        ! Trace to bottom and top edges using upwind PPM
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

              ! plus state on face j
              drhom  = flatn(i,j)*(rho  - Im(i,j,2,1,QRHO))
              dum    = flatn(i,j)*(u    - Im(i,j,2,1,QU))
              dvm    = flatn(i,j)*(v    - Im(i,j,2,1,QV))
              dpm    = flatn(i,j)*(p    - Im(i,j,2,1,QPRES))
              drhoem = flatn(i,j)*(rhoe - Im(i,j,2,1,QREINT))

              drho  = flatn(i,j)*(rho  - Im(i,j,2,2,QRHO))
              du    = flatn(i,j)*(u    - Im(i,j,2,2,QU))
              dv    = flatn(i,j)*(v    - Im(i,j,2,2,QV))
              dp    = flatn(i,j)*(p    - Im(i,j,2,2,QPRES))
              drhoe = flatn(i,j)*(rhoe - Im(i,j,2,2,QREINT))

              drhop  = flatn(i,j)*(rho  - Im(i,j,2,3,QRHO))
              dup    = flatn(i,j)*(u    - Im(i,j,2,3,QU))
              dvp    = flatn(i,j)*(v    - Im(i,j,2,3,QV))
              dpp    = flatn(i,j)*(p    - Im(i,j,2,3,QPRES))
              drhoep = flatn(i,j)*(rhoe - Im(i,j,2,3,QREINT))

              alpham = 0.5d0*(dpm/(rho*cc) - dvm)*rho/cc
              alphap = 0.5d0*(dpp/(rho*cc) + dvp)*rho/cc
              alpha0r = drho - dp/csq
              alpha0e = drhoe - dp*enth
              alpha0u = du

              if (v-cc .gt. 0.d0) then
                 amright = 0.d0
              else if (v-cc .lt. 0.d0) then
                 amright = -alpham
              else
                 amright = -0.5d0*alpham
              endif
              if (v+cc .gt. 0.d0) then
                 apright = 0.d0
              else if (v+cc .lt. 0.d0) then
                 apright = -alphap
              else
                 apright = -0.5d0*alphap
              endif
              if (v .gt. 0.d0) then
                 azrright = 0.d0
                 azeright = 0.d0
                 azu1rght = 0.d0
              else if (v .lt. 0.d0) then
                 azrright = -alpha0r
                 azeright = -alpha0e
                 azu1rght = -alpha0u
              else
                 azrright = -0.5d0*alpha0r
                 azeright = -0.5d0*alpha0e
                 azu1rght = -0.5d0*alpha0u
              endif

              if (j .ge. ilo2) then
                 qyp(i,j,QRHO) = rho + apright + amright + azrright
                 qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))
                 qyp(i,j,QV) = v + (apright - amright)*cc/rho
                 qyp(i,j,QU) = u + azu1rght
                 qyp(i,j,QPRES) = p + (apright + amright)*csq
                 qyp(i,j,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
              end if

              ! minus state on face j+1
              drhom  = flatn(i,j)*(rho  - Ip(i,j,2,1,QRHO))
              dum    = flatn(i,j)*(u    - Ip(i,j,2,1,QU))
              dvm    = flatn(i,j)*(v    - Ip(i,j,2,1,QV))
              dpm    = flatn(i,j)*(p    - Ip(i,j,2,1,QPRES))
              drhoem = flatn(i,j)*(rhoe - Ip(i,j,2,1,QREINT))

              drho  = flatn(i,j)*(rho  - Ip(i,j,2,2,QRHO))
              du    = flatn(i,j)*(u    - Ip(i,j,2,2,QU))
              dv    = flatn(i,j)*(v    - Ip(i,j,2,2,QV))
              dp    = flatn(i,j)*(p    - Ip(i,j,2,2,QPRES))
              drhoe = flatn(i,j)*(rhoe - Ip(i,j,2,2,QREINT))

              drhop  = flatn(i,j)*(rho  - Ip(i,j,2,3,QRHO))
              dup    = flatn(i,j)*(u    - Ip(i,j,2,3,QU))
              dvp    = flatn(i,j)*(v    - Ip(i,j,2,3,QV))
              dpp    = flatn(i,j)*(p    - Ip(i,j,2,3,QPRES))
              drhoep = flatn(i,j)*(rhoe - Ip(i,j,2,3,QREINT))

              alpham = 0.5d0*(dpm/(rho*cc) - dvm)*rho/cc
              alphap = 0.5d0*(dpp/(rho*cc) + dvp)*rho/cc
              alpha0r = drho - dp/csq
              alpha0e = drhoe - dp*enth
              alpha0u = du

              if (v-cc .gt. 0.d0) then
                 amleft = -alpham
              else if (v-cc .lt. 0.d0) then
                 amleft = 0.d0
              else
                 amleft = -0.5d0*alpham
              endif
              if (v+cc .gt. 0.d0) then
                 apleft = -alphap
              else if (v+cc .lt. 0.d0) then
                 apleft = 0.d0
              else
                 apleft = -0.5d0*alphap
              endif
              if (v .gt. 0.d0) then
                 azrleft = -alpha0r
                 azeleft = -alpha0e
                 azu1left = -alpha0u
              else if (v .lt. 0.d0) then
                 azrleft = 0.d0
                 azeleft = 0.d0
                 azu1left = 0.d0
              else
                 azrleft = -0.5d0*alpha0r
                 azeleft = -0.5d0*alpha0e
                 azu1left = -0.5d0*alpha0u
              endif

              if (j .le. ihi2) then
                 qym(i,j+1,QRHO) = rho + apleft + amleft + azrleft
                 qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))
                 qym(i,j+1,QV) = v + (apleft - amleft)*cc/rho
                 qym(i,j+1,QU) = u + azu1left
                 qym(i,j+1,QPRES) = p + (apleft + amleft)*csq
                 qym(i,j+1,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
              end if

           end do
        end do

        ! Now do the passively advected quantities
        do iadv = 1, nadv
           n = QFA + iadv - 1
           do i = ilo1-1, ihi1+1

              ! plus state on face j
              do j = ilo2, ihi2+1
                 v = q(i,j,QV)
                 if (v .gt. 0.d0) then
                    qyp(i,j,n) = q(i,j,n)
                 else if (v .lt. 0.d0) then
                    qyp(i,j,n) = q(i,j,n) + flatn(i,j)*(Im(i,j,2,2,n) - q(i,j,n))
                 else
                    qyp(i,j,n) = q(i,j,n) + 0.5d0*flatn(i,j)*(Im(i,j,2,2,n) - q(i,j,n))
                 endif
              enddo

              ! minus state on face j+1
              do j = ilo2-1, ihi2
                 v = q(i,j,QV)
                 if (v .gt. 0.d0) then
                    qym(i,j+1,n) = q(i,j,n) + flatn(i,j)*(Ip(i,j,2,2,n) - q(i,j,n))
                 else if (v .lt. 0.d0) then
                    qym(i,j+1,n) = q(i,j,n)
                 else
                    qym(i,j+1,n) = q(i,j,n) + 0.5d0*flatn(i,j)*(Ip(i,j,2,2,n) - q(i,j,n))
                 endif
              enddo

           enddo
        enddo

        do ispec = 1, nspec
           ns = QFS + ispec - 1
           do i = ilo1-1, ihi1+1

              ! plus state on face j
              do j = ilo2, ihi2+1
                 v = q(i,j,QV)
                 if (v .gt. 0.d0) then
                    qyp(i,j,ns) = q(i,j,ns)
                 else if (v .lt. 0.d0) then
                    qyp(i,j,ns) = q(i,j,ns) + flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
                 else
                    qyp(i,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
                 endif
              enddo

              ! minus state on face j+1
              do j = ilo2-1, ihi2
                 v = q(i,j,QV)
                 if (v .gt. 0.d0) then
                    qym(i,j+1,ns) = q(i,j,ns) + flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
                 else if (v .lt. 0.d0) then
                    qym(i,j+1,ns) = q(i,j,ns)
                 else
                    qym(i,j+1,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
                 endif
              enddo

           enddo
        enddo

        do iaux = 1, naux
           ns = QFX + iaux - 1
           do i = ilo1-1, ihi1+1

              ! plus state on face j
              do j = ilo2, ihi2+1
                 v = q(i,j,QV)
                 if (v .gt. 0.d0) then
                    qyp(i,j,ns) = q(i,j,ns)
                 else if (v .lt. 0.d0) then
                    qyp(i,j,ns) = q(i,j,ns) + flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
                 else
                    qyp(i,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
                 endif
              enddo

              ! minus state on face j+1
              do j = ilo2-1, ihi2
                 v = q(i,j,QV)
                 if (v .gt. 0.d0) then
                    qym(i,j+1,ns) = q(i,j,ns) + flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
                 else if (v .lt. 0.d0) then
                    qym(i,j+1,ns) = q(i,j,ns)
                 else
                    qym(i,j+1,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
                 endif
              enddo

           enddo
        enddo

        deallocate(Ip,Im)

      end subroutine trace_ppm
