module trace_ppm_rad_module

  implicit none

contains

  subroutine trace_ppm_rad(lam, lam_l1, lam_l2, lam_h1, lam_h2, &   
                           q,c,cg,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                           dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                           qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                           ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

    use network, only : nspec, naux
    use meth_params_module, only : QRHO, QU, QV, &
         QREINT, QPRES, QFA, QFS, QFX, &
         small_dens, ppm_type, &
         npassive, qpass_map
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot
    use rad_params_module, only : ngroups
    use ppm_module, only : ppm
    
    implicit none

    integer ilo1,ilo2,ihi1,ihi2
    integer lam_l1,lam_l2,lam_h1,lam_h2
    integer qd_l1,qd_l2,qd_h1,qd_h2
    integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
    integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
    
    double precision dx, dy, dt
    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision    cg(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
    double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
    
    ! Local variables
    integer i, j, g
    integer n, ns, ipassive
    
    double precision dtdx, dtdy
    
    double precision, dimension(0:ngroups-1) :: er, der, alphar, sourcer, qrtmp,hr
    double precision, dimension(0:ngroups-1) :: lam0, lamp, lamm
    double precision cc, csq, rho, u, v, p, ptot, rhoe, enth, cgassq
    double precision dum, dvm, dptotm
    double precision drho, du, dv, drhoe, dptot
    double precision dup, dvp, dptotp
    
    double precision alpham, alphap, alpha0, alphae, alphau, alphav
    double precision sourcr,sourcp,source,courn,eta,dlogatmp
    
    double precision rhoe_g, h_g, alphae_g, drhoe_g
    
    double precision, allocatable :: Ip(:,:,:,:,:)
    double precision, allocatable :: Im(:,:,:,:,:)
    
    double precision :: er_foo
    
    dtdx = dt/dx
    dtdy = dt/dy
    
    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
       call bl_error("Error:: ppm_2d.f90 :: trace_ppm")
    end if
    
    allocate(Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QRADVAR))
    allocate(Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QRADVAR))
    
    ! Compute Ip and Im
    do n=1,QRADVAR
       call ppm(q(:,:,n),qd_l1,qd_l2,qd_h1,qd_h2, &
                q(:,:,QU:), c, qd_l1,qd_l2,qd_h1,qd_h2,&
                flatn, &
                Ip(:,:,:,:,n),Im(:,:,:,:,n), &
                ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    end do
    
    ! Trace to left and right edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1
          
          do g=0, ngroups-1
             lam0(g) = lam(i,j,g)
             lamp(g) = lam(i,j,g)
             lamm(g) = lam(i,j,g)
          end do
          
          cgassq = cg(i,j)**2
          cc = c(i,j)
          csq = cc**2
          
          rho = q(i,j,QRHO)
          u = q(i,j,QU)
          v = q(i,j,QV)
          p = q(i,j,QPRES)
          rhoe_g = q(i,j,QREINT)
          h_g = (p+rhoe_g) / rho
          er(:) = q(i,j,qrad:qradhi)
          hr(:) = (lam0+1.d0)*er/rho
          ptot = q(i,j,qptot)
          rhoe = q(i,j,qreitot)
          enth = ( (rhoe+ptot)/rho )/csq
          
          ! plus state on face i
          dum    = flatn(i,j)*(u    - Im(i,j,1,1,QU))
          dptotm = flatn(i,j)*(ptot - Im(i,j,1,1,qptot))
          
          drho    = flatn(i,j)*(rho    - Im(i,j,1,2,QRHO))
          dv      = flatn(i,j)*(v      - Im(i,j,1,2,QV))
          dptot   = flatn(i,j)*(ptot   - Im(i,j,1,2,qptot))
          drhoe   = flatn(i,j)*(rhoe   - Im(i,j,1,2,qreitot))
          drhoe_g = flatn(i,j)*(rhoe_g - Im(i,j,1,2,QREINT))
          der(:)  = flatn(i,j)*(er(:)  - Im(i,j,1,2,qrad:qradhi))
          
          dup    = flatn(i,j)*(u    - Im(i,j,1,3,QU))
          dptotp = flatn(i,j)*(ptot - Im(i,j,1,3,qptot))
          
          alpham = 0.5d0*(dptotm/(rho*cc) - dum)*rho/cc
          alphap = 0.5d0*(dptotp/(rho*cc) + dup)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphav = dv
          alphar(:) = der(:) - dptot/csq*hr
          
          if (u-cc .gt. 0.d0) then
             alpham = 0.d0
          else if (u-cc .lt. 0.d0) then
             alpham = -alpham
          else
             alpham = -0.5d0*alpham
          endif
          if (u+cc .gt. 0.d0) then
             alphap = 0.d0
          else if (u+cc .lt. 0.d0) then
             alphap = -alphap
          else
             alphap = -0.5d0*alphap
          endif
          if (u .gt. 0.d0) then
             alpha0 = 0.d0
             alphav = 0.d0
             alphae = 0.d0 
             alphae_g = 0.d0 
             alphar(:) = 0.d0
          else if (u .lt. 0.d0) then
             alpha0 = -alpha0
             alphav = -alphav
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else
             alpha0 = -0.5d0*alpha0
             alphav = -0.5d0*alphav
             alphae = -0.5d0*alphae
             alphae_g = -0.5d0*alphae_g
             alphar(:) = -0.5d0*alphar(:)
          endif
          
          if (i .ge. ilo1) then
             qxp(i,j,QRHO) = rho + alphap + alpham + alpha0
             qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))
             qxp(i,j,QU) = u + (alphap - alpham)*cc/rho
             qxp(i,j,QV) = v + alphav 
             qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
             qxp(i,j,qrad:qradhi) = qrtmp
             qxp(i,j,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
             qxp(i,j,QPRES) = p + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))
             qxp(i,j,qptot) = ptot + (alphap + alpham)*csq 
             qxp(i,j,qreitot) = qxp(i,j,QREINT) + sum(qrtmp)
             
             do g=0, ngroups-1
                if (qxp(i,j,qrad+g) < 0.d0) then
                   er_foo = - qxp(i,j,qrad+g)
                   qxp(i,j,qrad+g) = 0.d0
                   qxp(i,j,qptot) = qxp(i,j,qptot) + lamp(g) * er_foo
                   qxp(i,j,qreitot) = qxp(i,j,qreitot) + er_foo
                end if
             end do
             
             if (qxp(i,j,QPRES) < 0.d0) then
                qxp(i,j,QPRES) = p
             end if
          end if
          
          ! minus state on face i+1
          dum    = flatn(i,j)*(u    - Ip(i,j,1,1,QU))
          dptotm = flatn(i,j)*(ptot - Ip(i,j,1,1,qptot))
          
          drho    = flatn(i,j)*(rho    - Ip(i,j,1,2,QRHO))
          dv      = flatn(i,j)*(v      - Ip(i,j,1,2,QV))
          dptot   = flatn(i,j)*(ptot   - Ip(i,j,1,2,qptot))
          drhoe   = flatn(i,j)*(rhoe   - Ip(i,j,1,2,qreitot))
          drhoe_g = flatn(i,j)*(rhoe_g - Ip(i,j,1,2,QREINT))
          der(:)  = flatn(i,j)*(er(:)  - Ip(i,j,1,2,qrad:qradhi))
          
          dup    = flatn(i,j)*(u    - Ip(i,j,1,3,QU))
          dptotp = flatn(i,j)*(ptot - Ip(i,j,1,3,qptot))
          
          alpham = 0.5d0*(dptotm/(rho*cc) - dum)*rho/cc
          alphap = 0.5d0*(dptotp/(rho*cc) + dup)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphav = dv
          alphar(:) = der(:)- dptot/csq*hr
          
          if (u-cc .gt. 0.d0) then
             alpham = -alpham
          else if (u-cc .lt. 0.d0) then
             alpham = 0.d0
          else
             alpham = -0.5d0*alpham
          endif
          if (u+cc .gt. 0.d0) then
             alphap = -alphap
          else if (u+cc .lt. 0.d0) then
             alphap = 0.d0
          else
             alphap = -0.5d0*alphap
          endif
          if (u .gt. 0.d0) then
             alpha0 = -alpha0
             alphav = -alphav
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else if (u .lt. 0.d0) then
             alpha0 = 0.d0
             alphav = 0.d0
             alphae = 0.d0
             alphae_g = 0.d0
             alphar(:) = 0.d0
          else
             alpha0 = -0.5d0*alpha0
             alphav = -0.5d0*alphav
             alphae = -0.5d0*alphae
             alphae_g = -0.5d0*alphae_g
             alphar(:) = -0.5d0*alphar(:)
          endif
          
          if (i .le. ihi1) then
             qxm(i+1,j,QRHO) = rho + alphap + alpham + alpha0
             qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
             qxm(i+1,j,QU) = u + (alphap - alpham)*cc/rho
             qxm(i+1,j,QV) = v + alphav
             qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
             qxm(i+1,j,qrad:qradhi) = qrtmp
             qxm(i+1,j,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
             qxm(i+1,j,QPRES) = p + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))
             qxm(i+1,j,qptot) = ptot + (alphap + alpham)*csq
             qxm(i+1,j,qreitot) = qxm(i+1,j,QREINT) + sum(qrtmp)
             
             do g=0, ngroups-1
                if (qxm(i+1,j,qrad+g) < 0.d0) then
                   er_foo = - qxm(i+1,j,qrad+g)
                   qxm(i+1,j,qrad+g) = 0.d0
                   qxm(i+1,j,qptot) = qxm(i+1,j,qptot) + lamm(g) * er_foo
                   qxm(i+1,j,qreitot) = qxm(i+1,j,qreitot) + er_foo
                end if
             end do
             
             if (qxm(i+1,j,QPRES) < 0.d0) then
                qxm(i+1,j,QPRES) = p
             end if
          end if
          
          if(dloga(i,j).ne.0)then
             courn = dtdx*(cc+abs(u))
             eta = (1.d0-courn)/(cc*dt*abs(dloga(i,j)))
             dlogatmp = min(eta,1.d0)*dloga(i,j)
             sourcr = -0.5d0*dt*rho*dlogatmp*u
             sourcp = sourcr*cgassq
             source = sourcr*h_g
             sourcer(:) = -0.5d0*dt*dlogatmp*u*(lam0(:)+1.d0)*er(:)
             if (i .le. ihi1) then
                qxm(i+1,j,QRHO) = qxm(i+1,j,QRHO) + sourcr
                qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
                qxm(i+1,j,QPRES) = qxm(i+1,j,QPRES) + sourcp
                qxm(i+1,j,QREINT) = qxm(i+1,j,QREINT) + source
                qxm(i+1,j,qrad:qradhi) = qxm(i+1,j,qrad:qradhi) + sourcer(:)
                !              qxm(i+1,j,qptot ) = sum(lamm(:)*qxm(i+1,j,qrad:qradhi)) + qxm(i+1,j,QPRES)
                qxm(i+1,j,qptot) = qxm(i+1,j,qptot) + sum(lamm(:)*sourcer(:)) + sourcp
                qxm(i+1,j,qreitot) = sum(qxm(i+1,j,qrad:qradhi))  + qxm(i+1,j,QREINT)
             end if
             if (i .ge. ilo1) then
                qxp(i,j,QRHO) = qxp(i,j,QRHO) + sourcr
                qxp(i,j,QRHO) = max(qxp(i,j,QRHO),small_dens)
                qxp(i,j,QPRES) = qxp(i,j,QPRES) + sourcp
                qxp(i,j,QREINT) = qxp(i,j,QREINT) + source
                qxp(i,j,qrad:qradhi) = qxp(i,j,qrad:qradhi) + sourcer(:)
                !              qxp(i,j,qptot ) = sum(lamp(:)*qxp(i,j,qrad:qradhi)) + qxp(i,j,QPRES)
                qxp(i,j,qptot) = qxp(i,j,qptot) + sum(lamp(:)*sourcer(:)) + sourcp
                qxp(i,j,qreitot) = sum(qxp(i,j,qrad:qradhi))  + qxp(i,j,QREINT)
             end if
          endif
          
       end do
    end do
    
    ! Now do the passively advected quantities
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
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


    
    ! Trace to bottom and top edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1
          
          do g=0, ngroups-1
             lam0(g) = lam(i,j,g)
             lamp(g) = lam(i,j,g)
             lamm(g) = lam(i,j,g)
          end do
          
          cgassq = cg(i,j)**2
          cc = c(i,j)
          csq = cc**2
          
          rho = q(i,j,QRHO)
          u = q(i,j,QU)
          v = q(i,j,QV)
          p = q(i,j,QPRES)
          rhoe_g = q(i,j,QREINT)
          h_g = (p+rhoe_g) / rho
          er(:) = q(i,j,qrad:qradhi)
          hr(:) = (lam0+1.d0)*er/rho
          ptot = q(i,j,qptot)
          rhoe = q(i,j,qreitot)
          enth = ( (rhoe+ptot)/rho )/csq
          
          ! plus state on face j
          dvm    = flatn(i,j)*(v    - Im(i,j,2,1,QV))
          dptotm = flatn(i,j)*(ptot - Im(i,j,2,1,qptot))
          
          drho    = flatn(i,j)*(rho    - Im(i,j,2,2,QRHO))
          du      = flatn(i,j)*(u      - Im(i,j,2,2,QU))
          dptot   = flatn(i,j)*(ptot   - Im(i,j,2,2,qptot))
          drhoe   = flatn(i,j)*(rhoe   - Im(i,j,2,2,qreitot))
          drhoe_g = flatn(i,j)*(rhoe_g - Im(i,j,2,2,QREINT))
          der(:)  = flatn(i,j)*(er(:)  - Im(i,j,2,2,qrad:qradhi))
          
          dvp    = flatn(i,j)*(v    - Im(i,j,2,3,QV))
          dptotp = flatn(i,j)*(ptot - Im(i,j,2,3,qptot))
          
          alpham = 0.5d0*(dptotm/(rho*cc) - dvm)*rho/cc
          alphap = 0.5d0*(dptotp/(rho*cc) + dvp)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphau = du
          alphar(:) = der(:) - dptot/csq*hr
          
          if (v-cc .gt. 0.d0) then
             alpham = 0.d0
          else if (v-cc .lt. 0.d0) then 
             alpham = -alpham
          else
             alpham = -0.5d0*alpham
          endif
          if (v+cc .gt. 0.d0) then
             alphap = 0.d0
          else if (v+cc .lt. 0.d0) then
             alphap = -alphap
          else
             alphap = -0.5d0*alphap
          endif
          if (v .gt. 0.d0) then
             alpha0 = 0.d0
             alphau = 0.d0
             alphae = 0.d0
             alphae_g = 0.d0 
             alphar(:) = 0.d0
          else if (v .lt. 0.d0) then
             alpha0 = -alpha0
             alphau = -alphau
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else
             alpha0 = -0.5d0*alpha0
             alphau = -0.5d0*alphau
             alphae = -0.5d0*alphae
             alphae_g = -0.5d0*alphae_g
             alphar(:) = -0.5d0*alphar(:)
          endif
          
          if (j .ge. ilo2) then
             qyp(i,j,QRHO) = rho + alphap + alpham + alpha0
             qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))
             qyp(i,j,QV) = v + (alphap - alpham)*cc/rho
             qyp(i,j,QU) = u + alphau
             qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
             qyp(i,j,qrad:qradhi) = qrtmp
             qyp(i,j,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
             qyp(i,j,QPRES) = p + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))
             qyp(i,j,qptot) = ptot + (alphap + alpham)*csq 
             qyp(i,j,qreitot) = qyp(i,j,QREINT) + sum(qrtmp)
             
             do g=0, ngroups-1
                if (qyp(i,j,qrad+g) < 0.d0) then
                   er_foo = - qyp(i,j,qrad+g)
                   qyp(i,j,qrad+g) = 0.d0
                   qyp(i,j,qptot) = qyp(i,j,qptot) + lamp(g) * er_foo
                   qyp(i,j,qreitot) = qyp(i,j,qreitot) + er_foo
                end if
             end do
             
             if (qyp(i,j,QPRES) < 0.d0) then
                qyp(i,j,QPRES) = p
             end if
          end if
          
          ! minus state on face j+1
          dvm    = flatn(i,j)*(v    - Ip(i,j,2,1,QV))
          dptotm = flatn(i,j)*(ptot - Ip(i,j,2,1,qptot))
          
          drho    = flatn(i,j)*(rho    - Ip(i,j,2,2,QRHO))
          du      = flatn(i,j)*(u      - Ip(i,j,2,2,QU))
          dptot   = flatn(i,j)*(ptot   - Ip(i,j,2,2,qptot))
          drhoe   = flatn(i,j)*(rhoe   - Ip(i,j,2,2,qreitot))
          drhoe_g = flatn(i,j)*(rhoe_g - Ip(i,j,2,2,QREINT))
          der(:)  = flatn(i,j)*(er(:)  - Ip(i,j,2,2,qrad:qradhi))
          
          dvp    = flatn(i,j)*(v    - Ip(i,j,2,3,QV))
          dptotp = flatn(i,j)*(ptot - Ip(i,j,2,3,qptot))
          
          alpham = 0.5d0*(dptotm/(rho*cc) - dvm)*rho/cc
          alphap = 0.5d0*(dptotp/(rho*cc) + dvp)*rho/cc
          alpha0 = drho - dptot/csq 
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphau = du
          alphar(:) = der(:)- dptot/csq*hr
          
          if (v-cc .gt. 0.d0) then
             alpham = -alpham
          else if (v-cc .lt. 0.d0) then
             alpham = 0.d0
          else
             alpham = -0.5d0*alpham
          endif
          if (v+cc .gt. 0.d0) then
             alphap = -alphap
          else if (v+cc .lt. 0.d0) then
             alphap = 0.d0
          else
             alphap = -0.5d0*alphap
          endif
          if (v .gt. 0.d0) then
             alpha0 = -alpha0
             alphau = -alphau
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else if (v .lt. 0.d0) then
             alpha0 = 0.d0
             alphau = 0.d0
             alphae = 0.d0
             alphae_g = 0.d0
             alphar(:) = 0.d0
          else
             alpha0 = -0.5d0*alpha0
             alphau = -0.5d0*alphau
             alphae = -0.5d0*alphae
             alphae_g = -0.5d0*alphae_g
             alphar(:) = -0.5d0*alphar(:)
          endif
          
          if (j .le. ihi2) then
             qym(i,j+1,QRHO) = rho + alphap + alpham + alpha0
             qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))
             qym(i,j+1,QV) = v + (alphap - alpham)*cc/rho
             qym(i,j+1,QU) = u + alphau
             qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
             qym(i,j+1,qrad:qradhi) = qrtmp
             qym(i,j+1,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
             qym(i,j+1,QPRES) = p + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))
             qym(i,j+1,qptot) = ptot + (alphap + alpham)*csq
             qym(i,j+1,qreitot) = qym(i,j+1,QREINT) + sum(qrtmp)
             
             do g=0, ngroups-1
                if (qym(i,j+1,qrad+g) < 0.d0) then
                   er_foo = - qym(i,j+1,qrad+g)
                   qym(i,j+1,qrad+g) = 0.d0
                   qym(i,j+1,qptot) = qym(i,j+1,qptot) + lamm(g) * er_foo
                   qym(i,j+1,qreitot) = qym(i,j+1,qreitot) + er_foo
                end if
             end do
             
             if (qym(i,j+1,QPRES) < 0.d0) then
                qym(i,j+1,QPRES) = p
             end if
          end if
          
       end do
    end do
    
    ! Now do the passively advected quantities
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
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
    
    deallocate(Ip,Im)
    
  end subroutine trace_ppm_rad

end module trace_ppm_rad_module
