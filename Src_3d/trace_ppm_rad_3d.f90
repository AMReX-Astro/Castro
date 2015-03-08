module trace_ppm_rad_module

  implicit none

contains
  
  subroutine tracexy_ppm_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       q,c,cg,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
       Ip,Im, &
       qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
       ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QPRES, QFA, QFS, QFX, nadv, small_dens, &
         ppm_type
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot
    use rad_params_module, only : ngroups

    implicit none

    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer kc,k3d

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision    cg(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)
    double precision   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)

    double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QRADVAR)
    double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QRADVAR)
    double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QRADVAR)
    double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QRADVAR)
    double precision dx, dy, dt

    ! Local variables
    integer i, j, g
    integer n, iadv
    integer ns, ispec, iaux

    double precision, dimension(0:ngroups-1) :: er, der, alphar, qrtmp,hr
    double precision, dimension(0:ngroups-1) :: lam0, lamp, lamm
    double precision cc, csq, rho, u, v, w, p, ptot, rhoe, enth, cgassq
    double precision dum, dvm, dptotm
    double precision drho, du, dv, dw, drhoe, dptot
    double precision dup, dvp, dptotp

    double precision alpham, alphap, alpha0, alphae, alphau, alphav, alphaw

    double precision rhoe_g, h_g, alphae_g, drhoe_g

    double precision :: er_foo

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call bl_error("Error:: RadHydro_3d.f90 :: tracexy_ppm_rad")
    end if

    ! Trace to left and right edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          do g=0, ngroups-1
             lam0(g) = lam(i,j,k3d,g)
             lamp(g) = lam(i,j,k3d,g)
             lamm(g) = lam(i,j,k3d,g)
          end do

          cgassq = cg(i,j,k3d)**2
          cc = c(i,j,k3d)
          csq = cc**2

          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          h_g = (p+rhoe_g) / rho
          er(:) = q(i,j,k3d,qrad:qradhi)
          hr(:) = (lam0+1.d0)*er/rho
          ptot = q(i,j,k3d,qptot)
          rhoe = q(i,j,k3d,qreitot)
          enth = ( (rhoe+ptot)/rho )/csq

          ! plus state on face i
          dum    = flatn(i,j,k3d)*(u    - Im(i,j,kc,1,1,QU))
          dptotm = flatn(i,j,k3d)*(ptot - Im(i,j,kc,1,1,qptot))

          drho    = flatn(i,j,k3d)*(rho    - Im(i,j,kc,1,2,QRHO))
          dv      = flatn(i,j,k3d)*(v      - Im(i,j,kc,1,2,QV))
          dw      = flatn(i,j,k3d)*(w      - Im(i,j,kc,1,2,QW))
          dptot   = flatn(i,j,k3d)*(ptot   - Im(i,j,kc,1,2,qptot))
          drhoe   = flatn(i,j,k3d)*(rhoe   - Im(i,j,kc,1,2,qreitot))
          drhoe_g = flatn(i,j,k3d)*(rhoe_g - Im(i,j,kc,1,2,QREINT))
          der(:)  = flatn(i,j,k3d)*(er(:)  - Im(i,j,kc,1,2,qrad:qradhi))

          dup    = flatn(i,j,k3d)*(u    - Im(i,j,kc,1,3,QU))
          dptotp = flatn(i,j,k3d)*(ptot - Im(i,j,kc,1,3,qptot))

          alpham = 0.5d0*(dptotm/(rho*cc) - dum)*rho/cc
          alphap = 0.5d0*(dptotp/(rho*cc) + dup)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphav = dv
          alphaw = dw
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
             alphaw = 0.d0
             alphae = 0.d0
             alphae_g = 0.d0
             alphar(:) = 0.d0
          else if (u .lt. 0.d0) then
             alpha0 = -alpha0
             alphav = -alphav
             alphaw = -alphaw
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else
             alpha0 = -0.5d0*alpha0
             alphav = -0.5d0*alphav
             alphaw = -0.5d0*alphaw
             alphae = -0.5d0*alphae
             alphae_g = -0.5d0*alphae_g
             alphar(:) = -0.5d0*alphar(:)
          endif

          if (i .ge. ilo1) then
             qxp(i,j,kc,QRHO) = rho + alphap + alpham + alpha0
             qxp(i,j,kc,QRHO) = max(small_dens,qxp(i,j,kc,QRHO))
             qxp(i,j,kc,QU) = u + (alphap - alpham)*cc/rho
             qxp(i,j,kc,QV) = v + alphav
             qxp(i,j,kc,QW) = w + alphaw
             qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
             qxp(i,j,kc,qrad:qradhi) = qrtmp
             qxp(i,j,kc,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
             qxp(i,j,kc,QPRES) = p + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))
             qxp(i,j,kc,qptot) = ptot + (alphap + alpham)*csq 
             qxp(i,j,kc,qreitot) = qxp(i,j,kc,QREINT) + sum(qrtmp)

             do g=0,ngroups-1
                if (qxp(i,j,kc,qrad+g) < 0.d0) then
                   er_foo = - qxp(i,j,kc,qrad+g)
                   qxp(i,j,kc,qrad+g) = 0.d0
                   qxp(i,j,kc,qptot) = qxp(i,j,kc,qptot) + lamp(g) * er_foo
                   qxp(i,j,kc,qreitot) = qxp(i,j,kc,qreitot) + er_foo
                end if
             end do

             if (qxp(i,j,kc,QPRES) < 0.d0) then
                qxp(i,j,kc,QPRES) = p
             end if
          end if

          ! minus state on face i+1
          dum    = flatn(i,j,k3d)*(u    - Ip(i,j,kc,1,1,QU))
          dptotm = flatn(i,j,k3d)*(ptot - Ip(i,j,kc,1,1,qptot))

          drho    = flatn(i,j,k3d)*(rho    - Ip(i,j,kc,1,2,QRHO))
          dv      = flatn(i,j,k3d)*(v      - Ip(i,j,kc,1,2,QV))
          dw      = flatn(i,j,k3d)*(w      - Ip(i,j,kc,1,2,QW))
          dptot   = flatn(i,j,k3d)*(ptot   - Ip(i,j,kc,1,2,qptot))
          drhoe   = flatn(i,j,k3d)*(rhoe   - Ip(i,j,kc,1,2,qreitot))
          drhoe_g = flatn(i,j,k3d)*(rhoe_g - Ip(i,j,kc,1,2,QREINT))
          der(:)  = flatn(i,j,k3d)*(er(:)  - Ip(i,j,kc,1,2,qrad:qradhi))

          dup    = flatn(i,j,k3d)*(u    - Ip(i,j,kc,1,3,QU))
          dptotp = flatn(i,j,k3d)*(ptot - Ip(i,j,kc,1,3,qptot))

          alpham = 0.5d0*(dptotm/(rho*cc) - dum)*rho/cc
          alphap = 0.5d0*(dptotp/(rho*cc) + dup)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphav = dv
          alphaw = dw
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
             alphaw = -alphaw
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else if (u .lt. 0.d0) then
             alpha0 = 0.d0
             alphav = 0.d0
             alphaw = 0.d0
             alphae = 0.d0
             alphae_g = 0.d0
             alphar(:) = 0.d0
          else
             alpha0 = -0.5d0*alpha0
             alphav = -0.5d0*alphav
             alphaw = -0.5d0*alphaw
             alphae = -0.5d0*alphae
             alphae_g = -0.5d0*alphae_g
             alphar(:) = -0.5d0*alphar(:)
          endif

          if (i .le. ihi1) then
             qxm(i+1,j,kc,QRHO) = rho + alphap + alpham + alpha0
             qxm(i+1,j,kc,QRHO) = max(qxm(i+1,j,kc,QRHO),small_dens)
             qxm(i+1,j,kc,QU) = u + (alphap - alpham)*cc/rho
             qxm(i+1,j,kc,QV) = v + alphav
             qxm(i+1,j,kc,QW) = w + alphaw
             qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
             qxm(i+1,j,kc,qrad:qradhi) = qrtmp
             qxm(i+1,j,kc,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
             qxm(i+1,j,kc,QPRES) = p + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))
             qxm(i+1,j,kc,qptot) = ptot + (alphap + alpham)*csq
             qxm(i+1,j,kc,qreitot) = qxm(i+1,j,kc,QREINT) + sum(qrtmp)

             do g=0,ngroups-1
                if (qxm(i+1,j,kc,qrad+g) < 0.d0) then
                   er_foo = - qxm(i+1,j,kc,qrad+g)
                   qxm(i+1,j,kc,qrad+g) = 0.d0
                   qxm(i+1,j,kc,qptot) = qxm(i+1,j,kc,qptot) + lamm(g) * er_foo
                   qxm(i+1,j,kc,qreitot) = qxm(i+1,j,kc,qreitot) + er_foo
                end if
             end do

             if (qxm(i+1,j,kc,QPRES) < 0.d0) then
                qxm(i+1,j,kc,QPRES) = p
             end if
          end if
       end do
    end do

    ! Now do the passively advected quantities
    do iadv = 1, nadv
       n = QFA + iadv - 1
       do j = ilo2-1, ihi2+1

          ! plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,QU)
             if (u .gt. 0.d0) then
                qxp(i,j,kc,n) = q(i,j,k3d,n)
             else if (u .lt. 0.d0) then
                qxp(i,j,kc,n) = q(i,j,k3d,n) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else
                qxp(i,j,kc,n) = q(i,j,k3d,n) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo

          ! minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)
             if (u .gt. 0.d0) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) &
                     + flatn(i,j,k3d)*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else if (u .lt. 0.d0) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n)
             else
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) &
                     + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo

       enddo
    enddo

    do ispec = 1, nspec
       ns = QFS + ispec - 1

       do j = ilo2-1, ihi2+1

          ! plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,QU)
             if (u .gt. 0.d0) then
                qxp(i,j,kc,ns) = q(i,j,k3d,ns)
             else if (u .lt. 0.d0) then
                qxp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             else
                qxp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

          ! minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)
             if (u .gt. 0.d0) then
                qxm(i+1,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Ip(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             else if (u .lt. 0.d0) then
                qxm(i+1,j,kc,ns) = q(i,j,k3d,ns)
             else
                qxm(i+1,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

       enddo
    enddo

    do iaux = 1, naux
       ns = QFX + iaux - 1
       do j = ilo2-1, ihi2+1

          ! plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,QU)
             if (u .gt. 0.d0) then
                qxp(i,j,kc,ns) = q(i,j,k3d,ns)
             else if (u .lt. 0.d0) then
                qxp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             else
                qxp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

          ! minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)
             if (u .gt. 0.d0) then
                qxm(i+1,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Ip(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             else if (u .lt. 0.d0) then
                qxm(i+1,j,kc,ns) = q(i,j,k3d,ns)
             else
                qxm(i+1,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

       enddo
    enddo

    ! Trace to bottom and top edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          do g=0, ngroups-1
             lam0(g) = lam(i,j,k3d,g)
             lamp(g) = lam(i,j,k3d,g)
             lamm(g) = lam(i,j,k3d,g)
          end do

          cgassq = cg(i,j,k3d)**2
          cc = c(i,j,k3d)
          csq = cc**2

          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          h_g = (p+rhoe_g) / rho
          er(:)= q(i,j,k3d,qrad:qradhi)
          hr(:) = (lam0+1.d0)*er/rho
          ptot = q(i,j,k3d,qptot)
          rhoe = q(i,j,k3d,qreitot)
          enth = ( (rhoe+ptot)/rho )/csq

          ! plus state on face j
          dvm    = flatn(i,j,k3d)*(v    - Im(i,j,kc,2,1,QV))
          dptotm = flatn(i,j,k3d)*(ptot - Im(i,j,kc,2,1,qptot))

          drho    = flatn(i,j,k3d)*(rho    - Im(i,j,kc,2,2,QRHO))
          du      = flatn(i,j,k3d)*(u      - Im(i,j,kc,2,2,QU))
          dw      = flatn(i,j,k3d)*(w      - Im(i,j,kc,2,2,QW))
          dptot   = flatn(i,j,k3d)*(ptot   - Im(i,j,kc,2,2,qptot))
          drhoe   = flatn(i,j,k3d)*(rhoe   - Im(i,j,kc,2,2,qreitot))
          drhoe_g = flatn(i,j,k3d)*(rhoe_g - Im(i,j,kc,2,2,QREINT))
          der(:)  = flatn(i,j,k3d)*(er(:)  - Im(i,j,kc,2,2,qrad:qradhi))

          dvp    = flatn(i,j,k3d)*(v    - Im(i,j,kc,2,3,QV))
          dptotp = flatn(i,j,k3d)*(ptot - Im(i,j,kc,2,3,qptot))

          alpham = 0.5d0*(dptotm/(rho*cc) - dvm)*rho/cc
          alphap = 0.5d0*(dptotp/(rho*cc) + dvp)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphau = du
          alphaw = dw
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
             alphaw = 0.d0
             alphae = 0.d0
             alphae_g = 0.d0
             alphar(:) = 0.d0
          else if (v .lt. 0.d0) then
             alpha0 = -alpha0
             alphau = -alphau
             alphaw = -alphaw
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else
             alpha0 = -0.5d0*alpha0
             alphau = -0.5d0*alphau
             alphaw = -0.5d0*alphaw
             alphae = -0.5d0*alphae
             alphae_g = -0.5d0*alphae_g
             alphar(:) = -0.5d0*alphar(:)
          endif

          if (j .ge. ilo2) then
             qyp(i,j,kc,QRHO) = rho + alphap + alpham + alpha0
             qyp(i,j,kc,QRHO) = max(small_dens, qyp(i,j,kc,QRHO))
             qyp(i,j,kc,QV) = v + (alphap - alpham)*cc/rho
             qyp(i,j,kc,QU) = u + alphau
             qyp(i,j,kc,QW) = w + alphaw
             qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
             qyp(i,j,kc,qrad:qradhi) = qrtmp
             qyp(i,j,kc,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
             qyp(i,j,kc,QPRES) = p + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))
             qyp(i,j,kc,qptot) = ptot + (alphap + alpham)*csq 
             qyp(i,j,kc,qreitot) = qyp(i,j,kc,QREINT) + sum(qrtmp)

             do g=0,ngroups-1
                if (qyp(i,j,kc,qrad+g) < 0.d0) then
                   er_foo = - qyp(i,j,kc,qrad+g)
                   qyp(i,j,kc,qrad+g) = 0.d0
                   qyp(i,j,kc,qptot) = qyp(i,j,kc,qptot) + lamp(g) * er_foo
                   qyp(i,j,kc,qreitot) = qyp(i,j,kc,qreitot) + er_foo
                end if
             end do

             if (qyp(i,j,kc,QPRES) < 0.d0) then
                qyp(i,j,kc,QPRES) = p
             end if
          end if

          ! minus state on face j+1
          dvm    = flatn(i,j,k3d)*(v    - Ip(i,j,kc,2,1,QV))
          dptotm = flatn(i,j,k3d)*(ptot - Ip(i,j,kc,2,1,qptot))

          drho    = flatn(i,j,k3d)*(rho    - Ip(i,j,kc,2,2,QRHO))
          du      = flatn(i,j,k3d)*(u      - Ip(i,j,kc,2,2,QU))
          dw      = flatn(i,j,k3d)*(w      - Ip(i,j,kc,2,2,QW))
          dptot   = flatn(i,j,k3d)*(ptot   - Ip(i,j,kc,2,2,qptot))
          drhoe   = flatn(i,j,k3d)*(rhoe   - Ip(i,j,kc,2,2,qreitot))
          drhoe_g = flatn(i,j,k3d)*(rhoe_g - Ip(i,j,kc,2,2,QREINT))
          der(:)  = flatn(i,j,k3d)*(er(:)  - Ip(i,j,kc,2,2,qrad:qradhi))

          dvp    = flatn(i,j,k3d)*(v    - Ip(i,j,kc,2,3,QV))
          dptotp = flatn(i,j,k3d)*(ptot - Ip(i,j,kc,2,3,qptot))

          alpham = 0.5d0*(dptotm/(rho*cc) - dvm)*rho/cc
          alphap = 0.5d0*(dptotp/(rho*cc) + dvp)*rho/cc
          alphae_g = drhoe_g - dptot/csq*h_g
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphau = du
          alphaw = dw
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
             alphaw = -alphaw
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else if (v .lt. 0.d0) then
             alpha0 = 0.d0
             alphau = 0.d0
             alphaw = 0.d0
             alphae = 0.d0
             alphae_g = 0.d0
             alphar(:) = 0.d0
          else
             alpha0 = -0.5d0*alpha0
             alphau = -0.5d0*alphau
             alphaw = -0.5d0*alphaw
             alphae = -0.5d0*alphae
             alphae_g = -0.5d0*alphae_g
             alphar(:) = -0.5d0*alphar(:)
          endif

          if (j .le. ihi2) then
             qym(i,j+1,kc,QRHO) = rho + alphap + alpham + alpha0
             qym(i,j+1,kc,QRHO) = max(small_dens, qym(i,j+1,kc,QRHO))
             qym(i,j+1,kc,QV) = v + (alphap - alpham)*cc/rho
             qym(i,j+1,kc,QU) = u + alphau
             qym(i,j+1,kc,QW) = w + alphaw
             qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
             qym(i,j+1,kc,qrad:qradhi) = qrtmp
             qym(i,j+1,kc,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
             qym(i,j+1,kc,QPRES) = p + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))
             qym(i,j+1,kc,qptot) = ptot + (alphap + alpham)*csq
             qym(i,j+1,kc,qreitot) = qym(i,j+1,kc,QREINT) + sum(qrtmp)

             do g=0,ngroups-1
                if (qym(i,j+1,kc,qrad+g) < 0.d0) then
                   er_foo = - qym(i,j+1,kc,qrad+g)
                   qym(i,j+1,kc,qrad+g) = 0.d0
                   qym(i,j+1,kc,qptot) = qym(i,j+1,kc,qptot) + lamm(g) * er_foo
                   qym(i,j+1,kc,qreitot) = qym(i,j+1,kc,qreitot) + er_foo
                end if
             end do

             if (qym(i,j+1,kc,QPRES) < 0.d0) then
                qym(i,j+1,kc,QPRES) = p
             end if
          end if
       end do
    end do

    ! Now do the passively advected quantities
    do iadv = 1, nadv
       n = QFA + iadv - 1
       do i = ilo1-1, ihi1+1

          ! plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,k3d,QV)
             if (v .gt. 0.d0) then
                qyp(i,j,kc,n) = q(i,j,k3d,n)
             else if (v .lt. 0.d0) then
                qyp(i,j,kc,n) = q(i,j,k3d,n) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else
                qyp(i,j,kc,n) = q(i,j,k3d,n) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo

          ! minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,k3d,QV)
             if (v .gt. 0.d0) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n) &
                     + flatn(i,j,k3d)*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else if (v .lt. 0.d0) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n)
             else
                qym(i,j+1,kc,n) = q(i,j,k3d,n) &
                     + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo

       enddo
    enddo

    do ispec = 1, nspec
       ns = QFS + ispec - 1
       do i = ilo1-1, ihi1+1

          ! plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,k3d,QV)
             if (v .gt. 0.d0) then
                qyp(i,j,kc,ns) = q(i,j,k3d,ns)
             else if (v .lt. 0.d0) then
                qyp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             else
                qyp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

          ! minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,k3d,QV)
             if (v .gt. 0.d0) then
                qym(i,j+1,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Ip(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             else if (v .lt. 0.d0) then
                qym(i,j+1,kc,ns) = q(i,j,k3d,ns)
             else
                qym(i,j+1,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

       enddo
    enddo

    do iaux = 1, naux
       ns = QFX + iaux - 1
       do i = ilo1-1, ihi1+1

          ! plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,k3d,QV)
             if (v .gt. 0.d0) then
                qyp(i,j,kc,ns) = q(i,j,k3d,ns)
             else if (v .lt. 0.d0) then
                qyp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             else
                qyp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

          ! minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,k3d,QV)
             if (v .gt. 0.d0) then
                qym(i,j+1,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Ip(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             else if (v .lt. 0.d0) then
                qym(i,j+1,kc,ns) = q(i,j,k3d,ns)
             else
                qym(i,j+1,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

       enddo
    enddo

  end subroutine tracexy_ppm_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine tracez_ppm_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       q,c,cg,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
       Ip,Im, &
       qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
       ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QPRES, QFA, QFS, QFX, nadv, small_dens, &
         ppm_type
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot
    use rad_params_module, only : ngroups

    implicit none

    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer km,kc,k3d

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision    cg(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)
    double precision   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)
    double precision qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QRADVAR)
    double precision qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QRADVAR)
    double precision dz, dt

    !     Local variables
    integer i, j, g
    integer n, iadv
    integer ns, ispec, iaux

    double precision, dimension(0:ngroups-1) :: er, der, alphar, qrtmp,hr
    double precision, dimension(0:ngroups-1) :: lam0, lamp, lamm
    double precision cc, csq, rho, u, v, w, p, ptot, rhoe, enth, cgassq
    double precision dwm, dptotm
    double precision drho, du, dv, drhoe, dptot
    double precision dwp, dptotp

    double precision alpham, alphap, alpha0, alphae, alphau, alphav

    double precision rhoe_g, h_g, alphae_g, drhoe_g

    double precision :: er_foo

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call bl_error("Error:: RadHydro_3d.f90 :: tracez_ppm_rad")
    end if

!!!!!!!!!!!!!!!
    ! PPM CODE
!!!!!!!!!!!!!!!

    ! Trace to left and right edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          do g=0, ngroups-1
             lam0(g) = lam(i,j,k3d,g)
             lamp(g) = lam(i,j,k3d,g)
             lamm(g) = lam(i,j,k3d,g)
          end do

          cgassq = cg(i,j,k3d)**2
          cc = c(i,j,k3d)
          csq = cc**2

          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          h_g = (p+rhoe_g) / rho
          er(:) = q(i,j,k3d,qrad:qradhi)
          hr(:) = (lam0+1.d0)*er/rho
          ptot = q(i,j,k3d,qptot)
          rhoe = q(i,j,k3d,qreitot)
          enth = ( (rhoe+ptot)/rho )/csq

          ! plus state on face kc
          dwm    = flatn(i,j,k3d)*(w    - Im(i,j,kc,3,1,QW))
          dptotm = flatn(i,j,k3d)*(ptot - Im(i,j,kc,3,1,qptot))

          drho    = flatn(i,j,k3d)*(rho    - Im(i,j,kc,3,2,QRHO))
          du      = flatn(i,j,k3d)*(u      - Im(i,j,kc,3,2,QU))
          dv      = flatn(i,j,k3d)*(v      - Im(i,j,kc,3,2,QV))
          dptot   = flatn(i,j,k3d)*(ptot   - Im(i,j,kc,3,2,qptot))
          drhoe   = flatn(i,j,k3d)*(rhoe   - Im(i,j,kc,3,2,qreitot))
          drhoe_g = flatn(i,j,k3d)*(rhoe_g - Im(i,j,kc,3,2,QREINT))
          der(:)  = flatn(i,j,k3d)*(er(:)  - Im(i,j,kc,3,2,qrad:qradhi))

          dwp    = flatn(i,j,k3d)*(w    - Im(i,j,kc,3,3,QW))
          dptotp = flatn(i,j,k3d)*(ptot - Im(i,j,kc,3,3,qptot))

          alpham = 0.5d0*(dptotm/(rho*cc) - dwm)*rho/cc
          alphap = 0.5d0*(dptotp/(rho*cc) + dwp)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphau = du
          alphav = dv
          alphar(:) = der(:) - dptot/csq*hr

          if (w-cc .gt. 0.d0) then
             alpham = 0.d0
          else if (w-cc .lt. 0.d0) then
             alpham = -alpham
          else
             alpham = -0.5d0*alpham
          endif
          if (w+cc .gt. 0.d0) then
             alphap = 0.d0
          else if (w+cc .lt. 0.d0) then
             alphap = -alphap
          else
             alphap = -0.5d0*alphap
          endif
          if (w .gt. 0.d0) then
             alpha0 = 0.d0
             alphau = 0.d0
             alphav = 0.d0
             alphae = 0.d0
             alphae_g = 0.d0
             alphar(:) = 0.d0
          else if (w .lt. 0.d0) then
             alpha0 = -alpha0
             alphau = -alphau
             alphav = -alphav
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else
             alpha0 = -0.5d0*alpha0
             alphau = -0.5d0*alphau
             alphav = -0.5d0*alphav
             alphae = -0.5d0*alphae
             alphae_g = -0.5d0*alphae_g
             alphar(:) = -0.5d0*alphar(:)
          endif

          qzp(i,j,kc,QRHO) = rho + alphap + alpham + alpha0
          qzp(i,j,kc,QRHO) = max(small_dens, qzp(i,j,kc,QRHO))
          qzp(i,j,kc,QW) = w + (alphap - alpham)*cc/rho 
          qzp(i,j,kc,QU) = u + alphau
          qzp(i,j,kc,QV) = v + alphav
          qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
          qzp(i,j,kc,qrad:qradhi) = qrtmp
          qzp(i,j,kc,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
          qzp(i,j,kc,QPRES) = p + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))
          qzp(i,j,kc,qptot) = ptot + (alphap + alpham)*csq 
          qzp(i,j,kc,qreitot) = qzp(i,j,kc,QREINT) + sum(qrtmp)

          do g=0,ngroups-1
             if (qzp(i,j,kc,qrad+g) < 0.d0) then
                er_foo = - qzp(i,j,kc,qrad+g)
                qzp(i,j,kc,qrad+g) = 0.d0
                qzp(i,j,kc,qptot) = qzp(i,j,kc,qptot) + lamp(g) * er_foo
                qzp(i,j,kc,qreitot) = qzp(i,j,kc,qreitot) + er_foo
             end if
          end do

          if (qzp(i,j,kc,QPRES) < 0.d0) then
             qzp(i,j,kc,QPRES) = p
          end if

          ! minus state on face kc
          ! note this is different from how we do 1D, 2D, and the
          ! x and y-faces in 3D, where the analogous thing would have
          ! been to find the minus state on face kc+1

          do g=0, ngroups-1
             lam0(g) = lam(i,j,k3d-1,g)
             lamp(g) = lam(i,j,k3d-1,g)
             lamm(g) = lam(i,j,k3d-1,g)
          end do

          cgassq = cg(i,j,k3d-1)**2
          cc = c(i,j,k3d-1)
          csq = cc**2

          rho = q(i,j,k3d-1,QRHO)
          u = q(i,j,k3d-1,QU)
          v = q(i,j,k3d-1,QV)
          w = q(i,j,k3d-1,QW)
          p = q(i,j,k3d-1,QPRES)
          rhoe_g = q(i,j,k3d-1,QREINT)
          h_g = (p+rhoe_g) / rho
          er(:) = q(i,j,k3d-1,qrad:qradhi)
          hr(:) = (lam0+1.d0)*er/rho
          ptot = q(i,j,k3d-1,qptot)
          rhoe = q(i,j,k3d-1,qreitot)
          enth = ( (rhoe+ptot)/rho )/csq

          dwm    = flatn(i,j,k3d-1)*(w    - Ip(i,j,km,3,1,QW))
          dptotm = flatn(i,j,k3d-1)*(ptot - Ip(i,j,km,3,1,qptot))

          drho    = flatn(i,j,k3d-1)*(rho    - Ip(i,j,km,3,2,QRHO))
          du      = flatn(i,j,k3d-1)*(u      - Ip(i,j,km,3,2,QU))
          dv      = flatn(i,j,k3d-1)*(v      - Ip(i,j,km,3,2,QV))
          dptot   = flatn(i,j,k3d-1)*(ptot   - Ip(i,j,km,3,2,qptot))
          drhoe   = flatn(i,j,k3d-1)*(rhoe   - Ip(i,j,km,3,2,qreitot))
          drhoe_g = flatn(i,j,k3d-1)*(rhoe_g - Ip(i,j,km,3,2,QREINT))
          der(:)  = flatn(i,j,k3d-1)*(er(:)  - Ip(i,j,km,3,2,qrad:qradhi))

          dwp    = flatn(i,j,k3d-1)*(w    - Ip(i,j,km,3,3,QW))
          dptotp = flatn(i,j,k3d-1)*(ptot - Ip(i,j,km,3,3,qptot))

          alpham = 0.5d0*(dptotm/(rho*cc) - dwm)*rho/cc
          alphap = 0.5d0*(dptotp/(rho*cc) + dwp)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphau = du
          alphav = dv
          alphar(:) = der(:) - dptot/csq*hr

          if (w-cc .gt. 0.d0) then
             alpham = -alpham
          else if (w-cc .lt. 0.d0) then
             alpham = 0.d0
          else
             alpham = -0.5d0*alpham
          endif
          if (w+cc .gt. 0.d0) then
             alphap = -alphap
          else if (w+cc .lt. 0.d0) then
             alphap = 0.d0
          else
             alphap = -0.5d0*alphap
          endif
          if (w .gt. 0.d0) then
             alpha0 = -alpha0
             alphau = -alphau
             alphav = -alphav
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else if (w .lt. 0.d0) then
             alpha0 = 0.d0
             alphau = 0.d0
             alphav = 0.d0
             alphae = 0.d0
             alphae_g = 0.d0
             alphar(:) = 0.d0
          else
             alpha0 = -0.5d0*alpha0
             alphau = -0.5d0*alphau
             alphav = -0.5d0*alphav
             alphae = -0.5d0*alphae
             alphae_g = -0.5d0*alphae_g
             alphar(:) = -0.5d0*alphar(:)
          endif

          qzm(i,j,kc,QRHO) = rho + alphap + alpham + alpha0
          qzm(i,j,kc,QRHO) = max(small_dens, qzm(i,j,kc,QRHO))
          qzm(i,j,kc,QW) = w + (alphap - alpham)*cc/rho 
          qzm(i,j,kc,QU) = u + alphau
          qzm(i,j,kc,QV) = v + alphav
          qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
          qzm(i,j,kc,qrad:qradhi) = qrtmp
          qzm(i,j,kc,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
          qzm(i,j,kc,QPRES) = p + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))
          qzm(i,j,kc,qptot) = ptot + (alphap + alpham)*csq
          qzm(i,j,kc,qreitot) = qzm(i,j,kc,QREINT) + sum(qrtmp)

          do g=0,ngroups-1
             if (qzm(i,j,kc,qrad+g) < 0.d0) then
                er_foo = - qzm(i,j,kc,qrad+g)
                qzm(i,j,kc,qrad+g) = 0.d0
                qzm(i,j,kc,qptot) = qzm(i,j,kc,qptot) + lamm(g) * er_foo
                qzm(i,j,kc,qreitot) = qzm(i,j,kc,qreitot) + er_foo
             end if
          end do

          if (qzm(i,j,kc,QPRES) < 0.d0) then
             qzm(i,j,kc,QPRES) = p
          end if
       end do
    end do

    ! Now do the passively advected quantities
    do iadv = 1, nadv
       n = QFA + iadv - 1
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! plus state on face kc
             w = q(i,j,k3d,QW)
             if (w .gt. 0.d0) then
                qzp(i,j,kc,n) = q(i,j,k3d,n)
             else if (w .lt. 0.d0) then
                qzp(i,j,kc,n) = q(i,j,k3d,n) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
             else
                qzp(i,j,kc,n) = q(i,j,k3d,n) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
             endif

             ! minus state on face k
             w = q(i,j,k3d-1,QW)
             if (w .gt. 0.d0) then
                qzm(i,j,kc,n) = q(i,j,k3d-1,n) &
                     + flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
             else if (w .lt. 0.d0) then
                qzm(i,j,kc,n) = q(i,j,k3d-1,n)
             else
                qzm(i,j,kc,n) = q(i,j,k3d-1,n) &
                     + 0.5d0*flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
             endif

          enddo
       enddo
    enddo

    do ispec = 1, nspec
       ns = QFS + ispec - 1
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! plus state on face kc
             w = q(i,j,k3d,QW)
             if (w .gt. 0.d0) then
                qzp(i,j,kc,ns) = q(i,j,k3d,ns)
             else if (w .lt. 0.d0) then
                qzp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,3,2,ns) - q(i,j,k3d,ns))
             else
                qzp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,3,2,ns) - q(i,j,k3d,ns))
             endif

             ! minus state on face k
             w = q(i,j,k3d-1,QW)
             if (w .gt. 0.d0) then
                qzm(i,j,kc,ns) = q(i,j,k3d-1,ns) &
                     + flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,ns) - q(i,j,k3d-1,ns))
             else if (w .lt. 0.d0) then
                qzm(i,j,kc,ns) = q(i,j,k3d-1,ns)
             else
                qzm(i,j,kc,ns) = q(i,j,k3d-1,ns) &
                     + 0.5d0*flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,ns) - q(i,j,k3d-1,ns))
             endif

          enddo
       enddo
    enddo

    do iaux = 1, naux
       ns = QFX + iaux - 1
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! plus state on face kc
             w = q(i,j,k3d,QW)
             if (w .gt. 0.d0) then
                qzp(i,j,kc,ns) = q(i,j,k3d,ns)
             else if (w .lt. 0.d0) then
                qzp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,3,2,ns) - q(i,j,k3d,ns))
             else
                qzp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,3,2,ns) - q(i,j,k3d,ns))
             endif

             ! minus state on face k
             w = q(i,j,k3d-1,QW)
             if (w .gt. 0.d0) then
                qzm(i,j,kc,ns) = q(i,j,k3d-1,ns) &
                     + flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,ns) - q(i,j,k3d-1,ns))
             else if (w .lt. 0.d0) then
                qzm(i,j,kc,ns) = q(i,j,k3d-1,ns)
             else
                qzm(i,j,kc,ns) = q(i,j,k3d-1,ns) &
                     + 0.5d0*flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,ns) - q(i,j,k3d-1,ns))
             endif

          enddo
       enddo
    enddo

  end subroutine tracez_ppm_rad

end module trace_ppm_rad_module
