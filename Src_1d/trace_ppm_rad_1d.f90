module trace_ppm_rad_module

  implicit none

contains

  subroutine trace_ppm_rad(lam, lam_l1, lam_h1, &
       q,dq,c,cg,flatn,qd_l1,qd_h1, &
       dloga,dloga_l1,dloga_h1, &
       srcQ,src_l1,src_h1,&
       grav,gv_l1,gv_h1, &
       qxm,qxp,qpd_l1,qpd_h1, &
       ilo,ihi,domlo,domhi,dx,dt)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QREINT, QPRES, &
                                   small_dens, ppm_type, fix_mass_flux, &
                                   npassive, qpass_map
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot
    use rad_params_module, only : ngroups
    use ppm_module, only : ppm

    implicit none

    integer ilo,ihi
    integer domlo(1),domhi(1)
    integer lam_l1, lam_h1
    integer    qd_l1,   qd_h1
    integer dloga_l1,dloga_h1
    integer   qpd_l1,  qpd_h1
    integer   src_l1,  src_h1
    integer    gv_l1,   gv_h1
    double precision dx, dt
    double precision lam(lam_l1:lam_h1, 0:ngroups-1)
    double precision     q( qd_l1: qd_h1,QRADVAR)
    double precision  srcQ(src_l1:src_h1,QVAR)
    double precision flatn(qd_l1:qd_h1)
    double precision     c(qd_l1:qd_h1)
    double precision    cg(qd_l1:qd_h1)
    double precision dloga(dloga_l1:dloga_h1)

    double precision   dq( qpd_l1: qpd_h1,QRADVAR)
    double precision  qxm( qpd_l1: qpd_h1,QRADVAR)
    double precision  qxp( qpd_l1: qpd_h1,QRADVAR)
    double precision grav(  gv_l1:  gv_h1)

    !     Local variables
    integer i, g
    integer n, ipassive

    double precision hdt,dtdx

    double precision, dimension(0:ngroups-1) :: er,der,alphar,sourcer,qrtmp,hr
    double precision, dimension(0:ngroups-1) :: lam0, lamp, lamm
    double precision cc, csq, rho, u, p, ptot, rhoe, enth, cgassq
    double precision dum, dptotm
    double precision drho, drhoe, dptot
    double precision dup, dptotp

    double precision alpham, alphap, alpha0, alphae
    double precision sourcr,sourcp,source,courn,eta,dlogatmp

    double precision rhoe_g, h_g, alphae_g, drhoe_g

    logical :: fix_mass_flux_lo, fix_mass_flux_hi

    double precision, allocatable :: Ip(:,:,:)
    double precision, allocatable :: Im(:,:,:)

    double precision :: er_foo

    fix_mass_flux_lo = (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) &
         .and. (ilo .eq. domlo(1))
    fix_mass_flux_hi = (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) &
         .and. (ihi .eq. domhi(1))

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
       call bl_error("Error:: RadHydro_1d.f90 :: trace_ppm_rad")
    end if

    hdt = 0.5d0 * dt
    dtdx = dt/dx

    allocate(Ip(ilo-1:ihi+1,3,QRADVAR))
    allocate(Im(ilo-1:ihi+1,3,QRADVAR))

    ! Compute Ip and Im
    do n=1,QRADVAR
       call ppm(q(:,n),qd_l1,qd_h1,q(:,QU),c, flatn, &
            Ip(:,:,n),Im(:,:,n),ilo,ihi,dx,dt)
    end do

    ! Trace to left and right edges using upwind PPM
    do i = ilo-1, ihi+1

       do g=0, ngroups-1
          lam0(g) = lam(i,g)
          lamp(g) = lam(i,g)
          lamm(g) = lam(i,g)
       end do

       cgassq = cg(i)**2
       cc = c(i)
       csq = cc**2

       rho = q(i,QRHO)
       u = q(i,QU)
       p = q(i,QPRES)
       rhoe_g = q(i,QREINT)
       h_g = (p+rhoe_g) / rho
       er(:) = q(i,qrad:qradhi)
       hr(:) = (lam0+1.d0)*er/rho
       ptot = q(i,qptot)
       rhoe = q(i,qreitot)
       enth = ( (rhoe+ptot)/rho )/csq

       ! plus state on face i
       dum    = flatn(i)*(u    - Im(i,1,QU))
       dptotm = flatn(i)*(ptot - Im(i,1,qptot))

       drho  = flatn(i)*(rho  - Im(i,2,QRHO))
       drhoe_g = flatn(i)*(rhoe_g  - Im(i,2,QREINT))
       drhoe = flatn(i)*(rhoe - Im(i,2,qreitot))
       dptot = flatn(i)*(ptot - Im(i,2,qptot))
       der(:)= flatn(i)*(er(:)- Im(i,2,qrad:qradhi))

       dup    = flatn(i)*(u    - Im(i,3,QU))
       dptotp = flatn(i)*(ptot - Im(i,3,qptot))

       alpham = 0.5d0*(dptotm/(rho*cc) - dum)*rho/cc
       alphap = 0.5d0*(dptotp/(rho*cc) + dup)*rho/cc
       alpha0 = drho - dptot/csq
       alphae = drhoe - dptot*enth
       alphae_g = drhoe_g - dptot/csq*h_g
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
          alphae = 0.d0
          alphae_g = 0.d0
          alphar(:) = 0.d0
       else if (u .lt. 0.d0) then
          alpha0 = -alpha0
          alphae = -alphae
          alphae_g = -alphae_g
          alphar(:) = -alphar(:)
       else
          alpha0 = -0.5d0*alpha0
          alphae = -0.5d0*alphae
          alphae_g = -0.5d0*alphae_g
          alphar(:) = -0.5d0*alphar(:)
       endif

       if (i .ge. ilo) then
          qxp(i,QRHO) = rho + alphap + alpham + alpha0
          qxp(i,QRHO) = max(small_dens,qxp(i,QRHO))
          qxp(i,QU) = u + (alphap - alpham)*cc/rho
          qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
          qxp(i,qrad:qradhi) = qrtmp
          qxp(i,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
          qxp(i,QPRES) = p + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))
          qxp(i,qptot) = ptot + (alphap + alpham)*csq
          qxp(i,qreitot) = qxp(i,QREINT) + sum(qrtmp)

          ! add non-gravitational source term
          qxp(i  ,QRHO  )  = qxp(i,QRHO   ) + hdt*srcQ(i,QRHO)
          qxp(i  ,QRHO  )  = max(small_dens,qxp(i,QRHO))
          qxp(i  ,QU    )  = qxp(i,QU     ) + hdt*srcQ(i,QU)
          qxp(i  ,QREINT)  = qxp(i,QREINT ) + hdt*srcQ(i,QREINT)
          qxp(i  ,QPRES )  = qxp(i,QPRES  ) + hdt*srcQ(i,QPRES)
          qxp(i  ,qptot )  = qxp(i,qptot  ) + hdt*srcQ(i,QPRES)
          qxp(i  ,qreitot) = qxp(i,qreitot) + hdt*srcQ(i,QREINT)

          ! add gravitational source term
          qxp(i  ,QU) = qxp(i,QU) + hdt*grav(i)

          do g=0, ngroups-1
             if (qxp(i,qrad+g) < 0.d0) then
                er_foo = - qxp(i,qrad+g)
                qxp(i,qrad+g) = 0.d0
                qxp(i,qptot) = qxp(i,qptot) + lamp(g) * er_foo
                qxp(i,qreitot) = qxp(i,qreitot) + er_foo
             end if
          end do

          if ( qxp(i,QPRES) < 0.d0 ) then
             qxp(i,QPRES) = p
          end if
       end if

       ! minus state on face i+1
       dum    = flatn(i)*(u    - Ip(i,1,QU))
       dptotm = flatn(i)*(ptot - Ip(i,1,qptot))

       drho  = flatn(i)*(rho  - Ip(i,2,QRHO))
       drhoe_g = flatn(i)*(rhoe_g - Ip(i,2,QREINT))
       drhoe = flatn(i)*(rhoe - Ip(i,2,qreitot))
       dptot = flatn(i)*(ptot - Ip(i,2,qptot))
       der(:)= flatn(i)*(er(:)- Ip(i,2,qrad:qradhi))

       dup    = flatn(i)*(u    - Ip(i,3,QU))
       dptotp = flatn(i)*(ptot - Ip(i,3,qptot))

       alpham = 0.5d0*(dptotm/(rho*cc) - dum)*rho/cc
       alphap = 0.5d0*(dptotp/(rho*cc) + dup)*rho/cc
       alpha0 = drho - dptot/csq
       alphae = drhoe - dptot*enth
       alphae_g = drhoe_g - dptot/csq*h_g
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
          alphae = -alphae
          alphae_g = -alphae_g
          alphar(:) = -alphar(:)
       else if (u .lt. 0.d0) then
          alpha0 = 0.d0
          alphae = 0.d0
          alphae_g = 0.d0
          alphar(:) = 0.d0
       else
          alpha0 = -0.5d0*alpha0
          alphae = -0.5d0*alphae
          alphae_g = -0.5d0*alphae_g
          alphar(:) = -0.5d0*alphar(:)
       endif

       if (i .le. ihi) then
          qxm(i+1,QRHO) = rho + alphap + alpham + alpha0
          qxm(i+1,QRHO) = max(qxm(i+1,QRHO),small_dens)
          qxm(i+1,QU) = u + (alphap - alpham)*cc/rho
          qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
          qxm(i+1,qrad:qradhi) = qrtmp
          qxm(i+1,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
          qxm(i+1,QPRES) = p + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))
          qxm(i+1,qptot) = ptot + (alphap + alpham)*csq
          qxm(i+1,qreitot) = qxm(i+1,QREINT) + sum(qrtmp)

          ! add non-gravitational source term
          qxm(i+1,QRHO   ) = qxm(i+1,QRHO   ) + hdt*srcQ(i,QRHO)
          qxm(i+1,QRHO   ) = max(small_dens, qxm(i+1,QRHO))
          qxm(i+1,QU     ) = qxm(i+1,QU     ) + hdt*srcQ(i,QU)
          qxm(i+1,QREINT ) = qxm(i+1,QREINT ) + hdt*srcQ(i,QREINT)
          qxm(i+1,QPRES  ) = qxm(i+1,QPRES  ) + hdt*srcQ(i,QPRES)
          qxm(i+1,qptot  ) = qxm(i+1,qptot  ) + hdt*srcQ(i,QPRES)
          qxm(i+1,qreitot) = qxm(i+1,qreitot) + hdt*srcQ(i,QREINT)

          ! add gravitational source term
          qxm(i+1,QU) = qxm(i+1,QU) + hdt*grav(i)

          do g=0, ngroups-1
             if (qxm(i+1,qrad+g) < 0.d0) then
                er_foo = - qxm(i+1,qrad+g)
                qxm(i+1,qrad+g) = 0.d0
                qxm(i+1,qptot) = qxm(i+1,qptot) + lamm(g) * er_foo
                qxm(i+1,qreitot) = qxm(i+1,qreitot) + er_foo
             end if
          end do

          if ( qxm(i+1,QPRES) < 0.d0 ) then
             qxm(i+1,QPRES) = p
          end if
       end if

       if(dloga(i).ne.0)then
          courn = dtdx*(cc+abs(u))
          eta = (1.d0-courn)/(cc*dt*abs(dloga(i)))
          dlogatmp = min(eta,1.d0)*dloga(i)
          sourcr = -0.5d0*dt*rho*dlogatmp*u
          sourcp = sourcr*cgassq
          source = sourcr*h_g
          sourcer(:) = -0.5d0*dt*dlogatmp*u*(lam0(:)+1.d0)*er(:)
          if (i .le. ihi) then
             qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + sourcr
             qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
             qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + sourcp
             qxm(i+1,QREINT) = qxm(i+1,QREINT) + source
             qxm(i+1,qrad:qradhi) = qxm(i+1,qrad:qradhi) + sourcer(:)
             !           qxm(i+1,qptot ) = sum(lamm(:)*qxm(i+1,qrad:qradhi)) + qxm(i+1,QPRES)
             qxm(i+1,qptot) = qxm(i+1,qptot) + sum(lamm(:)*sourcer(:)) + sourcp
             qxm(i+1,qreitot) = sum(qxm(i+1,qrad:qradhi))  + qxm(i+1,QREINT)
          end if
          if (i .ge. ilo) then
             qxp(i  ,QRHO  ) = qxp(i  ,QRHO  ) + sourcr
             qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
             qxp(i  ,QPRES ) = qxp(i  ,QPRES ) + sourcp
             qxp(i  ,QREINT) = qxp(i  ,QREINT) + source
             qxp(i  ,qrad:qradhi) = qxp(i  ,qrad:qradhi) + sourcer(:)
             !           qxp(i  ,qptot ) = sum(lamp(:)*qxp(i,qrad:qradhi)) + qxp(i,QPRES)
             qxp(i,qptot) = qxp(i,qptot) + sum(lamp(:)*sourcer(:)) + sourcp
             qxp(i  ,qreitot) = sum(qxp(i,qrad:qradhi))  + qxp(i,QREINT)
          end if
       endif

    end do

    ! Enforce constant mass flux rate if specified
    if (fix_mass_flux_lo) then
       qxm(ilo,QRHO   ) = q(domlo(1)-1,QRHO)
       qxm(ilo,QU     ) = q(domlo(1)-1,QU  )
       qxm(ilo,QPRES  ) = q(domlo(1)-1,QPRES)
       qxm(ilo,QREINT ) = q(domlo(1)-1,QREINT)
       qxm(ilo,qrad:qradhi) = q(domlo(1)-1,qrad:qradhi)
       qxm(ilo,qptot  ) = q(domlo(1)-1,qptot)
       qxm(ilo,qreitot) = q(domlo(1)-1,qreitot)
    end if

    ! Enforce constant mass flux rate if specified
    if (fix_mass_flux_hi) then
       qxp(ihi+1,QRHO   ) = q(domhi(1)+1,QRHO)
       qxp(ihi+1,QU     ) = q(domhi(1)+1,QU  )
       qxp(ihi+1,QPRES  ) = q(domhi(1)+1,QPRES)
       qxp(ihi+1,QREINT ) = q(domhi(1)+1,QREINT)
       qxp(ihi+1,qrad:qradhi) = q(domhi(1)+1,qrad:qradhi)
       qxp(ihi+1,qptot  ) = q(domhi(1)+1,qptot)
       qxp(ihi+1,qreitot) = q(domhi(1)+1,qreitot)
    end if

    ! Now do the passively advected quantities
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! plus state on face i
       do i = ilo, ihi+1
          u = q(i,QU)
          if (u .gt. 0.d0) then
             qxp(i,n) = q(i,n)
          else if (u .lt. 0.d0) then
             qxp(i,n) = q(i,n) + flatn(i)*(Im(i,2,n) - q(i,n))
          else
             qxp(i,n) = q(i,n) + 0.5d0*flatn(i)*(Im(i,2,n) - q(i,n))
          endif
       enddo
       if (fix_mass_flux_hi) qxp(ihi+1,n) = q(ihi+1,n)

       ! minus state on face i+1
       do i = ilo-1, ihi
          u = q(i,QU)
          if (u .gt. 0.d0) then
             qxm(i+1,n) = q(i,n) + flatn(i)*(Ip(i,2,n) - q(i,n))
          else if (u .lt. 0.d0) then
             qxm(i+1,n) = q(i,n)
          else
             qxm(i+1,n) = q(i,n) + 0.5d0*flatn(i)*(Ip(i,2,n) - q(i,n))
          endif
       enddo
       if (fix_mass_flux_lo) qxm(ilo,n) = q(ilo-1,n)

    enddo

    deallocate(Ip,Im)

  end subroutine trace_ppm_rad

end module trace_ppm_rad_module
