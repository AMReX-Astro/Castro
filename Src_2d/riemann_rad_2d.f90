module riemann_rad_module

  use bl_constants_module
  use bl_error_module

  implicit none

  private

  public cmpflx_rad

contains


  subroutine cmpflx_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2,&
                        qm,qp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                        rflx,rflx_l1,rflx_l2,rflx_h1,rflx_h2, &
                        pgd, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
                        ergd,ergd_l1,ergd_l2,ergd_h1,ergd_h2, &
                        lmgd,lmgd_l1,lmgd_l2,lmgd_h1,lmgd_h2, &
                        ugd, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                        vgd, &
                        gamc,gamcg,csml,c,qd_l1,qd_l2,qd_h1,qd_h2, &
                        idir,ilo,ihi,jlo,jhi)

    use meth_params_module, only : NVAR, hybrid_riemann, use_colglaz
    use radhydro_params_module, only : QRADVAR
    use rad_params_module, only : ngroups

    implicit none

    integer lam_l1, lam_l2, lam_h1, lam_h2
    integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer flx_l1,flx_l2,flx_h1,flx_h2
    integer rflx_l1,rflx_l2,rflx_h1,rflx_h2
    integer pgd_l1,pgd_l2,pgd_h1,pgd_h2
    integer ergd_l1,ergd_l2,ergd_h1,ergd_h2
    integer lmgd_l1,lmgd_l2,lmgd_h1,lmgd_h2
    integer ugd_l1,ugd_l2,ugd_h1,ugd_h2
    integer qd_l1,qd_l2,qd_h1,qd_h2
    integer idir,ilo,ihi,jlo,jhi

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
    double precision    qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision    qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision   flx( flx_l1: flx_h1, flx_l2: flx_h2,NVAR)
    double precision  rflx(rflx_l1:rflx_h1,rflx_l2:rflx_h2,0:ngroups-1)
    double precision  pgd( pgd_l1: pgd_h1, pgd_l2: pgd_h2)
    double precision ergd(ergd_l1:ergd_h1,ergd_l2:ergd_h2,0:ngroups-1)
    double precision lmgd(lmgd_l1:lmgd_h1,lmgd_l2:lmgd_h2,0:ngroups-1)
    double precision  ugd( ugd_l1: ugd_h1, ugd_l2: ugd_h2)
    double precision  vgd( ugd_l1: ugd_h1, ugd_l2: ugd_h2)
    double precision  gamc (qd_l1:qd_h1,qd_l2:qd_h2)
    double precision  gamcg(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision     c (qd_l1:qd_h1,qd_l2:qd_h2)
    double precision   csml(qd_l1:qd_h1,qd_l2:qd_h2)

    !     Local variables
    integer i, j

    double precision, allocatable :: smallc(:,:), cavg(:,:)
    double precision, allocatable :: gamcm(:,:), gamcp(:,:)
    double precision, allocatable :: gamcgm(:,:), gamcgp(:,:)

    allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate ( gamcgm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate ( gamcgp(ilo-1:ihi+1,jlo-1:jhi+1) )

    if (hybrid_riemann == 1) then
       call bl_error("ERROR: hybrid Riemann not supported for radiation")
    endif

    if (use_colglaz == 1) then
       call bl_error("ERROR: the Colella-Glaz Riemann solver is not supported for radiation")
    endif

    if(idir.eq.1) then
       do j = jlo, jhi
          do i = ilo, ihi+1
             smallc(i,j) = max( csml(i,j), csml(i-1,j) )
             cavg(i,j) = HALF*( c(i,j) + c(i-1,j) )
             gamcm(i,j) = gamc(i-1,j)
             gamcp(i,j) = gamc(i,j)
             gamcgm(i,j) = gamcg(i-1,j)
             gamcgp(i,j) = gamcg(i,j)
          enddo
       enddo
    else
       do j = jlo, jhi+1
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j), csml(i,j-1) )
             cavg(i,j) = HALF*( c(i,j) + c(i,j-1) )
             gamcm(i,j) = gamc(i,j-1)
             gamcp(i,j) = gamc(i,j)
             gamcgm(i,j) = gamcg(i,j-1)
             gamcgp(i,j) = gamcg(i,j)
          enddo
       enddo
    endif

    !     Solve Riemann problem (godunov state passed back, but only (u,p,Er) saved
    call riemannus_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2,&
                       qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                       gamcm, gamcp, gamcgm, gamcgp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                       flx,  flx_l1,  flx_l2,  flx_h1,  flx_h2, &
                       rflx, rflx_l1, rflx_l2, rflx_h1, rflx_h2, &
                       pgd,  pgd_l1,  pgd_l2,  pgd_h1,  pgd_h2, &
                       ergd, ergd_l1, ergd_l2, ergd_h1, ergd_h2, &
                       lmgd, lmgd_l1, lmgd_l2, lmgd_h1, lmgd_h2, &
                       ugd,  ugd_l1,  ugd_l2,  ugd_h1,  ugd_h2, &
                       vgd, &
                       idir, ilo, ihi, jlo, jhi)

    deallocate (smallc,cavg,gamcp,gamcm,gamcgp,gamcgm)

  end subroutine cmpflx_rad


  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine riemannus_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2, &
                           ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                           gamcl, gamcr, gamcgl, gamcgr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
                           uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
                           rflx, rflx_l1, rflx_l2, rflx_h1, rflx_h2, &
                           pgdnv,  pgd_l1,  pgd_l2,  pgd_h1,  pgd_h2, &
                           ergdnv, ergd_l1, ergd_l2, ergd_h1, ergd_h2, &
                           lmgdnv, lmgd_l1, lmgd_l2, lmgd_h1, lmgd_h2, &
                           ugdnv,  ugd_l1,  ugd_l2,  ugd_h1,  ugd_h2, &
                           vgdnv, &
                           idir, ilo1, ihi1, ilo2, ihi2)

    use network, only : nspec, naux
    use prob_params_module, only : physbc_lo,Symmetry
    use meth_params_module, only : NVAR, QRHO, QU, QV, QPRES, QREINT, &
         URHO, UMX, UMY, UEDEN, UEINT, small_dens, small_pres, &
         npassive, qpass_map, upass_map
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, fspace_type
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor

    implicit none

    double precision, parameter:: small = 1.d-8

    integer :: lam_l1, lam_l2, lam_h1, lam_h2
    integer :: qpd_l1, qpd_l2, qpd_h1, qpd_h2
    integer :: gd_l1, gd_l2, gd_h1, gd_h2
    integer :: uflx_l1, uflx_l2, uflx_h1, uflx_h2
    integer :: rflx_l1, rflx_l2, rflx_h1, rflx_h2
    integer :: pgd_l1, pgd_l2, pgd_h1, pgd_h2
    integer :: ergd_l1, ergd_l2, ergd_h1, ergd_h2
    integer :: lmgd_l1, lmgd_l2, lmgd_h1, lmgd_h2
    integer :: ugd_l1, ugd_l2, ugd_h1, ugd_h2
    integer :: idir, ilo1, ihi1, ilo2, ihi2
    integer :: ilo,ihi,jlo,jhi

    double precision :: lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision :: gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: gamcgl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: gamcgr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: rflx(rflx_l1:rflx_h1,rflx_l2:rflx_h2,0:ngroups-1)
    double precision ::  pgdnv( pgd_l1: pgd_h1, pgd_l2: pgd_h2)
    double precision :: ergdnv(ergd_l1:ergd_h1,ergd_l2:ergd_h2,0:ngroups-1)
    double precision :: lmgdnv(lmgd_l1:lmgd_h1,lmgd_l2:lmgd_h2,0:ngroups-1)
    double precision ::  ugdnv( ugd_l1: ugd_h1, ugd_l2: ugd_h2)
    double precision ::  vgdnv( ugd_l1: ugd_h1, ugd_l2: ugd_h2)

    double precision :: rgdnv, ustar
    double precision, dimension(0:ngroups-1) :: erl, err
    double precision :: rl, ul, vl, pl, rel, pl_g, rel_g, wl
    double precision :: rr, ur, vr, pr, rer, pr_g, rer_g, wr
    double precision :: rstar, cstar, pstar
    double precision :: ro, uo, po, reo, co, gamco, reo_g, po_g, co_g, gamco_g
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall, qavg
    double precision :: rhoetot, scr

    double precision :: regdnv_g, pgdnv_g, pgdnv_t
    double precision :: drho, estar_g, pstar_g
    double precision, dimension(0:ngroups-1) :: lambda, reo_r, po_r, estar_r, regdnv_r
    double precision :: eddf, f1
    integer :: i, j, g

    integer :: ipassive, n, nq

    !************************************************************
    !  set min/max based on normal direction
    if(idir.eq.1) then
       ilo = ilo1
       ihi = ihi1 + 1
       jlo = ilo2
       jhi = ihi2
    else
       ilo = ilo1
       ihi = ihi1
       jlo = ilo2
       jhi = ihi2+1
    endif

    !     Solve Riemann Problem
    !     NOTE: The calling routine will order velocity unknowns so that
    !     for the purposes of this routine, the normal component is always
    !     loaded in the QU slot.
    do j = jlo, jhi
       do i = ilo, ihi

          rl = ql(i,j,QRHO)

          !  pick left velocities based on direction
          if(idir.eq.1) then
             ul = ql(i,j,QU)
             vl = ql(i,j,QV)
          else
             ul = ql(i,j,QV)
             vl = ql(i,j,QU)
          endif

          pl = ql(i,j,qptot)
          rel = ql(i,j,qreitot)
          erl(:) = ql(i,j,qrad:qradhi)
          pl_g = ql(i,j,QPRES)
          rel_g = ql(i,j,QREINT)

          rr = qr(i,j,QRHO)

          !  pick right velocities based on direction
          if(idir.eq.1) then
             ur = qr(i,j,QU)
             vr = qr(i,j,QV)
          else
             ur = qr(i,j,QV)
             vr = qr(i,j,QU)
          endif

          pr  = qr(i,j,qptot)
          rer = qr(i,j,qreitot)
          err(:) = qr(i,j,qrad:qradhi)
          pr_g = qr(i,j,QPRES)
          rer_g = qr(i,j,QREINT)

          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = max(pstar,small_pres)
          ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

          if (ustar .gt. ZERO) then
             if (idir.eq.1) then
                lambda(:) = lam(i-1,j,:)
             else
                lambda(:) = lam(i,j-1,:)
             end if
             ro = rl
             uo = ul
             po = pl
             po_g = pl_g
             po_r(:) = erl(:) * lambda(:)
             reo = rel
             reo_r(:) = erl(:)
             reo_g = rel_g
             gamco = gamcl(i,j)
             gamco_g = gamcgl(i,j)

          else if (ustar .lt. ZERO) then
             lambda(:) = lam(i,j,:)
             ro = rr
             uo = ur
             po = pr
             po_g = pr_g
             po_r(:) = err(:) * lambda(:)
             reo = rer
             reo_r(:) = err(:)
             reo_g = rer_g
             gamco = gamcr(i,j)
             gamco_g = gamcgr(i,j)

          else
             if (idir.eq.1) then
                do g=0, ngroups-1
                   lambda(g) = 2.0d0*(lam(i-1,j,g)*lam(i,j,g))/(lam(i-1,j,g)+lam(i,j,g)+1.d-50)
                end do
             else
                do g=0, ngroups-1
                   lambda(g) = 2.0d0*(lam(i,j-1,g)*lam(i,j,g))/(lam(i,j-1,g)+lam(i,j,g)+1.d-50)
                end do
             end if
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             reo = HALF*(rel+rer)
             reo_r(:) = HALF*(erl(:)+err(:))
             reo_g = HALF*(rel_g+rer_g)
             po = HALF*(pl+pr)
             po_r(:) = lambda(:) * reo_r(:)
             po_g = HALF*(pr_g+pl_g)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
             gamco_g = HALF*(gamcgl(i,j)+gamcgr(i,j))
          endif
          ro = max(small_dens,ro)

          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)
          drho = (pstar - po)/co**2
          rstar = ro + drho
          rstar = max(small_dens,rstar)
          estar_g = reo_g + drho*(reo_g + po_g)/ro
          co_g = sqrt(abs(gamco_g*po_g/ro))
          co_g = max(csmall,co_g)
          pstar_g = po_g + drho*co_g**2
          pstar_g = max(pstar_g,small_pres)
          estar_r(:) = reo_r(:) + drho*(reo_r(:) + po_r(:))/ro
          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)

          sgnm = sign(ONE,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = HALF*(spin + spout)
          if (pstar-po .ge. ZERO) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin .eq. ZERO) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif
          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))

          if (ustar .gt. ZERO) then
             vgdnv(i,j) = vl
          else if (ustar .lt. ZERO) then
             vgdnv(i,j) = vr
          else
             vgdnv(i,j) = HALF*(vl+vr)
          endif

          rgdnv = frac*rstar + (ONE - frac)*ro
          ugdnv(i,j) = frac*ustar + (ONE - frac)*uo
          pgdnv_t = frac*pstar + (ONE - frac)*po
          pgdnv_g = frac*pstar_g + (ONE - frac)*po_g
          regdnv_g = frac*estar_g + (ONE - frac)*reo_g
          regdnv_r(:) = frac*estar_r(:) + (ONE - frac)*reo_r(:)

          if (spout .lt. ZERO) then
             rgdnv = ro
             ugdnv(i,j) = uo
             pgdnv_t = po
             pgdnv_g = po_g
             regdnv_g = reo_g
             regdnv_r(:) = reo_r(:)
          endif
          if (spin .ge. ZERO) then
             rgdnv = rstar
             ugdnv(i,j) = ustar
             pgdnv_t = pstar
             pgdnv_g = pstar_g
             regdnv_g = estar_g
             regdnv_r(:) = estar_r(:)
          endif

          ! Enforce that fluxes through a symmetry plane are hard zero.
          if (i.eq.0 .and. physbc_lo(1) .eq. Symmetry .and. idir .eq. 1) ugdnv(i,j) = ZERO
          if (j.eq.0 .and. physbc_lo(2) .eq. Symmetry .and. idir .eq. 2) ugdnv(i,j) = ZERO

          do g=0, ngroups-1
             ergdnv(i,j,g) = max(regdnv_r(g), ZERO)
          end do

          pgdnv(i,j) = pgdnv_g

          lmgdnv(i,j,:) = lambda

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,URHO) = rgdnv*ugdnv(i,j)

          ! note: here we do not include the pressure, since in 2-d,
          ! for some geometries, div{F} + grad{p} cannot be written
          ! in a flux difference form
          if(idir.eq.1) then
             uflx(i,j,UMX) = uflx(i,j,URHO)*ugdnv(i,j)
             uflx(i,j,UMY) = uflx(i,j,URHO)*vgdnv(i,j)
          else
             uflx(i,j,UMX) = uflx(i,j,URHO)*vgdnv(i,j)
             uflx(i,j,UMY) = uflx(i,j,URHO)*ugdnv(i,j)
          endif

          rhoetot = regdnv_g + HALF*rgdnv*(ugdnv(i,j)**2 + vgdnv(i,j)**2)
          uflx(i,j,UEDEN) = ugdnv(i,j)*(rhoetot + pgdnv_g)
          uflx(i,j,UEINT) = ugdnv(i,j)*regdnv_g

          if (fspace_type.eq.1) then
             do g=0,ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                rflx(i,j,g) = (ONE+f1) * ergdnv(i,j,g) * ugdnv(i,j)
             end do
          else ! type 2
             do g=0,ngroups-1
                rflx(i,j,g) = ergdnv(i,j,g) * ugdnv(i,j)
             end do
          end if

          do ipassive = 1, npassive
             n = upass_map(ipassive)
             nq = qpass_map(ipassive)

             if (ustar .gt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
             else
                qavg = HALF * (ql(i,j,nq) + qr(i,j,nq))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo

       end do
    end do

  end subroutine riemannus_rad

end module riemann_rad_module
