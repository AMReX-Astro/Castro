module riemann_rad_module

  use bl_constants_module

  implicit none

  private

  public cmpflx_rad

contains
  

  subroutine cmpflx_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
                        qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        flx,   flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3, &
                        rflx, rflx_l1,rflx_l2,rflx_l3,rflx_h1,rflx_h2,rflx_h3, &
                        ugdnv, v1gdnv, v2gdnv, pgdnv, ergdnv, lmgdnv, &
                        pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                        gamc,gamcg,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        idir,ilo,ihi,jlo,jhi,kc,kflux,k3d)
    
    use meth_params_module, only : QVAR, NVAR
    use radhydro_params_module, only : QRADVAR
    use rad_params_module, only : ngroups
    
    implicit none
    
    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer  flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3
    integer rflx_l1,rflx_l2,rflx_l3,rflx_h1,rflx_h2,rflx_h3
    integer pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer idir,ilo,ihi,jlo,jhi
    integer kc,kflux,k3d
    
    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QRADVAR)
    double precision qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QRADVAR)
    double precision  flx( flx_l1: flx_h1,  flx_l2: flx_h2,  flx_l3: flx_h3, NVAR)
    double precision rflx(rflx_l1:rflx_h1, rflx_l2:rflx_h2, rflx_l3:rflx_h3,0:ngroups-1)
    double precision  ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision v1gdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision v2gdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision  pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision ergdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3,0:ngroups-1)
    double precision lmgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3,0:ngroups-1)
    double precision gamc (qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision gamcg(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision csml (qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    
    ! Local variables
    
    integer i,j  
    
    double precision, allocatable :: smallc(:,:),cavg(:,:)
    double precision, allocatable :: gamcm(:,:),gamcp(:,:)
    double precision, allocatable :: gamcgm(:,:), gamcgp(:,:)
    
    allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate ( gamcgm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate ( gamcgp(ilo-1:ihi+1,jlo-1:jhi+1) )
    
    if(idir.eq.1) then
       do j = jlo, jhi
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j,k3d), csml(i-1,j,k3d) )
             cavg(i,j) = 0.5d0*( c(i,j,k3d) + c(i-1,j,k3d) )
             gamcm(i,j) = gamc(i-1,j,k3d)
             gamcp(i,j) = gamc(i,j,k3d)
             gamcgm(i,j) = gamcg(i-1,j,k3d)
             gamcgp(i,j) = gamcg(i,j,k3d)
          enddo
       enddo
    elseif(idir.eq.2) then
       do j = jlo, jhi
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j,k3d), csml(i,j-1,k3d) )
             cavg(i,j) = 0.5d0*( c(i,j,k3d) + c(i,j-1,k3d) )
             gamcm(i,j) = gamc(i,j-1,k3d)
             gamcp(i,j) = gamc(i,j,k3d)
             gamcgm(i,j) = gamcg(i,j-1,k3d)
             gamcgp(i,j) = gamcg(i,j,k3d)
          enddo
       enddo
    else
       do j = jlo, jhi
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j,k3d), csml(i,j,k3d-1) )
             cavg(i,j) = 0.5d0*( c(i,j,k3d) + c(i,j,k3d-1) )
             gamcm(i,j) = gamc(i,j,k3d-1)
             gamcp(i,j) = gamc(i,j,k3d)
             gamcgm(i,j) = gamcg(i,j,k3d-1)
             gamcgp(i,j) = gamcg(i,j,k3d)
          enddo
       enddo
    endif
    
    ! Solve Riemann problem
    call riemannus_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
                       qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                       gamcm,gamcp,gamcgm,gamcgp,cavg,smallc,ilo-1,jlo-1,ihi+1,jhi+1, &
                       flx,   flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3, &
                       rflx, rflx_l1,rflx_l2,rflx_l3,rflx_h1,rflx_h2,rflx_h3, &
                       ugdnv, v1gdnv, v2gdnv, pgdnv, ergdnv, lmgdnv, & 
                       pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                       idir,ilo,ihi,jlo,jhi,kc,kflux,k3d)
    
    deallocate(smallc,cavg,gamcm,gamcp,gamcgm,gamcgp)
  
  end subroutine cmpflx_rad

  
  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 
  
  subroutine riemannus_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
                           ql,qr,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                           gamcl,gamcr,gamcgl,gamcgr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                           uflx, uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3, &
                           rflx, rflx_l1,rflx_l2,rflx_l3,rflx_h1,rflx_h2,rflx_h3, &
                           ugdnv, v1gdnv, v2gdnv, pgdnv, ergdnv, lmgdnv, &
                           pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                           idir,ilo,ihi,jlo,jhi,kc,kflux,k3d)
    
    use network, only : nspec, naux
    use prob_params_module, only : physbc_lo,physbc_hi,Symmetry
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QFA, QFS, &
         QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX, &
         nadv, small_dens, small_pres
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, fspace_type
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
    
    implicit none
    
    double precision, parameter:: small = 1.d-8
    
    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer gd_l1,gd_l2,gd_h1,gd_h2
    integer uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3
    integer rflx_l1,rflx_l2,rflx_l3,rflx_h1,rflx_h2,rflx_h3
    integer pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
    integer idir,ilo,ihi,jlo,jhi
    integer kc,kflux,k3d
    
    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QRADVAR)
    double precision qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QRADVAR)
    double precision  gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision  gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision gamcgl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision gamcgr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision    cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,uflx_l3:uflx_h3,NVAR)
    double precision rflx(rflx_l1:rflx_h1,rflx_l2:rflx_h2,rflx_l3:rflx_h3,0:ngroups-1)
    double precision  ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision v1gdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision v2gdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision  pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision ergdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3,0:ngroups-1)
    double precision lmgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3,0:ngroups-1)
    
    ! Local variables
    
    integer i,j, g, n, nq
    integer iadv, ispec, iaux
    double precision :: qavg
    
    double precision rgdnv, ustar
    double precision, dimension(0:ngroups-1) :: erl, err
    double precision rl, ul, v1l, v2l, pl, rel, pl_g, rel_g, wl
    double precision rr, ur, v1r, v2r, pr, rer, pr_g, rer_g, wr
    double precision rstar, cstar, pstar
    double precision ro, uo, po, reo, co, gamco, reo_g, po_g, co_g, gamco_g
    double precision sgnm, spin, spout, ushock, frac
    double precision rhoetot, scr
    
    double precision wsmall, csmall
    
    double precision :: regdnv_g, pgdnv_g, pgdnv_t
    double precision :: drho, estar_g, pstar_g
    double precision, dimension(0:ngroups-1) :: lambda, laml, lamr, reo_r, po_r, estar_r, regdnv_r
    double precision :: eddf, f1
    
    do j = jlo, jhi
       do i = ilo, ihi
          
          if (idir.eq.1) then
             laml = lam(i-1,j,k3d,:)
          else if (idir.eq.2) then
             laml = lam(i,j-1,k3d,:)
          else
             laml = lam(i,j,k3d-1,:)
          end if
          lamr = lam(i,j,k3d,:)
          
          rl = ql(i,j,kc,QRHO)
          
          ! pick left velocities based on direction
          if(idir.eq.1) then
             ul  = ql(i,j,kc,QU)
             v1l = ql(i,j,kc,QV)
             v2l = ql(i,j,kc,QW)
          elseif(idir.eq.2) then
             ul  = ql(i,j,kc,QV)
             v1l = ql(i,j,kc,QU)
             v2l = ql(i,j,kc,QW)
          else
             ul  = ql(i,j,kc,QW)
             v1l = ql(i,j,kc,QU)
             v2l = ql(i,j,kc,QV)
          endif
          
          pl = ql(i,j,kc,qptot)
          rel = ql(i,j,kc,qreitot)
          erl(:) = ql(i,j,kc,qrad:qradhi)
          pl_g = ql(i,j,kc,QPRES)
          rel_g = ql(i,j,kc,QREINT)
          
          rr = qr(i,j,kc,QRHO)
          
          ! pick right velocities based on direction
          if(idir.eq.1) then
             ur  = qr(i,j,kc,QU)
             v1r = qr(i,j,kc,QV)
             v2r = qr(i,j,kc,QW)
          elseif(idir.eq.2) then
             ur  = qr(i,j,kc,QV)
             v1r = qr(i,j,kc,QU)
             v2r = qr(i,j,kc,QW)
          else
             ur  = qr(i,j,kc,QW)
             v1r = qr(i,j,kc,QU)
             v2r = qr(i,j,kc,QV)
          endif
          
          pr = qr(i,j,kc,qptot)
          rer = qr(i,j,kc,qreitot)
          err(:) = qr(i,j,kc,qrad:qradhi)
          pr_g = qr(i,j,kc,QPRES)
          rer_g = qr(i,j,kc,QREINT)
          
          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))
          
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)
          pstar = max(pstar,small_pres)
          
          if (ustar .gt. 0.d0) then
             lambda = laml
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
          else if (ustar .lt. 0.d0) then
             lambda = lamr
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
             do g=0, ngroups-1
                lambda(g) = 2.0d0*(laml(g)*lamr(g))/(laml(g)+lamr(g)+1.d-50)
             end do
             ro = 0.5d0*(rl+rr)
             uo = 0.5d0*(ul+ur)
             reo = 0.5d0*(rel+rer)
             reo_r(:) = 0.5d0*(erl(:)+err(:))
             reo_g = 0.5d0*(rel_g+rer_g)
             po = 0.5d0*(pl+pr)
             po_r(:) = lambda(:) * reo_r(:)
             po_g = 0.5*(pr_g+pl_g)
             gamco = 0.5d0*(gamcl(i,j)+gamcr(i,j))
             gamco_g = 0.5d0*(gamcgl(i,j)+gamcgr(i,j))
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
          estar_r = reo_r(:) + drho*(reo_r(:) + po_r(:))/ro
          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)
          
          sgnm = sign(1.d0,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = 0.5d0*(spin + spout)
          if (pstar-po .ge. 0.d0) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin .eq. 0.d0) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif
          frac = (1.d0 + (spout + spin)/scr)*0.5d0
          frac = max(0.d0,min(1.d0,frac))
          
          if (ustar .gt. 0.d0) then
             v1gdnv(i,j,kc) = v1l
             v2gdnv(i,j,kc) = v2l
          else if (ustar .lt. 0.d0) then
             v1gdnv(i,j,kc) = v1r
             v2gdnv(i,j,kc) = v2r
          else
             v1gdnv(i,j,kc) = 0.5d0*(v1l+v1r)
             v2gdnv(i,j,kc) = 0.5d0*(v2l+v2r)
          endif
          
          rgdnv = frac*rstar + (1.d0 - frac)*ro
          ugdnv(i,j,kc) = frac*ustar + (1.d0 - frac)*uo
          pgdnv_t       = frac*pstar + (1.d0 - frac)*po
          pgdnv_g       = frac*pstar_g + (1.d0 - frac)*po_g
          regdnv_g = frac*estar_g + (1.d0 - frac)*reo_g
          regdnv_r(:) = frac*estar_r(:) + (1.d0 - frac)*reo_r(:)
          
          if (spout .lt. 0.d0) then
             rgdnv = ro
             ugdnv(i,j,kc) = uo
             pgdnv_t = po
             pgdnv_g = po_g
             regdnv_g = reo_g
             regdnv_r = reo_r(:)
          endif
          if (spin .ge. 0.d0) then
             rgdnv = rstar
             ugdnv(i,j,kc) = ustar
             pgdnv_t = pstar
             pgdnv_g = pstar_g
             regdnv_g = estar_g
             regdnv_r = estar_r(:)
          endif
          
          ! Enforce that fluxes through a symmetry plane are hard zero.
          if (i    .eq.0 .and. physbc_lo(1) .eq. Symmetry .and. idir .eq. 1) &
               ugdnv(i,j,kc) = 0.d0
          if (j    .eq.0 .and. physbc_lo(2) .eq. Symmetry .and. idir .eq. 2) &
               ugdnv(i,j,kc) = 0.d0
          if (k3d  .eq.0 .and. physbc_lo(3) .eq. Symmetry .and. idir .eq. 3) &
               ugdnv(i,j,kc) = 0.d0
          
          do g=0, ngroups-1
             ergdnv(i,j,kc,g) = max(regdnv_r(g), 0.d0)
          end do
          
          pgdnv(i,j,kc) = pgdnv_g
          
          lmgdnv(i,j,kc,:) = lambda
          
          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = rgdnv*ugdnv(i,j,kc)
          
          if(idir.eq.1) then
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv_g
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v1gdnv(i,j,kc)
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv(i,j,kc)
          elseif(idir.eq.2) then
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv(i,j,kc)
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv_g
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv(i,j,kc)
          else
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv(i,j,kc)
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v2gdnv(i,j,kc)
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv_g
          endif
          
          rhoetot = regdnv_g + 0.5d0*rgdnv*(ugdnv(i,j,kc)**2 + v1gdnv(i,j,kc)**2 + v2gdnv(i,j,kc)**2)
          
          uflx(i,j,kflux,UEDEN) = ugdnv(i,j,kc)*(rhoetot + pgdnv_g)
          
          uflx(i,j,kflux,UEINT) = ugdnv(i,j,kc)*regdnv_g
          
          if (fspace_type.eq.1) then
             do g=0,ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                rflx(i,j,kflux,g) = (1.d0+f1) * ergdnv(i,j,kc,g) * ugdnv(i,j,kc)
             end do
          else ! type 2
             do g=0,ngroups-1
                rflx(i,j,kflux,g) = ergdnv(i,j,kc,g) * ugdnv(i,j,kc)
             end do
          end if
          
          do iadv = 1, nadv
             n  = UFA + iadv - 1
             nq = QFA + iadv - 1
             if (ustar .gt. 0.d0) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (ustar .lt. 0.d0) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = 0.5d0 * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo
          
          do ispec = 1, nspec
             n  = UFS + ispec - 1
             nq = QFS + ispec - 1
             if (ustar .gt. 0.d0) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (ustar .lt. 0.d0) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = 0.5d0 * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo
          
          do iaux = 1, naux
             n  = UFX + iaux - 1
             nq = QFX + iaux - 1
             if (ustar .gt. 0.d0) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (ustar .lt. 0.d0) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = 0.5d0 * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo
          
       enddo
    enddo
    
  end subroutine riemannus_rad

end module riemann_rad_module
