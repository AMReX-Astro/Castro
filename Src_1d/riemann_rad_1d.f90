module riemann_rad_module

  use bl_constants_module

  implicit none

  private

  public cmpflx_rad

contains
  
  subroutine cmpflx_rad(lo,hi,domlo,domhi, &
                        lam, lam_l1, lam_h1, &       
                        qm,qp,qpd_l1,qpd_h1, &
                        flx,flx_l1,flx_h1, &
                        rflx,rflx_l1,rflx_h1, &
                        pgdnv,pg_l1,pg_h1, &
                        ergdnv,erg_l1,erg_h1, &
                        lamgdnv,lg_l1,lg_h1, &
                        ugdnv,ug_l1,ug_h1, &
                        gamc,gamcg,csml,c,qd_l1,qd_h1,ilo,ihi)
    
    use meth_params_module, only : NVAR
    use radhydro_params_module, only : QRADVAR
    use rad_params_module, only : ngroups
    
    implicit none
    
    integer lo(1),hi(1)
    integer domlo(1),domhi(1)
    integer ilo,ihi
    integer lam_l1, lam_h1
    integer qpd_l1,qpd_h1
    integer flx_l1, flx_h1
    integer rflx_l1, rflx_h1
    integer  pg_l1, pg_h1
    integer erg_l1,erg_h1
    integer  lg_l1, lg_h1
    integer  ug_l1, ug_h1
    integer  qd_l1,  qd_h1
    double precision lam(lam_l1:lam_h1, 0:ngroups-1)
    double precision    qm(qpd_l1:qpd_h1, QRADVAR)
    double precision    qp(qpd_l1:qpd_h1, QRADVAR)
    double precision   flx(flx_l1:flx_h1, NVAR)
    double precision rflx(rflx_l1:rflx_h1, 0:ngroups-1)
    double precision  pgdnv( pg_l1: pg_h1)
    double precision ergdnv(erg_l1:erg_h1, 0:ngroups-1)
    double precision lamgdnv( lg_l1: lg_h1, 0:ngroups-1)
    double precision ugdnv( ug_l1: ug_h1)
    double precision  gamc( qd_l1: qd_h1)
    double precision gamcg( qd_l1: qd_h1)
    double precision     c( qd_l1: qd_h1)
    double precision  csml( qd_l1: qd_h1)
    
    !     Local variables
    integer i
    double precision, allocatable :: smallc(:),cavg(:),gamcp(:), gamcm(:), gamcgp(:), gamcgm(:)
    
    allocate ( smallc(ilo:ihi+1) )
    allocate ( cavg(ilo:ihi+1) )
    allocate ( gamcp(ilo:ihi+1) )
    allocate ( gamcm(ilo:ihi+1) )
    allocate (gamcgp(ilo:ihi+1) )
    allocate (gamcgm(ilo:ihi+1) )
    
    do i = ilo, ihi+1 
       smallc(i) = max( csml(i), csml(i-1) )
       cavg(i) = 0.5d0*( c(i) + c(i-1) )
       gamcgm(i) = gamcg(i-1)
       gamcgp(i) = gamcg(i)
       gamcm (i) = gamc (i-1)
       gamcp (i) = gamc (i)
    enddo
    
    !     Solve Riemann problem (gdnv state passed back, but only (u,p,Er) saved)
    call riemannus_rad(lam, lam_l1, lam_h1, &       
                       qm, qp,qpd_l1,qpd_h1, smallc, cavg, &
                       gamcm, gamcp, gamcgm, gamcgp, &
                       flx, flx_l1, flx_h1, &
                       rflx, rflx_l1, rflx_h1, &
                       pgdnv, pg_l1, pg_h1, &
                       ergdnv,erg_l1,erg_h1, &
                       lamgdnv, lg_l1, lg_h1, &
                       ugdnv, ug_l1, ug_h1, ilo, ihi, domlo, domhi )
    
    deallocate (smallc,cavg,gamcp,gamcm,gamcgp,gamcgm)
  
  end subroutine cmpflx_rad


  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 
  
  subroutine riemannus_rad(lam, lam_l1, lam_h1, & 
                           ql,qr,qpd_l1,qpd_h1,smallc,cav, &
                           gamcl,gamcr,gamcgl,gamcgr,uflx,uflx_l1,uflx_h1,&
                           rflx,rflx_l1,rflx_h1, &
                           pgdnv,pg_l1,pg_h1,  &
                           ergdnv,erg_l1,erg_h1,  &
                           lamgdnv,lg_l1,lg_h1,  &
                           ugdnv,ug_l1,ug_h1, &
                           ilo,ihi,domlo,domhi)
    
    use network, only : nspec
    use meth_params_module, only : NVAR, QRHO, QU, QPRES, QREINT, &
         URHO, UMX, UEDEN, UEINT, small_dens, small_pres, &
         fix_mass_flux, npassive, qpass_map, upass_map
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow, Symmetry
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, fspace_type
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
    
    implicit none
    
    double precision, parameter:: small = 1.d-8
    
    integer ilo,ihi
    integer domlo(1),domhi(1)
    integer lam_l1, lam_h1
    integer  qpd_l1,  qpd_h1
    integer   pg_l1,   pg_h1
    integer  erg_l1,  erg_h1
    integer   lg_l1,   lg_h1
    integer   ug_l1,   ug_h1
    integer uflx_l1, uflx_h1
    integer rflx_l1, rflx_h1
    double precision lam(lam_l1:lam_h1, 0:ngroups-1)
    double precision ql(qpd_l1:qpd_h1, QRADVAR)
    double precision qr(qpd_l1:qpd_h1, QRADVAR)
    double precision    cav(ilo:ihi+1),smallc(ilo:ihi+1)
    double precision  gamcl(ilo:ihi+1), gamcr(ilo:ihi+1)
    double precision gamcgl(ilo:ihi+1),gamcgr(ilo:ihi+1)
    double precision  uflx(uflx_l1:uflx_h1, NVAR)
    double precision  rflx(rflx_l1:rflx_h1, 0:ngroups-1)
    double precision  pgdnv( pg_l1: pg_h1)
    double precision ergdnv(erg_l1:erg_h1, 0:ngroups-1)
    double precision lamgdnv( lg_l1: lg_h1, 0:ngroups-1)
    double precision ugdnv( ug_l1: ug_h1)
    
    double precision rgdnv, ustar
    double precision, dimension(0:ngroups-1) :: erl, err
    double precision rl, ul, pl, rel, pl_g, rel_g, wl
    double precision rr, ur, pr, rer, pr_g, rer_g, wr
    double precision rstar, cstar, pstar
    double precision ro, uo, po, reo, co, gamco, reo_g, po_g, co_g, gamco_g
    double precision sgnm, spin, spout, ushock, frac
    double precision rhoetot, scr

    double precision qavg
    double precision wsmall, csmall
    integer n, nq
    integer k, ipassive
    logical :: fix_mass_flux_lo, fix_mass_flux_hi
    
    double precision :: regdnv_g, pgdnv_g, pgdnv_t
    double precision :: drho, estar_g, pstar_g
    double precision, dimension(0:ngroups-1) :: lambda, reo_r, po_r, estar_r, regdnv_r
    double precision :: eddf, f1
    integer :: g
    
    !     Solve Riemann Problem
    
    fix_mass_flux_lo = (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) .and. (ilo .eq. domlo(1))
    fix_mass_flux_hi = (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) .and. (ihi .eq. domhi(1))
    
    do k = ilo, ihi+1
       
       rl  = ql(k,QRHO)
       ul  = ql(k,QU)
       pl  = ql(k,qptot)
       rel = ql(k,qreitot)
       erl(:) = ql(k,qrad:qradhi)
       pl_g = ql(k,QPRES)
       rel_g = ql(k,QREINT)
       
       rr  = qr(k,QRHO)
       ur  = qr(k,QU)
       pr  = qr(k,qptot)
       rer = qr(k,qreitot)
       err(:) = qr(k,qrad:qradhi)
       pr_g = qr(k,QPRES)
       rer_g = qr(k,QREINT) 
       
       csmall = smallc(k)
       wsmall = small_dens*csmall
       wl = max(wsmall,sqrt(abs(gamcl(k)*pl*rl)))
       wr = max(wsmall,sqrt(abs(gamcr(k)*pr*rr)))
       pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
       pstar = max(pstar,small_pres)
       ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)
       
       if (ustar .gt. 0.d0) then
          lambda(:) = lam(k-1,:)
          ro = rl
          uo = ul
          po = pl
          po_g = pl_g
          po_r(:) = erl(:) * lambda(:)
          reo = rel
          reo_r(:) = erl(:)
          reo_g = rel_g
          gamco = gamcl(k)
          gamco_g = gamcgl(k)
       else if (ustar .lt. 0.d0) then
          lambda(:) = lam(k,:)
          ro = rr
          uo = ur
          po = pr
          po_g = pr_g
          po_r(:) = err(:) * lambda(:)
          reo = rer
          reo_r(:) = err(:)
          reo_g = rer_g
          gamco = gamcr(k)
          gamco_g = gamcgr(k)
       else
          
          do g=0, ngroups-1
             lambda(g) = 0.5d0*(lam(k-1,g)+lam(k,g))
          end do
          
          ro = 0.5d0*(rl+rr)
          uo = 0.5d0*(ul+ur)
          reo = 0.5d0*(rel+rer)
          reo_r(:) = 0.5d0*(erl(:)+err(:))
          reo_g = 0.5d0*(rel_g+rer_g)
          po = 0.5d0*(pl+pr)
          po_r(:) = lambda(:) * reo_r(:)
          po_g = 0.5*(pr_g+pl_g)
          gamco = 0.5d0*(gamcl(k)+gamcr(k))
          gamco_g = 0.5d0*(gamcgl(k)+gamcgr(k))
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
       
       sgnm = sign(1.d0,ustar)
       spout = co - sgnm*uo
       spin = cstar - sgnm*ustar
       ushock = 0.5d0*(spin + spout)
       if (pstar-po .ge. 0.d0) then
          spin = ushock
          spout = ushock
       endif
       if (spout-spin .eq. 0.d0) then
          scr = small*cav(k)  
       else
          scr = spout-spin
       endif
       frac = (1.d0 + (spout + spin)/scr)*0.5d0
       frac = max(0.d0,min(1.d0,frac))
       
       rgdnv = frac*rstar + (1.d0 - frac)*ro
       ugdnv(k) = frac*ustar + (1.d0 - frac)*uo
       pgdnv_t = frac*pstar + (1.d0 - frac)*po
       pgdnv_g = frac*pstar_g + (1.d0 - frac)*po_g
       regdnv_g = frac*estar_g + (1.d0 - frac)*reo_g
       regdnv_r(:) = frac*estar_r(:) + (1.d0 - frac)*reo_r(:)
       
       if (spout .lt. 0.d0) then
          rgdnv = ro
          ugdnv(k) = uo
          pgdnv_t = po
          pgdnv_g = po_g
          regdnv_g = reo_g
          regdnv_r(:) = reo_r(:)
       endif
       if (spin .ge. 0.d0) then
          rgdnv = rstar
          ugdnv(k) = ustar
          pgdnv_t = pstar
          pgdnv_g = pstar_g
          regdnv_g = estar_g
          regdnv_r(:) = estar_r(:)
       endif
       
       if (k.eq.0 .and. physbc_lo(1) .eq. Symmetry) ugdnv(k) = 0.d0
       
       if (fix_mass_flux_lo .and. k.eq.domlo(1) .and. ugdnv(k) .ge. 0.d0) then
          rgdnv    = ql(k,QRHO)
          ugdnv(k) = ql(k,QU)
          regdnv_g = rel_g
          regdnv_r(:) = erl(:)
       end if
       if (fix_mass_flux_hi .and. k.eq.domhi(1)+1 .and. ugdnv(k) .le. 0.d0) then
          rgdnv    = qr(k,QRHO)
          ugdnv(k) = qr(k,QU)
          regdnv_g = rer_g
          regdnv_r(:) = err(:)
       end if
       
       do g=0, ngroups-1
          ergdnv(k,g) = max(regdnv_r(g), 0.d0)
       end do
       
       pgdnv(k) = pgdnv_g
       
       lamgdnv(k,:) = lambda
       
       ! Compute fluxes, order as conserved state (not q)
       uflx(k,URHO) = rgdnv*ugdnv(k)
       uflx(k,UMX) = uflx(k,URHO)*ugdnv(k) 
       
       rhoetot = regdnv_g + 0.5d0*rgdnv*ugdnv(k)**2 
       uflx(k,UEDEN) = ugdnv(k)*(rhoetot + pgdnv_g)
       uflx(k,UEINT) = ugdnv(k)*regdnv_g
       
       if (fspace_type.eq.1) then
          do g=0,ngroups-1
             eddf = Edd_factor(lambda(g))
             f1 = 0.5d0*(1.d0-eddf)
             rflx(k,g) = (1.d0+f1) * ergdnv(k,g) * ugdnv(k)
          end do
       else ! type 2
          do g=0,ngroups-1
             rflx(k,g) = ergdnv(k,g) * ugdnv(k)
          end do
       end if
       
       do ipassive = 1, npassive
          n = upass_map(ipassive)
          nq = qpass_map(ipassive)
          
          if (ustar .gt. ZERO) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else if (ustar .lt. ZERO) then
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          else
             qavg = HALF * (ql(k,nq) + qr(k,nq))
             uflx(k,n) = uflx(k,URHO)*qavg
          endif
       enddo
       
    enddo
  end subroutine riemannus_rad

end module riemann_rad_module

