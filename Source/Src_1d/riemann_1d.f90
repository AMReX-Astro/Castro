module riemann_module

  implicit none

  private

  public cmpflx

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine cmpflx(lo,hi,domlo,domhi, &
                    qm,qp,qpd_l1,qpd_h1, &
                    flx,flx_l1,flx_h1, &
                    pgdnv,pg_l1,pg_h1, &
                    ugdnv,ug_l1,ug_h1, &
                    gamc,csml,c,qd_l1,qd_h1,ilo,ihi)
    
    use meth_params_module, only : QVAR, NVAR
    
    implicit none
    integer lo(1),hi(1)
    integer domlo(1),domhi(1)
    integer ilo,ihi
    integer qpd_l1,qpd_h1
    integer flx_l1, flx_h1
    integer  pg_l1, pg_h1
    integer  ug_l1, ug_h1
    integer  qd_l1,  qd_h1
    double precision    qm(qpd_l1:qpd_h1, QVAR)
    double precision    qp(qpd_l1:qpd_h1, QVAR)
    double precision   flx(flx_l1:flx_h1, NVAR)
    double precision pgdnv( pg_l1: pg_h1)
    double precision ugdnv( ug_l1: ug_h1)
    double precision  gamc( qd_l1: qd_h1)
    double precision     c( qd_l1: qd_h1)
    double precision  csml( qd_l1: qd_h1)
    
    ! Local variables
    integer i
    double precision, allocatable :: smallc(:),cavg(:),gamcp(:), gamcm(:)
    
    allocate ( smallc(ilo:ihi+1) )
    allocate ( cavg(ilo:ihi+1) )
    allocate ( gamcp(ilo:ihi+1) )
    allocate ( gamcm(ilo:ihi+1) )
    
    do i = ilo, ihi+1 
       smallc(i) = max( csml(i), csml(i-1) )
       cavg(i) = 0.5d0*( c(i) + c(i-1) )
       gamcm(i) = gamc(i-1)
       gamcp(i) = gamc(i)
    enddo
    
    ! Solve Riemann problem (gdnv state passed back, but only (u,p) saved)
    call riemannus(qm, qp,qpd_l1,qpd_h1, smallc, cavg, &
                   gamcm, gamcp, flx, flx_l1, flx_h1, &
                   pgdnv, pg_l1, pg_h1, &
                   ugdnv, ug_l1, ug_h1, ilo, ihi, domlo, domhi )
    
    deallocate (smallc,cavg,gamcp,gamcm)
    
  end subroutine cmpflx


! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine riemannus(ql,qr,qpd_l1,qpd_h1,smallc,cav, &
                       gamcl,gamcr,uflx,uflx_l1,uflx_h1,&
                       pgdnv,pg_l1,pg_h1,ugdnv,ug_l1,ug_h1, &
                       ilo,ihi,domlo,domhi)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QPRES, QREINT, QFA, QFS, QFX, &
                                   URHO, UMX, UEDEN, UEINT, &
                                   UFA, UFS, UFX, nadv, small_dens, small_pres, &
                                   fix_mass_flux
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow, Symmetry

    implicit none

    double precision, parameter:: small = 1.d-8

    integer ilo,ihi
    integer domlo(1),domhi(1)
    integer  qpd_l1,  qpd_h1
    integer   pg_l1,   pg_h1
    integer   ug_l1,   ug_h1
    integer uflx_l1, uflx_h1
    double precision ql(qpd_l1:qpd_h1, QVAR)
    double precision qr(qpd_l1:qpd_h1, QVAR)
    double precision   cav(ilo:ihi+1), smallc(ilo:ihi+1)
    double precision gamcl(ilo:ihi+1), gamcr(ilo:ihi+1)
    double precision  uflx(uflx_l1:uflx_h1, NVAR)
    double precision pgdnv( pg_l1: pg_h1)
    double precision ugdnv( ug_l1: ug_h1)
    
    double precision rgdnv, regdnv, ustar
    double precision rl, ul, pl, rel
    double precision rr, ur, pr, rer
    double precision wl, wr, rhoetot, scr
    double precision rstar, cstar, estar, pstar
    double precision ro, uo, po, reo, co, gamco, entho
    double precision sgnm, spin, spout, ushock, frac
    
    double precision wsmall, csmall
    integer iadv, n, nq
    integer k,ispec, iaux
    logical :: fix_mass_flux_lo, fix_mass_flux_hi
    
    ! Solve Riemann Problem

    fix_mass_flux_lo = (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) .and. (ilo .eq. domlo(1))
    fix_mass_flux_hi = (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) .and. (ihi .eq. domhi(1))

    do k = ilo, ihi+1
       rl  = ql(k,QRHO)
       ul  = ql(k,QU)
       pl  = ql(k,QPRES)
       rel = ql(k,QREINT)
       rr  = qr(k,QRHO)
       ur  = qr(k,QU)
       pr  = qr(k,QPRES)
       rer = qr(k,QREINT)
       
       csmall = smallc(k)
       wsmall = small_dens*csmall
       wl = max(wsmall,sqrt(abs(gamcl(k)*pl*rl)))
       wr = max(wsmall,sqrt(abs(gamcr(k)*pr*rr)))
       pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
       pstar = max(pstar,small_pres)
       ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)
       
       if (ustar .gt. 0.d0) then
          ro = rl
          uo = ul
          po = pl
          reo = rel
          gamco = gamcl(k)
       else if (ustar .lt. 0.d0) then
          ro = rr
          uo = ur
          po = pr
          reo = rer
          gamco = gamcr(k)
       else
          ro = 0.5d0*(rl+rr)
          uo = 0.5d0*(ul+ur)
          po = 0.5d0*(pl+pr)
          reo = 0.5d0*(rel+rer)
          gamco = 0.5d0*(gamcl(k)+gamcr(k))
       endif
       ro = max(small_dens,ro)
       
       co = sqrt(abs(gamco*po/ro))
       co = max(csmall,co)
       entho = (reo/ro + po/ro)/co**2
       rstar = ro + (pstar - po)/co**2
       rstar = max(small_dens,rstar)
       estar = reo + (pstar - po)*entho
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
       pgdnv(k) = frac*pstar + (1.d0 - frac)*po
       regdnv = frac*estar + (1.d0 - frac)*reo
       
       if (spout .lt. 0.d0) then
          rgdnv = ro
          ugdnv(k) = uo
          pgdnv(k) = po
          regdnv = reo
       endif
       if (spin .ge. 0.d0) then
          rgdnv = rstar
          ugdnv(k) = ustar
          pgdnv(k) = pstar
          regdnv = estar
       endif
       
       if (k.eq.0 .and. physbc_lo(1) .eq. Symmetry) ugdnv(k) = 0.d0
       
       if (fix_mass_flux_lo .and. k.eq.domlo(1) .and. ugdnv(k) .ge. 0.d0) then
          rgdnv    = ql(k,QRHO)
          ugdnv(k) = ql(k,QU)
          regdnv   = ql(k,QREINT)
       end if
       if (fix_mass_flux_hi .and. k.eq.domhi(1)+1 .and. ugdnv(k) .le. 0.d0) then
          rgdnv    = qr(k,QRHO)
          ugdnv(k) = qr(k,QU)
          regdnv   = qr(k,QREINT)
       end if

       ! Compute fluxes, order as conserved state (not q)
       uflx(k,URHO) = rgdnv*ugdnv(k)
       uflx(k,UMX) = uflx(k,URHO)*ugdnv(k) 
       
       rhoetot = regdnv + 0.5d0*rgdnv*ugdnv(k)**2 
       uflx(k,UEDEN) = ugdnv(k)*(rhoetot + pgdnv(k))
       uflx(k,UEINT) = ugdnv(k)*regdnv
       
       do iadv = 1, nadv
          n = UFA + iadv - 1
          nq = QFA + iadv - 1
          if (ustar .ge. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          endif
       enddo
       
       do ispec = 1, nspec
          n  = UFS + ispec - 1
          nq = QFS + ispec - 1
          if (ustar .ge. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          endif
       enddo
       
       do iaux = 1, naux
          n  = UFX + iaux - 1
          nq = QFX + iaux - 1
          if (ustar .ge. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          endif
       enddo
       
    enddo
  end subroutine riemannus

end module riemann_module
