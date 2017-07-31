module riemann_module

  use bl_types
  use bl_constants_module

  use meth_params_module, only : NQ, NQAUX, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, &
                                 QFS, QFX, &
                                 NGDNV, GDU, GDPRES, QGAMC, QC, QCSML, &
#ifdef RADIATION
                                 GDERADS, GDLAMS, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UEDEN, UEINT, &
                                 small_temp, small_dens, small_pres, &
                                 npassive, upass_map, qpass_map, &
                                 cg_maxiter, cg_tol, cg_blend, &
                                 fix_mass_flux, &
                                 riemann_solver, ppm_temp_fix, hybrid_riemann, &
                                 allow_negative_energy

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public cmpflx

  real(rt), parameter :: smallu = 1.e-12_rt

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(lo, hi, domlo, domhi, &
                    qm, qp, qpd_lo, qpd_hi, &
                    flx, flx_lo, flx_hi, &
                    qint, qg_lo, qg_hi, &
#ifdef RADIATION
                    rflx, rflx_lo, rflx_hi, &
#endif
                    qaux, qa_lo, qa_hi, &
                    ilo, ihi)


    use eos_type_module, only: eos_input_re, eos_input_rt, eos_t
    use eos_module, only: eos
    use network, only: nspec, naux
#ifdef RADIATION
    use rad_params_module, only : ngroups
#endif
    use actual_riemann_module, only : riemanncg

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer lo(1),hi(1)
    integer domlo(1),domhi(1)
    integer ilo,ihi
    integer qpd_lo(3), qpd_hi(3)
    integer flx_lo(3), flx_hi(3)
    integer  qg_lo(3), qg_hi(3)
    integer  qa_lo(3), qa_hi(3)

    real(rt)            qm(qpd_lo(1):qpd_hi(1), NQ)
    real(rt)            qp(qpd_lo(1):qpd_hi(1), NQ)

    real(rt)           flx(flx_lo(1):flx_hi(1), NVAR)
    real(rt)          qint( qg_lo(1): qg_hi(1), NGDNV)
    real(rt)          qaux( qa_lo(1): qa_hi(1), NQAUX)

#ifdef RADIATION
    integer rflx_lo(3), rflx_hi(3)
    real(rt)         rflx(rflx_lo(1):rflx_hi(1), 0:ngroups-1)
#endif

    ! Local variables
    integer i
    real(rt)        , allocatable :: smallc(:), cavg(:), gamcp(:), gamcm(:)
#ifdef RADIATION
    real(rt)        , allocatable :: gamcgp(:), gamcgm(:), lam(:,:)
#endif

    type (eos_t) :: eos_state

    allocate ( smallc(ilo:ihi+1) )
    allocate ( cavg(ilo:ihi+1) )
    allocate ( gamcp(ilo:ihi+1) )
    allocate ( gamcm(ilo:ihi+1) )
#ifdef RADIATION
    allocate (gamcgp(ilo:ihi+1) )
    allocate (gamcgm(ilo:ihi+1) )
    allocate (lam(ilo-1:ihi+2,0:ngroups-1) )
#endif

    do i = ilo, ihi+1
       smallc(i) = max( qaux(i,QCSML), qaux(i-1,QCSML) )
       cavg(i) = HALF*( qaux(i,QC) + qaux(i-1,QC) )
       gamcm(i) = qaux(i-1,QGAMC)
       gamcp(i) = qaux(i,QGAMC)
#ifdef RADIATION
       gamcgm (i) = qaux(i-1,QGAMCG)
       gamcgp (i) = qaux(i,QGAMCG)
#endif
    enddo

#ifdef RADIATION
    do i = ilo-1, ihi+2
       lam(i,:) = qaux(i,QLAMS:QLAMS+ngroups-1)
    enddo
#endif


    if (ppm_temp_fix == 2) then
       ! recompute the thermodynamics on the interface to make it
       ! all consistent -- THIS PROBABLY DOESN"T WORK WITH RADIATION

       ! we want to take the edge states of rho, p, and X, and get
       ! new values for gamc and (rho e) on the edges that are
       ! thermodynamically consistent.

       do i = ilo, ihi+1

          ! this is an initial guess for iterations, since we
          ! can't be certain that temp is on interfaces
          eos_state%T = 10000.0e0_rt

          ! minus state
          eos_state % rho = qm(i,QRHO)
          eos_state % p   = qm(i,QPRES)
          eos_state % e   = qm(i,QREINT)/qm(i,QRHO)
          eos_state % xn  = qm(i,QFS:QFS-1+nspec)
          eos_state % aux = qm(i,QFX:QFX-1+naux)

          ! Protect against negative energies

          if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
             eos_state % T = small_temp
             call eos(eos_input_rt, eos_state)
          else
             call eos(eos_input_re, eos_state)
          endif

          qm(i,QREINT) = qm(i,QRHO)*eos_state%e
          qm(i,QPRES) = eos_state%p
          gamcm(i) = eos_state%gam1


          ! plus state
          eos_state % rho = qp(i,QRHO)
          eos_state % p   = qp(i,QPRES)
          eos_state % e   = qp(i,QREINT)/qp(i,QRHO)
          eos_state % xn  = qp(i,QFS:QFS-1+nspec)
          eos_state % aux = qp(i,QFX:QFX-1+naux)
          
          ! Protect against negative energies

          if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
             eos_state % T = small_temp
             call eos(eos_input_rt, eos_state)
          else
             call eos(eos_input_re, eos_state)
          endif
          
          qp(i,QREINT) = qp(i,QRHO)*eos_state%e
          qp(i,QPRES) = eos_state%p
          gamcp(i) = eos_state%gam1
          
       enddo

    endif


    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    if (riemann_solver == 0) then
       ! Colella, Glaz, & Ferguson
       call riemannus(qm,  qp, qpd_lo, qpd_hi, &
                      smallc, cavg, &
                      gamcm, gamcp, &
                      flx, flx_lo, flx_hi, &
                      qint, qg_lo, qg_hi, &
#ifdef RADIATION
                      lam, gamcgm, gamcgp, &
                      rflx, rflx_lo, rflx_hi, &
#endif
                      ilo, ihi+1, domlo, domhi )

    elseif (riemann_solver == 1) then
       ! Colella & Glaz
       call riemanncg(qm, qp, qpd_lo, qpd_hi, &
                      gamcm, gamcp, cavg, smallc, [ilo, 0, 0], [ihi+1, 0, 0], &
                      flx, flx_lo, flx_hi, &
                      qint, qg_lo, qg_hi, &
                      1, ilo, ihi+1, 0, 0, 0, 0, 0, &
                      [domlo(1), 0, 0], [domhi(1), 0, 0])
    else
       call bl_error("ERROR: HLLC not support in 1-d yet")
    endif

    deallocate (smallc,cavg,gamcm,gamcp)
#ifdef RADIATION
    deallocate(gamcgm,gamcgp,lam)
#endif

  end subroutine cmpflx


! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemannus(ql, qr, qpd_lo, qpd_hi, &
                       smallc, cav, &
                       gamcl, gamcr, &
                       uflx, uflx_lo, uflx_hi, &
                       qint, qg_lo, qg_hi, &
#ifdef RADIATION
                       lam, gamcgl, gamcgr, &
                       rflx, rflx_lo, rflx_hi, &
#endif
                       ilo, ihi, domlo, domhi)

    use prob_params_module, only : physbc_lo, physbc_hi, Outflow, Symmetry
#ifdef RADIATION
    use meth_params_module, only : qrad, qradhi, qptot, qreitot, fspace_type
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif

    use amrex_fort_module, only : rt => amrex_real
    real(rt)        , parameter:: small = 1.e-8_rt

    integer ilo,ihi
    integer domlo(1),domhi(1)
    integer  qpd_lo(3), qpd_hi(3)
    integer   qg_lo(3), qg_hi(3)
    integer uflx_lo(3), uflx_hi(3)

#ifdef RADIATION
    integer rflx_lo(3), rflx_hi(3)
#endif

    real(rt)         ql(qpd_lo(1):qpd_hi(1), NQ)
    real(rt)         qr(qpd_lo(1):qpd_hi(1), NQ)

    real(rt)           cav(ilo:ihi), smallc(ilo:ihi)
    real(rt)         gamcl(ilo:ihi), gamcr(ilo:ihi)
    real(rt)          uflx(uflx_lo(1):uflx_hi(1), NVAR)
    real(rt)          qint( qg_lo(1): qg_hi(1), NGDNV)

#ifdef RADIATION
    real(rt)         lam(ilo-1:ihi+1, 0:ngroups-1)
    real(rt)         gamcgl(ilo:ihi),gamcgr(ilo:ihi)
    real(rt)          rflx(rflx_lo(1):rflx_hi(1), 0:ngroups-1)

    real(rt)        , dimension(0:ngroups-1) :: erl, err
    real(rt)         :: regdnv_g, pgdnv_g, pgdnv_t
    real(rt)         :: estar_g, pstar_g
    real(rt)        , dimension(0:ngroups-1) :: lambda, reo_r, po_r, estar_r, regdnv_r
    real(rt)         :: eddf, f1
    real(rt)         :: co_g, gamco_g, pl_g, po_g, pr_g, rel_g, reo_g, rer_g

    integer :: g
#endif

    real(rt)         rgdnv, regdnv, ustar, v1gdnv, v2gdnv
    real(rt)         rl, ul, v1l, v2l, pl, rel
    real(rt)         rr, ur, v1r, v2r, pr, rer
    real(rt)         wl, wr, rhoetot, scr
    real(rt)         rstar, cstar, estar, pstar, drho
    real(rt)         ro, uo, po, reo, co, gamco, entho
    real(rt)         sgnm, spin, spout, ushock, frac

    real(rt)         wsmall, csmall
    integer ipassive, n, nqp
    integer k
    logical :: fix_mass_flux_lo, fix_mass_flux_hi

    ! Solve Riemann Problem

    fix_mass_flux_lo = (fix_mass_flux == 1) .and. &
                       (physbc_lo(1) == Outflow) .and. &
                       (ilo == domlo(1))

    fix_mass_flux_hi = (fix_mass_flux == 1) .and. &
                       (physbc_hi(1) == Outflow) .and. &
                       (ihi == domhi(1)+1)

    do k = ilo, ihi
       rl  = ql(k,QRHO)
       ul  = ql(k,QU)
       v1l = ql(k,QV)
       v2l = ql(k,QW)

#ifdef RADIATION
       pl  = ql(k,QPTOT)
       rel = ql(k,QREITOT)

       erl(:) = ql(k,qrad:qradhi)
       pl_g = ql(k,QPRES)
       rel_g = ql(k,QREINT)
#else
       pl  = ql(k,QPRES)
       rel = ql(k,QREINT)
#endif

       rr  = qr(k,QRHO)
       ur  = qr(k,QU)
       v1r  = qr(k,QV)
       v2r  = qr(k,QW)

#ifdef RADIATION
       pr  = qr(k,QPTOT)
       rer = qr(k,QREITOT)

       err(:) = qr(k,qrad:qradhi)
       pr_g = qr(k,QPRES)
       rer_g = qr(k,QREINT)
#else
       pr  = qr(k,QPRES)
       rer = qr(k,QREINT)
#endif

       csmall = smallc(k)
       wsmall = small_dens*csmall
       wl = max(wsmall,sqrt(abs(gamcl(k)*pl*rl)))
       wr = max(wsmall,sqrt(abs(gamcr(k)*pr*rr)))

       pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
       ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

       pstar = max(pstar,small_pres)
       ! for symmetry preservation, if ustar is really small, then we
       ! set it to zero
       if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
          ustar = ZERO
       endif

       if (ustar > ZERO) then
          ro = rl
          uo = ul
          po = pl
          reo = rel
          gamco = gamcl(k)

#ifdef RADIATION
          lambda(:) = lam(k-1,:)
          po_g = pl_g
          po_r(:) = erl(:) * lambda(:)
          reo_r(:) = erl(:)
          reo_g = rel_g
          gamco_g = gamcgl(k)
#endif

       else if (ustar < ZERO) then
          ro = rr
          uo = ur
          po = pr
          reo = rer
          gamco = gamcr(k)

#ifdef RADIATION
          lambda(:) = lam(k,:)
          po_g = pr_g
          po_r(:) = err(:) * lambda(:)
          reo_r(:) = err(:)
          reo_g = rer_g
          gamco_g = gamcgr(k)
#endif

       else
          ro = HALF*(rl+rr)
          uo = HALF*(ul+ur)
          po = HALF*(pl+pr)
          reo = HALF*(rel+rer)
          gamco = HALF*(gamcl(k)+gamcr(k))

#ifdef RADIATION
          do g=0, ngroups-1
             lambda(g) = 0.5e0_rt*(lam(k-1,g)+lam(k,g))
          end do
          reo_r(:) = 0.5e0_rt*(erl(:)+err(:))
          reo_g = 0.5e0_rt*(rel_g+rer_g)
          po_r(:) = lambda(:) * reo_r(:)
          gamco_g = 0.5e0_rt*(gamcgl(k)+gamcgr(k))
          po_g = 0.5*(pr_g+pl_g)
#endif
       endif

       ro = max(small_dens,ro)

       co = sqrt(abs(gamco*po/ro))
       co = max(csmall,co)

       drho = (pstar - po)/co**2
       rstar = ro + drho
       rstar = max(small_dens,rstar)

#ifdef RADIATION
       estar_g = reo_g + drho*(reo_g + po_g)/ro

       co_g = sqrt(abs(gamco_g*po_g/ro))
       co_g = max(csmall,co_g)

       pstar_g = po_g + drho*co_g**2
       pstar_g = max(pstar_g,small_pres)

       estar_r(:) = reo_r(:) + drho*(reo_r(:) + po_r(:))/ro
#else
       entho = (reo/ro + po/ro)/co**2
       estar = reo + (pstar - po)*entho
#endif

       cstar = sqrt(abs(gamco*pstar/rstar))
       cstar = max(cstar,csmall)

       sgnm = sign(ONE, ustar)
       spout = co - sgnm*uo
       spin = cstar - sgnm*ustar
       ushock = HALF*(spin + spout)

       if (pstar-po >= ZERO) then
          spin = ushock
          spout = ushock
       endif

       if (spout-spin == ZERO) then
          scr = small*cav(k)
       else
          scr = spout-spin
       endif

       frac = (ONE + (spout + spin)/scr)*HALF
       frac = max(ZERO,min(ONE,frac))

       ! transverse velocities just go along for the ride
       if (ustar > ZERO) then
          v1gdnv = v1l
          v2gdnv = v2l
       else if (ustar < ZERO) then
          v1gdnv = v1r
          v2gdnv = v2r
       else
          v1gdnv = HALF*(v1l+v1r)
          v2gdnv = HALF*(v2l+v2r)
       endif

       rgdnv = frac*rstar + (ONE - frac)*ro
       qint(k,GDU) = frac*ustar + (ONE - frac)*uo
#ifdef RADIATION
       pgdnv_t = frac*pstar + (1.e0_rt - frac)*po
       pgdnv_g = frac*pstar_g + (1.e0_rt - frac)*po_g
       regdnv_g = frac*estar_g + (1.e0_rt - frac)*reo_g
       regdnv_r(:) = frac*estar_r(:) + (1.e0_rt - frac)*reo_r(:)
#else
       qint(k,GDPRES) = frac*pstar + (ONE - frac)*po
       regdnv = frac*estar + (ONE - frac)*reo
#endif

       if (spout < ZERO) then
          rgdnv = ro
          qint(k,GDU) = uo
#ifdef RADIATION
          pgdnv_t = po
          pgdnv_g = po_g
          regdnv_g = reo_g
          regdnv_r(:) = reo_r(:)
#else
          qint(k,GDPRES) = po
          regdnv = reo
#endif
       endif

       if (spin >= ZERO) then
          rgdnv = rstar
          qint(k,GDU) = ustar
#ifdef RADIATION
          pgdnv_t = pstar
          pgdnv_g = pstar_g
          regdnv_g = estar_g
          regdnv_r(:) = estar_r(:)
#else
          qint(k,GDPRES) = pstar
          regdnv = estar
#endif
       endif

       if (k == 0 .and. physbc_lo(1) == Symmetry) qint(k,GDU) = ZERO

       if (fix_mass_flux_lo .and. k == domlo(1) .and. qint(k,GDU) >= ZERO) then
          rgdnv    = ql(k,QRHO)
          qint(k,GDU) = ql(k,QU)
#ifdef RADIATION
          regdnv_g = rel_g
          regdnv_r(:) = erl(:)
#else
          regdnv   = ql(k,QREINT)
#endif
       end if

       if (fix_mass_flux_hi .and. k == domhi(1)+1 .and. qint(k,GDU) <= ZERO) then
          rgdnv    = qr(k,QRHO)
          qint(k,GDU) = qr(k,QU)
#ifdef RADIATION
          regdnv_g = rer_g
          regdnv_r(:) = err(:)
#else
          regdnv   = qr(k,QREINT)
#endif
       end if

#ifdef RADIATION
       do g=0, ngroups-1
          qint(k,GDERADS+g) = max(regdnv_r(g), ZERO)
       end do

       qint(k,GDPRES) = pgdnv_g

       qint(k,GDLAMS:GDLAMS-1+ngroups) = lambda(:)
#endif

       ! Compute fluxes, order as conserved state (not q)

       ! Note: currently in 1-d, we do not include p in the momentum flux
       ! this is to allow for the spherical gradient
       uflx(k,URHO) = rgdnv*qint(k,GDU)
       uflx(k,UMX) = uflx(k,URHO)*qint(k,GDU)

#ifdef RADIATION
       rhoetot = regdnv_g + HALF*rgdnv*(qint(k,GDU)**2 + v1gdnv**2 + v2gdnv**2)
       uflx(k,UEDEN) = qint(k,GDU)*(rhoetot + pgdnv_g)
       uflx(k,UEINT) = qint(k,GDU)*regdnv_g

       if (fspace_type.eq.1) then
          do g = 0, ngroups-1
             eddf = Edd_factor(lambda(g))
             f1 = 0.5e0_rt*(1.e0_rt-eddf)
             rflx(k,g) = (1.e0_rt+f1) * qint(k,GDERADS+g) * qint(k,GDU)
          end do
       else ! type 2
          do g = 0, ngroups-1
             rflx(k,g) = qint(k,GDERADS+g) * qint(k,GDU)
          end do
       end if
#else
       rhoetot = regdnv + HALF*rgdnv*(qint(k,GDU)**2 + v1gdnv**2 + v2gdnv**2)
       uflx(k,UEDEN) = qint(k,GDU)*(rhoetot + qint(k,GDPRES))
       uflx(k,UEINT) = qint(k,GDU)*regdnv
#endif

       ! advected quantities -- only the contact matters
       ! note: this includes the y,z-velocity flux
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nqp = qpass_map(ipassive)

          if (ustar >= ZERO) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nqp)
          else
             uflx(k,n) = uflx(k,URHO)*qr(k,nqp)
          endif

       enddo

    enddo
  end subroutine riemannus

end module riemann_module
