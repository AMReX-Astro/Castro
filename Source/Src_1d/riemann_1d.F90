module riemann_module

  use bl_types
  use bl_constants_module

  implicit none

  private

  public cmpflx

  real (kind=dp_t), parameter :: smallu = 1.e-12_dp_t

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(lo,hi,domlo,domhi, &
                    qm,qp,qpd_l1,qpd_h1, &
                    flx,flx_l1,flx_h1, &
                    pgdnv,pg_l1,pg_h1, &
                    ugdnv,ug_l1,ug_h1, &
#ifdef RADIATION
                    lam, lam_l1, lam_h1, &
                    rflx,rflx_l1,rflx_h1, &
                    ergdnv,erg_l1,erg_h1, &
                    lamgdnv,lg_l1,lg_h1, &
                    gamcg, &
#endif
                    gamc,csml,c,qd_l1,qd_h1,ilo,ihi)

    use meth_params_module, only : QVAR, NVAR, &
                                   riemann_solver
#ifdef RADIATION
    use radhydro_params_module, only : QRADVAR
    use rad_params_module, only : ngroups
#endif
    
    integer lo(1),hi(1)
    integer domlo(1),domhi(1)
    integer ilo,ihi
    integer qpd_l1,qpd_h1
    integer flx_l1, flx_h1
    integer  pg_l1, pg_h1
    integer  ug_l1, ug_h1
    integer  qd_l1,  qd_h1

#ifdef RADIATION
    integer lam_l1, lam_h1
    integer rflx_l1, rflx_h1
    integer erg_l1,erg_h1
    integer  lg_l1, lg_h1
#endif

#ifdef RADIATION    
    double precision    qm(qpd_l1:qpd_h1, QRADVAR)
    double precision    qp(qpd_l1:qpd_h1, QRADVAR)
#else
    double precision    qm(qpd_l1:qpd_h1, QVAR)
    double precision    qp(qpd_l1:qpd_h1, QVAR)
#endif
    double precision   flx(flx_l1:flx_h1, NVAR)
    double precision pgdnv( pg_l1: pg_h1)
    double precision ugdnv( ug_l1: ug_h1)
    double precision  gamc( qd_l1: qd_h1)
    double precision     c( qd_l1: qd_h1)
    double precision  csml( qd_l1: qd_h1)

#ifdef RADIATION    
    double precision lam(lam_l1:lam_h1, 0:ngroups-1)
    double precision rflx(rflx_l1:rflx_h1, 0:ngroups-1)
    double precision ergdnv(erg_l1:erg_h1, 0:ngroups-1)
    double precision lamgdnv( lg_l1: lg_h1, 0:ngroups-1)
    double precision gamcg( qd_l1: qd_h1)
#endif
    
    ! Local variables
    integer i
    double precision, allocatable :: smallc(:),cavg(:),gamcp(:), gamcm(:)
#ifdef RADIATION
    double precision, allocatable :: gamcgp(:), gamcgm(:)
#endif
    
    allocate ( smallc(ilo:ihi+1) )
    allocate ( cavg(ilo:ihi+1) )
    allocate ( gamcp(ilo:ihi+1) )
    allocate ( gamcm(ilo:ihi+1) )
#ifdef RADIATION
    allocate (gamcgp(ilo:ihi+1) )
    allocate (gamcgm(ilo:ihi+1) )
#endif
    
    do i = ilo, ihi+1
       smallc(i) = max( csml(i), csml(i-1) )
       cavg(i) = HALF*( c(i) + c(i-1) )
       gamcm(i) = gamc(i-1)
       gamcp(i) = gamc(i)
#ifdef RADIATION       
       gamcgm (i) = gamcg(i-1)
       gamcgp (i) = gamcg(i)
#endif
    enddo

    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    if (riemann_solver == 0) then
       ! Colella, Glaz, & Ferguson
       call riemannus(qm, qp,qpd_l1,qpd_h1, &
                      smallc, cavg, &
                      gamcm, gamcp, &
                      flx, flx_l1, flx_h1, &
                      pgdnv, pg_l1, pg_h1, &
                      ugdnv, ug_l1, ug_h1, &
#ifdef RADIATION
                      lam, lam_l1, lam_h1, &
                      gamcgm, gamcgp, &
                      rflx, rflx_l1, rflx_h1, &
                      ergdnv,erg_l1,erg_h1, &
                      lamgdnv, lg_l1, lg_h1, &
#endif
                      ilo, ihi, domlo, domhi )

    elseif (riemann_solver == 1) then
       ! Colella & Glaz
       call riemanncg(qm, qp,qpd_l1,qpd_h1, smallc, cavg, &
                      gamcm, gamcp, flx, flx_l1, flx_h1, &
                      pgdnv, pg_l1, pg_h1, &
                      ugdnv, ug_l1, ug_h1, ilo, ihi)
    else
       call bl_error("ERROR: HLLC not support in 1-d yet")
    endif

    deallocate (smallc,cavg,gamcm,gamcp)

  end subroutine cmpflx


! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemanncg(ql,qr,qpd_l1,qpd_h1,smallc,cav, &
                       gamcl,gamcr,uflx,uflx_l1,uflx_h1, &
                       pgdnv, pg_l1, pg_h1, &
                       ugdnv, ug_l1, ug_h1,ilo,ihi)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use bl_error_module
    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QFS, &
                                   QFX, URHO, UMX, UEDEN, UEINT, &
                                   UFS, UFX, &
                                   npassive, upass_map, qpass_map, small_dens, small_pres, small_temp, &
                                   cg_maxiter, cg_tol, cg_blend
    use riemann_util_module

    double precision, parameter :: small = 1.d-8

    integer :: qpd_l1, qpd_h1
    integer :: uflx_l1, uflx_h1
    integer :: ug_l1, ug_h1
    integer :: pg_l1, pg_h1


    double precision :: ql(qpd_l1:qpd_h1, QVAR+3)
    double precision :: qr(qpd_l1:qpd_h1, QVAR+3)
    double precision ::  gamcl(ilo:ihi+1)
    double precision ::  gamcr(ilo:ihi+1)
    double precision ::    cav(ilo:ihi+1)
    double precision :: smallc(ilo:ihi+1)
    double precision :: uflx(uflx_l1:uflx_h1, NVAR)
    double precision :: ugdnv(ug_l1:ug_h1)
    double precision :: pgdnv(pg_l1:pg_h1)

    integer :: ilo,ihi
    integer :: n, nq, ipassive

    double precision :: rgdnv,ustar,gamgdnv,v1gdnv,v2gdnv
    double precision :: rl, ul, v1l, v2l, pl, rel
    double precision :: rr, ur, v1r, v2r, pr, rer
    double precision :: wl, wr, rhoetot
    double precision :: rstar, cstar, pstar
    double precision :: ro, uo, po, co, gamco
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg

    double precision :: gcl, gcr
    double precision :: clsq, clsql, clsqr, wlsq, wosq, wrsq, wo
    double precision :: zm, zp
    double precision :: denom, dpditer, dpjmp
    double precision :: gamc_bar, game_bar
    double precision :: gamel, gamer, gameo, gamstar, gmin, gmax, gdot

    integer :: iter, iter_max
    double precision :: tol
    double precision :: err

    logical :: converged

    double precision :: pstnm1
    double precision :: taul, taur, tauo
    double precision :: ustarm, ustarp, ustnm1, ustnp1
    double precision :: pstarl, pstarc, pstaru, pfuncc, pfuncu

    double precision, parameter :: weakwv = 1.d-3

    double precision, allocatable :: pstar_hist(:), pstar_hist_extra(:)

    type (eos_t) :: eos_state

    integer :: k

    if (cg_blend .eq. 2 .and. cg_maxiter < 5) then

       call bl_error("Error: need cg_maxiter >= 5 to do a bisection search on secant iteration failure.")

    endif

    tol = cg_tol
    iter_max = cg_maxiter

    allocate (pstar_hist(iter_max))
    allocate (pstar_hist_extra(2*iter_max))

    do k = ilo, ihi+1

       ! left state
       rl  = max( ql(k,QRHO), small_dens)

       ul  = ql(k,QU)
       v1l = ql(k,QV)
       v2l = ql(k,QW)

       pl  = ql(k,QPRES)
       rel = ql(k,QREINT)
       gcl = gamcl(k)

       ! sometimes we come in here with negative energy or pressure
       ! note: reset both in either case, to remain thermo
       ! consistent
       if (rel <= ZERO .or. pl <= small_pres) then
          print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres
          eos_state % T   = small_temp
          eos_state % rho = rl
          eos_state % xn  = ql(k,QFS:QFS-1+nspec)
          eos_state % aux = ql(k,QFX:QFX-1+naux)

          call eos(eos_input_rt, eos_state)

          rel = rl*eos_state%e
          pl  = eos_state%p
          gcl = eos_state%gam1
       endif

       ! right state
       rr  = max( qr(k,QRHO), small_dens)

       ur  = qr(k,QU)
       v1r  = qr(k,QV)
       v2r  = qr(k,QW)

       pr  = qr(k,QPRES)
       rer = qr(k,QREINT)
       gcr = gamcr(k)

       if (rer <= ZERO .or. pr <= small_pres) then
          print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres
          eos_state % T   = small_temp
          eos_state % rho = rr
          eos_state % xn  = qr(k,QFS:QFS-1+nspec)
          eos_state % aux = qr(k,QFX:QFX-1+naux)

          call eos(eos_input_rt, eos_state)

          rer = rr*eos_state%e
          pr  = eos_state%p
          gcr = eos_state%gam1
       endif

       ! common quantities
       taul = ONE/rl
       taur = ONE/rr

       ! lagrangian sound speeds
       clsql = gcl*pl*rl
       clsqr = gcr*pr*rr


       ! Note: in the original Colella & Glaz paper, they predicted
       ! gamma_e to the interfaces using a special (non-hyperbolic)
       ! evolution equation.  In Castro, we instead bring (rho e)
       ! to the edges, so we construct the necessary gamma_e here from
       ! what we have on the interfaces.
       gamel = pl/rel + ONE
       gamer = pr/rer + ONE

       ! these should consider a wider average of the cell-centered
       ! gammas
       gmin = min(gamel, gamer, ONE, 4.d0/3.d0)
       gmax = max(gamel, gamer, TWO, 5.d0/3.d0)

       game_bar = HALF*(gamel + gamer)
       gamc_bar = HALF*(gcl + gcr)

       gdot = TWO*(ONE - game_bar/gamc_bar)*(game_bar - ONE)

       csmall = smallc(k)
       wsmall = small_dens*csmall
       wl = max(wsmall,sqrt(abs(clsql)))
       wr = max(wsmall,sqrt(abs(clsqr)))

       ! make an initial guess for pstar -- this is a two-shock
       ! approximation
       !  pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
       pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)
       pstar = max(pstar,small_pres)

       ! get the shock speeds -- this computes W_s from CG Eq. 34
       call wsqge(pl,taul,gamel,gdot,  &
                  gamstar,pstar,wlsq,clsql,gmin,gmax)

       call wsqge(pr,taur,gamer,gdot,  &
                  gamstar,pstar,wrsq,clsqr,gmin,gmax)

       pstnm1 = pstar

       wl = sqrt(wlsq)
       wr = sqrt(wrsq)

       ! R-H jump conditions give ustar across each wave -- these should
       ! be equal when we are done iterating
       ustarp = ul - (pstar-pl)/wl
       ustarm = ur + (pstar-pr)/wr

       ! revise our pstar guess
       !pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
       pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)
       pstar = max(pstar,small_pres)

       ! secant iteration
       converged = .false.
       iter = 1
       do while ((iter <= iter_max .and. .not. converged) .or. iter <= 2)

          call wsqge(pl,taul,gamel,gdot,  &
                     gamstar,pstar,wlsq,clsql,gmin,gmax)

          call wsqge(pr,taur,gamer,gdot,  &
                     gamstar,pstar,wrsq,clsqr,gmin,gmax)

          wl = ONE / sqrt(wlsq)
          wr = ONE / sqrt(wrsq)
          ustnm1 = ustarm
          ustnp1 = ustarp
          ustarm = ur-(pr-pstar)*wr
          ustarp = ul+(pl-pstar)*wl

          dpditer=abs(pstnm1-pstar)
          zp=abs(ustarp-ustnp1)
          if(zp-weakwv*cav(k) <= ZERO)then
             zp = dpditer*wl
          endif

          zm=abs(ustarm-ustnm1)
          if(zm-weakwv*cav(k) <= ZERO)then
             zm = dpditer*wr
          endif

          ! the new pstar is found via CG Eq. 18
          denom=dpditer/max(zp+zm,small*(cav(k)))
          pstnm1 = pstar
          pstar=pstar-denom*(ustarm-ustarp)
          pstar=max(pstar,small_pres)

          err = abs(pstar - pstnm1)
          if (err < tol*pstar) converged = .true.

          pstar_hist(iter) = pstar

          iter = iter + 1

       enddo

       ! If we failed to converge using the secant iteration, we can either
       ! stop here; or, revert to the original two-shock estimate for pstar;
       ! or do a bisection root find using the bounds established by the most
       ! recent iterations.

       if (.not. converged) then

          if (cg_blend .eq. 0) then

             print *, 'pstar history: '
             do iter = 1, iter_max
                print *, iter, pstar_hist(iter)
             enddo

             print *, ' '
             print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
             print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr

             call bl_error("ERROR: non-convergence in the Riemann solver")

          else if (cg_blend .eq. 1) then

             pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)

          else if (cg_blend .eq. 2) then

             pstarl = minval(pstar_hist(iter_max-5:iter_max))
             pstaru = maxval(pstar_hist(iter_max-5:iter_max))

             iter = 1

             do while (iter <= 2 * iter_max .and. .not. converged)

                pstarc = HALF * (pstaru + pstarl)

                pstar_hist_extra(iter) = pstarc

                call wsqge(pl,taul,gamel,gdot,  &
                     gamstar,pstaru,wlsq,clsql,gmin,gmax)

                call wsqge(pr,taur,gamer,gdot,  &
                     gamstar,pstaru,wrsq,clsqr,gmin,gmax)

                wl = ONE / sqrt(wlsq)
                wr = ONE / sqrt(wrsq)

                ustarm = ur-(pr-pstar)*wr
                ustarp = ul+(pl-pstar)*wl

                pfuncc = ustarp - ustarm

                if ( HALF * (pstaru - pstarl) < cg_tol * pstarc ) then
                   converged = .true.
                endif

                iter = iter + 1

                call wsqge(pl,taul,gamel,gdot,  &
                     gamstar,pstaru,wlsq,clsql,gmin,gmax)

                call wsqge(pr,taur,gamer,gdot,  &
                     gamstar,pstaru,wrsq,clsqr,gmin,gmax)

                wl = ONE / sqrt(wlsq)
                wr = ONE / sqrt(wrsq)

                ustarm = ur-(pr-pstar)*wr
                ustarp = ul+(pl-pstar)*wl

                pfuncu = ustarp - ustarm

                if (pfuncc * pfuncu < ZERO) then

                   pstarl = pstarc

                else

                   pstaru = pstarc

                endif

             enddo

             if (.not. converged) then

                print *, 'pstar history: '
                do iter = 1, iter_max
                   print *, iter, pstar_hist(iter)
                enddo
                do iter = 1, 2 * iter_max
                   print *, iter + iter_max, pstar_hist_extra(iter)
                enddo

                print *, ' '
                print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
                print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr

                call bl_error("ERROR: non-convergence in the Riemann solver")

             endif

          else

             call bl_error("ERROR: unrecognized cg_blend option.")

          endif

       endif

       ! we converged! construct the single ustar for the region
       ! between the left and right waves, using the updated wave speeds
       ustarm = ur-(pr-pstar)*wr  ! careful -- here wl, wr are 1/W
       ustarp = ul+(pl-pstar)*wl

       ustar = HALF* ( ustarp + ustarm )

       ! for symmetry preservation, if ustar is really small, then we
       ! set it to zero
       if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
          ustar = ZERO
       endif

       ! sample the solution -- here we look first at the direction
       ! that the contact is moving.  This tells us if we need to
       ! worry about the L/L* states or the R*/R states.
       if (ustar > ZERO) then
          ro = rl
          uo = ul
          po = pl
          tauo = taul
          !reo = rel
          gamco = gcl
          gameo = gamel

       else if (ustar < ZERO) then
          ro = rr
          uo = ur
          po = pr
          tauo = taur
          !reo = rer
          gamco = gcr
          gameo = gamer

       else
          ro = HALF*(rl+rr)
          uo = HALF*(ul+ur)
          po = HALF*(pl+pr)
          tauo = HALF*(taul+taur)

          gamco = HALF*(gcl+gcr)
          gameo = HALF*(gamel + gamer)
       endif

       ! use tau = 1/rho as the independent variable here
       ro = max(small_dens, ONE/tauo)
       tauo = ONE/ro

       co = sqrt(abs(gamco*po/ro))
       co = max(csmall,co)
       clsq = (co*ro)**2

       ! now that we know which state (left or right) we need to worry
       ! about, get the value of gamstar and wosq across the wave we
       ! are dealing with.
       call wsqge(po,tauo,gameo,gdot,   &
                  gamstar,pstar,wosq,clsq,gmin,gmax)

       sgnm = sign(ONE,ustar)

       wo = sqrt(wosq)
       dpjmp = pstar - po

       ! is this max really necessary?
       !rstar=max(1.-ro*dpjmp/wosq,(gameo-1.)/(gameo+1.))
       rstar=ONE-ro*dpjmp/wosq
       rstar=ro/rstar
       rstar = max(small_dens,rstar)

       !entho = (reo/ro + po/ro)/co**2
       !estar = reo + (pstar - po)*entho

       cstar = sqrt(abs(gamco*pstar/rstar))
       cstar = max(cstar,csmall)

       spout = co - sgnm*uo
       spin = cstar - sgnm*ustar

       !ushock = HALF*(spin + spout)
       ushock = wo/ro - sgnm*uo

       if (pstar-po >= ZERO) then
          spin = ushock
          spout = ushock
       endif

       ! if (spout-spin .eq. ZERO) then
       !   scr = small*cav(k)
       !  else
       !     scr = spout-spin
       !  endif
       !  frac = (ONE + (spout + spin)/scr)*HALF
       !  frac = max(ZERO,min(ONE,frac))

       frac = HALF*(ONE + (spin + spout)/max(spout-spin,spin+spout, small*cav(k)))

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

       ! linearly interpolate between the star and normal state -- this covers the
       ! case where we are inside the rarefaction fan.
       rgdnv = frac*rstar + (ONE - frac)*ro
       ugdnv(k) = frac*ustar + (ONE - frac)*uo
       pgdnv(k) = frac*pstar + (ONE - frac)*po
       gamgdnv =  frac*gamstar + (ONE-frac)*gameo

       ! now handle the cases where instead we are fully in the
       ! star or fully in the original (l/r) state
       if (spout < ZERO) then
          rgdnv = ro
          ugdnv(k) = uo
          pgdnv(k) = po
          gamgdnv = gameo
       endif

       if (spin >= ZERO) then
          rgdnv = rstar
          ugdnv(k) = ustar
          pgdnv(k) = pstar
          gamgdnv = gamstar
       endif

       pgdnv(k) = max(pgdnv(k),small_pres)

       ! Compute fluxes, order as conserved state (not q)
       uflx(k,URHO) = rgdnv*ugdnv(k)

       ! note: here we do not include the pressure, since in 1-d,
       ! for some geometries, div{F} + grad{p} cannot be written
       ! in a flux difference form
       uflx(k,UMX) = uflx(k,URHO)*ugdnv(k)

       ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
       rhoetot = pgdnv(k)/(gamgdnv - ONE) + &
            HALF*rgdnv*(ugdnv(k)**2 + v1gdnv**2 + v2gdnv**2)

       uflx(k,UEDEN) = ugdnv(k)*(rhoetot + pgdnv(k))
       uflx(k,UEINT) = ugdnv(k)*pgdnv(k)/(gamgdnv - ONE)

       ! passively advected quantities -- only the contact matters
       ! note: this also includes the y- and z-velocity flux
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nq = qpass_map(ipassive)

          if (ustar > ZERO) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else if (ustar < ZERO) then
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          else
             qavg = HALF * (ql(k,nq) + qr(k,nq))
             uflx(k,n) = uflx(k,URHO)*qavg
          endif

       enddo

    enddo

    deallocate(pstar_hist)
    deallocate(pstar_hist_extra)

  end subroutine riemanncg

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemannus(ql,qr,qpd_l1,qpd_h1, &
                       smallc,cav, &
                       gamcl,gamcr, &
                       uflx,uflx_l1,uflx_h1,&
                       pgdnv,pg_l1,pg_h1, &
                       ugdnv,ug_l1,ug_h1, &
#ifdef RADIATION
                       lam, lam_l1, lam_h1, &
                       gamcgl, gamcgr, &
                       rflx,rflx_l1,rflx_h1, &
                       ergdnv,erg_l1,erg_h1,  &
                       lamgdnv,lg_l1,lg_h1,  &
#endif                       
                       ilo,ihi,domlo,domhi)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, &
                                   URHO, UMX, UEDEN, UEINT, small_dens, small_pres, &
                                   npassive, upass_map, qpass_map, &
                                   fix_mass_flux
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow, Symmetry
#ifdef RADIATION
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, fspace_type
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif

    double precision, parameter:: small = 1.d-8

    integer ilo,ihi
    integer domlo(1),domhi(1)
    integer  qpd_l1,  qpd_h1
    integer   pg_l1,   pg_h1
    integer   ug_l1,   ug_h1
    integer uflx_l1, uflx_h1

#ifdef RADIATION
    integer rflx_l1, rflx_h1
    integer lam_l1, lam_h1
    integer  erg_l1,  erg_h1
    integer   lg_l1,   lg_h1    
#endif

#ifdef RADIATION    
    double precision ql(qpd_l1:qpd_h1, QRADVAR)
    double precision qr(qpd_l1:qpd_h1, QRADVAR)
#else
    double precision ql(qpd_l1:qpd_h1, QVAR)
    double precision qr(qpd_l1:qpd_h1, QVAR)
#endif
    double precision   cav(ilo:ihi+1), smallc(ilo:ihi+1)
    double precision gamcl(ilo:ihi+1), gamcr(ilo:ihi+1)
    double precision  uflx(uflx_l1:uflx_h1, NVAR)
    double precision pgdnv( pg_l1: pg_h1)
    double precision ugdnv( ug_l1: ug_h1)

#ifdef RADIATION
    double precision lam(lam_l1:lam_h1, 0:ngroups-1)
    double precision gamcgl(ilo:ihi+1),gamcgr(ilo:ihi+1)
    double precision  rflx(rflx_l1:rflx_h1, 0:ngroups-1)
    double precision ergdnv(erg_l1:erg_h1, 0:ngroups-1)
    double precision lamgdnv( lg_l1: lg_h1, 0:ngroups-1)    

    double precision, dimension(0:ngroups-1) :: erl, err
    double precision :: regdnv_g, pgdnv_g, pgdnv_t
    double precision :: estar_g, pstar_g
    double precision, dimension(0:ngroups-1) :: lambda, reo_r, po_r, estar_r, regdnv_r    
    double precision :: eddf, f1
    double precision :: co_g, gamco_g, pl_g, po_g, pr_g, rel_g, reo_g, rer_g    

    integer :: g
#endif
    
    double precision rgdnv, regdnv, ustar, v1gdnv, v2gdnv
    double precision rl, ul, v1l, v2l, pl, rel
    double precision rr, ur, v1r, v2r, pr, rer
    double precision wl, wr, rhoetot, scr
    double precision rstar, cstar, estar, pstar, drho
    double precision ro, uo, po, reo, co, gamco, entho
    double precision sgnm, spin, spout, ushock, frac

    double precision wsmall, csmall
    integer ipassive, n, nq
    integer k
    logical :: fix_mass_flux_lo, fix_mass_flux_hi

    ! Solve Riemann Problem

    fix_mass_flux_lo = (fix_mass_flux == 1) .and. &
                       (physbc_lo(1) == Outflow) .and. &
                       (ilo == domlo(1))

    fix_mass_flux_hi = (fix_mass_flux == 1) .and. &
                       (physbc_hi(1) == Outflow) .and. &
                       (ihi == domhi(1))

    do k = ilo, ihi+1
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
             lambda(g) = 0.5d0*(lam(k-1,g)+lam(k,g))
          end do
          reo_r(:) = 0.5d0*(erl(:)+err(:))
          reo_g = 0.5d0*(rel_g+rer_g)
          po_r(:) = lambda(:) * reo_r(:)
          gamco_g = 0.5d0*(gamcgl(k)+gamcgr(k))
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
       ugdnv(k) = frac*ustar + (ONE - frac)*uo
#ifdef RADIATION
       pgdnv_t = frac*pstar + (1.d0 - frac)*po
       pgdnv_g = frac*pstar_g + (1.d0 - frac)*po_g
       regdnv_g = frac*estar_g + (1.d0 - frac)*reo_g
       regdnv_r(:) = frac*estar_r(:) + (1.d0 - frac)*reo_r(:)
#else
       pgdnv(k) = frac*pstar + (ONE - frac)*po
       regdnv = frac*estar + (ONE - frac)*reo
#endif
       
       if (spout < ZERO) then
          rgdnv = ro
          ugdnv(k) = uo
#ifdef RADIATION
          pgdnv_t = po
          pgdnv_g = po_g
          regdnv_g = reo_g
          regdnv_r(:) = reo_r(:)
#else
          pgdnv(k) = po
          regdnv = reo
#endif
       endif

       if (spin >= ZERO) then
          rgdnv = rstar
          ugdnv(k) = ustar
#ifdef RADIATION
          pgdnv_t = pstar
          pgdnv_g = pstar_g
          regdnv_g = estar_g
          regdnv_r(:) = estar_r(:)
#else
          pgdnv(k) = pstar
          regdnv = estar
#endif
       endif

       if (k == 0 .and. physbc_lo(1) == Symmetry) ugdnv(k) = ZERO

       if (fix_mass_flux_lo .and. k == domlo(1) .and. ugdnv(k) >= ZERO) then
          rgdnv    = ql(k,QRHO)
          ugdnv(k) = ql(k,QU)
#ifdef RADIATION
          regdnv_g = rel_g
          regdnv_r(:) = erl(:)
#else
          regdnv   = ql(k,QREINT)
#endif
       end if

       if (fix_mass_flux_hi .and. k == domhi(1)+1 .and. ugdnv(k) <= ZERO) then
          rgdnv    = qr(k,QRHO)
          ugdnv(k) = qr(k,QU)
#ifdef RADIATION
          regdnv_g = rer_g
          regdnv_r(:) = err(:)
#else
          regdnv   = qr(k,QREINT)
#endif
       end if

#ifdef RADIATION
       do g=0, ngroups-1
          ergdnv(k,g) = max(regdnv_r(g), 0.d0)
       end do
       
       pgdnv(k) = pgdnv_g
       
       lamgdnv(k,:) = lambda
#endif       
       
       ! Compute fluxes, order as conserved state (not q)

       ! Note: currently in 1-d, we do not include p in the momentum flux
       ! this is to allow for the spherical gradient
       uflx(k,URHO) = rgdnv*ugdnv(k)
       uflx(k,UMX) = uflx(k,URHO)*ugdnv(k)

#ifdef RADIATION
       rhoetot = regdnv_g + HALF*rgdnv*(ugdnv(k)**2 + v1gdnv**2 + v2gdnv**2)
       uflx(k,UEDEN) = ugdnv(k)*(rhoetot + pgdnv_g)
       uflx(k,UEINT) = ugdnv(k)*regdnv_g
       
       if (fspace_type.eq.1) then
          do g = 0, ngroups-1
             eddf = Edd_factor(lambda(g))
             f1 = 0.5d0*(1.d0-eddf)
             rflx(k,g) = (1.d0+f1) * ergdnv(k,g) * ugdnv(k)
          end do
       else ! type 2
          do g = 0, ngroups-1
             rflx(k,g) = ergdnv(k,g) * ugdnv(k)
          end do
       end if
#else       
       rhoetot = regdnv + HALF*rgdnv*(ugdnv(k)**2 + v1gdnv**2 + v2gdnv**2)
       uflx(k,UEDEN) = ugdnv(k)*(rhoetot + pgdnv(k))
       uflx(k,UEINT) = ugdnv(k)*regdnv
#endif
       
       ! advected quantities -- only the contact matters
       ! note: this includes the y,z-velocity flux       
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nq = qpass_map(ipassive)

          if (ustar >= ZERO) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          endif

       enddo

    enddo
  end subroutine riemannus

end module riemann_module
