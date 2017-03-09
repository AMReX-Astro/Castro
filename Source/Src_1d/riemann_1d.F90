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

  use amrex_fort_module, only : rt => c_real
  implicit none

  private

  public cmpflx, shock

  real(rt)        , parameter :: smallu = 1.e-12_rt

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(lo, hi, domlo, domhi, &
                    qm, qp, qpd_l1, qpd_h1, &
                    flx, flx_l1, flx_h1, &
                    qint, qg_l1, qg_h1, &
#ifdef RADIATION
                    rflx, rflx_l1, rflx_h1, &
#endif
                    qaux, qa_l1, qa_h1, ilo, ihi)


    use eos_type_module
    use eos_module
    use network, only: nspec, naux

#ifdef RADIATION
    use rad_params_module, only : ngroups
#endif

    use amrex_fort_module, only : rt => c_real
    integer lo(1),hi(1)
    integer domlo(1),domhi(1)
    integer ilo,ihi
    integer qpd_l1,qpd_h1
    integer flx_l1, flx_h1
    integer  qg_l1,  qg_h1
    integer  qa_l1,  qa_h1

    real(rt)            qm(qpd_l1:qpd_h1, NQ)
    real(rt)            qp(qpd_l1:qpd_h1, NQ)

    real(rt)           flx(flx_l1:flx_h1, NVAR)
    real(rt)          qint( qg_l1: qg_h1, NGDNV)
    real(rt)          qaux( qa_l1: qa_h1, NQAUX)

#ifdef RADIATION
    integer rflx_l1, rflx_h1
    real(rt)         rflx(rflx_l1:rflx_h1, 0:ngroups-1)
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
       call riemannus(qm, qp,qpd_l1,qpd_h1, &
                      smallc, cavg, &
                      gamcm, gamcp, &
                      flx, flx_l1, flx_h1, &
                      qint, qg_l1, qg_h1, &
#ifdef RADIATION
                      lam, gamcgm, gamcgp, &
                      rflx, rflx_l1, rflx_h1, &
#endif
                      ilo, ihi, domlo, domhi )

    elseif (riemann_solver == 1) then
       ! Colella & Glaz
       call riemanncg(qm, qp,qpd_l1,qpd_h1, smallc, cavg, &
                      gamcm, gamcp, flx, flx_l1, flx_h1, &
                      qint, qg_l1, qg_h1, &
                      ilo, ihi)
    else
       call bl_error("ERROR: HLLC not support in 1-d yet")
    endif

    deallocate (smallc,cavg,gamcm,gamcp)
#ifdef RADIATION
    deallocate(gamcgm,gamcgp,lam)
#endif

  end subroutine cmpflx


  subroutine shock(q,qd_l1,qd_h1, &
                   shk,s_l1,s_h1, &
                   ilo1,ihi1,dx)

    use prob_params_module, only : coord_type

    use amrex_fort_module, only : rt => c_real
    integer, intent(in) :: qd_l1, qd_h1
    integer, intent(in) :: s_l1, s_h1
    integer, intent(in) :: ilo1, ihi1
    real(rt)        , intent(in) :: dx
    real(rt)        , intent(in) :: q(qd_l1:qd_h1,NQ)
    real(rt)        , intent(inout) :: shk(s_l1:s_h1)

    integer :: i

    real(rt)         :: divU
    real(rt)         :: px_pre, px_post
    real(rt)         :: e_x, d
    real(rt)         :: p_pre, p_post, pjump

    real(rt)         :: rc, rm, rp

    real(rt)        , parameter :: small = 1.e-10_rt
    real(rt)        , parameter :: eps = 0.33e0_rt

    ! This is a basic multi-dimensional shock detection algorithm.
    ! This implementation follows Flash, which in turn follows
    ! AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
    !
    ! The spirit of this follows the shock detection in Colella &
    ! Woodward (1984)

    do i = ilo1-1, ihi1+1

       ! construct div{U}
       if (coord_type == 0) then
          divU = HALF*(q(i+1,QU) - q(i-1,QU))/dx

       else if (coord_type == 1) then
          ! cylindrical r
          rc = dble(i + HALF)*dx
          rm = dble(i - 1 + HALF)*dx
          rp = dble(i + 1 + HALF)*dx

          divU = HALF*(rp*q(i+1,QU) - rm*q(i-1,QU))/(rc*dx)
       else if (coord_type == 2) then
          rc = dble(i + HALF)*dx
          rm = dble(i - 1 + HALF)*dx
          rp = dble(i + 1 + HALF)*dx

          divU = HALF*(rp**2*q(i+1,QU) - rm**2*q(i-1,QU))/(rc**2*dx)
       else
          call bl_error("ERROR: invalid coord_type in shock")
       endif

       ! find the pre- and post-shock pressures in each direction
       if (q(i+1,QPRES) - q(i-1,QPRES) < ZERO) then
          px_pre  = q(i+1,QPRES)
          px_post = q(i-1,QPRES)
       else
          px_pre  = q(i-1,QPRES)
          px_post = q(i+1,QPRES)
       endif

       ! project the pressures onto the shock direction (trivial in 1-d)
       p_pre  = px_pre
       p_post = px_post

       ! test for compression + pressure jump to flag a shock
       if (p_pre == ZERO) then
          ! this can arise if e_x = e_y = 0 (U = 0)
          pjump = ZERO
       else
          pjump = eps - (p_post - p_pre)/p_pre
       endif

       if (pjump < ZERO .and. divU < ZERO) then
          shk(i) = ONE
       else
          shk(i) = ZERO
       endif

    enddo

  end subroutine shock




! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemanncg(ql,qr,qpd_l1,qpd_h1,smallc,cav, &
                       gamcl,gamcr,uflx,uflx_l1,uflx_h1, &
                       qint,qg_l1,qg_h1, &
                       ilo,ihi)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use bl_error_module
    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use riemann_util_module

    use amrex_fort_module, only : rt => c_real
    real(rt)        , parameter :: small = 1.e-8_rt

    integer :: qpd_l1, qpd_h1
    integer :: uflx_l1, uflx_h1
    integer :: qg_l1, qg_h1

    real(rt)         :: ql(qpd_l1:qpd_h1, NQ)
    real(rt)         :: qr(qpd_l1:qpd_h1, NQ)
    real(rt)         ::  gamcl(ilo:ihi+1)
    real(rt)         ::  gamcr(ilo:ihi+1)
    real(rt)         ::    cav(ilo:ihi+1)
    real(rt)         :: smallc(ilo:ihi+1)
    real(rt)         :: uflx(uflx_l1:uflx_h1, NVAR)
    real(rt)         ::  qint(qg_l1:qg_h1, NGDNV)

    integer :: ilo,ihi
    integer :: n, nqp, ipassive

    real(rt)         :: rgdnv,ustar,gamgdnv,v1gdnv,v2gdnv
    real(rt)         :: rl, ul, v1l, v2l, pl, rel
    real(rt)         :: rr, ur, v1r, v2r, pr, rer
    real(rt)         :: wl, wr, rhoetot
    real(rt)         :: rstar, cstar, pstar
    real(rt)         :: ro, uo, po, co, gamco
    real(rt)         :: sgnm, spin, spout, ushock, frac
    real(rt)         :: wsmall, csmall,qavg

    real(rt)         :: gcl, gcr
    real(rt)         :: clsq, clsql, clsqr, wlsq, wosq, wrsq, wo
    real(rt)         :: zm, zp
    real(rt)         :: denom, dpditer, dpjmp
    real(rt)         :: gamc_bar, game_bar
    real(rt)         :: gamel, gamer, gameo, gamstar, gmin, gmax, gdot

    integer :: iter, iter_max
    real(rt)         :: tol
    real(rt)         :: err

    logical :: converged

    real(rt)         :: pstnm1
    real(rt)         :: taul, taur, tauo
    real(rt)         :: ustarm, ustarp, ustnm1, ustnp1
    real(rt)         :: pstarl, pstarc, pstaru, pfuncc, pfuncu

    real(rt)        , parameter :: weakwv = 1.e-3_rt

    real(rt)        , allocatable :: pstar_hist(:), pstar_hist_extra(:)

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
       gmin = min(gamel, gamer, ONE, 4.e0_rt/3.e0_rt)
       gmax = max(gamel, gamer, TWO, 5.e0_rt/3.e0_rt)

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
       qint(k,GDU) = frac*ustar + (ONE - frac)*uo
       qint(k,GDPRES) = frac*pstar + (ONE - frac)*po
       gamgdnv =  frac*gamstar + (ONE-frac)*gameo

       ! now handle the cases where instead we are fully in the
       ! star or fully in the original (l/r) state
       if (spout < ZERO) then
          rgdnv = ro
          qint(k,GDU) = uo
          qint(k,GDPRES) = po
          gamgdnv = gameo
       endif

       if (spin >= ZERO) then
          rgdnv = rstar
          qint(k,GDU) = ustar
          qint(k,GDPRES) = pstar
          gamgdnv = gamstar
       endif

       qint(k,GDPRES) = max(qint(k,GDPRES),small_pres)

       ! Compute fluxes, order as conserved state (not q)
       uflx(k,URHO) = rgdnv*qint(k,GDU)

       ! note: here we do not include the pressure, since in 1-d,
       ! for some geometries, div{F} + grad{p} cannot be written
       ! in a flux difference form
       uflx(k,UMX) = uflx(k,URHO)*qint(k,GDU)

       ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
       rhoetot = qint(k,GDPRES)/(gamgdnv - ONE) + &
            HALF*rgdnv*(qint(k,GDU)**2 + v1gdnv**2 + v2gdnv**2)

       uflx(k,UEDEN) = qint(k,GDU)*(rhoetot + qint(k,GDPRES))
       uflx(k,UEINT) = qint(k,GDU)*qint(k,GDPRES)/(gamgdnv - ONE)

       ! passively advected quantities -- only the contact matters
       ! note: this also includes the y- and z-velocity flux
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nqp = qpass_map(ipassive)

          if (ustar > ZERO) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nqp)
          else if (ustar < ZERO) then
             uflx(k,n) = uflx(k,URHO)*qr(k,nqp)
          else
             qavg = HALF * (ql(k,nqp) + qr(k,nqp))
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
                       qint,qg_l1,qg_h1, &
#ifdef RADIATION
                       lam, gamcgl, gamcgr, &
                       rflx,rflx_l1,rflx_h1, &
#endif
                       ilo,ihi,domlo,domhi)

    use prob_params_module, only : physbc_lo, physbc_hi, Outflow, Symmetry
#ifdef RADIATION
    use meth_params_module, only : qrad, qradhi, qptot, qreitot, fspace_type
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif

    use amrex_fort_module, only : rt => c_real
    real(rt)        , parameter:: small = 1.e-8_rt

    integer ilo,ihi
    integer domlo(1),domhi(1)
    integer  qpd_l1,  qpd_h1
    integer   qg_l1,   qg_h1
    integer uflx_l1, uflx_h1

#ifdef RADIATION
    integer rflx_l1, rflx_h1
#endif

    real(rt)         ql(qpd_l1:qpd_h1, NQ)
    real(rt)         qr(qpd_l1:qpd_h1, NQ)

    real(rt)           cav(ilo:ihi+1), smallc(ilo:ihi+1)
    real(rt)         gamcl(ilo:ihi+1), gamcr(ilo:ihi+1)
    real(rt)          uflx(uflx_l1:uflx_h1, NVAR)
    real(rt)          qint( qg_l1: qg_h1, NGDNV)

#ifdef RADIATION
    real(rt)         lam(ilo-1:ihi+2, 0:ngroups-1)
    real(rt)         gamcgl(ilo:ihi+1),gamcgr(ilo:ihi+1)
    real(rt)          rflx(rflx_l1:rflx_h1, 0:ngroups-1)

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
