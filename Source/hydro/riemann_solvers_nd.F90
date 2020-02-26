module riemann_solvers_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use meth_params_module, only : NQ, NVAR, NQAUX, &
                                 URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UFX, UTEMP, &
                                 QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QFS, QFX, &
                                 QC, QGAMC, QGC, &
                                 NGDNV, GDRHO, GDPRES, &

#ifdef RADIATION
                                 qrad, qptot, qreitot, &
                                 GDERADS, QGAMCG, QLAMS, QREITOT, &
#endif
                                 npassive, upass_map, qpass_map, &
                                 small_dens, small_pres, small_temp, &
                                 use_eos_in_riemann, use_reconstructed_gamma1, &
                                 hybrid_riemann, riemann_solver

  use riemann_util_module

  implicit none

  real(rt), parameter :: smallu = 1.e-12_rt
  real(rt), parameter :: small = 1.e-8_rt

  private :: smallu, small

contains

  subroutine riemanncg(ql, ql_lo, ql_hi, &
                       qr, qr_lo, qr_hi, &
                       qaux, qa_lo, qa_hi, &
                       qint, q_lo, q_hi, &
                       idir, lo, hi, &
                       domlo, domhi)
    ! this implements the approximate Riemann solver of Colella & Glaz
    ! (1985)
    !
    ! this version is dimension agnostic -- for 1- and 2-d, set kc,
    ! kflux, and k3d to 0

    use castro_error_module
#ifndef AMREX_USE_CUDA
    use amrex_mempool_module, only : bl_allocate, bl_deallocate
#endif
    use prob_params_module, only : physbc_lo, physbc_hi, &
         Symmetry, SlipWall, NoSlipWall
    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use meth_params_module, only : cg_maxiter, cg_tol, cg_blend
#ifndef AMREX_USE_CUDA
    use riemann_util_module, only : pstar_bisection
#endif
    use riemann_util_module, only : wsqge

    implicit none

    integer, intent(in) :: ql_lo(3), ql_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(ql_lo(1):ql_hi(1),ql_lo(2):ql_hi(2),ql_lo(3):ql_hi(3),NQ)
    real(rt), intent(in) :: qr(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

    integer :: i, j, k
    integer :: n, nqp, ipassive

    real(rt) :: ustar
    real(rt) :: rl, ul, v1l, v2l, pl, rel
    real(rt) :: rr, ur, v1r, v2r, pr, rer
    real(rt) :: wl, wr

    real(rt) :: rstar, cstar, pstar
    real(rt) :: ro, uo, po, co, gamco
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall, qavg
    real(rt) :: cavg

    real(rt) :: gcl, gcr, game_int
    real(rt) :: clsq, clsql, clsqr, wlsq, wosq, wrsq, wo
    real(rt) :: zm, zp
    real(rt) :: denom, dpditer, dpjmp
    real(rt) :: gamc_bar, game_bar
    real(rt) :: gamel, gamer, gameo, gamstar, gmin, gmax, gdot

    integer :: iter, iter_max
    real(rt) :: tol
    real(rt) :: err

    logical :: converged

    real(rt) :: pstar_old
    real(rt) :: taul, taur, tauo
    real(rt) :: ustar_r, ustar_l, ustar_r_old, ustar_l_old
    real(rt) :: pstarl, pstaru

    real(rt), parameter :: weakwv = 1.e-3_rt

#ifndef AMREX_USE_CUDA
    real(rt), pointer :: pstar_hist(:), pstar_hist_extra(:)
#endif

    type (eos_t) :: eos_state

    real(rt) :: u_adv

    integer :: iu, iv1, iv2, im1, im2, im3, sx, sy, sz

    !$gpu

#ifndef AMREX_USE_CUDA
    if (cg_blend == 2 .and. cg_maxiter < 5) then

       call castro_error("Error: need cg_maxiter >= 5 to do a bisection search on secant iteration failure.")

    endif
#endif

    if (idir == 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       im1 = UMX
       im2 = UMY
       im3 = UMZ
       sx = 1
       sy = 0
       sz = 0
    else if (idir == 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       im1 = UMY
       im2 = UMX
       im3 = UMZ
       sx = 0
       sy = 1
       sz = 0
    else
       iu = QW
       iv1 = QU
       iv2 = QV
       im1 = UMZ
       im2 = UMX
       im3 = UMY
       sx = 0
       sy = 0
       sz = 1
    end if

    tol = cg_tol
    iter_max = cg_maxiter

#ifndef AMREX_USE_CUDA
    call bl_allocate(pstar_hist, 1,iter_max)
    call bl_allocate(pstar_hist_extra, 1,2*iter_max)
#endif

    ! left state
    rl = max(ql(QRHO), small_dens)
    pl  = ql(QPRES)
    rel = ql(QREINT)

    ! pick left velocities based on direction
    ul  = ql(iu)
    v1l = ql(iv1)
    v2l = ql(iv2)

    ! sometime we come in here with negative energy or pressure
    ! note: reset both in either case, to remain thermo
    ! consistent
    if (rel <= ZERO .or. pl < small_pres) then
#ifndef AMREX_USE_CUDA
       print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres
#endif

       eos_state % T   = small_temp
       eos_state % rho = rl
       eos_state % xn  = ql(QFS:QFS-1+nspec)
       eos_state % aux = ql(QFX:QFX-1+naux)

       call eos(eos_input_rt, eos_state)

       rel = rl*eos_state % e
       pl  = eos_state % p
       gcl = eos_state % gam1
    end if

    ! right state
    rr = max(qr(QRHO), small_dens)
    pr  = qr(QPRES)
    rer = qr(QREINT)

    ! pick right velocities based on direction
    ur  = qr(iu)
    v1r = qr(iv1)
    v2r = qr(iv2)

    if (rer <= ZERO .or. pr < small_pres) then
#ifndef AMREX_USE_CUDA
       print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres
#endif

       eos_state % T   = small_temp
       eos_state % rho = rr
       eos_state % xn  = qr(QFS:QFS-1+nspec)
       eos_state % aux = qr(QFX:QFX-1+naux)

       call eos(eos_input_rt, eos_state)

       rer = rr*eos_state % e
       pr  = eos_state % p
       gcr = eos_state % gam1
    end if


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
    gmin = min(gamel, gamer, ONE, FOUR3RD)
    gmax = max(gamel, gamer, TWO, FIVE3RD)

    game_bar = HALF*(gamel + gamer)
    gamc_bar = HALF*(gcl + gcr)

    gdot = TWO*(ONE - game_bar/gamc_bar)*(game_bar - ONE)

    wsmall = small_dens*csmall
    wl = max(wsmall,sqrt(abs(clsql)))
    wr = max(wsmall,sqrt(abs(clsqr)))

    ! make an initial guess for pstar -- this is a two-shock
    ! approximation
    !pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
    pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)
    pstar = max(pstar, small_pres)

    ! get the shock speeds -- this computes W_s from CG Eq. 34
    call wsqge(pl, taul, gamel, gdot,  &
               gamstar, pstar, wlsq, clsql, gmin, gmax)

    call wsqge(pr, taur, gamer, gdot,  &
               gamstar, pstar, wrsq, clsqr, gmin, gmax)

    pstar_old = pstar

    wl = sqrt(wlsq)
    wr = sqrt(wrsq)

    ! R-H jump conditions give ustar across each wave -- these
    ! should be equal when we are done iterating.  Our notation
    ! here is a little funny, comparing to CG, ustar_l = u*_L and
    ! ustar_r = u*_R.
    ustar_l = ul - (pstar-pl)/wl
    ustar_r = ur + (pstar-pr)/wr

    ! revise our pstar guess
    !pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
    pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)
    pstar = max(pstar, small_pres)

    ! secant iteration
    converged = .false.
    iter = 1
    do while ((iter <= iter_max .and. .not. converged) .or. iter <= 2)

       call wsqge(pl, taul, gamel, gdot,  &
                  gamstar, pstar, wlsq, clsql, gmin, gmax)

       call wsqge(pr, taur, gamer, gdot,  &
                  gamstar, pstar, wrsq, clsqr, gmin, gmax)

       ! NOTE: these are really the inverses of the wave speeds!
       wl = ONE / sqrt(wlsq)
       wr = ONE / sqrt(wrsq)

       ustar_r_old = ustar_r
       ustar_l_old = ustar_l

       ustar_r = ur - (pr-pstar)*wr
       ustar_l = ul + (pl-pstar)*wl

       dpditer = abs(pstar_old-pstar)

       ! Here we are going to do the Secant iteration version in
       ! CG.  Note that what we call zp and zm here are not
       ! actually the Z_p = |dp*/du*_p| defined in CG, by rather
       ! simply |du*_p| (or something that looks like dp/Z!).
       zp = abs(ustar_l - ustar_l_old)
       if (zp - weakwv*cavg <= ZERO) then
          zp = dpditer*wl
       endif

       zm = abs(ustar_r - ustar_r_old)
       if (zm - weakwv*cavg <= ZERO) then
          zm = dpditer*wr
       endif

       ! the new pstar is found via CG Eq. 18
       denom = dpditer/max(zp+zm, small*cavg)
       pstar_old = pstar
       pstar = pstar - denom*(ustar_r - ustar_l)
       pstar = max(pstar, small_pres)

       err = abs(pstar - pstar_old)
       if (err < tol*pstar) converged = .true.

#ifndef AMREX_USE_CUDA
       pstar_hist(iter) = pstar
#endif

       iter = iter + 1

    end do

    ! If we failed to converge using the secant iteration, we
    ! can either stop here; or, revert to the original
    ! two-shock estimate for pstar; or do a bisection root
    ! find using the bounds established by the most recent
    ! iterations.

    if (.not. converged) then

       if (cg_blend == 0) then

#ifndef AMREX_USE_CUDA
          print *, 'pstar history: '
          do iter = 1, iter_max
             print *, iter, pstar_hist(iter)
          end do

          print *, ' '
          print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
          print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr
          print *, 'cavg, smallc:', cavg, csmall
          call castro_error("ERROR: non-convergence in the Riemann solver")
#endif
       else if (cg_blend == 1) then

          pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)

       else if (cg_blend == 2) then

          ! we don't store the history if we are in CUDA, so
          ! we can't do this
#ifndef AMREX_USE_CUDA
          ! first try to find a reasonable bounds
          pstarl = minval(pstar_hist(iter_max-5:iter_max))
          pstaru = maxval(pstar_hist(iter_max-5:iter_max))

          call pstar_bisection(pstarl, pstaru, &
                               ul, pl, taul, gamel, clsql, &
                               ur, pr, taur, gamer, clsqr, &
                               gdot, gmin, gmax, &
                               pstar, gamstar, converged, pstar_hist_extra)

          if (.not. converged) then

             print *, 'pstar history: '
             do iter = 1, iter_max
                print *, iter, pstar_hist(iter)
             end do
             do iter = 1, 2 * iter_max
                print *, iter + iter_max, pstar_hist_extra(iter)
             end do

             print *, ' '
             print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
             print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr
             print *, 'cavg, smallc:', cavg, csmall
             call castro_error("ERROR: non-convergence in the Riemann solver")

          end if
#endif
       else

#ifndef AMREX_USE_CUDA
          call castro_error("ERROR: unrecognized cg_blend option.")
#endif
       end if

    end if

    ! we converged!  construct the single ustar for the region
    ! between the left and right waves, using the updated wave speeds
    ustar_r = ur - (pr-pstar)*wr  ! careful -- here wl, wr are 1/W
    ustar_l = ul + (pl-pstar)*wl

    ustar = HALF* (ustar_l + ustar_r)

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
       gamco = gcl
       gameo = gamel

    else if (ustar < ZERO) then
       ro = rr
       uo = ur
       po = pr
       tauo = taur
       gamco = gcr
       gameo = gamer

    else
       ro = HALF*(rl+rr)
       uo = HALF*(ul+ur)
       po = HALF*(pl+pr)
       tauo = HALF*(taul+taur)
       gamco = HALF*(gcl+gcr)
       gameo = HALF*(gamel + gamer)
    end if

    ! use tau = 1/rho as the independent variable here
    ro = max(small_dens, ONE/tauo)
    tauo = ONE/ro

    co = sqrt(abs(gamco*po*tauo))
    co = max(csmall, co)
    clsq = (co*ro)**2

    ! now that we know which state (left or right) we need to worry
    ! about, get the value of gamstar and wosq across the wave we
    ! are dealing with.
    call wsqge(po, tauo, gameo, gdot,   &
               gamstar, pstar, wosq, clsq, gmin, gmax)

    sgnm = sign(ONE, ustar)

    wo = sqrt(wosq)
    dpjmp = pstar - po

    ! is this max really necessary?
    !rstar=max(ONE-ro*dpjmp/wosq, (gameo-ONE)/(gameo+ONE))
    rstar = ONE - ro*dpjmp/wosq
    rstar = ro/rstar
    rstar = max(small_dens, rstar)

    cstar = sqrt(abs(gamco*pstar/rstar))
    cstar = max(cstar, csmall)

    spout = co - sgnm*uo
    spin = cstar - sgnm*ustar

    !ushock = HALF*(spin + spout)
    ushock = wo*tauo - sgnm*uo

    if (pstar-po >= ZERO) then
       spin = ushock
       spout = ushock
    endif

    frac = HALF*(ONE + (spin + spout)/max(spout-spin, spin+spout, small*cavg))

    ! the transverse velocity states only depend on the
    ! direction that the contact moves
    if (ustar > ZERO) then
       qint(iv1) = v1l
       qint(iv2) = v2l
    else if (ustar < ZERO) then
       qint(iv1) = v1r
       qint(iv2) = v2r
    else
       qint(iv1) = HALF*(v1l+v1r)
       qint(iv2) = HALF*(v2l+v2r)
    endif

    ! linearly interpolate between the star and normal state -- this covers the
    ! case where we are inside the rarefaction fan.
    qint(QRHO) = frac*rstar + (ONE - frac)*ro
    qint(iu) = frac*ustar + (ONE - frac)*uo
    qint(QPRES) = frac*pstar + (ONE - frac)*po
    game_int = frac*gamstar + (ONE-frac)*gameo

    ! now handle the cases where instead we are fully in the
    ! star or fully in the original (l/r) state
    if (spout < ZERO) then
       qint(QRHO) = ro
       qint(iu) = uo
       qint(QPRES) = po
       game_int = gameo
    end if

    if (spin >= ZERO) then
       qint(QRHO) = rstar
       qint(iu) = ustar
       qint(QPRES) = pstar
       game_int = gamstar
    end if

    qint(QPRES) = max(qint(QPRES), small_pres)

    u_adv = qint(iu)

    ! Enforce that fluxes through a symmetry plane or wall are hard zero.
    u_adv = u_adv * bnd_fac

    ! Compute fluxes, order as conserved state (not q)
    qint(iu) = u_adv

    ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
    qint(QREINT) = qint(QPRES)/(game_int - ONE)

    ! advected quantities -- only the contact matters
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       if (ustar > ZERO) then
          qint(nqp) = ql(nqp)
       else if (ustar < ZERO) then
          qint(nqp) = qr(nqp)
       else
          qavg = HALF * (ql(nqp) + qr(nqp))
          qint(nqp) = qavg
       end if
    end do

#ifndef AMREX_USE_CUDA
    call bl_deallocate(pstar_hist)
    call bl_deallocate(pstar_hist_extra)
#endif

  end subroutine riemanncg


  subroutine riemannus(ql, ql_lo, ql_hi, &
                       qr, qr_lo, qr_hi, &
                       qaux, qa_lo, qa_hi, &
                       qint, q_lo, q_hi, &
#ifdef RADIATION
                       lambda_int, l_lo, l_hi, &
#endif
                       idir, compute_gammas, lo, hi, &
                       domlo, domhi)
    ! Colella, Glaz, and Ferguson solver
    !
    ! this is a 2-shock solver that uses a very simple approximation for the
    ! star state, and carries an auxiliary jump condition for (rho e) to
    ! deal with a real gas

    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall
    use eos_type_module, only : eos_t, eos_input_rp
    use eos_module, only : eos
    use network, only : nspec, naux
    use meth_params_module, only: T_guess

    implicit none

    integer, intent(in) :: ql_lo(3), ql_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, lo(3), hi(3)
    integer, intent(in) :: domlo(3),domhi(3)

#ifdef RADIATION
    integer, intent(in) :: l_lo(3), l_hi(3)
#endif

    integer, intent(in) :: compute_gammas
    real(rt), intent(in) :: ql(ql_lo(1):ql_hi(1),ql_lo(2):ql_hi(2),ql_lo(3):ql_hi(3),NQ)
    real(rt), intent(in) :: qr(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
#ifdef RADIATION
    real(rt), intent(inout) :: lambda_int(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3),0:ngroups-1)
#endif

    integer :: i, j, k
    integer :: nqp, ipassive

    real(rt) :: regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel
    real(rt) :: rr, ur, v1r, v2r, pr, rer
    real(rt) :: wl, wr, scr
    real(rt) :: rstar, cstar, estar, pstar, ustar
    real(rt) :: ro, uo, po, reo, co, gamco, entho, drho
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall
    real(rt) :: cavg, gamcl, gamcr

#ifdef RADIATION
    real(rt), dimension(0:ngroups-1) :: erl, err
    real(rt) :: reo_g, po_g, co_g, gamco_g
    real(rt) :: pl_g, rel_g, pr_g, rer_g
    real(rt) :: regdnv_g, pgdnv_g, pgdnv_t
    real(rt) :: estar_g, pstar_g
    real(rt), dimension(0:ngroups-1) :: lambda, laml, lamr, reo_r, po_r, estar_r, regdnv_r
    integer :: g
    real(rt) :: gamcgl, gamcgr
#endif

    real(rt) :: u_adv

    integer :: iu, iv1, iv2, im1, im2, im3
    real(rt) :: wwinv, roinv, co2inv
    real(rt) :: fp, fm

    type(eos_t) :: eos_state
    real(rt), dimension(nspec) :: xn

    !$gpu

    ! set integer pointers for the normal and transverse velocity and
    ! momentum

    if (idir == 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       im1 = UMX
       im2 = UMY
       im3 = UMZ
    else if (idir == 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       im1 = UMY
       im2 = UMX
       im3 = UMZ
    else
       iu = QW
       iv1 = QU
       iv2 = QV
       im1 = UMZ
       im2 = UMX
       im3 = UMY
    end if


    rl = max(ql(QRHO), small_dens)

    ! pick left velocities based on direction
    ul  = ql(iu)
    v1l = ql(iv1)
    v2l = ql(iv2)

#ifdef RADIATION
    pl = ql(qptot)
    rel = ql(qreitot)
    erl(:) = ql(qrad:qrad-1+ngroups)
    pl_g = ql(QPRES)
    rel_g = ql(QREINT)
#else
    pl  = max(ql(QPRES), small_pres)
    rel = ql(QREINT)
#endif

    rr = max(qr(QRHO), small_dens)

    ! pick right velocities based on direction
    ur  = qr(iu)
    v1r = qr(iv1)
    v2r = qr(iv2)

#ifdef RADIATION
    pr = qr(qptot)
    rer = qr(qreitot)
    err(:) = qr(qrad:qrad-1+ngroups)
    pr_g = qr(QPRES)
    rer_g = qr(QREINT)
#else
    pr  = max(qr(QPRES), small_pres)
    rer = qr(QREINT)
#endif

    ! ------------------------------------------------------------------
    ! estimate the star state: pstar, ustar
    ! ------------------------------------------------------------------

    wsmall = small_dens*csmall

    ! this is Castro I: Eq. 33
    wl = max(wsmall, sqrt(abs(gamcl*pl*rl)))
    wr = max(wsmall, sqrt(abs(gamcr*pr*rr)))

    wwinv = ONE/(wl + wr)
    pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
    ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

    pstar = max(pstar, small_pres)

    ! for symmetry preservation, if ustar is really small, then we
    ! set it to zero
    if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
       ustar = ZERO
    end if

    ! ------------------------------------------------------------------
    ! look at the contact to determine which region we are in
    ! ------------------------------------------------------------------

    ! this just determines which of the left or right states is still
    ! in play.  We still need to look at the other wave to determine
    ! if the star state or this state is on the interface.
    sgnm = sign(ONE, ustar)
    if (ustar == ZERO) sgnm = ZERO

    fp = HALF*(ONE + sgnm)
    fm = HALF*(ONE - sgnm)

    ro = fp*rl + fm*rr
    uo = fp*ul + fm*ur
    po = fp*pl + fm*pr
    reo = fp*rel + fm*rer
    gamco = fp*gamcl + fm*gamcr
#ifdef RADIATION
    lambda = fp*laml + fm*lamr

    if (ustar == 0) then
       ! harmonic average
       do g=0, ngroups-1
          lambda(g) = 2.0e0_rt*(laml(g)*lamr(g))/(laml(g)+lamr(g)+1.e-50_rt)
       end do
    end if

    po_g = fp*pl_g + fm*pr_g
    reo_r(:) = fp*erl(:) + fm*err(:)
    po_r(:) = lambda(:)*reo_r(:)
    reo_g = fp*rel_g + fm*rer_g
    gamco_g = fp*gamcgl + fm*gamcgr
#endif

    ro = max(small_dens, ro)

    roinv = ONE/ro

    co = sqrt(abs(gamco*po*roinv))
    co = max(csmall,co)
    co2inv = ONE/(co*co)

    ! we can already deal with the transverse velocities -- they
    ! only jump across the contact
    qint(iv1) = fp*v1l + fm*v1r
    qint(iv2) = fp*v2l + fm*v2r

    ! ------------------------------------------------------------------
    ! compute the rest of the star state
    ! ------------------------------------------------------------------

    drho = (pstar - po)*co2inv
    rstar = ro + drho
    rstar = max(small_dens, rstar)

#ifdef RADIATION
    estar_g = reo_g + drho*(reo_g + po_g)*roinv
    co_g = sqrt(abs(gamco_g*po_g*roinv))
    co_g = max(csmall, co_g)
    pstar_g = po_g + drho*co_g**2
    pstar_g = max(pstar_g, small_pres)
    estar_r = reo_r(:) + drho*(reo_r(:) + po_r(:))*roinv
#else
    entho = (reo + po)*roinv*co2inv
    estar = reo + (pstar - po)*entho
#endif
    cstar = sqrt(abs(gamco*pstar/rstar))
    cstar = max(cstar, csmall)

    ! ------------------------------------------------------------------
    ! finish sampling the solution
    ! ------------------------------------------------------------------

    ! look at the remaining wave to determine if the star state or the
    ! 'o' state above is on the interface


    ! the values of u +/- c on either side of the non-contact
    ! wave
    spout = co - sgnm*uo
    spin = cstar - sgnm*ustar

    ! a simple estimate of the shock speed
    ushock = HALF*(spin + spout)

    if (pstar-po > ZERO) then
       spin = ushock
       spout = ushock
    endif

    if (spout-spin == ZERO) then
       scr = small*cavg
    else
       scr = spout-spin
    endif

    ! interpolate for the case that we are in a rarefaction
    frac = (ONE + (spout + spin)/scr)*HALF
    frac = max(ZERO, min(ONE, frac))

    qint(QRHO) = frac*rstar + (ONE - frac)*ro
    qint(iu) = frac*ustar + (ONE - frac)*uo

#ifdef RADIATION
    pgdnv_t = frac*pstar + (ONE - frac)*po
    pgdnv_g = frac*pstar_g + (ONE - frac)*po_g
    regdnv_g = frac*estar_g + (ONE - frac)*reo_g
    regdnv_r(:) = frac*estar_r(:) + (ONE - frac)*reo_r(:)
#else
    qint(QPRES) = frac*pstar + (ONE - frac)*po
    regdnv = frac*estar + (ONE - frac)*reo
#endif

    ! as it stands now, we set things assuming that the rarefaction
    ! spans the interface.  We overwrite that here depending on the
    ! wave speeds

    ! look at the speeds on either side of the remaining wave
    ! to determine which region we are in
    if (spout < ZERO) then
       ! the l or r state is on the interface
       qint(QRHO) = ro
       qint(iu) = uo
#ifdef RADIATION
       pgdnv_t = po
       pgdnv_g = po_g
       regdnv_g = reo_g
       regdnv_r = reo_r(:)
#else
       qint(PRES) = po
       regdnv = reo
#endif
    end if

    if (spin >= ZERO) then
       ! the star state is on the interface
       qint(QRHO) = rstar
       qint(iu  ) = ustar
#ifdef RADIATION
       pgdnv_t = pstar
       pgdnv_g = pstar_g
       regdnv_g = estar_g
       regdnv_r = estar_r(:)
#else
       qint(QPRES) = pstar
       regdnv = estar
#endif
    end if


#ifdef RADIATION
    do g=0, ngroups-1
       qint(QRAD+g) = max(regdnv_r(g), 0.e0_rt)
    end do

    qint(QPRES) = pgdnv_g
    qint(QPTOT) = pgdnv_t
    qint(QREINT) = regdnv_g
    qint(QREITOT) = sum(regdnv_r(:)) + regdnv_g

    lambda_int(:) = lambda(:)

#else
    qint(QPRES) = max(qint(QPRES),small_pres)
    qint(QREINT) = regdnv
#endif


    ! we are potentially thermodynamically inconsistent, fix that
    ! here
    if (use_eos_in_riemann == 1) then
       ! we need to know the species -- they only jump across
       ! the contact
       xn(:) = fp*ql(QFS:QFS-1+nspec) + fm*qr(QFS:QFS-1+nspec)

       eos_state % rho = qint(QRHO)
       eos_state % p = qint(QPRES)
       eos_state % xn(:) = xn(:)
       eos_state % T = T_guess

       call eos(eos_input_rp, eos_state)

       qint(QREINT) = eos_state % rho * eos_state % e
    end if

    ! Enforce that fluxes through a symmetry plane or wall are hard zero.
    u_adv = qint(i,j,k,iu)

    u_adv = u_adv * bnd_fac

    qint(iu) = u_adv

    ! passively advected quantities
    do ipassive = 1, npassive
       nqp = qpass_map(ipassive)
       qint(nqp) = fp*ql(nqp) + fm*qr(nqp)
    end do

  end subroutine riemannus


  subroutine HLLC(ql, ql_lo, ql_hi, &
                  qr, qr_lo, qr_hi, &
                  qaux, qa_lo, qa_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  qint, q_lo, q_hi, &
                  idir, lo, hi, &
                  domlo, domhi)
    ! this is an implementation of the HLLC solver described in Toro's
    ! book.  it uses the simplest estimate of the wave speeds, since
    ! those should work for a general EOS.  We also initially do the
    ! CGF Riemann construction to get pstar and ustar, since we'll
    ! need to know the pressure and velocity on the interface for the
    ! pdV term in the internal energy update.

    use prob_params_module, only : physbc_lo, physbc_hi, &
         Symmetry, SlipWall, NoSlipWall

    implicit none

    integer, intent(in) :: ql_lo(3), ql_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(ql_lo(1):ql_hi(1),ql_lo(2):ql_hi(2),ql_lo(3):ql_hi(3),NQ)
    real(rt), intent(in) :: qr(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

    integer :: i, j, k

    real(rt) :: rgdnv, regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel
    real(rt) :: rr, ur, v1r, v2r, pr, rer
    real(rt) :: wl, wr, scr
    real(rt) :: rstar, cstar, estar, pstar, ustar
    real(rt) :: ro, uo, po, reo, co, gamco, entho
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall
    real(rt) :: cavg, gamcl, gamcr

    integer :: iu, iv1, iv2, sx, sy, sz
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac
    real(rt) :: wwinv, roinv, co2inv

    real(rt) :: U_hllc_state(nvar), U_state(nvar), F_state(nvar)
    real(rt) :: S_l, S_r, S_c

    real(rt) :: q_zone(NQ)

    !$gpu

    if (idir == 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       sx = 1
       sy = 0
       sz = 0
    else if (idir == 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       sx = 0
       sy = 1
       sz = 0
    else
       iu = QW
       iv1 = QU
       iv2 = QV
       sx = 0
       sy = 0
       sz = 1
    end if

    special_bnd_lo = (physbc_lo(idir) == Symmetry &
         .or.         physbc_lo(idir) == SlipWall &
         .or.         physbc_lo(idir) == NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) == Symmetry &
         .or.         physbc_hi(idir) == SlipWall &
         .or.         physbc_hi(idir) == NoSlipWall)

    if (idir == 1) then
       special_bnd_lo_x = special_bnd_lo
       special_bnd_hi_x = special_bnd_hi
    else
       special_bnd_lo_x = .false.
       special_bnd_hi_x = .false.
    end if

    do k = lo(3), hi(3)

       bnd_fac_z = 1
       if (idir == 3) then
          if ( k == domlo(3)   .and. special_bnd_lo .or. &
               k == domhi(3)+1 .and. special_bnd_hi ) then
             bnd_fac_z = 0
          end if
       end if

       do j = lo(2), hi(2)

          bnd_fac_y = 1
          if (idir == 2) then
             if ( j == domlo(2)   .and. special_bnd_lo .or. &
                  j == domhi(2)+1 .and. special_bnd_hi ) then
                bnd_fac_y = 0
             end if
          end if

          !dir$ ivdep
          do i = lo(1), hi(1)

             rl = max(ql(i,j,k,QRHO), small_dens)

             ! pick left velocities based on direction
             ul  = ql(i,j,k,iu)
             v1l = ql(i,j,k,iv1)
             v2l = ql(i,j,k,iv2)

             pl  = max(ql(i,j,k,QPRES), small_pres)
             rel = ql(i,j,k,QREINT)

             rr = max(qr(i,j,k,QRHO), small_dens)

             ! pick right velocities based on direction
             ur  = qr(i,j,k,iu)
             v1r = qr(i,j,k,iv1)
             v2r = qr(i,j,k,iv2)

             pr  = max(qr(i,j,k,QPRES), small_pres)
             rer = qr(i,j,k,QREINT)

             ! now we essentially do the CGF solver to get p and u on the
             ! interface, but we won't use these in any flux construction.
             csmall = max( small, max(small * qaux(i,j,k,QC) , small * qaux(i-sx,j-sy,k-sz,QC)) )
             cavg = HALF*(qaux(i,j,k,QC) + qaux(i-sx,j-sy,k-sz,QC))

             if (use_reconstructed_gamma1 == 1) then
                gamcl = ql(i,j,k,QGC)
                gamcr = qr(i,j,k,QGC)
             else
                gamcl = qaux(i-sx,j-sy,k-sz,QGAMC)
                gamcr = qaux(i,j,k,QGAMC)
             endif

             wsmall = small_dens*csmall
             wl = max(wsmall, sqrt(abs(gamcl*pl*rl)))
             wr = max(wsmall, sqrt(abs(gamcr*pr*rr)))

             wwinv = ONE/(wl + wr)
             pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
             ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

             pstar = max(pstar, small_pres)
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
                gamco = gamcl

             else if (ustar < ZERO) then
                ro = rr
                uo = ur
                po = pr
                reo = rer
                gamco = gamcr
             else
                ro = HALF*(rl + rr)
                uo = HALF*(ul + ur)
                po = HALF*(pl + pr)
                reo = HALF*(rel + rer)
                gamco = HALF*(gamcl + gamcr)
             endif
             ro = max(small_dens, ro)

             roinv = ONE/ro
             co = sqrt(abs(gamco*po*roinv))
             co = max(csmall, co)
             co2inv = ONE/(co*co)

             rstar = ro + (pstar - po)*co2inv
             rstar = max(small_dens, rstar)

             entho = (reo + po)*co2inv * roinv
             estar = reo + (pstar - po)*entho

             cstar = sqrt(abs(gamco*pstar/rstar))
             cstar = max(cstar, csmall)

             sgnm = sign(ONE, ustar)
             spout = co - sgnm*uo
             spin = cstar - sgnm*ustar
             ushock = HALF*(spin + spout)

             if (pstar-po > ZERO) then
                spin = ushock
                spout = ushock
             endif
             if (spout-spin == ZERO) then
                scr = small*cavg
             else
                scr = spout-spin
             endif
             frac = (ONE + (spout + spin)/scr)*HALF
             frac = max(ZERO, min(ONE, frac))

             rgdnv = frac*rstar + (ONE - frac)*ro
             regdnv = frac*estar + (ONE - frac)*reo

             qint(i,j,k,iu) = frac*ustar + (ONE - frac)*uo
             qint(i,j,k,QPRES) = frac*pstar + (ONE - frac)*po


             ! now we do the HLLC construction


             ! Enforce that the fluxes through a symmetry plane or wall are hard zero.
             if ( special_bnd_lo_x .and. i== domlo(1) .or. &
                  special_bnd_hi_x .and. i== domhi(1)+1 ) then
                bnd_fac_x = 0
             else
                bnd_fac_x = 1
             end if

             bnd_fac = bnd_fac_x*bnd_fac_y*bnd_fac_z

             ! use the simplest estimates of the wave speeds
             S_l = min(ul - sqrt(gamcl*pl/rl), ur - sqrt(gamcr*pr/rr))
             S_r = max(ul + sqrt(gamcl*pl/rl), ur + sqrt(gamcr*pr/rr))

             ! estimate of the contact speed -- this is Toro Eq. 10.8
             S_c = (pr - pl + rl*ul*(S_l - ul) - rr*ur*(S_r - ur))/ &
                  (rl*(S_l - ul) - rr*(S_r - ur))

             if (S_r <= ZERO) then
                ! R region
                q_zone(:) = qr(i,j,k,:)
                call cons_state(q_zone, U_state)
                call compute_flux(idir, bnd_fac, U_state, pr, F_state)

             else if (S_r > ZERO .and. S_c <= ZERO) then
                ! R* region
                q_zone(:) = qr(i,j,k,:)
                call cons_state(q_zone, U_state)
                call compute_flux(idir, bnd_fac, U_state, pr, F_state)

                call HLLC_state(idir, S_r, S_c, q_zone, U_hllc_state)

                ! correct the flux
                F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

             else if (S_c > ZERO .and. S_l < ZERO) then
                ! L* region
                q_zone(:) = ql(i,j,k,:)
                call cons_state(q_zone, U_state)
                call compute_flux(idir, bnd_fac, U_state, pl, F_state)

                call HLLC_state(idir, S_l, S_c, q_zone, U_hllc_state)

                ! correct the flux
                F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

             else
                ! L region
                q_zone(:) = ql(i,j,k,:)
                call cons_state(q_zone, U_state)
                call compute_flux(idir, bnd_fac, U_state, pl, F_state)

             endif

             uflx(i,j,k,:) = F_state(:)
          end do
       end do
    end do

  end subroutine HLLC

end module riemann_solvers_module
