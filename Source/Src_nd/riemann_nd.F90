module actual_riemann_module

  use amrex_fort_module, only : rt => amrex_real
  use bl_constants_module
  use meth_params_module, only : NQ, NVAR, &
                                 URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UFX, &
                                 QRHO, QU, QV, QW, &
                                 QPRES, QGAME, QREINT, QFS, QFX, &
                                 NGDNV, GDRHO, GDPRES, GDGAME, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, fspace_type, &
                                 GDERADS, GDLAMS, QGAMCG, QLAMS, &
#endif
                                 npassive, upass_map, qpass_map, &
                                 small_dens, small_pres, small_temp
  use riemann_util_module

#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
  
  implicit none

  private

  public :: riemanncg, riemannus, hllc

  real(rt), parameter :: smallu = 1.e-12_rt
  real(rt), parameter :: small = 1.e-8_rt

contains

  subroutine riemanncg(ql, qr, qpd_lo, qpd_hi, &
                       gamcl, gamcr, cav, smallc, gd_lo, gd_hi, &
                       uflx, uflx_lo, uflx_hi, &
                       qint, q_lo, q_hi, &
                       idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, &
                       domlo, domhi)

    ! this implements the approximate Riemann solver of Colella & Glaz
    ! (1985)

    ! this version is dimension agnostic -- for 1- and 2-d, set kc,
    ! kflux, and k3d to 0

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall, &
                                   mom_flux_has_p
    use network, only : nspec, naux
    use eos_type_module
    use eos_module
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : compute_hybrid_flux
#endif
    use meth_params_module, only : cg_maxiter, cg_tol, cg_blend
    use riemann_util_module, only : wsqge, pstar_bisection

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: gd_lo(2), gd_hi(2)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, ilo, ihi, jlo, jhi
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) ::  gamcl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) ::  gamcr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) ::    cav(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) :: smallc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

    ! Note:
    !
    !  k3d: the k corresponding to the full 3d array -- it should be
    !       used for print statements or tests against domlo, domhi,
    !       etc
    !
    !  kc: the k corresponding to the 2-wide slab of k-planes, so in
    !      this routine it takes values only of 1 or 2
    !
    !  kflux: used for indexing the uflx array -- in the initial calls
    !         to cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, 
    !         kflux = kc, but in later calls, when uflx = {flux1,flux2,flux3},
    !         kflux = k3d

    integer :: i,j,kc,kflux,k3d
    integer :: n, nqp, ipassive

    real(rt)         :: ustar,gamgdnv
    real(rt)         :: rl, ul, v1l, v2l, pl, rel
    real(rt)         :: rr, ur, v1r, v2r, pr, rer
    real(rt)         :: wl, wr, rhoetot
   !real(rt)         :: scr
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

    real(rt)         :: pstar_old
    real(rt)         :: taul, taur, tauo
    real(rt)         :: ustar_r, ustar_l, ustar_r_old, ustar_l_old
    real(rt)         :: pstarl, pstarc, pstaru, pfuncc, pfuncu

    real(rt)        , parameter :: weakwv = 1.e-3_rt

    real(rt)        , pointer :: pstar_hist(:), pstar_hist_extra(:)

    type (eos_t) :: eos_state

    real(rt)        , pointer :: us1d(:)

#ifdef ROTATION
    real(rt)         :: vel(3)
#endif

    real(rt)         :: u_adv

    integer :: iu, iv1, iv2, im1, im2, im3
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    real(rt)         :: bnd_fac_x, bnd_fac_y, bnd_fac_z

    if (cg_blend == 2 .and. cg_maxiter < 5) then

       call bl_error("Error: need cg_maxiter >= 5 to do a bisection search on secant iteration failure.")

    endif

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

    ! do we want to force the flux to zero at the boundary?
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

    bnd_fac_z = ONE
    if (idir==3) then
       if ( k3d == domlo(3)   .and. special_bnd_lo .or. &
            k3d == domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac_z = ZERO
       end if
    end if

    tol = cg_tol
    iter_max = cg_maxiter

    call bl_allocate(pstar_hist, 1,iter_max)
    call bl_allocate(pstar_hist_extra, 1,2*iter_max)
    call bl_allocate(us1d, ilo,ihi)

    do j = jlo, jhi

       bnd_fac_y = ONE
       if (idir == 2) then
          if ( j == domlo(2)   .and. special_bnd_lo .or. &
               j == domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac_y = ZERO
          end if
       end if

       do i = ilo, ihi

          ! left state
          rl = max(ql(i,j,kc,QRHO), small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

          pl  = ql(i,j,kc,QPRES)
          rel = ql(i,j,kc,QREINT)
          gcl = gamcl(i,j)

          ! sometime we come in here with negative energy or pressure
          ! note: reset both in either case, to remain thermo
          ! consistent
          if (rel <= ZERO .or. pl < small_pres) then
             print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres

             eos_state % T   = small_temp
             eos_state % rho = rl
             eos_state % xn  = ql(i,j,kc,QFS:QFS-1+nspec)
             eos_state % aux = ql(i,j,kc,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

             rel = rl*eos_state % e
             pl  = eos_state % p
             gcl = eos_state % gam1
          endif

          ! right state
          rr = max(qr(i,j,kc,QRHO), small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

          pr  = qr(i,j,kc,QPRES)
          rer = qr(i,j,kc,QREINT)
          gcr = gamcr(i,j)

          if (rer <= ZERO .or. pr < small_pres) then
             print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres

             eos_state % T   = small_temp
             eos_state % rho = rr
             eos_state % xn  = qr(i,j,kc,QFS:QFS-1+nspec)
             eos_state % aux = qr(i,j,kc,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

             rer = rr*eos_state % e
             pr  = eos_state % p
             gcr = eos_state % gam1
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
          gmin = min(gamel, gamer, ONE, FOUR3RD)
          gmax = max(gamel, gamer, TWO, FIVE3RD)

          game_bar = HALF*(gamel + gamer)
          gamc_bar = HALF*(gcl + gcr)

          gdot = TWO*(ONE - game_bar/gamc_bar)*(game_bar - ONE)

          csmall = smallc(i,j)
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
             if (zp-weakwv*cav(i,j) <= ZERO) then
                zp = dpditer*wl
             endif

             zm = abs(ustar_r - ustar_r_old)
             if (zm-weakwv*cav(i,j) <= ZERO) then
                zm = dpditer*wr
             endif

             ! the new pstar is found via CG Eq. 18
             denom = dpditer/max(zp+zm, small*(cav(i,j)))
             pstar_old = pstar
             pstar = pstar - denom*(ustar_r - ustar_l)
             pstar = max(pstar, small_pres)

             err = abs(pstar - pstar_old)
             if (err < tol*pstar) converged = .true.

             pstar_hist(iter) = pstar

             iter = iter + 1

          enddo

          ! If we failed to converge using the secant iteration, we can either
          ! stop here; or, revert to the original two-shock estimate for pstar;
          ! or do a bisection root find using the bounds established by the most
          ! recent iterations.

          if (.not. converged) then

             if (cg_blend == 0) then

                print *, 'pstar history: '
                do iter = 1, iter_max
                   print *, iter, pstar_hist(iter)
                enddo

                print *, ' '
                print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
                print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr
                print *, 'cav, smallc:', cav(i,j), csmall
                call bl_error("ERROR: non-convergence in the Riemann solver")

             else if (cg_blend == 1) then

                pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)

             else if (cg_blend == 2) then

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
                   enddo
                   do iter = 1, 2 * iter_max
                      print *, iter + iter_max, pstar_hist_extra(iter)
                   enddo

                   print *, ' '
                   print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
                   print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr
                   print *, 'cav, smallc:', cav(i,j), csmall
                   call bl_error("ERROR: non-convergence in the Riemann solver")

                endif

             else

                call bl_error("ERROR: unrecognized cg_blend option.")

             endif

          endif

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
          endif

          ! use tau = 1/rho as the independent variable here
          ro = max(small_dens, ONE/tauo)
          tauo = ONE/ro

          co = sqrt(abs(gamco*po/ro))
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
          ushock = wo/ro - sgnm*uo

          if (pstar-po >= ZERO) then
             spin = ushock
             spout = ushock
          endif

          !if (spout-spin == ZERO) then
          !   scr = small*cav(i,j)
          !else
          !   scr = spout-spin
          !endif
          !frac = (ONE + (spout + spin)/scr)*HALF
          !frac = max(ZERO,min(ONE,frac))

          frac = HALF*(ONE + (spin + spout)/max(spout-spin, spin+spout, small*cav(i,j)))

          ! the transverse velocity states only depend on the
          ! direction that the contact moves
          if (ustar > ZERO) then
             qint(i,j,kc,iv1) = v1l
             qint(i,j,kc,iv2) = v2l
          else if (ustar < ZERO) then
             qint(i,j,kc,iv1) = v1r
             qint(i,j,kc,iv2) = v2r
          else
             qint(i,j,kc,iv1) = HALF*(v1l+v1r)
             qint(i,j,kc,iv2) = HALF*(v2l+v2r)
          endif

          ! linearly interpolate between the star and normal state -- this covers the
          ! case where we are inside the rarefaction fan.
          qint(i,j,kc,GDRHO ) = frac*rstar + (ONE - frac)*ro
          qint(i,j,kc,iu   ) = frac*ustar + (ONE - frac)*uo
          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          gamgdnv =  frac*gamstar + (ONE-frac)*gameo

          ! now handle the cases where instead we are fully in the
          ! star or fully in the original (l/r) state
          if (spout < ZERO) then
             qint(i,j,kc,GDRHO ) = ro
             qint(i,j,kc,iu   ) = uo
             qint(i,j,kc,GDPRES) = po
             gamgdnv = gameo
          endif

          if (spin >= ZERO) then
             qint(i,j,kc,GDRHO ) = rstar
             qint(i,j,kc,iu   ) = ustar
             qint(i,j,kc,GDPRES) = pstar
             gamgdnv = gamstar
          endif

          qint(i,j,kc,GDGAME) = gamgdnv

          qint(i,j,kc,GDPRES) = max(qint(i,j,kc,GDPRES), small_pres)

          u_adv = qint(i,j,kc,iu)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i ==  domlo(1) .or. &
               special_bnd_hi_x .and. i ==  domhi(1)+1 ) then
             bnd_fac_x = ZERO
          else
             bnd_fac_x = ONE
          end if
          u_adv = u_adv * bnd_fac_x*bnd_fac_y*bnd_fac_z

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = qint(i,j,kc,GDRHO)*u_adv

          uflx(i,j,kflux,im1) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iu)
          if (mom_flux_has_p(idir) %  comp(im1)) then
             uflx(i,j,kflux,im1) = uflx(i,j,kflux,im1) + qint(i,j,kc,GDPRES)
          endif
          uflx(i,j,kflux,im2) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv1)
          uflx(i,j,kflux,im3) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv2)

#ifdef HYBRID_MOMENTUM
          call compute_hybrid_flux(qint(i,j,kc,:), uflx(i,j,kflux,:), idir, [i, j, k3d])
#endif

          ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
          rhoetot = qint(i,j,kc,GDPRES)/(gamgdnv - ONE) + &
               HALF*qint(i,j,kc,GDRHO)*(qint(i,j,kc,iu)**2 + qint(i,j,kc,iv1)**2 + qint(i,j,kc,iv2)**2)

          uflx(i,j,kflux,UEDEN) = u_adv*(rhoetot + qint(i,j,kc,GDPRES))
          uflx(i,j,kflux,UEINT) = u_adv*qint(i,j,kc,GDPRES)/(gamgdnv - ONE)

          us1d(i) = ustar
       end do

       ! advected quantities -- only the contact matters
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nqp = qpass_map(ipassive)

          do i = ilo, ihi
             if (us1d(i) > ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nqp)
             else if (us1d(i) < ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nqp)
             else
                qavg = HALF * (ql(i,j,kc,nqp) + qr(i,j,kc,nqp))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo

       enddo
    enddo

    call bl_deallocate(pstar_hist)
    call bl_deallocate(pstar_hist_extra)
    call bl_deallocate(us1d)

  end subroutine riemanncg


  subroutine riemannus(ql, qr, qpd_lo, qpd_hi, &
                       gamcl, gamcr, cav, smallc, gd_lo, gd_hi, &
                       uflx, uflx_lo, uflx_hi, &
                       qint, q_lo, q_hi, &
#ifdef RADIATION
                       lam, lam_lo, lam_hi, &
                       gamcgl, gamcgr, &
                       rflx, rflx_lo, rflx_hi, &
#endif
                       idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, &
                       domlo, domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall, &
                                   mom_flux_has_p
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : compute_hybrid_flux
#endif

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: gd_lo(2), gd_hi(2)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir,ilo,ihi,jlo,jhi
    integer, intent(in) :: domlo(3),domhi(3)

#ifdef RADIATION
    integer, intent(in) :: lam_lo(3), lam_hi(3)
    integer, intent(in) :: rflx_lo(3),rflx_hi(3)
#endif

    real(rt), intent(in) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    real(rt), intent(in) :: gamcl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) :: gamcr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) :: cav(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) :: smallc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

#ifdef RADIATION
    real(rt), intent(in) :: lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    real(rt), intent(in) :: gamcgl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) :: gamcgr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(inout) :: rflx(rflx_lo(1):rflx_hi(1),rflx_lo(2):rflx_hi(2),rflx_lo(3):rflx_hi(3),0:ngroups-1)
#endif

    ! Note:
    !
    !  k3d: the k corresponding to the full 3d array -- it should be
    !       used for print statements or tests against domlo, domhi,
    !       etc
    !
    !  kc: the k corresponding to the 2-wide slab of k-planes, so in
    !      this routine it takes values only of 1 or 2
    !
    !  kflux: used for indexing the uflx array -- in the initial calls
    !         to cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, 
    !         kflux = kc, but in later calls, when uflx = {flux1,flux2,flux3},
    !         kflux = k3d

    integer :: i,j,kc,kflux,k3d
    integer :: n, nqp, ipassive

    real(rt)         :: regdnv
    real(rt)         :: rl, ul, v1l, v2l, pl, rel
    real(rt)         :: rr, ur, v1r, v2r, pr, rer
    real(rt)         :: wl, wr, rhoetot, scr
    real(rt)         :: rstar, cstar, estar, pstar, ustar
    real(rt)         :: ro, uo, po, reo, co, gamco, entho, drho
    real(rt)         :: sgnm, spin, spout, ushock, frac
    real(rt)         :: wsmall, csmall,qavg

#ifdef RADIATION
    real(rt)        , dimension(0:ngroups-1) :: erl, err
    real(rt)         :: reo_g, po_g, co_g, gamco_g
    real(rt)         :: pl_g, rel_g, pr_g, rer_g
    real(rt)         :: regdnv_g, pgdnv_g, pgdnv_t
    real(rt)         :: estar_g, pstar_g
    real(rt)        , dimension(0:ngroups-1) :: lambda, laml, lamr, reo_r, po_r, estar_r, regdnv_r
    real(rt)         :: eddf, f1
    integer :: g
#endif

    real(rt)        , pointer :: us1d(:)

#ifdef ROTATION
    real(rt)         :: vel(3)
#endif

    real(rt)         :: u_adv

    integer :: iu, iv1, iv2, im1, im2, im3
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    real(rt)         :: bnd_fac_x, bnd_fac_y, bnd_fac_z
    real(rt)         :: wwinv, roinv, co2inv

    call bl_allocate(us1d,ilo,ihi)

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

    bnd_fac_z = ONE
    if (idir == 3) then
       if ( k3d == domlo(3)   .and. special_bnd_lo .or. &
            k3d == domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac_z = ZERO
       end if
    end if

    do j = jlo, jhi

       bnd_fac_y = ONE
       if (idir == 2) then
          if ( j == domlo(2)   .and. special_bnd_lo .or. &
               j == domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac_y = ZERO
          end if
       end if

       !dir$ ivdep
       do i = ilo, ihi

#ifdef RADIATION
          if (idir == 1) then
             laml = lam(i-1,j,k3d,:)
          elseif (idir == 2) then
             laml = lam(i,j-1,k3d,:)
          else
             laml = lam(i,j,k3d-1,:)
          end if
          lamr = lam(i,j,k3d,:)
#endif

          rl = max(ql(i,j,kc,QRHO), small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

#ifdef RADIATION
          pl = ql(i,j,kc,qptot)
          rel = ql(i,j,kc,qreitot)
          erl(:) = ql(i,j,kc,qrad:qradhi)
          pl_g = ql(i,j,kc,QPRES)
          rel_g = ql(i,j,kc,QREINT)
#else
          pl  = max(ql(i,j,kc,QPRES ), small_pres)
          rel =     ql(i,j,kc,QREINT)
#endif

          rr = max(qr(i,j,kc,QRHO), small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

#ifdef RADIATION
          pr = qr(i,j,kc,qptot)
          rer = qr(i,j,kc,qreitot)
          err(:) = qr(i,j,kc,qrad:qradhi)
          pr_g = qr(i,j,kc,QPRES)
          rer_g = qr(i,j,kc,QREINT)
#else
          pr  = max(qr(i,j,kc,QPRES), small_pres)
          rer =     qr(i,j,kc,QREINT)
#endif

          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

          wwinv = ONE/(wl + wr)
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
          ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

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
             gamco = gamcl(i,j)
#ifdef RADIATION
             lambda = laml
             po_g = pl_g
             po_r(:) = erl(:) * lambda(:)
             reo_r(:) = erl(:)
             reo_g = rel_g
             gamco_g = gamcgl(i,j)
#endif

          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr(i,j)
#ifdef RADIATION
             lambda = lamr
             po_g = pr_g
             po_r(:) = err(:) * lambda(:)
             reo_r(:) = err(:)
             reo_g = rer_g
             gamco_g = gamcgr(i,j)
#endif

          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             reo = HALF*(rel+rer)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
#ifdef RADIATION
             do g=0, ngroups-1
                lambda(g) = 2.0e0_rt*(laml(g)*lamr(g))/(laml(g)+lamr(g)+1.e-50_rt)
             end do

             reo_r(:) = 0.5e0_rt*(erl(:)+err(:))
             reo_g = 0.5e0_rt*(rel_g+rer_g)
             po_r(:) = lambda(:) * reo_r(:)
             gamco_g = 0.5e0_rt*(gamcgl(i,j)+gamcgr(i,j))
             po_g = 0.5*(pr_g+pl_g)
#endif

          endif

          ro = max(small_dens, ro)

          roinv = ONE/ro

          co = sqrt(abs(gamco*po*roinv))
          co = max(csmall,co)
          co2inv = ONE/(co*co)

          drho = (pstar - po)*co2inv
          rstar = ro + drho
          rstar = max(small_dens, rstar)

#ifdef RADIATION
          estar_g = reo_g + drho*(reo_g + po_g)/ro
          co_g = sqrt(abs(gamco_g*po_g/ro))
          co_g = max(csmall, co_g)
          pstar_g = po_g + drho*co_g**2
          pstar_g = max(pstar_g, small_pres)
          estar_r = reo_r(:) + drho*(reo_r(:) + po_r(:))/ro
#else
          entho = (reo + po)*roinv*co2inv
          estar = reo + (pstar - po)*entho
#endif
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
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif

          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO, min(ONE, frac))

          if (ustar > ZERO) then
             qint(i,j,kc,iv1) = v1l
             qint(i,j,kc,iv2) = v2l

          else if (ustar < ZERO) then
             qint(i,j,kc,iv1) = v1r
             qint(i,j,kc,iv2) = v2r

          else
             qint(i,j,kc,iv1) = HALF*(v1l+v1r)
             qint(i,j,kc,iv2) = HALF*(v2l+v2r)
          endif
          qint(i,j,kc,GDRHO) = frac*rstar + (ONE - frac)*ro
          qint(i,j,kc,iu  ) = frac*ustar + (ONE - frac)*uo

#ifdef RADIATION
          pgdnv_t = frac*pstar + (1.e0_rt - frac)*po
          pgdnv_g = frac*pstar_g + (1.e0_rt - frac)*po_g
          regdnv_g = frac*estar_g + (1.e0_rt - frac)*reo_g
          regdnv_r(:) = frac*estar_r(:) + (1.e0_rt - frac)*reo_r(:)
#else
          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          regdnv = frac*estar + (ONE - frac)*reo
#endif

          if (spout < ZERO) then
             qint(i,j,kc,GDRHO) = ro
             qint(i,j,kc,iu  ) = uo
#ifdef RADIATION
             pgdnv_t = po
             pgdnv_g = po_g
             regdnv_g = reo_g
             regdnv_r = reo_r(:)
#else
             qint(i,j,kc,GDPRES) = po
             regdnv = reo
#endif
          endif

          if (spin >= ZERO) then
             qint(i,j,kc,GDRHO) = rstar
             qint(i,j,kc,iu  ) = ustar
#ifdef RADIATION
             pgdnv_t = pstar
             pgdnv_g = pstar_g
             regdnv_g = estar_g
             regdnv_r = estar_r(:)
#else
             qint(i,j,kc,GDPRES) = pstar
             regdnv = estar
#endif
          endif


#ifdef RADIATION
          do g=0, ngroups-1
             qint(i,j,kc,GDERADS+g) = max(regdnv_r(g), 0.e0_rt)
          end do

          qint(i,j,kc,GDPRES) = pgdnv_g
          qint(i,j,kc,GDLAMS:GDLAMS-1+ngroups) = lambda(:)

          qint(i,j,kc,GDGAME) = pgdnv_g/regdnv_g + ONE
#else
          qint(i,j,kc,GDGAME) = qint(i,j,kc,GDPRES)/regdnv + ONE
          qint(i,j,kc,GDPRES) = max(qint(i,j,kc,GDPRES),small_pres)
#endif

          u_adv = qint(i,j,kc,iu)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i==domlo(1) .or. &
               special_bnd_hi_x .and. i==domhi(1)+1 ) then
             bnd_fac_x = ZERO
          else
             bnd_fac_x = ONE
          end if
          u_adv = u_adv * bnd_fac_x*bnd_fac_y*bnd_fac_z


          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = qint(i,j,kc,GDRHO)*u_adv

          uflx(i,j,kflux,im1) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iu )
          if (mom_flux_has_p(idir) %  comp(im1)) then
             uflx(i,j,kflux,im1) = uflx(i,j,kflux,im1) + qint(i,j,kc,GDPRES)
          endif
          uflx(i,j,kflux,im2) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv1)
          uflx(i,j,kflux,im3) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv2)

#ifdef HYBRID_MOMENTUM
          call compute_hybrid_flux(qint(i,j,kc,:), uflx(i,j,kflux,:), idir, [i, j, k3d])
#endif

#ifdef RADIATION
          rhoetot = regdnv_g + HALF*qint(i,j,kc,GDRHO)*(qint(i,j,kc,iu)**2 + qint(i,j,kc,iv1)**2 + qint(i,j,kc,iv2)**2)

          uflx(i,j,kflux,UEDEN) = u_adv*(rhoetot + pgdnv_g)

          uflx(i,j,kflux,UEINT) = u_adv*regdnv_g
#else
          rhoetot = regdnv + HALF*qint(i,j,kc,GDRHO)*(qint(i,j,kc,iu)**2 + qint(i,j,kc,iv1)**2 + qint(i,j,kc,iv2)**2)

          uflx(i,j,kflux,UEDEN) = u_adv*(rhoetot + qint(i,j,kc,GDPRES))
          uflx(i,j,kflux,UEINT) = u_adv*regdnv
#endif

          ! store this for vectorization
          us1d(i) = ustar

#ifdef RADIATION
          if (fspace_type==1) then
             do g=0,ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5e0_rt*(1.e0_rt-eddf)
                rflx(i,j,kflux,g) = (1.e0_rt+f1) * qint(i,j,kc,GDERADS+g) * u_adv
             end do
          else ! type 2
             do g=0,ngroups-1
                rflx(i,j,kflux,g) = qint(i,j,kc,GDERADS+g) * u_adv
             end do
          end if
#endif
       end do

       ! passively advected quantities
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nqp = qpass_map(ipassive)

          !dir$ ivdep
          do i = ilo, ihi
             if (us1d(i) > ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nqp)

             else if (us1d(i) < ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nqp)

             else
                qavg = HALF * (ql(i,j,kc,nqp) + qr(i,j,kc,nqp))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo

       enddo
    enddo

    call bl_deallocate(us1d)

  end subroutine riemannus


  subroutine HLLC(ql, qr, qpd_lo, qpd_hi, &
                  gamcl, gamcr, cav, smallc, gd_lo, gd_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  qint, q_lo, q_hi, &
                  idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, &
                  domlo, domhi)


    ! this is an implementation of the HLLC solver described in Toro's
    ! book.  it uses the simplest estimate of the wave speeds, since
    ! those should work for a general EOS.  We also initially do the
    ! CGF Riemann construction to get pstar and ustar, since we'll
    ! need to know the pressure and velocity on the interface for the
    ! pdV term in the internal energy update.

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: gd_lo(2), gd_hi(2)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, ilo, ihi, jlo, jhi
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) ::  gamcl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) ::  gamcr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) ::    cav(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) :: smallc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)
    integer, intent(in) :: kc, kflux, k3d

    ! Note:
    !
    !  k3d: the k corresponding to the full 3d array -- it should be
    !       used for print statements or tests against domlo, domhi,
    !       etc
    !
    !  kc: the k corresponding to the 2-wide slab of k-planes, so in
    !      this routine it takes values only of 1 or 2
    !
    !  kflux: used for indexing the uflx array -- in the initial calls
    !         to cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, 
    !         kflux = kc, but in later calls, when uflx = {flux1,flux2,flux3},
    !         kflux = k3d

    integer :: i, j

    real(rt) :: rgdnv, regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel
    real(rt) :: rr, ur, v1r, v2r, pr, rer
    real(rt) :: wl, wr, scr
    real(rt) :: rstar, cstar, estar, pstar, ustar
    real(rt) :: ro, uo, po, reo, co, gamco, entho
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall

    integer :: iu, iv1, iv2
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac
    real(rt) :: wwinv, roinv, co2inv

    real(rt) :: U_hllc_state(nvar), U_state(nvar), F_state(nvar)
    real(rt) :: S_l, S_r, S_c

    if (idir == 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
    else if (idir == 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
    else
       iu = QW
       iv1 = QU
       iv2 = QV
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

    bnd_fac_z = 1
    if (idir == 3) then
       if ( k3d == domlo(3)   .and. special_bnd_lo .or. &
            k3d == domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac_z = 0
       end if
    end if

    do j = jlo, jhi

       bnd_fac_y = 1
       if (idir == 2) then
          if ( j == domlo(2)   .and. special_bnd_lo .or. &
               j == domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac_y = 0
          end if
       end if

       !dir$ ivdep
       do i = ilo, ihi

          rl = max(ql(i,j,kc,QRHO), small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

          pl  = max(ql(i,j,kc,QPRES ), small_pres)
          rel = ql(i,j,kc,QREINT)

          rr = max(qr(i,j,kc,QRHO), small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

          pr  = max(qr(i,j,kc,QPRES), small_pres)
          rer = qr(i,j,kc,QREINT)

          ! now we essentially do the CGF solver to get p and u on the
          ! interface, but we won't use these in any flux construction.
          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall, sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall, sqrt(abs(gamcr(i,j)*pr*rr)))

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
             gamco = gamcl(i,j)

          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr(i,j)
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             reo = HALF*(rel+rer)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
          endif
          ro = max(small_dens, ro)

          roinv = ONE/ro
          co = sqrt(abs(gamco*po*roinv))
          co = max(csmall, co)
          co2inv = ONE/(co*co)

          rstar = ro + (pstar - po)*co2inv
          rstar = max(small_dens, rstar)

          entho = (reo + po)*co2inv/ro
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
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif
          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO, min(ONE, frac))

          rgdnv = frac*rstar + (ONE - frac)*ro
          regdnv = frac*estar + (ONE - frac)*reo

          qint(i,j,kc,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          qint(i,j,kc,GDGAME) = qint(i,j,kc,GDPRES)/regdnv + ONE


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
          S_l = min(ul - sqrt(gamcl(i,j)*pl/rl), ur - sqrt(gamcr(i,j)*pr/rr))
          S_r = max(ul + sqrt(gamcl(i,j)*pl/rl), ur + sqrt(gamcr(i,j)*pr/rr))

          ! estimate of the contact speed -- this is Toro Eq. 10.8
          S_c = (pr - pl + rl*ul*(S_l - ul) - rr*ur*(S_r - ur))/ &
               (rl*(S_l - ul) - rr*(S_r - ur))

          if (S_r <= ZERO) then
             ! R region
             call cons_state(qr(i,j,kc,:), U_state)
             call compute_flux(idir, bnd_fac, U_state, pr, F_state)

          else if (S_r > ZERO .and. S_c <= ZERO) then
             ! R* region
             call cons_state(qr(i,j,kc,:), U_state)
             call compute_flux(idir, bnd_fac, U_state, pr, F_state)

             call HLLC_state(idir, S_r, S_c, qr(i,j,kc,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

          else if (S_c > ZERO .and. S_l < ZERO) then
             ! L* region
             call cons_state(ql(i,j,kc,:), U_state)
             call compute_flux(idir, bnd_fac, U_state, pl, F_state)

             call HLLC_state(idir, S_l, S_c, ql(i,j,kc,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

          else
             ! L region
             call cons_state(ql(i,j,kc,:), U_state)
             call compute_flux(idir, bnd_fac, U_state, pl, F_state)

          endif

          uflx(i,j,kflux,:) = F_state(:)
       enddo
    enddo

  end subroutine HLLC


end module actual_riemann_module
