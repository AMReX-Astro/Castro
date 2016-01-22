module riemann_module

  use bl_types
  use bl_constants_module
  use riemann_util_module

  use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QFS, &
                                 QFX, URHO, UMX, UMY, UEDEN, UEINT, &
                                 GDPRES, GDGAME, GDERADS, GDLAMS, ngdnv, &
                                 small_dens, small_pres, small_temp, &
                                 cg_maxiter, cg_tol, &
                                 npassive, upass_map, qpass_map, &
                                 riemann_solver, ppm_temp_fix, hybrid_riemann, &
                                 allow_negative_energy, use_colglaz

#ifdef RADIATION
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, fspace_type
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif

  implicit none

  private

  public cmpflx, shock

  real (kind=dp_t), parameter :: smallu = 1.e-12_dp_t

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(qm,qp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                    flx,flx_l1,flx_l2,flx_h1,flx_h2, &
                    qint, qg_l1,qg_l2,qg_h1,qg_h2, &
#ifdef RADIATION
                    lam,lam_l1,lam_l2,lam_h1,lam_h2, &
                    rflx,rflx_l1,rflx_l2,rflx_h1,rflx_h2, &
                    gamcg, &
#endif
                    gamc,csml,c,qd_l1,qd_l2,qd_h1,qd_h2, &
                    shk,s_l1,s_l2,s_h1,s_h2, &
                    idir,ilo,ihi,jlo,jhi,domlo,domhi)

    use eos_type_module
    use eos_module

#ifdef RADIATION
    integer, intent(in) :: lam_l1,lam_l2,lam_h1,lam_h2
    integer, intent(in) :: rflx_l1,rflx_l2,rflx_h1,rflx_h2
#endif

    integer, intent(in) :: qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer, intent(in) :: flx_l1,flx_l2,flx_h1,flx_h2
    integer, intent(in) :: qg_l1,qg_l2,qg_h1,qg_h2
    integer, intent(in) :: qd_l1,qd_l2,qd_h1,qd_h2

    integer, intent(in) :: s_l1,s_l2,s_h1,s_h2
    integer, intent(in) :: idir,ilo,ihi,jlo,jhi
    integer, intent(in) :: domlo(2),domhi(2)

#ifdef RADIATION
    double precision, intent(inout) :: lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
    double precision, intent(inout) :: rflx(rflx_l1:rflx_h1,rflx_l2:rflx_h2,0:ngroups-1)
#endif
    double precision, intent(inout) :: qint(qg_l1:qg_h1,qg_l2:qg_h2,ngdnv)

    double precision, intent(inout) ::  qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision, intent(inout) ::  qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision, intent(inout) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,NVAR)

    double precision, intent(in) :: gamc(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision, intent(in) ::    c(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision, intent(in) :: csml(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision, intent(in) ::  shk( s_l1: s_h1, s_l2: s_h2)

    ! Local variables
    integer i, j

    double precision, allocatable :: smallc(:,:), cavg(:,:)
    double precision, allocatable :: gamcm(:,:), gamcp(:,:)
#ifdef RADIATION    
    double precision, allocatable :: gamcgm(:,:), gamcgp(:,:)
#endif
    
    integer :: imin, imax, jmin, jmax
    integer :: is_shock
    double precision :: cl, cr
    type (eos_t) :: eos_state

    allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )
#ifdef RADIATION   
    allocate ( gamcgm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate ( gamcgp(ilo-1:ihi+1,jlo-1:jhi+1) )
#endif

#ifdef RADIATION
    if (hybrid_riemann == 1) then
       call bl_error("ERROR: hybrid Riemann not supported for radiation")
    endif

    if (use_colglaz == 1) then
       call bl_error("ERROR: the Colella-Glaz Riemann solver is not supported for radiation")
    endif

    if (riemann_solver > 0) then
       call bl_error("ERROR: only the CGF Riemann solver is supported for radiation")
    endif
#endif

    if (idir == 1) then
       do j = jlo, jhi
          do i = ilo, ihi+1
             smallc(i,j) = max( csml(i,j), csml(i-1,j) )
             cavg(i,j) = HALF*( c(i,j) + c(i-1,j) )
             gamcm(i,j) = gamc(i-1,j)
             gamcp(i,j) = gamc(i,j)
#ifdef RADIATION
             gamcgm(i,j) = gamcg(i-1,j)
             gamcgp(i,j) = gamcg(i,j)
#endif
          enddo
       enddo

    else
       do j = jlo, jhi+1
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j), csml(i,j-1) )
             cavg(i,j) = HALF*( c(i,j) + c(i,j-1) )
             gamcm(i,j) = gamc(i,j-1)
             gamcp(i,j) = gamc(i,j)
#ifdef RADIATION
             gamcgm(i,j) = gamcg(i,j-1)
             gamcgp(i,j) = gamcg(i,j)
#endif
          enddo
       enddo
    endif

    if (ppm_temp_fix == 2) then
       ! recompute the thermodynamics on the interface to make it
       ! all consistent -- THIS PROBABLY DOESN"T WORK WITH RADIATION

       ! we want to take the edge states of rho, p, and X, and get
       ! new values for gamc and (rho e) on the edges that are
       ! thermodynamically consistent.

       if (idir == 1) then
          imin = ilo
          imax = ihi+1
          jmin = jlo
          jmax = jhi
       else
          imin = ilo
          imax = ihi
          jmin = jlo
          jmax = jhi+1
       endif

       do j = jmin, jmax
          do i = imin, imax

             ! this is an initial guess for iterations, since we
             ! can't be certain that temp is on interfaces
             eos_state%T = 10000.0d0

             ! minus state
             eos_state % rho = qm(i,j,QRHO)
             eos_state % p   = qm(i,j,QPRES)
             eos_state % e   = qm(i,j,QREINT)/qm(i,j,QRHO)
             eos_state % xn  = qm(i,j,QFS:QFS-1+nspec)
             eos_state % aux = qm(i,j,QFX:QFX-1+naux)

             ! Protect against negative energies

             if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
                eos_state % T = small_temp
                call eos(eos_input_rt, eos_state)
             else
                call eos(eos_input_re, eos_state)
             endif

             qm(i,j,QREINT) = qm(i,j,QRHO)*eos_state%e
             qm(i,j,QPRES) = eos_state%p
             gamcm(i,j) = eos_state%gam1


             ! plus state
             eos_state % rho = qp(i,j,QRHO)
             eos_state % p   = qp(i,j,QPRES)
             eos_state % e   = qp(i,j,QREINT)/qp(i,j,QRHO)
             eos_state % xn  = qp(i,j,QFS:QFS-1+nspec)
             eos_state % aux = qp(i,j,QFX:QFX-1+naux)

             ! Protect against negative energies

             if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
                eos_state % T = small_temp
                call eos(eos_input_rt, eos_state)
             else
                call eos(eos_input_re, eos_state)
             endif

             qp(i,j,QREINT) = qp(i,j,QRHO)*eos_state%e
             qp(i,j,QPRES) = eos_state%p
             gamcp(i,j) = eos_state%gam1

          enddo
       enddo

    endif

    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    if (riemann_solver == 0) then
       ! Colella, Glaz, & Ferguson solver
       call riemannus(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                      gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                      flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                      qint, qg_l1, qg_l2, qg_h1, qg_h2, &
#ifdef RADIATION
                      lam,lam_l1,lam_l2,lam_h1,lam_h2,&
                      gamcgm, gamcgp, &
                      rflx, rflx_l1, rflx_l2, rflx_h1, rflx_h2, &
#endif
                      idir, ilo, ihi, jlo, jhi, domlo, domhi)

    elseif (riemann_solver == 1) then
       ! Colella & Glaz solver
       call riemanncg(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                      gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                      flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                      qint, qg_l1, qg_l2, qg_h1, qg_h2, &
                      idir, ilo, ihi, jlo, jhi, domlo, domhi)

    elseif (riemann_solver == 2) then
       ! HLLC
       call HLLC(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                 gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                 flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                 qint, qg_l1, qg_l2, qg_h1, qg_h2, &
                 idir, ilo, ihi, jlo, jhi, domlo, domhi)
    else
       call bl_error("ERROR: invalid value of riemann_solver")
    endif

    if (hybrid_riemann == 1) then
       ! correct the fluxes using an HLL scheme if we are in a shock
       ! and doing the hybrid approach
       if (idir == 1) then
          imin = ilo
          imax = ihi+1
          jmin = jlo
          jmax = jhi
       else
          imin = ilo
          imax = ihi
          jmin = jlo
          jmax = jhi+1
       endif

       do j = jmin, jmax
          do i = imin, imax

             if (idir == 1) then
                is_shock = shk(i-1,j) + shk(i,j)
             else
                is_shock = shk(i,j-1) + shk(i,j)
             endif

             if (is_shock >= 1) then

                if (idir == 1) then
                   cl = c(i-1,j)
                   cr = c(i,j)
                else
                   cl = c(i,j-1)
                   cr = c(i,j)
                endif

                call HLL(qm(i,j,:), qp(i,j,:), cl, cr, &
                         idir, 2, flx(i,j,:))

             endif

          enddo
       enddo

    endif

    deallocate(smallc,cavg,gamcm,gamcp)
#ifdef RADIATION
    deallocate(gamcgm,gamcgp)
#endif

  end subroutine cmpflx


  subroutine shock(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                   shk,s_l1,s_l2,s_h1,s_h2, &
                   ilo1,ilo2,ihi1,ihi2,dx,dy)

    use prob_params_module, only : coord_type

    integer, intent(in) :: qd_l1, qd_l2, qd_h1, qd_h2
    integer, intent(in) :: s_l1, s_l2, s_h1, s_h2
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    double precision, intent(in) :: dx, dy
    double precision, intent(in) :: q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision, intent(inout) :: shk(s_l1:s_h1,s_l2:s_h2)

    integer :: i, j

    double precision :: divU
    double precision :: px_pre, px_post, py_pre, py_post
    double precision :: e_x, e_y, d
    double precision :: p_pre, p_post, pjump

    double precision, parameter :: small = 1.d-10
    double precision, parameter :: eps = 0.33d0

    ! This is a basic multi-dimensional shock detection algorithm.
    ! This implementation follows Flash, which in turn follows
    ! AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
    !
    ! The spirit of this follows the shock detection in Colella &
    ! Woodward (1984)

    if (coord_type /= 0) then
       call bl_error("ERROR: invalid geometry in shock()")
    endif

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          ! construct div{U}
          divU = HALF*(q(i+1,j,QU) - q(i-1,j,QU))/dx + &
                 HALF*(q(i,j+1,QV) - q(i,j-1,QV))/dy

          ! find the pre- and post-shock pressures in each direction
          if (q(i+1,j,QPRES) - q(i-1,j,QPRES) < ZERO) then
             px_pre  = q(i+1,j,QPRES)
             px_post = q(i-1,j,QPRES)
          else
             px_pre  = q(i-1,j,QPRES)
             px_post = q(i+1,j,QPRES)
          endif

          if (q(i,j+1,QPRES) - q(i,j-1,QPRES) < ZERO) then
             py_pre  = q(i,j+1,QPRES)
             py_post = q(i,j-1,QPRES)
          else
             py_pre  = q(i,j-1,QPRES)
             py_post = q(i,j+1,QPRES)
          endif

          ! use compression to create unit vectors for the shock direction
          e_x = (q(i+1,j,QU) - q(i-1,j,QU))**2
          e_y = (q(i,j+1,QV) - q(i,j-1,QV))**2
          d = ONE/(e_x + e_y + small)

          e_x = e_x*d
          e_y = e_y*d

          ! project the pressures onto the shock direction
          p_pre  = e_x*px_pre + e_y*py_pre
          p_post = e_x*px_post + e_y*py_post

          ! test for compression + pressure jump to flag a shock
          pjump = eps - (p_post - p_pre)/p_pre

          if (pjump < ZERO .and. divU < ZERO) then
             shk(i,j) = ONE
          else
             shk(i,j) = ZERO
          endif

       enddo
    enddo

  end subroutine shock


! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemanncg(ql,qr,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       gamcl,gamcr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                       uflx,uflx_l1,uflx_l2,uflx_h1,uflx_h2, &
                       qint, qg_l1,qg_l2,qg_h1,qg_h2, &
                       idir,ilo1,ihi1,ilo2,ihi2,domlo,domhi)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use bl_error_module
    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use prob_params_module, only : coord_type

    double precision, parameter:: small = 1.d-8

    integer :: qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer :: gd_l1,gd_l2,gd_h1,gd_h2
    integer :: uflx_l1,uflx_l2,uflx_h1,uflx_h2
    integer :: qg_l1,qg_l2,qg_h1,qg_h2
    integer :: idir,ilo1,ihi1,ilo2,ihi2
    integer :: domlo(2),domhi(2)

    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision ::  gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::  gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::    cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: qint(qg_l1:qg_h1,qg_l2:qg_h2,ngdnv)

    integer :: i,j,ilo,jlo,ihi,jhi, ipassive
    integer :: n, nq

    double precision :: rgdnv,vgdnv,wgdnv,ustar,gamgdnv
    double precision :: rl, ul, vl, v2l, pl, rel
    double precision :: rr, ur, vr, v2r, pr, rer
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

    double precision, parameter :: weakwv = 1.d-3

    double precision, allocatable :: pstar_hist(:)

    type (eos_t) :: eos_state

    integer :: iu, iv1, iv2

    tol = cg_tol
    iter_max = cg_maxiter

    !  set min/max based on normal direction
    if (idir == 1) then
       ilo = ilo1
       ihi = ihi1 + 1
       jlo = ilo2
       jhi = ihi2

       iu = QU
       iv1 = QV
       iv2 = QW
    else
       ilo = ilo1
       ihi = ihi1
       jlo = ilo2
       jhi = ihi2+1

       iu = QV
       iv1 = QU
       iv2 = QW
    endif

    allocate (pstar_hist(iter_max))

    do j = jlo, jhi
       do i = ilo, ihi

          ! left state
          rl = max(ql(i,j,QRHO),small_dens)

          ! pick left velocities based on direction
          ul = ql(i,j,iu)
          vl = ql(i,j,iv1)
          v2l = ql(i,j,iv2)

          pl  = ql(i,j,QPRES )
          rel = ql(i,j,QREINT)
          gcl = gamcl(i,j)

          ! sometimes we come in here with negative energy or pressure
          ! note: reset both in either case, to remain thermo
          ! consistent
          if (rel <= ZERO .or. pl <= small_pres) then
             print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres
             eos_state % T   = small_temp
             eos_state % rho = rl
             eos_state % xn  = ql(i,j,QFS:QFS-1+nspec)
             eos_state % aux = ql(i,j,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

             rel = rl*eos_state%e
             pl  = eos_state%p
             gcl = eos_state%gam1
          endif

          ! right state
          rr = max(qr(i,j,QRHO),small_dens)

          ! pick right velocities based on direction
          if (idir == 1) then
             ur = qr(i,j,QU)
             vr = qr(i,j,QV)
             v2r = qr(i,j,QW)
          else
             ur = qr(i,j,QV)
             vr = qr(i,j,QU)
             v2r = qr(i,j,QW)
          endif

          pr  = qr(i,j,QPRES)
          rer = qr(i,j,QREINT)
          gcr = gamcr(i,j)

          if (rer <= ZERO .or. pr <= small_pres) then
             print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres
             eos_state % T   = small_temp
             eos_state % rho = rr
             eos_state % xn  = qr(i,j,QFS:QFS-1+nspec)
             eos_state % aux = qr(i,j,QFX:QFX-1+naux)

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
             if(zp-weakwv*cav(i,j) <= ZERO)then
                zp = dpditer*wl
             endif

             zm=abs(ustarm-ustnm1)
             if(zm-weakwv*cav(i,j) <= ZERO)then
                zm = dpditer*wr
             endif

             ! the new pstar is found via CG Eq. 18
             denom=dpditer/max(zp+zm,small*(cav(i,j)))
             pstnm1 = pstar
             pstar=pstar-denom*(ustarm-ustarp)
             pstar=max(pstar,small_pres)

             err = abs(pstar - pstnm1)
             if (err < tol*pstar) converged = .true.

             pstar_hist(iter) = pstar

             iter = iter + 1

          enddo

          if (.not. converged) then
             print *, 'pstar history: '
             do iter = 1, iter_max
                print *, iter, pstar_hist(iter)
             enddo

             print *, ' '
             print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
             print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr
             call bl_error("ERROR: non-convergence in the Riemann solver")
          endif


          ! we converged!  construct the single ustar for the region
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
          if (ustar .gt. ZERO) then
             ro = rl
             uo = ul
             po = pl
             tauo = taul
             !reo = rel
             gamco = gcl
             gameo = gamel

          else if (ustar .lt. ZERO) then
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
             !reo = HALF*(rel+rer)
             gamco = HALF*(gcl+gcr)
             gameo = HALF*(gamel + gamer)
          endif

          ! use tau = 1/rho as the independent variable here
          ro = max(small_dens,ONE/tauo)
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
          !rstar=max(ONE-ro*dpjmp/wosq, (gameo-ONE)/(gameo+ONE))
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

          if (pstar-po .ge. ZERO) then
             spin = ushock
             spout = ushock
          endif
          ! if (spout-spin .eq. ZERO) then
          !    scr = small*cav(i,j)
          ! else
          !    scr = spout-spin
          ! endif
          ! frac = (ONE + (spout + spin)/scr)*HALF
          ! frac = max(ZERO,min(ONE,frac))

          frac = HALF*(ONE + (spin + spout)/max(spout-spin,spin+spout, small*cav(i,j)))

          ! the transverse velocity states only depend on the
          ! direction that the contact moves
          if (ustar .gt. ZERO) then
             vgdnv = vl
             wgdnv = v2l
          else if (ustar .lt. ZERO) then
             vgdnv = vr
             wgdnv = v2r
          else
             vgdnv = HALF*(vl+vr)
             wgdnv = HALF*(v2l+v2r)
          endif

          ! linearly interpolate between the star and normal state -- this covers the
          ! case where we are inside the rarefaction fan.
          rgdnv = frac*rstar + (ONE - frac)*ro
          qint(i,j,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,GDPRES) = frac*pstar + (ONE - frac)*po
          gamgdnv =  frac*gamstar + (ONE-frac)*gameo

          ! now handle the cases where instead we are fully in the
          ! star or fully in the original (l/r) state
          if (spout .lt. ZERO) then
             rgdnv = ro
             qint(i,j,iu) = uo
             qint(i,j,GDPRES) = po
             gamgdnv = gameo
          endif
          if (spin .ge. ZERO) then
             rgdnv = rstar
             qint(i,j,iu) = ustar
             qint(i,j,GDPRES) = pstar
             gamgdnv = gamstar
          endif

          qint(i,j,GDGAME) = gamgdnv

          qint(i,j,GDPRES) = max(qint(i,j,GDPRES),small_pres)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          qint(i,j,iu) = bc_test(idir, i, j, domlo, domhi) * qint(i,j,iu)

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,URHO) = rgdnv*qint(i,j,iu)

          ! note: for axisymmetric geometries, we do not include the
          ! pressure in the r-direction, since div{F} + grad{p} cannot
          ! be written in a flux difference form
          if (idir == 1) then
             uflx(i,j,UMX) = uflx(i,j,URHO)*qint(i,j,iu)
             uflx(i,j,UMY) = uflx(i,j,URHO)*vgdnv
             if (coord_type == 0) then
                uflx(i,j,UMX) = uflx(i,j,UMX) + qint(i,j,GDPRES)
             endif
          else
             uflx(i,j,UMX) = uflx(i,j,URHO)*vgdnv
             uflx(i,j,UMY) = uflx(i,j,URHO)*qint(i,j,iu) + qint(i,j,GDPRES)
          endif

          ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
          rhoetot = qint(i,j,GDPRES)/(gamgdnv - ONE) + &
               HALF*rgdnv*(qint(i,j,iu)**2 + vgdnv**2 + wgdnv**2)

          uflx(i,j,UEDEN) = qint(i,j,iu)*(rhoetot + qint(i,j,GDPRES))
          uflx(i,j,UEINT) = qint(i,j,iu)*qint(i,j,GDPRES)/(gamgdnv - ONE)

          ! advected quantities -- only the contact matters
          ! note: this includes the z-velocity flux
          do ipassive = 1, npassive
             n  = upass_map(ipassive)
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

       enddo
    enddo

    deallocate(pstar_hist)

  end subroutine riemanncg

! :::
! ::: ------------------------------------------------------------------
! :::


  subroutine riemannus(ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                       gamcl, gamcr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
                       uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
                       qint, qg_l1, qg_l2, qg_h1, qg_h2, &
#ifdef RADIATION
                       lam,lam_l1,lam_l2,lam_h1,lam_h2, &
                       gamcgl,gamcgr, &
                       rflx, rflx_l1, rflx_l2, rflx_h1, rflx_h2, &
#endif
                       idir, ilo1, ihi1, ilo2, ihi2, domlo, domhi)

    use prob_params_module, only : coord_type

    double precision, parameter:: small = 1.d-8

    integer :: qpd_l1, qpd_l2, qpd_h1, qpd_h2
    integer :: gd_l1, gd_l2, gd_h1, gd_h2
    integer :: uflx_l1, uflx_l2, uflx_h1, uflx_h2
    integer :: qg_l1, qg_l2, qg_h1, qg_h2
#ifdef RADIATION
    integer :: lam_l1, lam_l2, lam_h1, lam_h2
    integer :: rflx_l1, rflx_l2, rflx_h1, rflx_h2
#endif
    integer :: idir, ilo1, ihi1, ilo2, ihi2
    integer :: domlo(2),domhi(2)

#ifdef RADIATION
    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
#else
    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
#endif

    double precision :: gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: qint(qg_l1:qg_h1,qg_l2:qg_h2,ngdnv)
#ifdef RADIATION
    double precision :: lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
    double precision :: gamcgl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: gamcgr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: rflx(rflx_l1:rflx_h1,rflx_l2:rflx_h2,0:ngroups-1)
#endif

    integer :: ilo,ihi,jlo,jhi
    integer :: n, nq
    integer :: i, j, ipassive

#ifdef RADIATION
    integer :: g
#endif

    double precision :: rgd, vgd, wgd, regd, ustar
#ifdef RADIATION
    double precision, dimension(0:ngroups-1) :: erl, err
#endif
    double precision :: rl, ul, vl, v2l, pl, rel
    double precision :: rr, ur, vr, v2r, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, estar, pstar
    double precision :: ro, uo, po, reo, co, gamco, entho, drho
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg

#ifdef RADIATION
    double precision :: regdnv_g, pgdnv_g, pgdnv_t
    double precision :: estar_g, pstar_g
    double precision, dimension(0:ngroups-1) :: lambda, reo_r, po_r, estar_r, regdnv_r
    double precision :: eddf, f1
    double precision :: co_g, gamco_g, pl_g, po_g, pr_g, rel_g, reo_g, rer_g
#endif

    integer :: iu, iv1, iv2

    !************************************************************
    !  set min/max based on normal direction
    if (idir == 1) then
       ilo = ilo1
       ihi = ihi1 + 1
       jlo = ilo2
       jhi = ihi2

       iu = QU
       iv1 = QV
       iv2 = QW
    else
       ilo = ilo1
       ihi = ihi1
       jlo = ilo2
       jhi = ihi2+1

       iu = QV
       iv1 = QU
       iv2 = QW
    endif

    do j = jlo, jhi
       do i = ilo, ihi

          rl = ql(i,j,QRHO)

          !  pick left velocities based on direction
          ul = ql(i,j,iu)
          vl = ql(i,j,iv1)
          v2l = ql(i,j,iv2)

#ifdef RADIATION
          pl = ql(i,j,QPTOT)
          rel = ql(i,j,QREITOT)
#else
          pl = ql(i,j,QPRES)
          rel = ql(i,j,QREINT)
#endif

#ifdef RADIATION
          erl(:) = ql(i,j,qrad:qradhi)
          pl_g = ql(i,j,QPRES)
          rel_g = ql(i,j,QREINT)
#endif

          rr = qr(i,j,QRHO)

          !  pick right velocities based on direction
          if (idir == 1) then
             ur = qr(i,j,QU)
             vr = qr(i,j,QV)
             v2r = qr(i,j,QW)
          else
             ur = qr(i,j,QV)
             vr = qr(i,j,QU)
             v2r = qr(i,j,QW)
          endif

#ifdef RADIATION
          pr = qr(i,j,QPTOT)
          rer = qr(i,j,QREITOT)
#else
          pr = qr(i,j,QPRES)
          rer = qr(i,j,QREINT)
#endif

#ifdef RADIATION
          err(:) = qr(i,j,qrad:qradhi)
          pr_g = qr(i,j,QPRES)
          rer_g = qr(i,j,QREINT)
#endif

          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

          pstar = max(pstar,small_pres)
          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

          if (ustar .gt. ZERO) then
#ifdef RADIATION
             if (idir == 1) then
                lambda(:) = lam(i-1,j,:)
             else
                lambda(:) = lam(i,j-1,:)
             end if
#endif

             ro = rl
             uo = ul
             po = pl

#ifdef RADIATION
             po_g = pl_g
             po_r(:) = erl(:) * lambda(:)
#endif

             reo = rel
             gamco = gamcl(i,j)

#ifdef RADIATION
             reo_r(:) = erl(:)
             reo_g = rel_g
             gamco_g = gamcgl(i,j)
#endif

          else if (ustar .lt. ZERO) then
#ifdef RADIATION
             lambda(:) = lam(i,j,:)
#endif

             ro = rr
             uo = ur
             po = pr

#ifdef RADIATION
             po_g = pr_g
             po_r(:) = err(:) * lambda(:)
#endif

             reo = rer
             gamco = gamcr(i,j)

#ifdef RADIATION
             reo_r(:) = err(:)
             reo_g = rer_g
             gamco_g = gamcgr(i,j)
#endif

          else
#ifdef RADIATION
             if (idir == 1) then
                do g = 0, ngroups-1
                   lambda(g) = 2.0d0*(lam(i-1,j,g)*lam(i,j,g))/ &
                        (lam(i-1,j,g)+lam(i,j,g)+1.d-50)
                end do
             else
                do g = 0, ngroups-1
                   lambda(g) = 2.0d0*(lam(i,j-1,g)*lam(i,j,g))/ &
                        (lam(i,j-1,g)+lam(i,j,g)+1.d-50)
                end do
             end if
#endif

             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)

#ifdef RADIATION
             po_r(:) = lambda(:) * reo_r(:)
             po_g = HALF*(pr_g+pl_g)
#endif

             reo = HALF*(rel+rer)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))

#ifdef RADIATION
             reo_r(:) = HALF*(erl(:)+err(:))
             reo_g = HALF*(rel_g+rer_g)
             gamco_g = HALF*(gamcgl(i,j)+gamcgr(i,j))
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

          sgnm = sign(ONE,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar

          ushock = HALF*(spin + spout)

          if (pstar-po >= ZERO) then
             spin = ushock
             spout = ushock
          endif

          if (spout-spin == ZERO) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif

          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))

          if (ustar .gt. ZERO) then
             vgd = vl
             wgd = v2l
          else if (ustar .lt. ZERO) then
             vgd = vr
             wgd = v2r
          else
             vgd = HALF*(vl+vr)
             wgd = HALF*(v2l+v2r)
          endif

          rgd = frac*rstar + (ONE - frac)*ro

          qint(i,j,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,iv1) = vgd
          qint(i,j,iv2) = wgd

#ifdef RADIATION
          pgdnv_t = frac*pstar + (ONE - frac)*po
          pgdnv_g = frac*pstar_g + (ONE - frac)*po_g
          regdnv_g = frac*estar_g + (ONE - frac)*reo_g
          regdnv_r(:) = frac*estar_r(:) + (ONE - frac)*reo_r(:)
#else
          qint(i,j,GDPRES) = frac*pstar + (ONE - frac)*po
          regd = frac*estar + (ONE - frac)*reo
#endif

          if (spout < ZERO) then
             rgd = ro
             qint(i,j,iu) = uo
#ifdef RADIATION
             pgdnv_t = po
             pgdnv_g = po_g
             regdnv_g = reo_g
             regdnv_r(:) = reo_r(:)
#else
             qint(i,j,GDPRES) = po
             regd = reo
#endif
          endif

          if (spin >= ZERO) then
             rgd = rstar
             qint(i,j,iu) = ustar
#ifdef RADIATION
             pgdnv_t = pstar
             pgdnv_g = pstar_g
             regdnv_g = estar_g
             regdnv_r(:) = estar_r(:)
#else
             qint(i,j,GDPRES) = pstar
             regd = estar
#endif
          endif

          ! not sure what this should be for radiation?
#ifdef RADIATION
          qint(i,j,GDGAME) = pgdnv_g/regdnv_g + ONE
#else
          qint(i,j,GDGAME) = qint(i,j,GDPRES)/regd + ONE
#endif

          ! enforce that the fluxes through a symmetry plane or wall are zero
          qint(i,j,iu) = bc_test(idir, i, j, domlo, domhi) * qint(i,j,iu)

#ifdef RADIATION
          do g=0, ngroups-1
             qint(i,j,GDERADS+g) = max(regdnv_r(g), ZERO)
          end do

          qint(i,j,GDPRES) = pgdnv_g

          qint(i,j,GDLAMS:GDLAMS-1+ngroups) = lambda(:)
#endif

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,URHO) = rgd*qint(i,j,iu)

          ! note: for axisymmetric geometries, we do not include the
          ! pressure in the r-direction, since div{F} + grad{p} cannot
          ! be written in a flux difference form
          if (idir == 1) then
             uflx(i,j,UMX) = uflx(i,j,URHO)*qint(i,j,iu)
             uflx(i,j,UMY) = uflx(i,j,URHO)*vgd
             if (coord_type == 0) then
                uflx(i,j,UMX) = uflx(i,j,UMX) + qint(i,j,GDPRES)
             endif
          else
             uflx(i,j,UMX) = uflx(i,j,URHO)*vgd
             uflx(i,j,UMY) = uflx(i,j,URHO)*qint(i,j,iu) + qint(i,j,GDPRES)
          endif

#ifdef RADIATION
          rhoetot = regdnv_g + HALF*rgd*(qint(i,j,iu)**2 + vgd**2 + wgd**2)

          uflx(i,j,UEDEN) = qint(i,j,iu)*(rhoetot + pgdnv_g)
          uflx(i,j,UEINT) = qint(i,j,iu)*regdnv_g
#else
          rhoetot = regd + HALF*rgd*(qint(i,j,iu)**2 + vgd**2 + wgd**2)

          uflx(i,j,UEDEN) = qint(i,j,iu)*(rhoetot + qint(i,j,GDPRES))
          uflx(i,j,UEINT) = qint(i,j,iu)*regd
#endif

#ifdef RADIATION
          if (fspace_type == 1) then
             do g = 0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                rflx(i,j,g) = (ONE+f1) * qint(i,j,GDERADS+g) * qint(i,j,iu)
             end do
          else ! type 2
             do g = 0, ngroups-1
                rflx(i,j,g) = qint(i,j,GDERADS+g) * qint(i,j,iu)
             end do
          end if
#endif

          do ipassive = 1, npassive
             n  = upass_map(ipassive)
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

       enddo
    enddo
  end subroutine riemannus


! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine HLLC(ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                  gamcl, gamcr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
                  uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
                  qint, qg_l1, qg_l2, qg_h1, qg_h2, &
                  idir, ilo1, ihi1, ilo2, ihi2, domlo, domhi)

    ! this is an implementation of the HLLC solver described in Toro's
    ! book.  it uses the simplest estimate of the wave speeds, since
    ! those should work for a general EOS.  We also initially do the
    ! CGF Riemann construction to get pstar and ustar, since we'll need
    ! to know the pressure and velocity on the interface for the grad p
    ! term in momentum and for an internal energy update

    double precision, parameter:: small = 1.d-8

    integer :: qpd_l1, qpd_l2, qpd_h1, qpd_h2
    integer :: gd_l1, gd_l2, gd_h1, gd_h2
    integer :: uflx_l1, uflx_l2, uflx_h1, uflx_h2
    integer :: qg_l1, qg_l2, qg_h1, qg_h2
    integer :: idir, ilo1, ihi1, ilo2, ihi2
    integer :: domlo(2),domhi(2)

    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: qint(qg_l1:qg_h1,qg_l2:qg_h2,ngdnv)

    integer :: ilo,ihi,jlo,jhi
    integer :: n, nq
    integer :: i, j
    integer :: bnd_fac
    
    double precision :: rgd, vgd, regd, ustar
    double precision :: rl, ul, pl, rel
    double precision :: rr, ur, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, pstar
    double precision :: ro, uo, po, co, gamco
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall, qavg

    double precision :: U_hllc_state(nvar), U_state(nvar), F_state(nvar)
    double precision :: S_l, S_r, S_c

    integer :: iu, iv1, iv2

    !  set min/max based on normal direction
    if (idir == 1) then
       ilo = ilo1
       ihi = ihi1 + 1
       jlo = ilo2
       jhi = ihi2

       iu = QU
       iv1 = QV
       iv2 = QW
    else
       ilo = ilo1
       ihi = ihi1
       jlo = ilo2
       jhi = ihi2+1

       iu = QV
       iv1 = QU
       iv2 = QW
    endif

    do j = jlo, jhi
       do i = ilo, ihi

          rl = ql(i,j,QRHO)

          ! pick left velocities based on direction
          ! ul is always normal to the interface
          if (idir == 1) then
             ul  = ql(i,j,QU)
          else
             ul  = ql(i,j,QV)
          endif

          pl = ql(i,j,QPRES)
          rel = ql(i,j,QREINT)

          rr = qr(i,j,QRHO)

          ! pick right velocities based on direction
          ! ur is always normal to the interface
          if (idir == 1) then
             ur  = qr(i,j,QU)
          else
             ur  = qr(i,j,QV)
          endif

          pr = qr(i,j,QPRES)
          rer = qr(i,j,QREINT)

          ! now we essentially do the CGF solver to get p and u on the
          ! interface, but we won't use these in any flux construction.
          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

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
             gamco = gamcl(i,j)

          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             gamco = gamcr(i,j)

          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
          endif

          ro = max(small_dens,ro)

          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)

          rstar = ro + (pstar - po)/co**2
          rstar = max(small_dens,rstar)

          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)

          sgnm = sign(ONE,ustar)
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
          frac = max(ZERO,min(ONE,frac))

          qint(i,j,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,GDPRES) = frac*pstar + (ONE - frac)*po

          ! TODO
          !gegdnv(i,j) = pgdnv(i,j)/regd + ONE

          ! now we do the HLLC construction

          bnd_fac = bc_test(idir, i, j, domlo, domhi)
          
          ! use the simplest estimates of the wave speeds
          S_l = min(ul - sqrt(gamcl(i,j)*pl/rl), ur - sqrt(gamcr(i,j)*pr/rr))
          S_r = max(ul + sqrt(gamcl(i,j)*pl/rl), ur + sqrt(gamcr(i,j)*pr/rr))

          ! estimate of the contact speed -- this is Toro Eq. 10.8
          S_c = (pr - pl + rl*ul*(S_l - ul) - rr*ur*(S_r - ur))/ &
             (rl*(S_l - ul) - rr*(S_r - ur))

          if (S_r <= ZERO) then
             ! R region
             call cons_state(qr(i,j,:), U_state)
             call compute_flux(idir, 2, bnd_fac, U_state, pr, F_state)

          else if (S_r > ZERO .and. S_c <= ZERO) then
             ! R* region
             call cons_state(qr(i,j,:), U_state)
             call compute_flux(idir, 2, bnd_fac, U_state, pr, F_state)

             call HLLC_state(idir, S_r, S_c, qr(i,j,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

          else if (S_c > ZERO .and. S_l < ZERO) then
             ! L* region
             call cons_state(ql(i,j,:), U_state)
             call compute_flux(idir, 2, bnd_fac, U_state, pl, F_state)

             call HLLC_state(idir, S_l, S_c, ql(i,j,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

          else
             ! L region
             call cons_state(ql(i,j,:), U_state)
             call compute_flux(idir, 2, bnd_fac, U_state, pl, F_state)

          endif


          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          ! and store the fluxes
          uflx(i,j,:) = F_state(:)

       enddo
    enddo
  end subroutine HLLC

end module riemann_module
