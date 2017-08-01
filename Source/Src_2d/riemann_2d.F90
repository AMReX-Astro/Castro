module riemann_module

  use bl_types
  use bl_constants_module
  use riemann_util_module
  use meth_params_module, only : NQ, NQAUX, NVAR, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QFS, &
                                 QFX, URHO, UMX, UMY, UEDEN, UEINT, &
                                 GDPRES, GDGAME, QGAMC, QC, QCSML, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, fspace_type, &
                                 GDERADS, GDLAMS, QGAMCG, QLAMS, &
#endif
                                 NGDNV, small_dens, small_pres, small_temp, &
                                 cg_maxiter, cg_tol, cg_blend, &
                                 npassive, upass_map, qpass_map, &
                                 riemann_solver, ppm_temp_fix, hybrid_riemann, &
                                 allow_negative_energy

#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public cmpflx

  real(rt), parameter :: smallu = 1.e-12_rt

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(qm, qp, qpd_lo, qpd_hi, &
                    flx, flx_lo, flx_hi, &
                    qint, qg_lo, qg_hi, &
#ifdef RADIATION
                    rflx, rflx_lo, rflx_hi, &
#endif
                    qaux, qa_lo, qa_hi, &
                    shk, s_lo, s_hi, &
                    idir, ilo, ihi, jlo, jhi, domlo, domhi)

    use actual_riemann_module
    use eos_type_module, only: eos_input_re, eos_input_rt, eos_t
    use eos_module, only: eos
    use network, only: nspec, naux
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: qg_lo(3), qg_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: idir,ilo,ihi,jlo,jhi
    integer, intent(in) :: domlo(2),domhi(2)

#ifdef RADIATION
    integer, intent(in) :: rflx_lo(3), rflx_hi(3)
    real(rt)        , intent(inout) :: rflx(rflx_lo(1):rflx_hi(1),rflx_lo(2):rflx_hi(2),0:ngroups-1)
#endif

    real(rt)        , intent(inout) :: qint(qg_lo(1):qg_hi(1),qg_lo(2):qg_hi(2),NGDNV)

    real(rt)        , intent(inout) ::  qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),NQ)
    real(rt)        , intent(inout) ::  qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),NQ)
    real(rt)        , intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),NVAR)

    real(rt)        , intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)
    real(rt)        , intent(in) ::  shk( s_lo(1): s_hi(1), s_lo(2): s_hi(2))

    ! Local variables
    integer i, j

    real(rt)        , allocatable :: smallc(:,:), cavg(:,:)
    real(rt)        , allocatable :: gamcm(:,:), gamcp(:,:)
#ifdef RADIATION
    real(rt)        , allocatable :: gamcgm(:,:), gamcgp(:,:), lam(:,:,:)
#endif

    ! these will refer to the zone interfaces that we solve the
    ! Riemann problem across
    integer :: imin, imax, jmin, jmax

    integer :: is_shock
    real(rt)         :: cl, cr
    type (eos_t) :: eos_state

    allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )
#ifdef RADIATION
    allocate ( gamcgm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate ( gamcgp(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (    lam(ilo-1:ihi+1,jlo-1:jhi+1,0:ngroups-1) )
#endif

#ifdef RADIATION
    if (hybrid_riemann == 1) then
       call bl_error("ERROR: hybrid Riemann not supported for radiation")
    endif

    if (riemann_solver > 0) then
       call bl_error("ERROR: only the CGF Riemann solver is supported for radiation")
    endif
#endif

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

    if (idir == 1) then
       do j = jmin, jmax
          do i = imin, imax
             smallc(i,j) = max( qaux(i,j,QCSML), qaux(i-1,j,QCSML) )
             cavg(i,j) = HALF*( qaux(i,j,QC) + qaux(i-1,j,QC) )
             gamcm(i,j) = qaux(i-1,j,QGAMC)
             gamcp(i,j) = qaux(i,j,QGAMC)
#ifdef RADIATION
             gamcgm(i,j) = qaux(i-1,j,QGAMCG)
             gamcgp(i,j) = qaux(i,j,QGAMCG)
#endif
          enddo
       enddo

    else
       do j = jmin, jmax
          do i = imin, imax
             smallc(i,j) = max( qaux(i,j,QCSML), qaux(i,j-1,QCSML) )
             cavg(i,j) = HALF*( qaux(i,j,QC) + qaux(i,j-1,QC) )
             gamcm(i,j) = qaux(i,j-1,QGAMC)
             gamcp(i,j) = qaux(i,j,QGAMC)
#ifdef RADIATION
             gamcgm(i,j) = qaux(i,j-1,QGAMCG)
             gamcgp(i,j) = qaux(i,j,QGAMCG)
#endif
          enddo
       enddo
    endif

#ifdef RADIATION
    do j = jlo-1, jhi+1
       do i = ilo-1, ihi+1
          lam(i,j,:) = qaux(i,j,QLAMS:QLAMS+ngroups-1)
       enddo
    enddo
#endif

    if (ppm_temp_fix == 2) then
       ! recompute the thermodynamics on the interface to make it
       ! all consistent -- THIS PROBABLY DOESN"T WORK WITH RADIATION

       ! we want to take the edge states of rho, p, and X, and get
       ! new values for gamc and (rho e) on the edges that are
       ! thermodynamically consistent.

       do j = jmin, jmax
          do i = imin, imax

             ! this is an initial guess for iterations, since we
             ! can't be certain that temp is on interfaces
             eos_state%T = 10000.0e0_rt

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
       call riemannus(qm, qp, qpd_lo, qpd_hi, &
                      gamcm, gamcp, cavg, smallc, [ilo-1, jlo-1, 0], [ihi+1, jhi+1, 0], &
                      flx, flx_lo, flx_hi, &
                      qint, qg_lo, qg_hi, &
#ifdef RADIATION
                      lam, [ilo-1, jlo-1, 0], [ihi+1, jhi+1, 0], &
                      gamcgm, gamcgp, &
                      rflx, rflx_lo, rflx_hi, &
#endif
                      idir, imin, imax, jmin, jmax, 0, 0, 0, &
                      [domlo(1), domlo(2), 0], [domhi(1), domhi(2), 0])

    elseif (riemann_solver == 1) then
       ! Colella & Glaz solver
       call riemanncg(qm, qp, qpd_lo, qpd_hi, &
                      gamcm, gamcp, cavg, smallc, [ilo-1, jlo-1, 0], [ihi+1, jhi+1, 0], &
                      flx, flx_lo, flx_hi, &
                      qint, qg_lo, qg_hi, &
                      idir, imin, imax, jmin, jmax, 0, 0, 0, &
                      [domlo(1), domlo(2), 0], [domhi(1), domhi(2), 0])

    elseif (riemann_solver == 2) then
       ! HLLC
       call HLLC(qm, qp, qpd_lo, qpd_hi, &
                 gamcm, gamcp, cavg, smallc, [ilo-1, jlo-1, 0], [ihi+1, jhi+1, 0], &
                 flx, flx_lo, flx_hi, &
                 qint, qg_lo, qg_hi, &
                 idir, imin, imax, jmin, jmax, domlo, domhi)
    else
       call bl_error("ERROR: invalid value of riemann_solver")
    endif

    if (hybrid_riemann == 1) then
       ! correct the fluxes using an HLL scheme if we are in a shock
       ! and doing the hybrid approach

       do j = jmin, jmax
          do i = imin, imax

             if (idir == 1) then
                is_shock = shk(i-1,j) + shk(i,j)
             else
                is_shock = shk(i,j-1) + shk(i,j)
             endif

             if (is_shock >= 1) then

                if (idir == 1) then
                   cl = qaux(i-1,j,QC)
                   cr = qaux(i,j,QC)
                else
                   cl = qaux(i,j-1,QC)
                   cr = qaux(i,j,QC)
                endif

                call HLL(qm(i,j,:), qp(i,j,:), cl, cr, &
                         idir, flx(i,j,:))

             endif

          enddo
       enddo

    endif

    deallocate(smallc,cavg,gamcm,gamcp)
#ifdef RADIATION
    deallocate(gamcgm,gamcgp,lam)
#endif

  end subroutine cmpflx


! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine HLLC(ql, qr, qpd_lo, qpd_hi, &
                  gamcl, gamcr, cav, smallc, gd_lo, gd_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  qint, qg_lo, qg_hi, &
                  idir, ilo, ihi, jlo, jhi, domlo, domhi)

    ! this is an implementation of the HLLC solver described in Toro's
    ! book.  it uses the simplest estimate of the wave speeds, since
    ! those should work for a general EOS.  We also initially do the
    ! CGF Riemann construction to get pstar and ustar, since we'll need
    ! to know the pressure and velocity on the interface for the grad p
    ! term in momentum and for an internal energy update

    use amrex_fort_module, only : rt => amrex_real
    real(rt)        , parameter:: small = 1.e-8_rt

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: gd_lo(3), gd_hi(3)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: qg_lo(3), qg_hi(3)
    integer, intent(in) :: idir 
    integer, intent(in) :: ilo, ihi, jlo, jhi
    integer, intent(in) :: domlo(2),domhi(2)

    real(rt), intent(in) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),NQ)

    real(rt), intent(in) :: gamcl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) :: gamcr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) :: cav(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(in) :: smallc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),NVAR)
    real(rt), intent(inout) :: qint(qg_lo(1):qg_hi(1),qg_lo(2):qg_hi(2),NGDNV)

    integer :: i, j
    integer :: bnd_fac

    !real(rt)         :: regd
    real(rt) :: ustar
    real(rt) :: rl, ul, pl, rel
    real(rt) :: rr, ur, pr, rer
    real(rt) :: wl, wr, scr
    real(rt) :: rstar, cstar, pstar
    real(rt) :: ro, uo, po, co, gamco
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall

    real(rt) :: U_hllc_state(nvar), U_state(nvar), F_state(nvar)
    real(rt) :: S_l, S_r, S_c

    integer :: iu, iv1, iv2

    !  set min/max based on normal direction
    if (idir == 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
    else
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
             call compute_flux(idir, bnd_fac, U_state, pr, F_state)

          else if (S_r > ZERO .and. S_c <= ZERO) then
             ! R* region
             call cons_state(qr(i,j,:), U_state)
             call compute_flux(idir, bnd_fac, U_state, pr, F_state)

             call HLLC_state(idir, S_r, S_c, qr(i,j,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

          else if (S_c > ZERO .and. S_l < ZERO) then
             ! L* region
             call cons_state(ql(i,j,:), U_state)
             call compute_flux(idir, bnd_fac, U_state, pl, F_state)

             call HLLC_state(idir, S_l, S_c, ql(i,j,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

          else
             ! L region
             call cons_state(ql(i,j,:), U_state)
             call compute_flux(idir, bnd_fac, U_state, pl, F_state)

          endif


          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          ! and store the fluxes
          uflx(i,j,:) = F_state(:)

       enddo
    enddo
  end subroutine HLLC

end module riemann_module
