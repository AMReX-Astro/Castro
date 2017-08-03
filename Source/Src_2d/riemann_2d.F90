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
    real(rt), intent(inout) :: rflx(rflx_lo(1):rflx_hi(1),rflx_lo(2):rflx_hi(2),0:ngroups-1)
#endif

    real(rt), intent(inout) :: qint(qg_lo(1):qg_hi(1),qg_lo(2):qg_hi(2),NGDNV)

    real(rt), intent(inout) ::  qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),NQ)
    real(rt), intent(inout) ::  qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),NQ)
    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),NVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)
    real(rt), intent(in) ::  shk( s_lo(1): s_hi(1), s_lo(2): s_hi(2))

    ! Local variables
    integer i, j

    ! these will refer to the zone interfaces that we solve the
    ! Riemann problem across
    integer :: imin, imax, jmin, jmax

    integer :: is_shock
    real(rt)         :: cl, cr
    type (eos_t) :: eos_state

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
             !gamcm(i,j) = eos_state%gam1


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
             !gamcp(i,j) = eos_state%gam1

          enddo
       enddo

    endif

    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    if (riemann_solver == 0) then
       ! Colella, Glaz, & Ferguson solver
       call riemannus(qm, qp, qpd_lo, qpd_hi, &
                      qaux, qa_lo, qa_hi, &
                      flx, flx_lo, flx_hi, &
                      qint, qg_lo, qg_hi, &
#ifdef RADIATION
                      rflx, rflx_lo, rflx_hi, &
#endif
                      idir, imin, imax, jmin, jmax, 0, 0, 0, &
                      [domlo(1), domlo(2), 0], [domhi(1), domhi(2), 0])

    elseif (riemann_solver == 1) then
       ! Colella & Glaz solver
       call riemanncg(qm, qp, qpd_lo, qpd_hi, &
                      qaux, qa_lo, qa_hi, &
                      flx, flx_lo, flx_hi, &
                      qint, qg_lo, qg_hi, &
                      idir, imin, imax, jmin, jmax, 0, 0, 0, &
                      [domlo(1), domlo(2), 0], [domhi(1), domhi(2), 0])

    elseif (riemann_solver == 2) then
       ! HLLC
       call HLLC(qm, qp, qpd_lo, qpd_hi, &
                 qaux, qa_lo, qa_hi, &
                 flx, flx_lo, flx_hi, &
                 qint, qg_lo, qg_hi, &
                 idir, imin, imax, jmin, jmax, 0, 0, 0, &
                 [domlo(1), domlo(2), 0], [domhi(1), domhi(2), 0])

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

  end subroutine cmpflx

end module riemann_module
