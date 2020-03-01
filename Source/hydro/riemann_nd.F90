module riemann_module

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

  use riemann_solvers_module
  use riemann_util_module

#ifdef RADIATION
  use rad_params_module, only : ngroups
#endif

  implicit none

  real(rt), parameter :: small = 1.e-8_rt

  private :: small

contains

  subroutine cmpflx_plus_godunov(lo, hi, &
                                 qm, qm_lo, qm_hi, &
                                 qp, qp_lo, qp_hi, &
                                 flx, flx_lo, flx_hi, &
                                 qint, q_lo, q_hi, &
#ifdef RADIATION
                                 rflx, rflx_lo, rflx_hi, &
                                 lambda_int, li_lo, li_hi, &
#endif
                                 qgdnv, qg_lo, qg_hi, &
                                 qaux, qa_lo, qa_hi, &
                                 shk, s_lo, s_hi, &
                                 idir, domlo, domhi) bind(C, name="cmpflx_plus_godunov")

    use network, only: nspec, naux
    use castro_error_module
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    ! note: lo, hi necessarily the limits of the valid (no ghost
    ! cells) domain, but could be hi+1 in some dimensions.  We rely on
    ! the caller to specific the interfaces over which to solve the
    ! Riemann problems

    integer, intent(in) :: lo(3), hi(3)

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: qg_lo(3), qg_hi(3)

    integer, intent(in), value :: idir

    integer, intent(in) :: domlo(3),domhi(3)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)

    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

#ifdef RADIATION
    integer, intent(in) :: rflx_lo(3), rflx_hi(3)
    real(rt), intent(inout) :: rflx(rflx_lo(1):rflx_hi(1), rflx_lo(2):rflx_hi(2), rflx_lo(3):rflx_hi(3),0:ngroups-1)
    integer, intent(in) :: li_lo(3), li_hi(3)
    real(rt), intent(inout) :: lambda_int(li_lo(1),li_hi(1), li_lo(2):li_hi(2), li_lo(3):li_hi(3), 0:ngroups-1)
#endif

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) ::  shk(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    real(rt), intent(inout) :: qgdnv(qg_lo(1):qg_hi(1), qg_lo(2):qg_hi(2), qg_lo(3):qg_hi(3), NGDNV)

    integer :: i, j, k
    real(rt) :: cl, cr
    real(rt) :: ql_zone(NQ), qr_zone(NQ), flx_zone(NVAR)
    integer :: is_shock

    !$gpu

    ! do we want to force the flux to zero at the boundary?
    special_bnd_lo = (physbc_lo(idir) == Symmetry &
         .or.         physbc_lo(idir) == SlipWall &
         .or.         physbc_lo(idir) == NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) == Symmetry &
         .or.         physbc_hi(idir) == SlipWall &
         .or.         physbc_hi(idir) == NoSlipWall)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             call interface_input_states(i, j, k, &
                                         qm, qm_lo, qm_hi, &
                                         qp, qp_lo, qp_hi, &
                                         qaux, qa_lo, qa_hi, &
                                         idir, compute_gammas, &
                                         qleft, qright, &
                                         gamcl, gamcr, &
#ifdef RADIATION
                                         laml, lamr, &
                                         gamcgl, gamcgr, &
#endif
                                         csmall, cavg, bnd_fac, &
                                         special_bnd_lo, special_bnd_hi, &
                                         domlo, domhi)

             if (riemann_solver == 0) then
                ! Colella, Glaz, & Ferguson solver

                call riemannus(qleft, qright, &
                               gamcl, gamcr, &
#ifdef RADIATION
                               laml, lamr, &
                               gamcgl, gamcgr, &
#endif
                               csmall, cavg, &
                               bnd_fac, &
                               qint_zone, &
#ifdef RADIATION
                               lambda_int_zone, &
#endif
                               idir)

                qint(i,j,k,:) = qint_zone(:)
#ifdef RADIATION
                lambda_int(i,j,k,:) = lambda_int_zone(:)
#endif

                call compute_flux_q(lo, hi, &
                                    qint, q_lo, q_hi, &
                                    flx, flx_lo, flx_hi, &
#ifdef RADIATION
                                    lambda_int, q_lo, q_hi, &
                                    rflx, rflx_lo, rflx_hi, &
#endif
                                    idir, 0)


             else if (riemann_solver == 1) then
                ! Colella & Glaz solver

#ifndef RADIATION
                call riemanncg(qleft, qright, &
                               gamcl, gamcr, &
                               csmall, cavg, &
                               bnd_fac, &
                               qint_zone_zone, &
                               idir)

                qint(i,j,k,:) = qint_zone(:)
#ifdef RADIATION
                lambda_int(i,j,k,:) = lambda_int_zone(:)
#endif

                call compute_flux_q(lo, hi, &
                                    qint, q_lo, q_hi, &
                                    flx, flx_lo, flx_hi, &
#ifdef RADIATION
                                    lambda_int, q_lo, q_hi, &
                                    rflx, rflx_lo, rflx_hi, &
#endif
                                    idir, 0)


             else if (riemann_solver == 2) then

                ! HLLC
                call HLLC(qm, qm_lo, qm_hi, &
                          qp, qp_lo, qp_hi, &
                          qaux, qa_lo, qa_hi, &
                          flx, flx_lo, flx_hi, &
                          qint, q_lo, q_hi, &
                          idir, lo, hi, &
                          domlo, domhi)
#ifndef AMREX_USE_CUDA
             else
                call castro_error("ERROR: invalid value of riemann_solver")
#endif
             endif


             if (hybrid_riemann == 1) then
                ! correct the fluxes using an HLL scheme if we are in a shock
                ! and doing the hybrid approach

                is_shock = 0

                select case (idir)
                case (1)
                   is_shock = shk(i-1,j,k) + shk(i,j,k)
                case (2)
                   is_shock = shk(i,j-1,k) + shk(i,j,k)
                case (3)
                   is_shock = shk(i,j,k-1) + shk(i,j,k)
                end select

                if (is_shock >= 1) then

                   select case (idir)
                   case (1)
                      cl = qaux(i-1,j,k,QC)
                      cr = qaux(i,j,k,QC)
                   case (2)
                      cl = qaux(i,j-1,k,QC)
                      cr = qaux(i,j,k,QC)
                   case (3)
                      cl = qaux(i,j,k-1,QC)
                      cr = qaux(i,j,k,QC)
                   end select

                   call HLL(qleft, qright, cl, cr, idir, flx_zone)
                   flx(i,j,k,:) = flx_zone(:)
                end if

             end if

          end do
       end do
    end do

  end subroutine cmpflx_plus_godunov


  subroutine interface_input_states(i, j, k, &
                                    qm, qm_lo, qm_hi, &
                                    qp, qp_lo, qp_hi, &
                                    qaux, qa_lo, qa_hi, &
                                    idir, compute_gammas, &
                                    qleft, qright, &
                                    gamcl, gamcr, &
#ifdef RADIATION
                                    laml, lamr, &
                                    gamcgl, gamcgr, &
#endif
                                    csmall, cavg, bnd_fac, &
                                    special_bnd_lo, special_bnd_hi, &
                                    domlo, domhi) bind(C, name="riemann_state")

             ! Extract this zone's interface values

             qleft(:) = qm(i,j,k,:)
             qright(:) = qp(i,j,k,:)

#ifdef RADIATION
             if (idir == 1) then
                laml(:) = qaux(i-1,j,k,QLAMS:QLAMS+ngroups-1)
             else if (idir == 2) then
                laml(:) = qaux(i,j-1,k,QLAMS:QLAMS+ngroups-1)
             else
                laml(:) = qaux(i,j,k-1,QLAMS:QLAMS+ngroups-1)
             end if
             lamr(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)

#endif

             if (ppm_temp_fix == 2) then
                ! recompute the thermodynamics on the interface to make it
                ! all consistent

                ! we want to take the edge states of rho, e, and X, and get
                ! new values for p on the edges that are
                ! thermodynamically consistent.


                ! this is an initial guess for iterations, since we
                ! can't be certain what temp is on interfaces
                eos_state % T = T_guess

                ! minus state
                eos_state % rho = qleft(QRHO)
                eos_state % p   = qleft(QPRES)
                eos_state % e   = qleft(QREINT)/qleft(QRHO)
                eos_state % xn  = qleft(QFS:QFS+nspec-1)
                eos_state % aux = qleft(QFX:QFX+naux-1)

                call eos(eos_input_re, eos_state)

                qleft(QREINT) = eos_state % e * eos_state % rho
                qleft(QPRES)  = eos_state % p

                ! plus state
                eos_state % rho = qright(QRHO)
                eos_state % p   = qright(QPRES)
                eos_state % e   = qright(QREINT)/qright(QRHO)
                eos_state % xn  = qright(QFS:QFS+nspec-1)
                eos_state % aux = qright(QFX:QFX+naux-1)

                call eos(eos_input_re, eos_state)

                qright(QREINT) = eos_state % e * eos_state % rho
                qright(QPRES)  = eos_state % p

             end if

             if (idir == 1) then
                csmall = max( small, small * max(qaux(i,j,k,QC), qaux(i-1,j,k,QC)))
                cavg = HALF*(qaux(i,j,k,QC) + qaux(i-1,j,k,QC))
                gamcl = qaux(i-1,j,k,QGAMC)
#ifdef RADIATION
                gamcgl = qaux(i-1,j,k,QGAMCG)
#endif
             else if (idir == 2) then
                csmall = max( small, small * max(qaux(i,j,k,QC), qaux(i,j-1,k,QC)))
                cavg = HALF*(qaux(i,j,k,QC) + qaux(i,j-1,k,QC))
                gamcl = qaux(i,j-1,k,QGAMC)
#ifdef RADIATION
                gamcgl = qaux(i,j-1,k,QGAMCG)
#endif
             else
                csmall = max( small, small * max(qaux(i,j,k,QC), qaux(i,j,k-1,QC)))
                cavg = HALF*(qaux(i,j,k,QC) + qaux(i,j,k-1,QC))
                gamcl = qaux(i,j,k-1,QGAMC)
#ifdef RADIATION
                gamcgl = qaux(i,j,k-1,QGAMCG)
#endif
             end if
             gamcr = qaux(i,j,k,QGAMC)
#ifdef RADIATION
             gamcgr = qaux(i,j,k,QGAMCG)
#endif

             ! override the gammas if needed
#ifndef RADIATION
             if (use_reconstructed_gamma1 == 1) then
                gamcl = qleft(QGC)
                gamcr = qright(QGC)

             else if (compute_gammas == 1) then
                ! we come in with a good p, rho, and X on the interfaces
                ! -- use this to find the gamma used in the sound speed
                eos_state % p = qleft(QPRES)
                eos_state % rho = qleft(QRHO)
                eos_state % xn(:) = qleft(QFS:QFS-1+nspec)
                eos_state % T = T_guess ! initial guess
                eos_state % aux(:) = qleft(QFX:QFX-1+naux)

                call eos(eos_input_rp, eos_state)

                gamcl = eos_state % gam1

                eos_state % p = qright(QPRES)
                eos_state % rho = qright(QRHO)
                eos_state % xn(:) = qright(QFS:QFS-1+nspec)
                eos_state % T = T_guess ! initial guess
                eos_state % aux(:) = qright(QFX:QFX-1+naux)

                call eos(eos_input_rp, eos_state)

                gamcr = eos_state % gam1

             end if
#endif

             ! deal with hard walls
             bnd_fac = 1.0_rt

             if (idir == 1) then
                if ((i == domlo(1) .and. special_bnd_lo) .or. &
                    (i == domhi(1)+1 .and. special_bnd_hi)) then
                   bnd_fac = 0.0_rt
                end if
             else if (idir == 2) then
                if ((j == domlo(2) .and. special_bnd_lo) .or. &
                    (j == domhi(2)+1 .and. special_bnd_hi)) then
                   bnd_fac = 0.0_rt
                end if
             else
                if ((k == domlo(3) .and. special_bnd_lo) .or. &
                    (k == domhi(3)+1 .and. special_bnd_hi)) then
                   bnd_fac = 0.0_rt
                end if
             end if


  end subroutine interface_input_states


  subroutine riemann_state(lo, hi, &
                           qm, qm_lo, qm_hi, &
                           qp, qp_lo, qp_hi, &
                           qint, q_lo, q_hi, &
#ifdef RADIATION
                           lambda_int, l_lo, l_hi, &
#endif
                           qaux, qa_lo, qa_hi, &
                           idir, compute_gammas, &
                           domlo, domhi) bind(C, name="riemann_state")

    ! just compute the hydrodynamic state on the interfaces
    ! don't compute the fluxes

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re, eos_input_rp
    use network, only: nspec, naux
    use castro_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : hybrid_riemann, ppm_temp_fix, riemann_solver, &
                                   T_guess
    use prob_params_module, only : physbc_lo, physbc_hi, &
         Symmetry, SlipWall, NoSlipWall

    implicit none

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    integer, intent(in), value :: idir, compute_gammas
    ! note: lo, hi are not necessarily the limits of the valid (no
    ! ghost cells) domain, but could be hi+1 in some dimensions.  We
    ! rely on the caller to specific the interfaces over which to
    ! solve the Riemann problems
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)

    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
#ifdef RADIATION
    integer, intent(in) :: l_lo(3), l_hi(3)
    real(rt), intent(inout) :: lambda_int(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    ! local variables

    integer i, j, k

    type (eos_t) :: eos_state

    real(rt) :: qleft(NQ), qright(NQ)
    real(rt) :: qint_zone(NQ)
    real(rt) :: lambda_int_zone(0:ngroups-1)
    real(rt) :: bnd_fac
    real(rt) :: csmall, cavg
    real(rt) :: gamcl, gamcr
#ifdef RADIATION
    real(rt) :: gamcgl, gamcgr
#endif

    logical :: special_bnd_lo, special_bnd_hi

    !$gpu

#ifdef RADIATION
#ifndef AMREX_USE_CUDA
    if (hybrid_riemann == 1) then
       call castro_error("ERROR: hybrid Riemann not supported for radiation")
    endif

    if (riemann_solver > 0) then
       call castro_error("ERROR: only the CGF Riemann solver is supported for radiation")
    endif
#endif
#endif

#if AMREX_SPACEDIM == 1
#ifndef AMREX_USE_CUDA
    if (riemann_solver > 1) then
       call castro_error("ERROR: HLLC not implemented for 1-d")
    endif
#endif
#endif

    ! do we want to force the flux to zero at the boundary?
    special_bnd_lo = (physbc_lo(idir) == Symmetry &
         .or.         physbc_lo(idir) == SlipWall &
         .or.         physbc_lo(idir) == NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) == Symmetry &
         .or.         physbc_hi(idir) == SlipWall &
         .or.         physbc_hi(idir) == NoSlipWall)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             call interface_input_states(i, j, k, &
                                         qm, qm_lo, qm_hi, &
                                         qp, qp_lo, qp_hi, &
                                         qaux, qa_lo, qa_hi, &
                                         idir, compute_gammas, &
                                         qleft, qright, &
                                         gamcl, gamcr, &
#ifdef RADIATION
                                         laml, lamr, &
                                         gamcgl, gamcgr, &
#endif
                                         csmall, cavg, bnd_fac, &
                                         special_bnd_lo, special_bnd_hi, &
                                         domlo, domhi)

             if (riemann_solver == 0) then
                ! Colella, Glaz, & Ferguson solver

                call riemannus(qleft, qright, &
                               gamcl, gamcr, &
#ifdef RADIATION
                               laml, lamr, &
                               gamcgl, gamcgr, &
#endif
                               csmall, cavg, &
                               bnd_fac, &
                               qint_zone, &
#ifdef RADIATION
                               lambda_int_zone, &
#endif
                               idir)

             elseif (riemann_solver == 1) then
                ! Colella & Glaz solver

#ifndef RADIATION
                call riemanncg(qleft, qright, &
                               gamcl, gamcr, &
                               csmall, cavg, &
                               bnd_fac, &
                               qint_zone_zone, &
                               idir)

#else
#ifndef AMREX_USE_CUDA
                call castro_error("ERROR: CG solver does not support radiaiton")
#endif
#endif

#ifndef AMREX_USE_CUDA
             else
                call castro_error("ERROR: invalid value of riemann_solver")
#endif
             endif

             qint(i,j,k,:) = qint_zone(:)
#ifdef RADIATION
             lambda_int(i,j,k,:) = lambda_int_zone(:)
#endif

          end do
       end do
    end do

  end subroutine riemann_state


end module riemann_module
