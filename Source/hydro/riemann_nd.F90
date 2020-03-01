module riemann_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use meth_params_module, only : NQ, NVAR, NQAUX, &
                                 URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UFX, UTEMP, &
                                 QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QFS, QFX, &
                                 QC, QGAMC, QGC, &
                                 GDRHO, GDU, GDV, GDW, GDPRES, &
                                 NGDNV, GDRHO, GDPRES, &

#ifdef RADIATION
                                 qrad, qptot, qreitot, &
                                 GDERADS, QGAMCG, QLAMS, QREITOT, GDLAMS, &
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
#ifdef RADIATION
                                 rflx, rflx_lo, rflx_hi, &
#endif
                                 qgdnv, qg_lo, qg_hi, &
                                 qaux, qa_lo, qa_hi, &
                                 shk, s_lo, s_hi, &
                                 idir, domlo, domhi) bind(C, name="cmpflx_plus_godunov")

    use network, only: nspec, naux
    use castro_error_module
    use amrex_fort_module, only : rt => amrex_real
    use prob_params_module, only : physbc_lo, physbc_hi, &
         Symmetry, SlipWall, NoSlipWall

    implicit none

    ! note: lo, hi necessarily the limits of the valid (no ghost
    ! cells) domain, but could be hi+1 in some dimensions.  We rely on
    ! the caller to specific the interfaces over which to solve the
    ! Riemann problems

    integer, intent(in) :: lo(3), hi(3)

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: qg_lo(3), qg_hi(3)

    integer, intent(in), value :: idir

    integer, intent(in) :: domlo(3),domhi(3)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)

    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)

#ifdef RADIATION
    integer, intent(in) :: rflx_lo(3), rflx_hi(3)
    real(rt), intent(inout) :: rflx(rflx_lo(1):rflx_hi(1), rflx_lo(2):rflx_hi(2), rflx_lo(3):rflx_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) ::  shk(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    real(rt), intent(inout) :: qgdnv(qg_lo(1):qg_hi(1), qg_lo(2):qg_hi(2), qg_lo(3):qg_hi(3), NGDNV)

    integer :: i, j, k
    real(rt) :: cl, cr
    real(rt) :: qleft(NQ), qright(NQ), qint_zone(NQ), flx_zone(NVAR)
    real(rt) :: gamcl, gamcr
    real(rt) :: csmall, cavg, bnd_fac
    integer :: is_shock

#ifdef RADIATION
    real(rt) :: gamcgl, gamcgr
    real(rt) :: laml(0:ngroups-1), lamr(0:ngroups-1), lambda_int_zone(0:ngroups-1)
#endif

    logical :: special_bnd_lo, special_bnd_hi

    !$gpu

    ! do we want to force the flux to zero at the boundary?
    special_bnd_lo = (physbc_lo(idir) == Symmetry &
         .or.         physbc_lo(idir) == SlipWall &
         .or.         physbc_lo(idir) == NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) == Symmetry &
         .or.         physbc_hi(idir) == SlipWall &
         .or.         physbc_hi(idir) == NoSlipWall)

    if (riemann_solver < 2) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                call interface_input_states(i, j, k, &
                                            qm, qm_lo, qm_hi, &
                                            qp, qp_lo, qp_hi, &
                                            qaux, qa_lo, qa_hi, &
                                            idir, 0, &
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

                else if (riemann_solver == 1) then
                   ! Colella & Glaz solver

#ifndef RADIATION
                   call riemanncg(qleft, qright, &
                                  gamcl, gamcr, &
                                  csmall, cavg, &
                                  bnd_fac, &
                                  qint_zone, &
                                  idir)
#endif

                end if

                qgdnv(i,j,k,GDRHO) = qint_zone(QRHO)
                qgdnv(i,j,k,GDU) = qint_zone(QU)
                qgdnv(i,j,k,GDV) = qint_zone(QV)
                qgdnv(i,j,k,GDW) = qint_zone(QW)
                qgdnv(i,j,k,GDPRES) = qint_zone(QPRES)
#ifdef RADIATION
                qgdnv(i,j,k,GDLAMS:GDLAMS-1+ngroups) = lambda_int_zone(:)
                qgdnv(i,j,k,GDERADS:GDERADS-1+ngroups) = qint_zone(QRAD:QRAD-1+ngroups)
#endif

                call compute_flux_q_single(i, j, k, &
                                           qint_zone, &
                                           flx, flx_lo, flx_hi, &
#ifdef RADIATION
                                           lambda_int_zone, &
                                           rflx, rflx_lo, rflx_hi, &
#endif
                                           idir, 0)

             end do
          end do
       end do

    else if (riemann_solver == 2) then

       ! HLLC
       call HLLC(qm, qm_lo, qm_hi, &
                 qp, qp_lo, qp_hi, &
                 qaux, qa_lo, qa_hi, &
                 flx, flx_lo, flx_hi, &
                 qgdnv, qg_lo, qg_hi, &
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

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

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

             end do
          end do
       end do

    end if

  end subroutine cmpflx_plus_godunov

  subroutine compute_flux_q(lo, hi, &
                            qint, q_lo, q_hi, &
                            F, F_lo, F_hi, &
#ifdef RADIATION
                            lambda, l_lo, l_hi, &
                            rF, rF_lo, rF_hi, &
#endif
                            idir, enforce_eos) bind(C, name="compute_flux_q")

    ! given a primitive state, compute the flux in direction idir
    !

    use prob_params_module, only : mom_flux_has_p
    use meth_params_module, only : NQ, NVAR, NQAUX, &
                                   URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, &
#ifdef SHOCK_VAR
                                   USHK, &
#endif
                                   QRHO, QU, QV, QW, &
                                   QPRES, QREINT, &
                                   QGAMC, QFS, QFX, &
#ifdef HYBRID_MOMENTUM
                                   NGDNV, GDPRES, &
                                   GDRHO, GDU, GDV, GDW, &
#endif
#ifdef RADIATION
                                   QRAD, fspace_type, &
                                   GDERADS, GDLAMS, &
#endif
                                   npassive, upass_map, qpass_map, T_guess
#ifdef RADIATION
    use fluxlimiter_module, only: Edd_factor ! function
    use rad_params_module, only : ngroups
#endif
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : compute_hybrid_flux
#endif
    use eos_type_module, only : eos_t, eos_input_rp
    use eos_module, only : eos
    use network, only : nspec, naux

    implicit none

    integer, intent(in), value :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: F_lo(3), F_hi(3)
#ifdef RADIATION
    integer, intent(in) :: l_lo(3), l_hi(3)
    integer, intent(in) :: rF_lo(3), rF_hi(3)
#endif

    real(rt), intent(in) :: qint(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(out) :: F(F_lo(1):F_hi(1), F_lo(2):F_hi(2), F_lo(3):F_hi(3), NVAR)
#ifdef RADIATION
    real(rt), intent(in) :: lambda(l_lo(1):l_hi(1), l_lo(2):l_hi(2), l_lo(3):l_hi(3), 0:ngroups-1)
    real(rt), intent(out) :: rF(rF_lo(1):rF_hi(1), rF_lo(2):rF_hi(2), rF_lo(3):rF_hi(3), 0:ngroups-1)
#endif
    integer, intent(in), value :: enforce_eos
    integer, intent(in) :: lo(3), hi(3)

    integer :: i, j, k

    real(rt) :: qint_zone(NQ)
#ifdef RADIATION
    real(rt) :: lambda_zone(0:ngroups-1)
#endif

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             qint_zone(:) = qint(i,j,k,:)
#ifdef RADIATION
             lambda_zone(:) = lambda(i,j,k,:)
#endif
             call compute_flux_q_single(i, j, k, &
                                        qint_zone, &
                                        F, F_lo, F_hi, &
#ifdef RADIATION
                                        lambda_zone, &
                                        rF, rF_lo, rF_hi, &
#endif
                                        idir, enforce_eos)
          end do
       end do
    end do

  end subroutine compute_flux_q


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
    real(rt) :: bnd_fac
    real(rt) :: csmall, cavg
    real(rt) :: gamcl, gamcr
#ifdef RADIATION
    real(rt) :: laml(0:ngroups-1), lamr(0:ngroups-1), lambda_int_zone(0:ngroups-1)
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
                               qint_zone, &
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
