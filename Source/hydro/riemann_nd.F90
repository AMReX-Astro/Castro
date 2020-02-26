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

    use eos_module, only: eos
    use eos_type_module, only: eos_t
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

    ! Solve Riemann problem to get the fluxes
    if (riemann_solver == 0 .or. riemann_solver == 1) then
       ! Colella, Glaz, & Ferguson solver

       call riemann_state(lo, hi, &
                          qm, qm_lo, qm_hi, &
                          qp, qp_lo, qp_hi, &
                          qint, q_lo, q_hi, &
#ifdef RADIATION
                          lambda_int, q_lo, q_hi, &
#endif
                          qaux, qa_lo, qa_hi, &
                          idir, 0, &
                          domlo, domhi)

       call compute_flux_q(lo, hi, &
                           qint, q_lo, q_hi, &
                           flx, flx_lo, flx_hi, &
#ifdef RADIATION
                           lambda_int, q_lo, q_hi, &
                           rflx, rflx_lo, rflx_hi, &
#endif
                           idir, 0)

    elseif (riemann_solver == 2) then
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

                   ql_zone(:) = qm(i,j,k,:)
                   qr_zone(:) = qp(i,j,k,:)
                   call HLL(ql_zone, qr_zone, cl, cr, idir, flx_zone)
                   flx(i,j,k,:) = flx_zone(:)
                endif

             end do
          end do
       end do

    endif

    call ca_store_godunov_state(lo, hi, &
                                qint, q_lo, q_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnv, qg_lo, qg_hi)

  end subroutine cmpflx_plus_godunov


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
    use eos_type_module, only: eos_t, eos_input_re
    use network, only: nspec, naux
    use castro_error_module
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : hybrid_riemann, ppm_temp_fix, riemann_solver, &
                                   T_guess

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

    if (ppm_temp_fix == 2) then
       ! recompute the thermodynamics on the interface to make it
       ! all consistent

       ! we want to take the edge states of rho, e, and X, and get
       ! new values for p on the edges that are
       ! thermodynamically consistent.

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! this is an initial guess for iterations, since we
                ! can't be certain what temp is on interfaces
                eos_state % T = T_guess

                ! minus state
                eos_state % rho = qm(i,j,k,QRHO)
                eos_state % p   = qm(i,j,k,QPRES)
                eos_state % e   = qm(i,j,k,QREINT)/qm(i,j,k,QRHO)
                eos_state % xn  = qm(i,j,k,QFS:QFS+nspec-1)
                eos_state % aux = qm(i,j,k,QFX:QFX+naux-1)

                call eos(eos_input_re, eos_state)

                qm(i,j,k,QREINT) = eos_state % e * eos_state % rho
                qm(i,j,k,QPRES)  = eos_state % p
                !gamcm(i,j)        = eos_state % gam1

             end do
          end do
       end do

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! this is an initial guess for iterations, since we
                ! can't be certain what temp is on interfaces
                eos_state % T = T_guess

                ! plus state
                eos_state % rho = qp(i,j,k,QRHO)
                eos_state % p   = qp(i,j,k,QPRES)
                eos_state % e   = qp(i,j,k,QREINT)/qp(i,j,k,QRHO)
                eos_state % xn  = qp(i,j,k,QFS:QFS+nspec-1)
                eos_state % aux = qp(i,j,k,QFX:QFX+naux-1)

                call eos(eos_input_re, eos_state)

                qp(i,j,k,QREINT) = eos_state % e * eos_state % rho
                qp(i,j,k,QPRES)  = eos_state % p
                !gamcp(i,j)        = eos_state % gam1

             end do
          end do
       end do

    endif

    ! Solve Riemann problem
    if (riemann_solver == 0) then
       ! Colella, Glaz, & Ferguson solver

       call riemannus(qm, qm_lo, qm_hi, &
                      qp, qp_lo, qp_hi, &
                      qaux, qa_lo, qa_hi, &
                      qint, q_lo, q_hi, &
#ifdef RADIATION
                      lambda_int, q_lo, q_hi, &
#endif
                      idir, compute_gammas, lo, hi, &
                      domlo, domhi)

    elseif (riemann_solver == 1) then
       ! Colella & Glaz solver

#ifndef RADIATION
       call riemanncg(qm, qm_lo, qm_hi, &
                      qp, qp_lo, qp_hi, &
                      qaux, qa_lo, qa_hi, &
                      qint, q_lo, q_hi, &
                      idir, lo, hi, &
                      domlo, domhi)
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

  end subroutine riemann_state


end module riemann_module
