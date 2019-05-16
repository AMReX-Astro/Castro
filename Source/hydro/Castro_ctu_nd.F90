module ctu_module
  !
  ! advection routines in support of the CTU unsplit advection scheme

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

  subroutine ctu_ppm_states(lo, hi, &
                            vlo, vhi, &
                            q_core, qc_lo, qc_hi, &
                            q_pass, qp_lo, qp_hi, &
#ifdef RADIATION
                            q_rad, qr_lo, qr_hi, &
#endif
                            flatn, f_lo, f_hi, &
                            qaux, qa_lo, qa_hi, &
                            q_core_src, qcs_lo, qcs_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                            q_pass_src, qps_lo, qps_hi, &
#endif
                            shk, sk_lo, sk_hi, &
                            Ip_core, Icp_lo, Icp_hi, &
                            Im_core, Icm_lo, Icm_hi, &
                            Ip_pass, Ipp_lo, Ipp_hi, &
                            Im_pass, Ipm_lo, Ipm_hi, &
#ifdef RADIATION
                            Ip_rad, Irp_lo, Irp_hi, &
                            Im_rad, Irm_lo, Irm_hi, &
#endif
                            Ip_core_src, Icsp_lo, Icsp_hi, &
                            Im_core_src, Icsm_lo, Icsm_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                            Ip_pass_src, Ipsp_lo, Ipsp_hi, &
                            Im_pass_src, Ipsm_lo, Ipsm_hi, &
#endif
                            Ip_gc, Ipg_lo, Ipg_hi, &
                            Im_gc, Img_lo, Img_hi, &
                            sm, sm_lo, sm_hi, &
                            sp, sp_lo, sp_hi, &
                            qxm_core, qxcm_lo, qxcm_hi, &
                            qxp_core, qxcp_lo, qxcp_hi, &
                            qxm_pass, qxpm_lo, qxpm_hi, &
                            qxp_pass, qxpp_lo, qxpp_hi, &
#ifdef RADIATION
                            qxm_rad, qxrm_lo, qxrm_hi, &
                            qxp_rad, qxrp_lo, qxrp_hi, &
#endif
#if AMREX_SPACEDIM >= 2
                            qym_core, qycm_lo, qycm_hi, &
                            qyp_core, qycp_lo, qycp_hi, &
                            qym_pass, qypm_lo, qypm_hi, &
                            qyp_pass, qypp_lo, qypp_hi, &
#ifdef RADIATION
                            qym_rad, qyrm_lo, qyrm_hi, &
                            qyp_rad, qyrp_lo, qyrp_hi, &
#endif
#endif
#if AMREX_SPACEDIM == 3
                            qzm_core, qzcm_lo, qzcm_hi, &
                            qzp_core, qzcp_lo, qzcp_hi, &
                            qzm_pass, qzpm_lo, qzpm_hi, &
                            qzp_pass, qzpp_lo, qzpp_hi, &
#ifdef RADIATION
                            qzm_rad, qzrm_lo, qzrm_hi, &
                            qzp_rad, qzrp_lo, qzrp_hi, &
#endif
#endif
                            dx, dt, &
#if AMREX_SPACEDIM < 3
                            dloga, dloga_lo, dloga_hi, &
#endif
                            domlo, domhi) bind(C, name="ctu_ppm_states")
    ! Compute the normal interface states by reconstructing
    ! the primitive variables using the piecewise parabolic method
    ! and doing characteristic tracing.  We do not apply the
    ! transverse terms here.

    use meth_params_module, only : NQC_SRC, NQC, NQP, NVAR, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                                   NQP_SRC, &
#endif
#ifdef RADIATION
                                   NQR, &
#endif
                                   QFS, QFX, QTEMP, QREINT, &
                                   QC, QGAMC, NQAUX, QGAME, QREINT, &
                                   NGDNV, GDU, GDV, GDW, GDPRES, &
                                   ppm_predict_gammae, ppm_temp_fix, &
                                   hybrid_riemann
    use ppm_module, only : ca_ppm_reconstruct, ppm_int_profile, ppm_reconstruct_with_eos
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use trace_ppm_rad_module, only : trace_ppm_rad
#else
    use trace_ppm_module, only : trace_ppm
#endif
    use advection_util_module, only : ca_shock
    use prob_params_module, only : dg

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: qc_lo(3), qc_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: qcs_lo(3), qcs_hi(3)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    integer, intent(in) :: qps_lo(3), qps_hi(3)
#endif
    integer, intent(in) :: sk_lo(3), sk_hi(3)
    integer, intent(in) :: Icp_lo(3), Icp_hi(3)
    integer, intent(in) :: Icm_lo(3), Icm_hi(3)
    integer, intent(in) :: Ipp_lo(3), Ipp_hi(3)
    integer, intent(in) :: Ipm_lo(3), Ipm_hi(3)
#ifdef RADIATION
    integer, intent(in) :: Irp_lo(3), Irp_hi(3)
    integer, intent(in) :: Irm_lo(3), Irm_hi(3)
#endif
    integer, intent(in) :: Icsp_lo(3), Icsp_hi(3)
    integer, intent(in) :: Icsm_lo(3), Icsm_hi(3)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    integer, intent(in) :: Ipsp_lo(3), Ipsp_hi(3)
    integer, intent(in) :: Ipsm_lo(3), Ipsm_hi(3)
#endif
    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)
    integer, intent(in) :: sm_lo(3), sm_hi(3)
    integer, intent(in) :: sp_lo(3), sp_hi(3)

    integer, intent(in) :: qxcm_lo(3), qxcm_hi(3)
    integer, intent(in) :: qxcp_lo(3), qxcp_hi(3)
    integer, intent(in) :: qxpm_lo(3), qxpm_hi(3)
    integer, intent(in) :: qxpp_lo(3), qxpp_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qxrm_lo(3), qxrm_hi(3)
    integer, intent(in) :: qxrp_lo(3), qxrp_hi(3)
#endif
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: qycm_lo(3), qycm_hi(3)
    integer, intent(in) :: qycp_lo(3), qycp_hi(3)
    integer, intent(in) :: qypm_lo(3), qypm_hi(3)
    integer, intent(in) :: qypp_lo(3), qypp_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qyrm_lo(3), qyrm_hi(3)
    integer, intent(in) :: qyrp_lo(3), qyrp_hi(3)
#endif
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: qzcm_lo(3), qzcm_hi(3)
    integer, intent(in) :: qzcp_lo(3), qzcp_hi(3)
    integer, intent(in) :: qzpm_lo(3), qzpm_hi(3)
    integer, intent(in) :: qzpp_lo(3), qzpp_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qzrm_lo(3), qzrm_hi(3)
    integer, intent(in) :: qzrp_lo(3), qzrp_hi(3)
#endif
#endif
#if AMREX_SPACEDIM < 3
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    real(rt), intent(in) :: dx(3)   ! grid spacing in X, Y, Z direction
    real(rt), intent(in), value :: dt    ! time stepsize
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q_core(qc_lo(1):qc_hi(1),qc_lo(2):qc_hi(2),qc_lo(3):qc_hi(3),NQC)   ! input state, primitives
    real(rt), intent(in) :: q_pass(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQP)   ! input state, primitives
#ifdef RADIATION
    real(rt), intent(in) :: q_rad(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),NQR)   ! input state, primitives
#endif
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)   ! auxiliary hydro data
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))   ! flattening parameter
    real(rt), intent(in) ::  q_core_src(qcs_lo(1):qcs_hi(1),qcs_lo(2):qcs_hi(2),qcs_lo(3):qcs_hi(3),NQC_SRC)   ! primitive variable source
#ifdef PRIM_SPECIES_HAVE_SOURCES
    real(rt), intent(in) ::  q_prim_src(qps_lo(1):qps_hi(1),qps_lo(2):qps_hi(2),qps_lo(3):qps_hi(3),NQP_SRC)   ! primitive variable source
#endif
    real(rt), intent(inout) :: shk(sk_lo(1):sk_hi(1), sk_lo(2):sk_hi(2), sk_lo(3):sk_hi(3))
    real(rt), intent(inout) :: Ip_core(Icp_lo(1):Icp_hi(1),Icp_lo(2):Icp_hi(2),Icp_lo(3):Icp_hi(3),1:3,NQC)
    real(rt), intent(inout) :: Im_core(Icm_lo(1):Icm_hi(1),Icm_lo(2):Icm_hi(2),Icm_lo(3):Icm_hi(3),1:3,NQC)
    real(rt), intent(inout) :: Ip_pass(Ipp_lo(1):Ipp_hi(1),Ipp_lo(2):Ipp_hi(2),Ipp_lo(3):Ipp_hi(3),1:3,NQP)
    real(rt), intent(inout) :: Im_pass(Ipm_lo(1):Ipm_hi(1),Ipm_lo(2):Ipm_hi(2),Ipm_lo(3):Ipm_hi(3),1:3,NQP)

    real(rt), intent(inout) :: Ip_core_src(Icsp_lo(1):Icsp_hi(1),Icsp_lo(2):Icsp_hi(2),Icsp_lo(3):Icsp_hi(3),1:3,NQC_SRC)
    real(rt), intent(inout) :: Im_core_src(Icsm_lo(1):Icsm_hi(1),Icsm_lo(2):Icsm_hi(2),Icsm_lo(3):Icsm_hi(3),1:3,NQC_SRC)

#ifdef PRIM_SPECIES_HAVE_SOURCES
    real(rt), intent(inout) :: Ip_pass_src(Ipsp_lo(1):Ipsp_hi(1),Ipsp_lo(2):Ipsp_hi(2),Ipsp_lo(3):Ipsp_hi(3),1:3,NQP_SRC)
    real(rt), intent(inout) :: Im_pass_src(Ipsm_lo(1):Ipsm_hi(1),Ipsm_lo(2):Ipsm_hi(2),Ipsm_lo(3):Ipsm_hi(3),1:3,NQP_SRC)
#endif
    real(rt), intent(inout) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:3,1)
    real(rt), intent(inout) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:3,1)

    real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1), sm_lo(2):sm_hi(2), sm_lo(3):sm_hi(3))
    real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1), sp_lo(2):sp_hi(2), sp_lo(3):sp_hi(3))

    real(rt), intent(inout) :: qxm_core(qxcm_lo(1):qxcm_hi(1), qxcm_lo(2):qxcm_hi(2), qxcm_lo(3):qxcm_hi(3), NQC)
    real(rt), intent(inout) :: qxm_pass(qxpm_lo(1):qxpm_hi(1), qxpm_lo(2):qxpm_hi(2), qxpm_lo(3):qxpm_hi(3), NQP)
    real(rt), intent(inout) :: qxp_core(qxcp_lo(1):qxcp_hi(1), qxcp_lo(2):qxcp_hi(2), qxcp_lo(3):qxcp_hi(3), NQC)
    real(rt), intent(inout) :: qxp_pass(qxpp_lo(1):qxpp_hi(1), qxpp_lo(2):qxpp_hi(2), qxpp_lo(3):qxpp_hi(3), NQP)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: qym_core(qycm_lo(1):qycm_hi(1), qycm_lo(2):qycm_hi(2), qycm_lo(3):qycm_hi(3), NQC)
    real(rt), intent(inout) :: qym_pass(qypm_lo(1):qypm_hi(1), qypm_lo(2):qypm_hi(2), qypm_lo(3):qypm_hi(3), NQP)
    real(rt), intent(inout) :: qyp_core(qycp_lo(1):qycp_hi(1), qycp_lo(2):qycp_hi(2), qycp_lo(3):qycp_hi(3), NQC)
    real(rt), intent(inout) :: qyp_pass(qypp_lo(1):qypp_hi(1), qypp_lo(2):qypp_hi(2), qypp_lo(3):qypp_hi(3), NQP)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: qzm_core(qzcm_lo(1):qzcm_hi(1), qzcm_lo(2):qzcm_hi(2), qzcm_lo(3):qzcm_hi(3), NQC)
    real(rt), intent(inout) :: qzm_pass(qzpm_lo(1):qzpm_hi(1), qzpm_lo(2):qzpm_hi(2), qzpm_lo(3):qzpm_hi(3), NQP)
    real(rt), intent(inout) :: qzp_core(qzcp_lo(1):qzcp_hi(1), qzcp_lo(2):qzcp_hi(2), qzcp_lo(3):qzcp_hi(3), NQC)
    real(rt), intent(inout) :: qzp_pass(qzpp_lo(1):qzpp_hi(1), qzpp_lo(2):qzpp_hi(2), qzpp_lo(3):qzpp_hi(3), NQP)
#endif
#if AMREX_SPACEDIM < 3
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt) :: hdt
    integer :: i, j, k, n, idir

    logical :: source_nonzero(NQSRC)
    logical :: reconstruct_state(NQ)

    logical :: compute_shock

    !$gpu

    hdt = HALF*dt

    ! multidimensional shock detection

#ifdef SHOCK_VAR
    compute_shock = .true.
#else
    compute_shock = .false.
#endif

    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1 .or. compute_shock) then
       call ca_shock(lo, hi, &
                     q, qd_lo, qd_hi, &
                     shk, sk_lo, sk_hi, &
                     dx)
    else
       shk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ZERO
    endif

    ! we don't need to reconstruct all of the NQ state variables,
    ! depending on how we are tracing
    reconstruct_state(:) = .true.
    if (ppm_predict_gammae /= 1) then
       reconstruct_state(QGAME) = .false.
    else
       reconstruct_state(QREINT) = .false.
    endif
    if (ppm_temp_fix == 0 .or. ppm_temp_fix == 2) then
       reconstruct_state(QTEMP) = .false.
    endif

    ! preprocess the sources -- we don't want to trace under a source
    ! that is empty. This check only needs to be done over the tile
    ! we're working on, since the PPM reconstruction and integration
    ! done here is only local to this tile.

    do n = 1, NQSRC
       if (maxval(abs(srcQ(lo(1)-2:hi(1)+2,lo(2)-2*dg(2):hi(2)+2*dg(2),lo(3)-2*dg(3):hi(3)+2*dg(3),n))) == ZERO) then
          source_nonzero(n) = .false.
       else
          source_nonzero(n) = .true.
       endif
    enddo

    do idir = 1, AMREX_SPACEDIM

       ! Compute Ip and Im -- this does the parabolic reconstruction,
       ! limiting, and returns the integral of each profile under each
       ! wave to each interface
       do n = 1, NQ
          if (.not. reconstruct_state(n)) cycle

          ! core primitive variables
          call ca_ppm_reconstruct(lo, hi, 0, idir, &
                                  q_core, qc_lo, qc_hi, NQC, n, n, &
                                  flatn, f_lo, f_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  1, 1, 1)

          call ppm_int_profile(lo, hi, idir, &
                               q_core, qc_lo, qc_hi, NQC, n, &
                               q_core, qc_lo, qc_hi, &
                               qaux, qa_lo, qa_hi, &
                               sm, sm_lo, sm_hi, &
                               sp, sp_lo, sp_hi, &
                               Ip_core, Icp_lo, Icp_hi, &
                               Im_core, Icm_lo, Icm_hi, NQC, n, &
                               dx, dt)


          ! passives
          call ca_ppm_reconstruct(lo, hi, 0, idir, &
                                  q_pass, qp_lo, qp_hi, NQP, n, n, &
                                  flatn, f_lo, f_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  1, 1, 1)

          call ppm_int_profile(lo, hi, idir, &
                               q_pass, qp_lo, qp_hi, NQP, n, &
                               q_core, qc_lo, qc_hi, &
                               qaux, qa_lo, qa_hi, &
                               sm, sm_lo, sm_hi, &
                               sp, sp_lo, sp_hi, &
                               Ip_pass, Ipp_lo, Ipp_hi, &
                               Im_pass, Ipm_lo, Ipm_hi, NQP, n, &
                               dx, dt)

#ifdef RADIATION
          ! radiation
          call ca_ppm_reconstruct(lo, hi, 0, idir, &
                                  q_rad, qr_lo, qr_hi, NQR, n, n, &
                                  flatn, f_lo, f_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  1, 1, 1)

          call ppm_int_profile(lo, hi, idir, &
                               q_rad, qr_lo, qr_hi, NQR, n, &
                               q_core, qc_lo, qc_hi, &
                               qaux, qa_lo, qa_hi, &
                               sm, sm_lo, sm_hi, &
                               sp, sp_lo, sp_hi, &
                               Ip_rad, Irp_lo, Irp_hi, &
                               Im_rad, Irm_lo, Irm_hi, NQR, n, &
                               dx, dt)
       end do
#endif

       if (ppm_temp_fix /= 1) then
          call ca_ppm_reconstruct(lo, hi, 0, idir, &
                                  qaux, qa_lo, qa_hi, NQAUX, QGAMC, QGAMC, &
                                  flatn, f_lo, f_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  1, 1, 1)

          call ppm_int_profile(lo, hi, idir, &
                               qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               q_core, qc_lo, qc_hi, &
                               qaux, qa_lo, qa_hi, &
                               sm, sm_lo, sm_hi, &
                               sp, sp_lo, sp_hi, &
                               Ip_gc, Ipg_lo, Ipg_hi, &
                               Im_gc, Img_lo, Img_hi, 1, 1, &
                               dx, dt)
       else

          ! temperature-based PPM
          call ppm_reconstruct_with_eos(lo, hi, idir, &
                                        Ip_core, Ipc_lo, Ipc_hi, &
                                        Im_core, Imc_lo, Imc_hi, &
                                        Ip_pass, Ipp_lo, Ipp_hi, &
                                        Im_pass, Imp_lo, Imp_hi, &
                                        Ip_gc, Ipg_lo, Ipg_hi, &
                                        Im_gc, Img_lo, Img_hi)

       end if


       ! source terms
       do n = 1, NQSRC
          if (source_nonzero(n)) then
             call ca_ppm_reconstruct(lo, hi, 0, idir, &
                                     srcQ, src_lo, src_hi, NQSRC, n, n, &
                                     flatn, f_lo, f_hi, &
                                     sm, sm_lo, sm_hi, &
                                     sp, sp_lo, sp_hi, &
                                     1, 1, 1)

             call ppm_int_profile(lo, hi, idir, &
                                  srcQ, src_lo, src_hi, NQSRC, n, &
                                  q, qd_lo, qd_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  Ip_src, Ips_lo, Ips_hi, &
                                  Im_src, Ims_lo, Ims_hi, NQSRC, n, &
                                  dx, dt)
          else
             Ip_src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:,n) = ZERO
             Im_src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:,n) = ZERO
          endif

       enddo


       ! compute the interface states

#ifdef RADIATION
       if (idir == 1) then
          call trace_ppm_rad(lo, hi, &
                             1, &
                            q_core, qc_lo, qc_hi, &
                             q_pass, qp_lo, qp_hi, &
                             q_rad, qr_lo, qr_hi, &
                             qaux, qa_lo, qa_hi, &
                             Ip_core, Ip_lo, Ip_hi, &
                             Im_core, Im_lo, Im_hi, &
                             Ip_pass, Ipp_lo, Ipp_hi, &
                             Im_pass, Ipm_lo, Ipm_hi, &
                             Ip_rad, Irp_lo, Irp_hi, &
                             Im_rad, Irm_lo, Irm_hi, &
                             Ip_core_src, Icsp_lo, Icsp_hi, &
                             Im_core_src, Icsm_lo, Icsm_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                             Ip_pass_src, Ipsp_lo, Ipsp_hi, &
                             Im_pass_src, Ipsm_lo, Ipsm_hi, &
#endif
                             qxm_core, qxm_lo, qxm_hi, &
                             qxp_core, qxp_lo, qxp_hi, &
                             qxm_pass, qxpm_lo, qxpm_hi, &
                             qxp_pass, qxpp_lo, qxpp_hi, &
#if AMREX_SPACEDIM <= 2
                             dloga, dloga_lo, dloga_hi, &
#endif
                             vlo, vhi, domlo, domhi, &
                             dx, dt)

#if AMREX_SPACEDIM >= 2
       else if (idir == 2) then
          call trace_ppm_rad(lo, hi, &
                             2, &
                             q_core, qc_lo, qc_hi, &
                             q_pass, qp_lo, qp_hi, &
                             q_rad, qr_lo, qr_hi, &
                             qaux, qa_lo, qa_hi, &
                             Ip_core, Ip_lo, Ip_hi, &
                             Im_core, Im_lo, Im_hi, &
                             Ip_pass, Ipp_lo, Ipp_hi, &
                             Im_pass, Ipm_lo, Ipm_hi, &
                             Ip_rad, Irp_lo, Irp_hi, &
                             Im_rad, Irm_lo, Irm_hi, &
                             Ip_core_src, Icsp_lo, Icsp_hi, &
                             Im_core_src, Icsm_lo, Icsm_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                             Ip_pass_src, Ipsp_lo, Ipsp_hi, &
                             Im_pass_src, Ipsm_lo, Ipsm_hi, &
#endif
                             qym_core, qym_lo, qym_hi, &
                             qyp_core, qyp_lo, qyp_hi, &
                             qym_pass, qypm_lo, qypm_hi, &
                             qyp_pass, qypp_lo, qypp_hi, &
#if AMREX_SPACEDIM == 2
                             dloga, dloga_lo, dloga_hi, &
#endif
                             vlo, vhi, domlo, domhi, &
                             dx, dt)
#endif

#if AMREX_SPACEDIM == 3
       else
          call trace_ppm_rad(lo, hi, &
                             3, &
                             q_core, qc_lo, qc_hi, &
                             q_pass, qp_lo, qp_hi, &
                             q_rad, qr_lo, qr_hi, &
                             qaux, qa_lo, qa_hi, &
                             Ip_core, Ip_lo, Ip_hi, &
                             Im_core, Im_lo, Im_hi, &
                             Ip_pass, Ipp_lo, Ipp_hi, &
                             Im_pass, Ipm_lo, Ipm_hi, &
                             Ip_rad, Irp_lo, Irp_hi, &
                             Im_rad, Irm_lo, Irm_hi, &
                             Ip_core_src, Icsp_lo, Icsp_hi, &
                             Im_core_src, Icsm_lo, Icsm_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                             Ip_pass_src, Ipsp_lo, Ipsp_hi, &
                             Im_pass_src, Ipsm_lo, Ipsm_hi, &
#endif
                             qzm_core, qzm_lo, qzm_hi, &
                             qzp_core, qzp_lo, qzp_hi, &
                             qzm_pass, qzpm_lo, qzpm_hi, &
                             qzp_pass, qzpp_lo, qzpp_hi, &
                             vlo, vhi, domlo, domhi, &
                             dx, dt)
#endif
       endif
#else
       ! hydro (no radiation)
       if (idir == 1) then
          call trace_ppm(lo, hi, &
                         1, &
                         q_core, qc_lo, qc_hi, &
                         q_pass, qp_lo, qp_hi, &
                         qaux, qa_lo, qa_hi, &
                         Ip_core, Icp_lo, Icp_hi, &
                         Im_core, Icm_lo, Icm_hi, &
                         Ip_pass, Ipp_lo, Ipp_hi, &
                         Im_pass, Ipm_lo, Ipm_hi, &
                         Ip_core_src, Icsp_lo, Icsp_hi, &
                         Im_core_src, Icsm_lo, Icsm_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                         Ip_pass_src, Ipsp_lo, Ipsp_hi, &
                         Im_pass_src, Ipsm_lo, Ipsm_hi, &
#endif
                         Ip_gc, Ipg_lo, Ipg_hi, &
                         Im_gc, Img_lo, Img_hi, &
                         qxm_core, qxcm_lo, qxcm_hi, &
                         qxp_core, qxcp_lo, qxcp_hi, &
                         qxm_pass, qxpm_lo, qxpm_hi, &
                         qxp_pass, qxpp_lo, qxpp_hi, &
#if AMREX_SPACEDIM <= 2
                         dloga, dloga_lo, dloga_hi, &
#endif
                         vlo, vhi, domlo, domhi, &
                         dx, dt)

#if AMREX_SPACEDIM >= 2
       else if (idir == 2) then
          call trace_ppm(lo, hi, &
                         2, &
                         q_core, qc_lo, qc_hi, &
                         q_pass, qp_lo, qp_hi, &
                         qaux, qa_lo, qa_hi, &
                         Ip_core, Icp_lo, Icp_hi, &
                         Im_core, Icm_lo, Icm_hi, &
                         Ip_pass, Ipp_lo, Ipp_hi, &
                         Im_pass, Ipm_lo, Ipm_hi, &
                         Ip_core_src, Icsp_lo, Icsp_hi, &
                         Im_core_src, Icsm_lo, Icsm_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                         Ip_pass_src, Ipsp_lo, Ipsp_hi, &
                         Im_pass_src, Ipsm_lo, Ipsm_hi, &
#endif
                         Ip_gc, Ipg_lo, Ipg_hi, &
                         Im_gc, Img_lo, Img_hi, &

                         qym_core, qycm_lo, qycm_hi, &
                         qyp_core, qycp_lo, qycp_hi, &
                         qym_pass, qypm_lo, qypm_hi, &
                         qyp_pass, qypp_lo, qypp_hi, &
#if AMREX_SPACEDIM == 2
                         dloga, dloga_lo, dloga_hi, &
#endif
                         vlo, vhi, domlo, domhi, &
                         dx, dt)
#endif

#if AMREX_SPACEDIM == 3
       else
          call trace_ppm(lo, hi, &
                         3, &
                         q_core, qc_lo, qc_hi, &
                         q_pass, qp_lo, qp_hi, &
                         qaux, qa_lo, qa_hi, &
                         Ip_core, Icp_lo, Icp_hi, &
                         Im_core, Icm_lo, Icm_hi, &
                         Ip_pass, Ipp_lo, Ipp_hi, &
                         Im_pass, Ipm_lo, Ipm_hi, &
                         Ip_core_src, Icsp_lo, Icsp_hi, &
                         Im_core_src, Icsm_lo, Icsm_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                         Ip_pass_src, Ipsp_lo, Ipsp_hi, &
                         Im_pass_src, Ipsm_lo, Ipsm_hi, &
#endif
                         Ip_gc, Ipg_lo, Ipg_hi, &
                         Im_gc, Img_lo, Img_hi, &

                         qzm_core, qzcm_lo, qzcm_hi, &
                         qzp_core, qzcp_lo, qzcp_hi, &
                         qzm_pass, qzpm_lo, qzpm_hi, &
                         qzp_pass, qzpp_lo, qzpp_hi, &
                         vlo, vhi, domlo, domhi, &
                         dx, dt)
#endif
       end if
#endif

    end do

  end subroutine ctu_ppm_states


  subroutine ctu_plm_states(lo, hi, &
                            vlo, vhi, &
                            q_core, qc_lo, qc_hi, &
                            q_pass, qp_lo, qp_hi, &
                            flatn, f_lo, f_hi, &
                            qaux, qa_lo, qa_hi, &
                            q_core_src, qcs_lo, qcs_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                            q_pass_src, qps_lo, qps_hi, &
#endif
                            shk, sk_lo, sk_hi, &
                            dq_core, dqc_lo, dqc_hi, &
                            dq_pass, dqp_lo, dqp_hi, &
                            qxm_core, qxmc_lo, qxmc_hi, &
                            qxm_pass, qxmp_lo, qxmp_hi, &
                            qxp_core, qxpc_lo, qxpc_hi, &
                            qxp_pass, qxpp_lo, qxpp_hi, &
#if AMREX_SPACEDIM >= 2
                            qym_core, qymc_lo, qymc_hi, &
                            qym_pass, qymp_lo, qymp_hi, &
                            qyp_core, qypc_lo, qypc_hi, &
                            qyp_pass, qypp_lo, qypp_hi, &
#endif
#if AMREX_SPACEDIM == 3
                            qzm_core, qzmc_lo, qzmc_hi, &
                            qzm_pass, qzmp_lo, qzmp_hi, &
                            qzp_core, qzpc_lo, qzpc_hi, &
                            qzp_pass, qzpp_lo, qzpp_hi, &
#endif
                            dx, dt, &
#if AMREX_SPACEDIM < 3
                            dloga, dloga_lo, dloga_hi, &
#endif
                            domlo, domhi) bind(C, name="ctu_plm_states")
    ! Compute the normal interface states by reconstructing
    ! the primitive variables using piecewise linear slopes and doing
    ! characteristic tracing.  We do not apply the transverse terms here.
    !
    ! .. todo::
    !    we can get rid of the the different temporary q Godunov
    !    state arrays
    !

    use meth_params_module, only : NQSRC, NQ, NVAR, &
                                   QFS, QFX, QTEMP, QREINT, &
                                   QC, QGAMC, NQAUX, QGAME, QREINT, &
                                   NGDNV, GDU, GDV, GDW, GDPRES, &
                                   plm_iorder, use_pslope, hybrid_riemann
    use trace_plm_module, only : trace_plm
    use slope_module, only : uslope, pslope
    use advection_util_module, only : ca_shock
    use prob_params_module, only : dg

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: qc_lo(3), qc_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: sk_lo(3), sk_hi(3)
    integer, intent(in) :: dqc_lo(3), dqc_hi(3)
    integer, intent(in) :: dqp_lo(3), dqp_hi(3)
    integer, intent(in) :: qxmc_lo(3), qxmc_hi(3)
    integer, intent(in) :: qxmp_lo(3), qxmp_hi(3)
    integer, intent(in) :: qxpc_lo(3), qxpc_hi(3)
    integer, intent(in) :: qxpp_lo(3), qxpp_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: qymc_lo(3), qymc_hi(3)
    integer, intent(in) :: qymp_lo(3), qymp_hi(3)
    integer, intent(in) :: qypc_lo(3), qypc_hi(3)
    integer, intent(in) :: qypp_lo(3), qypp_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: qzmc_lo(3), qzmc_hi(3)
    integer, intent(in) :: qzmp_lo(3), qzmp_hi(3)
    integer, intent(in) :: qzpc_lo(3), qzpc_hi(3)
    integer, intent(in) :: qzpp_lo(3), qzpp_hi(3)
#endif
#if AMREX_SPACEDIM < 3
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    real(rt), intent(in) :: dx(3)   ! grid spacing in X, Y, Z direction
    real(rt), intent(in), value :: dt   ! time stepsize
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q_core(qc_lo(1):qc_hi(1),qc_lo(2):qc_hi(2),qc_lo(3):qc_hi(3),NQC)   ! input state, primitives
    real(rt), intent(in) :: q_pass(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQP)   ! input state, primitives
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)   ! auxiliary hydro data
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))   ! flattening parameter
    real(rt), intent(in) :: q_core_src(qcs_lo(1):qcs_hi(1),qcs_lo(2):qcs_hi(2),qcs_lo(3):qcs_hi(3),NQC_SRC)   ! primitive variable source
#ifdef PRIM_SPECIES_HAVE_SOURCES
    real(rt), intent(in) :: q_pass_src(qps_lo(1):qps_hi(1),qps_lo(2):qps_hi(2),qps_lo(3):qps_hi(3),NQP_SRC)   ! primitive variable source
#endif
    real(rt), intent(inout) :: shk(sk_lo(1):sk_hi(1), sk_lo(2):sk_hi(2), sk_lo(3):sk_hi(3))
    real(rt), intent(inout) :: dqc(dqc_lo(1):dqc_hi(1), dqc_lo(2):dqc_hi(2), dqc_lo(3):dqc_hi(3), NQC)
    real(rt), intent(inout) :: dqp(dqp_lo(1):dqp_hi(1), dqp_lo(2):dqp_hi(2), dqp_lo(3):dqp_hi(3), NQP)

    real(rt), intent(inout) :: qxm_core(qxmc_lo(1):qxmc_hi(1), qxmc_lo(2):qxmc_hi(2), qxmc_lo(3):qxmc_hi(3), NQC)
    real(rt), intent(inout) :: qxm_pass(qxmp_lo(1):qxmp_hi(1), qxmp_lo(2):qxmp_hi(2), qxmp_lo(3):qxmp_hi(3), NQP)
    real(rt), intent(inout) :: qxp_core(qxpc_lo(1):qxpc_hi(1), qxpc_lo(2):qxpc_hi(2), qxpc_lo(3):qxpc_hi(3), NQC)
    real(rt), intent(inout) :: qxp_pass(qxpp_lo(1):qxpp_hi(1), qxpp_lo(2):qxpp_hi(2), qxpp_lo(3):qxpp_hi(3), NQP)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: qym_core(qymc_lo(1):qymc_hi(1), qymc_lo(2):qymc_hi(2), qymc_lo(3):qymc_hi(3), NQC)
    real(rt), intent(inout) :: qym_pass(qymp_lo(1):qymp_hi(1), qymp_lo(2):qymp_hi(2), qymp_lo(3):qymp_hi(3), NQP)
    real(rt), intent(inout) :: qyp_core(qypc_lo(1):qypc_hi(1), qypc_lo(2):qypc_hi(2), qypc_lo(3):qypc_hi(3), NQC)
    real(rt), intent(inout) :: qyp_pass(qypp_lo(1):qypp_hi(1), qypp_lo(2):qypp_hi(2), qypp_lo(3):qypp_hi(3), NQP)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: qzm_core(qzmc_lo(1):qzmc_hi(1), qzmc_lo(2):qzmc_hi(2), qzmc_lo(3):qzmc_hi(3), NQC)
    real(rt), intent(inout) :: qzm_pass(qzmp_lo(1):qzmp_hi(1), qzmp_lo(2):qzmp_hi(2), qzmp_lo(3):qzmp_hi(3), NQP)
    real(rt), intent(inout) :: qzp_core(qzpc_lo(1):qzpc_hi(1), qzpc_lo(2):qzpc_hi(2), qzpc_lo(3):qzpc_hi(3), NQC)
    real(rt), intent(inout) :: qzp_pass(qzpp_lo(1):qzpp_hi(1), qzpp_lo(2):qzpp_hi(2), qzpp_lo(3):qzpp_hi(3), NQP)
#endif
#if AMREX_SPACEDIM < 3
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt) :: hdt
    integer :: i, j, k, n, idir

    logical :: reconstruct_state(NQC)

    logical :: compute_shock

    !$gpu

    hdt = HALF*dt

    ! multidimensional shock detection

#ifdef SHOCK_VAR
    compute_shock = .true.
#else
    compute_shock = .false.
#endif

    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1 .or. compute_shock) then
       call ca_shock(lo, hi, &
                     q_core, qc_lo, qc_hi, &
                     shk, sk_lo, sk_hi, &
                     dx)
    else
       shk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ZERO
    endif

    ! we don't need to reconstruct all of the NQ state variables,
    ! depending on how we are tracing
    reconstruct_state(:) = .true.
    reconstruct_state(QGAME) = .false.
    reconstruct_state(QTEMP) = .false.


#ifdef RADIATION
#ifndef AMREX_USE_CUDA
    call amrex_error("ppm_type <=0 is not supported in with radiation")
#endif
#endif

    ! Compute all slopes
    do idir = 1, AMREX_SPACEDIM

       do n = 1, NQC
          if (.not. reconstruct_state(n)) cycle
          call uslope(lo, hi, idir, &
                      q_core, qc_lo, qc_hi, n, NQC, &
                      flatn, f_lo, f_hi, &
                      dq_core, dqc_lo, dqc_hi)
       end do

       do n = 1, NQP
          call uslope(lo, hi, idir, &
                      q_pass, qp_lo, qp_hi, n, NQP, &
                      flatn, f_lo, f_hi, &
                      dp_pass, dqp_lo, dqp_hi)
       end do

       if (use_pslope == 1) then
          call pslope(lo, hi, idir, &
                      q_core, qc_lo, qc_hi, &
                      flatn, f_lo, f_hi, &
                      dq_core, dqc_lo, dqc_hi, &
                      q_core_src, qcs_lo, qcs_hi, &
                      dx)
       endif


       ! compute the interface states

       if (idir == 1) then
          call trace_plm(lo, hi, &
                         1, &
                         q_core, qc_lo, qc_hi, &
                         qaux, qa_lo, qa_hi, &
                         dq_core, dqc_lo, dqc_hi, &
                         qxm_core, qxmc_lo, qxmc_hi, &
                         qxp_core, qxpc_lo, qxpc_hi, &
#if AMREX_SPACEDIM < 3
                         dloga, dloga_lo, dloga_hi, &
#endif
                         q_core_src, qcs_lo, qcs_hi, &
                         vlo, vhi, domlo, domhi, &
                         dx, dt)

          call trace_plm_passive(lo, hi, &
                                 1, &
                                 q_pass, qp_lo, qp_hi, &
                                 q_core, qc_lo, qc_hi, &
                                 dq_pass, dqp_lo, dqp_hi, &
                                 qxm_pass, qxmp_lo, qxmp_hi, &
                                 qxp_pass, qxpp_lo, qxpp_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                                 q_pass_src, qps_lo, qps_hi, &
#endif
                                 vlo, vhi, domlo, domhi, &
                                 dx, dt)

#if AMREX_SPACEDIM >= 2
       else if (idir == 2) then
          call trace_plm(lo, hi, &
                         2, &
                         q_core, qc_lo, qc_hi, &
                         qaux, qa_lo, qa_hi, &
                         dq, dq_lo, dq_hi, &
                         qym, qym_lo, qym_hi, &
                         qyp, qyp_lo, qyp_hi, &
#if AMREX_SPACEDIM < 3
                         dloga, dloga_lo, dloga_hi, &
#endif
                         SrcQ, src_lo, src_hi, &
                         vlo, vhi, domlo, domhi, &
                         dx, dt)

          call trace_plm_passive(lo, hi, &
                                 2, &
                                 q_pass, qp_lo, qp_hi, &
                                 q_core, qc_lo, qc_hi, &
                                 dq_pass, dqp_lo, dqp_hi, &
                                 qym_pass, qymp_lo, qymp_hi, &
                                 qyp_pass, qypp_lo, qypp_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                                 q_pass_src, qps_lo, qps_hi, &
#endif
                                 vlo, vhi, domlo, domhi, &
                                 dx, dt)

#endif

#if AMREX_SPACEDIM == 3
       else
          call trace_plm(lo, hi, &
                         3, &
                         q_core, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         dq, dq_lo, dq_hi, &
                         qzm, qzm_lo, qzm_hi, &
                         qzp, qzp_lo, qzp_hi, &
                         SrcQ, src_lo, src_hi, &
                         vlo, vhi, domlo, domhi, &
                         dx, dt)

          call trace_plm_passive(lo, hi, &
                                 3, &
                                 q_pass, qp_lo, qp_hi, &
                                 q_core, qc_lo, qc_hi, &
                                 dq_pass, dqp_lo, dqp_hi, &
                                 qzm_pass, qzmp_lo, qzmp_hi, &
                                 qzp_pass, qzpp_lo, qzpp_hi, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                                 q_pass_src, qps_lo, qps_hi, &
#endif
                                 vlo, vhi, domlo, domhi, &
                                 dx, dt)

#endif
       end if

    end do

  end subroutine ctu_plm_states


  subroutine ctu_consup(lo, hi, &
                        uin, uin_lo, uin_hi, &
                        shk,  sk_lo, sk_hi, &
                        uout, uout_lo, uout_hi, &
                        update, updt_lo, updt_hi, &
                        flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
                        flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                        flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
                        Erin, Erin_lo, Erin_hi, &
                        Erout, Erout_lo, Erout_hi, &
                        radflux1, radflux1_lo, radflux1_hi, &
#if AMREX_SPACEDIM >= 2
                        radflux2, radflux2_lo, radflux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                        radflux3, radflux3_lo, radflux3_hi, &
#endif
                        nstep_fsp, &
#endif
                        qx, qx_lo, qx_hi, &
#if AMREX_SPACEDIM >= 2
                        qy, qy_lo, qy_hi, &
#endif
#if AMREX_SPACEDIM == 3
                        qz, qz_lo, qz_hi, &
#endif
                        area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                        area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                        area3, area3_lo, area3_hi, &
#endif
                        vol, vol_lo, vol_hi, &
                        pdivu, pdivu_lo, pdivu_hi, &
                        dx, dt) bind(C, name="ctu_consup")

    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, NGDNV, NQ, &
#ifdef RADIATION
                                   fspace_type, comoving, &
                                   GDU, GDV, GDW, GDLAMS, GDERADS, &
#endif
                                   GDPRES
    use advection_util_module, only : calc_pdivu
    use prob_params_module, only : mom_flux_has_p, center, dg
#ifdef RADIATION
    use rad_params_module, only : ngroups, nugroup, dlognu
    use radhydro_nd_module, only : advect_in_fspace
    use fluxlimiter_module, only : Edd_factor
#endif
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : add_hybrid_advection_source
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
    use amrex_constants_module, only : ZERO, ONE, TWO, FOURTH, HALF

    integer, intent(in) ::       lo(3),       hi(3)
    integer, intent(in) ::   uin_lo(3),   uin_hi(3)
    integer, intent(in) ::     q_lo(3),     q_hi(3)
    integer, intent(in) :: sk_lo(3), sk_hi(3)
    integer, intent(in) ::  uout_lo(3),  uout_hi(3)
    integer, intent(in) ::  updt_lo(3),  updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) ::    qy_lo(3),    qy_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)
    integer, intent(in) ::    qz_lo(3),    qz_hi(3)
#endif
    integer, intent(in) ::    qx_lo(3),    qx_hi(3)
    integer, intent(in) ::   vol_lo(3),   vol_hi(3)
    integer, intent(in) ::   pdivu_lo(3),   pdivu_hi(3)
#ifdef RADIATION
    integer, intent(in) :: Erout_lo(3), Erout_hi(3)
    integer, intent(in) :: Erin_lo(3), Erin_hi(3)
    integer, intent(in) :: radflux1_lo(3), radflux1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: radflux2_lo(3), radflux2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: radflux3_lo(3), radflux3_hi(3)
#endif
    integer, intent(inout) :: nstep_fsp
#endif

    real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(in) :: shk(sk_lo(1):sk_hi(1),sk_lo(2):sk_hi(2),sk_lo(3):sk_hi(3))
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),updt_lo(3):updt_hi(3),NVAR)

    real(rt), intent(in) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt), intent(in) ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)

#if AMREX_SPACEDIM >= 2
    real(rt), intent(in) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    real(rt), intent(in) ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
#endif

#if AMREX_SPACEDIM == 3
    real(rt), intent(in) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    real(rt), intent(in) ::    qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
#endif

    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: pdivu(pdivu_lo(1):pdivu_hi(1),pdivu_lo(2):pdivu_hi(2),pdivu_lo(3):pdivu_hi(3))
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

#ifdef RADIATION
    real(rt), intent(in) :: Erin(Erin_lo(1):Erin_hi(1),Erin_lo(2):Erin_hi(2),Erin_lo(3):Erin_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: Erout(Erout_lo(1):Erout_hi(1),Erout_lo(2):Erout_hi(2),Erout_lo(3):Erout_hi(3),0:ngroups-1)
    real(rt), intent(in) :: radflux1(radflux1_lo(1):radflux1_hi(1),radflux1_lo(2):radflux1_hi(2),radflux1_lo(3):radflux1_hi(3),0:ngroups-1)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(in) :: radflux2(radflux2_lo(1):radflux2_hi(1),radflux2_lo(2):radflux2_hi(2),radflux2_lo(3):radflux2_hi(3),0:ngroups-1)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(in) :: radflux3(radflux3_lo(1):radflux3_hi(1),radflux3_lo(2):radflux3_hi(2),radflux3_lo(3):radflux3_hi(3),0:ngroups-1)
#endif

#endif

    integer :: i, j, g, k, n
    integer :: domlo(3), domhi(3)
    real(rt) :: volInv

#ifdef RADIATION
    real(rt), dimension(0:ngroups-1) :: Erscale
    real(rt), dimension(0:ngroups-1) :: ustar, af
    real(rt) :: Eddf, Eddfxm, Eddfxp, Eddfym, Eddfyp, Eddfzm, Eddfzp
    real(rt) :: f1, f2, f1xm, f1xp, f1ym, f1yp, f1zm, f1zp
    real(rt) :: Gf1E(3)
    real(rt) :: ux, uy, uz, divu, lamc, Egdc
    real(rt) :: dudx(3), dudy(3), dudz(3), nhat(3), GnDotu(3), nnColonDotGu
    real(rt) :: dprdx, dprdy, dprdz, ek1, ek2, dek, dpdx
    real(rt) :: urho_new
    real(rt) :: umx_new1, umy_new1, umz_new1
    real(rt) :: umx_new2, umy_new2, umz_new2
#endif

    !$gpu

#ifdef RADIATION
    if (ngroups .gt. 1) then
       if (fspace_type .eq. 1) then
          Erscale = dlognu
       else
          Erscale = nugroup*dlognu
       end if
    end if
#endif

    call calc_pdivu(lo, hi, &
                    qx, qx_lo, qx_hi, &
                    area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                    qy, qy_lo, qy_hi, &
                    area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                    qz, qz_lo, qz_hi, &
                    area3, area3_lo, area3_hi, &
#endif
                    vol, vol_lo, vol_hi, &
                    dx, pdivu, pdivu_lo, pdivu_hi)


    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                volinv = ONE / vol(i,j,k)

                update(i,j,k,n) = update(i,j,k,n) + &
                     ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) &
#if AMREX_SPACEDIM >= 2
                     + flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) &
#endif
#if AMREX_SPACEDIM == 3
                     + flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) &
#endif
                     ) * volinv

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then
                   update(i,j,k,n) = update(i,j,k,n) - pdivu(i,j,k)
                endif

             enddo
          enddo
       enddo
    enddo

#ifdef SHOCK_VAR
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uout(i,j,k,USHK) = shk(i,j,k)
          end do
       end do
    end do
#endif

#ifdef HYBRID_MOMENTUM
    call add_hybrid_advection_source(lo, hi, dt, &
                                     update, uout_lo, uout_hi, &
                                     qx, qx_lo, qx_hi, &
                                     qy, qy_lo, qy_hi, &
                                     qz, qz_lo, qz_hi)
#endif


#ifndef RADIATION
    ! Add gradp term to momentum equation -- only for axisymmetric
    ! coords (and only for the radial flux).

    if (.not. mom_flux_has_p(1)%comp(UMX)) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                update(i,j,k,UMX) = update(i,j,k,UMX) - (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx(1)
                !update(i,j,UMY) = update(i,j,UMY) - (pgdy(i,j+1)-pgdy(i,j)) / dy
             enddo
          enddo
       enddo
    endif

#else
    ! radiation energy update.  For the moment, we actually update things
    ! fully here, instead of creating a source term for the update
    do g = 0, ngroups-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                Erout(i,j,k,g) = Erin(i,j,k,g) + dt * &
                     ( radflux1(i,j,k,g) * area1(i,j,k) - radflux1(i+1,j,k,g) * area1(i+1,j,k) &
#if AMREX_SPACEDIM >= 2
                     + radflux2(i,j,k,g) * area2(i,j,k) - radflux2(i,j+1,k,g) * area2(i,j+1,k) &
#endif
#if AMREX_SPACEDIM == 3
                     + radflux3(i,j,k,g) * area3(i,j,k) - radflux3(i,j,k+1,g) * area3(i,j,k+1) &
#endif
                     ) / vol(i,j,k)
             enddo
          enddo
       enddo
    enddo

    ! Add gradp term to momentum equation -- only for axisymmetry coords
    ! (and only for the radial flux);  also add the radiation pressure gradient
    ! to the momentum for all directions

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! pgdnv from the Riemann solver is only the gas contribution,
             ! not the radiation contribution.  Note that we've already included
             ! the gas pressure in the momentum flux for all Cartesian coordinate
             ! directions
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                dpdx = ( qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES))/ dx(1)
                update(i,j,k,UMX) = update(i,j,k,UMX) - dpdx
             endif

             ! radiation contribution -- this is sum{lambda E_r}
             dprdx = ZERO
             dprdy = ZERO
             dprdz = ZERO

             do g = 0, ngroups-1
#if AMREX_SPACEDIM == 1
                lamc = HALF*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 2
                lamc = FOURTH*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + &
                     qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 3
                lamc = (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + &
                     qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g) + &
                     qz(i,j,k,GDLAMS+g) + qz(i,j,k+1,GDLAMS+g) ) / 6.e0_rt
#endif

                dprdx = dprdx + lamc*(qx(i+1,j,k,GDERADS+g) - qx(i,j,k,GDERADS+g))/dx(1)
#if AMREX_SPACEDIM >= 2
                dprdy = dprdy + lamc*(qy(i,j+1,k,GDERADS+g) - qy(i,j,k,GDERADS+g))/dx(2)
#endif
#if AMREX_SPACEDIM == 3
                dprdz = dprdz + lamc*(qz(i,j,k+1,GDERADS+g) - qz(i,j,k,GDERADS+g))/dx(3)
#endif
             end do

             ! we now want to compute the change in the kinetic energy -- we
             ! base this off of uout, since that has the source terms in it.
             ! But note that this does not have the fluxes (since we are
             ! using update)

             ! note, we need to include the Z component here too (even
             ! for 1- and 2-d), since rotation might be in play

             urho_new = uout(i,j,k,URHO) + dt * update(i,j,k,URHO)

             ! this update includes the hydro fluxes and grad{p} from hydro
             umx_new1 = uout(i,j,k,UMX) + dt * update(i,j,k,UMX)
             umy_new1 = uout(i,j,k,UMY) + dt * update(i,j,k,UMY)
             umz_new1 = uout(i,j,k,UMZ) + dt * update(i,j,k,UMZ)

             ek1 = (umx_new1**2 + umy_new1**2 + umz_new1**2) / (TWO*urho_new)

             ! now add the radiation pressure gradient
             update(i,j,k,UMX) = update(i,j,k,UMX) - dprdx
             update(i,j,k,UMY) = update(i,j,k,UMY) - dprdy
             update(i,j,k,UMZ) = update(i,j,k,UMZ) - dprdz

             umx_new2 = umx_new1 - dt * dprdx
             umy_new2 = umy_new1 - dt * dprdy
             umz_new2 = umz_new1 - dt * dprdz

             ek2 = (umx_new2**2 + umy_new2**2 + umz_new2**2) / (TWO*urho_new)

             dek = ek2 - ek1

             ! the update is a source / dt, so scale accordingly
             update(i,j,k,UEDEN) = update(i,j,k,UEDEN) + dek/dt

             if (.not. comoving) then ! mixed-frame (single group only)
                Erout(i,j,k,0) = Erout(i,j,k,0) - dek
             end if

          end do
       end do
    end do

    ! Add radiation source terms to rho*u, rhoE, and Er
    if (comoving) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ux = HALF*(qx(i,j,k,GDU) + qx(i+1,j,k,GDU))
#if AMREX_SPACEDIM >= 2
                uy = HALF*(qy(i,j,k,GDV) + qy(i,j+dg(2),k,GDV))
#endif
#if AMREX_SPACEDIM == 3
                uz = HALF*(qz(i,j,k,GDW) + qz(i,j,k+dg(3),GDW))
#endif

                dudx(:) = ZERO
                dudy(:) = ZERO
                dudz(:) = ZERO

                dudx(1) = (qx(i+1,j,k,GDU) - qx(i,j,k,GDU))/dx(1)
                dudx(2) = (qx(i+1,j,k,GDV) - qx(i,j,k,GDV))/dx(1)
                dudx(3) = (qx(i+1,j,k,GDW) - qx(i,j,k,GDW))/dx(1)

#if AMREX_SPACEDIM >= 2
                dudy(1) = (qy(i,j+1,k,GDU) - qy(i,j,k,GDU))/dx(2)
                dudy(2) = (qy(i,j+1,k,GDV) - qy(i,j,k,GDV))/dx(2)
                dudy(3) = (qy(i,j+1,k,GDW) - qy(i,j,k,GDW))/dx(2)
#endif

#if AMREX_SPACEDIM == 3
                dudz(1) = (qz(i,j,k+1,GDU) - qz(i,j,k,GDU))/dx(3)
                dudz(2) = (qz(i,j,k+1,GDV) - qz(i,j,k,GDV))/dx(3)
                dudz(3) = (qz(i,j,k+1,GDW) - qz(i,j,k,GDW))/dx(3)
#endif

                divu = dudx(1) + dudy(2) + dudz(3)

                ! Note that for single group, fspace_type is always 1
                do g=0, ngroups-1

                   nhat(:) = ZERO

                   nhat(1) = (qx(i+1,j,k,GDERADS+g) - qx(i,j,k,GDERADS+g))/dx(1)
#if AMREX_SPACEDIM >= 2
                   nhat(2) = (qy(i,j+1,k,GDERADS+g) - qy(i,j,k,GDERADS+g))/dx(2)
#endif
#if AMREX_SPACEDIM == 3
                   nhat(3) = (qz(i,j,k+1,GDERADS+g) - qz(i,j,k,GDERADS+g))/dx(3)
#endif

                   GnDotu(1) = dot_product(nhat, dudx)
                   GnDotu(2) = dot_product(nhat, dudy)
                   GnDotu(3) = dot_product(nhat, dudz)

                   nnColonDotGu = dot_product(nhat, GnDotu) / (dot_product(nhat,nhat)+1.e-50_rt)

#if AMREX_SPACEDIM == 1
                   lamc = HALF*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 2
                   lamc = 0.25e0_rt*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + &
                        qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 3
                   lamc = (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + &
                        qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g) + &
                        qz(i,j,k,GDLAMS+g) + qz(i,j,k+1,GDLAMS+g) ) / 6.e0_rt
#endif

                   Eddf = Edd_factor(lamc)
                   f1 = (ONE-Eddf)*HALF
                   f2 = (3.e0_rt*Eddf-ONE)*HALF
                   af(g) = -(f1*divu + f2*nnColonDotGu)

                   if (fspace_type .eq. 1) then
                      Eddfxp = Edd_factor(qx(i+1,j  ,k  ,GDLAMS+g))
                      Eddfxm = Edd_factor(qx(i  ,j  ,k  ,GDLAMS+g))
#if AMREX_SPACEDIM >= 2
                      Eddfyp = Edd_factor(qy(i  ,j+1,k  ,GDLAMS+g))
                      Eddfym = Edd_factor(qy(i  ,j  ,k  ,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 3
                      Eddfzp = Edd_factor(qz(i  ,j  ,k+1,GDLAMS+g))
                      Eddfzm = Edd_factor(qz(i  ,j  ,k  ,GDLAMS+g))
#endif

                      f1xp = HALF*(ONE-Eddfxp)
                      f1xm = HALF*(ONE-Eddfxm)
#if AMREX_SPACEDIM >= 2
                      f1yp = HALF*(ONE-Eddfyp)
                      f1ym = HALF*(ONE-Eddfym)
#endif
#if AMREX_SPACEDIM == 3
                      f1zp = HALF*(ONE-Eddfzp)
                      f1zm = HALF*(ONE-Eddfzm)
#endif

                      Gf1E(1) = (f1xp*qx(i+1,j,k,GDERADS+g) - f1xm*qx(i,j,k,GDERADS+g)) / dx(1)
#if AMREX_SPACEDIM >= 2
                      Gf1E(2) = (f1yp*qy(i,j+1,k,GDERADS+g) - f1ym*qy(i,j,k,GDERADS+g)) / dx(2)
#endif
#if AMREX_SPACEDIM == 3
                      Gf1E(3) = (f1zp*qz(i,j,k+1,GDERADS+g) - f1zm*qz(i,j,k,GDERADS+g)) / dx(3)
#endif


#if AMREX_SPACEDIM == 1
                      Egdc = HALF*(qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g))
                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*ux*Gf1E(1) &
                           - dt*f2*Egdc*nnColonDotGu
#endif
#if AMREX_SPACEDIM == 2
                      Egdc = 0.25e0_rt*(qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g) + &
                           qy(i,j,k,GDERADS+g) + qy(i,j+1,k,GDERADS+g))
                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*(ux*Gf1E(1)+uy*Gf1E(2)) &
                           - dt*f2*Egdc*nnColonDotGu
#endif
#if AMREX_SPACEDIM == 3
                      Egdc = (qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g) &
                           +  qy(i,j,k,GDERADS+g) + qy(i,j+1,k,GDERADS+g) &
                           +  qz(i,j,k,GDERADS+g) + qz(i,j,k+1,GDERADS+g) ) / 6.e0_rt
                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*(ux*Gf1E(1)+uy*Gf1E(2)+uz*Gf1E(3)) &
                           - dt*f2*Egdc*nnColonDotGu
#endif

                   end if

                end do

                if (ngroups.gt.1) then
                   ustar = Erout(i,j,k,:) / Erscale
                   call advect_in_fspace(ustar, af, dt, nstep_fsp)
                   Erout(i,j,k,:) = ustar * Erscale
                end if
             enddo
          enddo
       enddo
    endif
#endif

  end subroutine ctu_consup

  subroutine ca_track_grid_losses(lo, hi, &
                                  flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
                                  flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                                  flux3, flux3_lo, flux3_hi, &
#endif
                                  mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                                  eden_lost, xang_lost, yang_lost, zang_lost) &
                                  bind(C, name="ca_track_grid_losses")


    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, NVAR
    use amrinfo_module, only : amr_level
    use prob_params_module, only : domlo_level, domhi_level, center
    use castro_util_module, only: position ! function
    use castro_util_module, only: linear_to_angular_momentum ! function
    use amrex_fort_module, only: amrex_reduce_add

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    real(rt), intent(in) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    real(rt), intent(in) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    real(rt), intent(in) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
#endif

    real(rt), intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    real(rt), intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

    real(rt) :: loc(3), ang_mom(3), flux(3)
    integer :: domlo(3), domhi(3)
    integer :: i, j, k

    !$gpu

    domlo = domlo_level(:,amr_level)
    domhi = domhi_level(:,amr_level)

#if AMREX_SPACEDIM == 3
    if (lo(3) .le. domlo(3) .and. hi(3) .ge. domlo(3)) then

       k = domlo(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccz=.false.) - center

             call amrex_reduce_add(mass_lost, -flux3(i,j,k,URHO))
             call amrex_reduce_add(xmom_lost, -flux3(i,j,k,UMX))
             call amrex_reduce_add(ymom_lost, -flux3(i,j,k,UMY))
             call amrex_reduce_add(zmom_lost, -flux3(i,j,k,UMZ))
             call amrex_reduce_add(eden_lost, -flux3(i,j,k,UEDEN))

             flux(:) = flux3(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call amrex_reduce_add(xang_lost, -ang_mom(1))
             call amrex_reduce_add(yang_lost, -ang_mom(2))
             call amrex_reduce_add(zang_lost, -ang_mom(3))

          enddo
       enddo

    endif

    if (lo(3) .le. domhi(3) .and. hi(3) .ge. domhi(3)) then

       k = domhi(3) + 1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccz=.false.) - center

             call amrex_reduce_add(mass_lost, flux3(i,j,k,URHO))
             call amrex_reduce_add(xmom_lost, flux3(i,j,k,UMX))
             call amrex_reduce_add(ymom_lost, flux3(i,j,k,UMY))
             call amrex_reduce_add(zmom_lost, flux3(i,j,k,UMZ))
             call amrex_reduce_add(eden_lost, flux3(i,j,k,UEDEN))

             flux(:) = flux3(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call amrex_reduce_add(xang_lost, ang_mom(1))
             call amrex_reduce_add(yang_lost, ang_mom(2))
             call amrex_reduce_add(zang_lost, ang_mom(3))

          enddo
       enddo

    endif
#endif

#if AMREX_SPACEDIM >= 2
    if (lo(2) .le. domlo(2) .and. hi(2) .ge. domlo(2)) then

       j = domlo(2)
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccy=.false.) - center

             call amrex_reduce_add(mass_lost, -flux2(i,j,k,URHO))
             call amrex_reduce_add(xmom_lost, -flux2(i,j,k,UMX))
             call amrex_reduce_add(ymom_lost, -flux2(i,j,k,UMY))
             call amrex_reduce_add(zmom_lost, -flux2(i,j,k,UMZ))
             call amrex_reduce_add(eden_lost, -flux2(i,j,k,UEDEN))

             flux(:) = flux2(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call amrex_reduce_add(xang_lost, -ang_mom(1))
             call amrex_reduce_add(yang_lost, -ang_mom(2))
             call amrex_reduce_add(zang_lost, -ang_mom(3))

          enddo
       enddo

    endif

    if (lo(2) .le. domhi(2) .and. hi(2) .ge. domhi(2)) then

       j = domhi(2) + 1
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccy=.false.) - center

             call amrex_reduce_add(mass_lost, flux2(i,j,k,URHO))
             call amrex_reduce_add(xmom_lost, flux2(i,j,k,UMX))
             call amrex_reduce_add(ymom_lost, flux2(i,j,k,UMY))
             call amrex_reduce_add(zmom_lost, flux2(i,j,k,UMZ))
             call amrex_reduce_add(eden_lost, flux2(i,j,k,UEDEN))

             flux(:) = flux2(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call amrex_reduce_add(xang_lost, ang_mom(1))
             call amrex_reduce_add(yang_lost, ang_mom(2))
             call amrex_reduce_add(zang_lost, ang_mom(3))

          enddo
       enddo

    endif
#endif

    if (lo(1) .le. domlo(1) .and. hi(1) .ge. domlo(1)) then

       i = domlo(1)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)

             loc = position(i,j,k,ccx=.false.) - center

             call amrex_reduce_add(mass_lost, -flux1(i,j,k,URHO))
             call amrex_reduce_add(xmom_lost, -flux1(i,j,k,UMX))
             call amrex_reduce_add(ymom_lost, -flux1(i,j,k,UMY))
             call amrex_reduce_add(zmom_lost, -flux1(i,j,k,UMZ))
             call amrex_reduce_add(eden_lost, -flux1(i,j,k,UEDEN))

             flux(:) = flux1(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call amrex_reduce_add(xang_lost, -ang_mom(1))
             call amrex_reduce_add(yang_lost, -ang_mom(2))
             call amrex_reduce_add(zang_lost, -ang_mom(3))

          enddo
       enddo

    endif

    if (lo(1) .le. domhi(1) .and. hi(1) .ge. domhi(1)) then

       i = domhi(1) + 1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)

             loc = position(i,j,k,ccx=.false.) - center

             call amrex_reduce_add(mass_lost, flux1(i,j,k,URHO))
             call amrex_reduce_add(xmom_lost, flux1(i,j,k,UMX))
             call amrex_reduce_add(ymom_lost, flux1(i,j,k,UMY))
             call amrex_reduce_add(zmom_lost, flux1(i,j,k,UMZ))
             call amrex_reduce_add(eden_lost, flux1(i,j,k,UEDEN))

             flux(:) = flux1(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call amrex_reduce_add(xang_lost, ang_mom(1))
             call amrex_reduce_add(yang_lost, ang_mom(2))
             call amrex_reduce_add(zang_lost, ang_mom(3))

          enddo
       enddo

    endif

  end subroutine ca_track_grid_losses

end module ctu_module
