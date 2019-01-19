! advection routines in support of the CTU unsplit advection scheme

module ctu_module

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains


  !> @brief Compute the normal interface states by reconstructing
  !! the primitive variables and doing characteristic tracing.  We
  !! do not apply the transverse terms here.
  !!
  !! @todo we can get rid of the the different temporary q Godunov
  !! state arrays
  !!
  !! @param[in] q            (const)  input state, primitives
  !! @param[in] qaux         (const)  auxiliary hydro data
  !! @param[in] flatn        (const)  flattening parameter
  !! @param[in] srcQ         (const)  primitive variable source
  !! @param[in] dx           (const)  grid spacing in X, Y, Z direction
  !! @param[in] dt           (const)  time stepsize
  !! @param[inout] flux1        (modify) flux in X direction on X edges
  !! @param[inout] flux2        (modify) flux in Y direction on Y edges
  !! @param[inout] flux3        (modify) flux in Z direction on Z edges
  !! @param[inout] q1           (modify) Godunov interface state in X
  !! @param[inout] q2           (modify) Godunov interface state in Y
  !! @param[inout] q3           (modify) Godunov interface state in Z
  !!
  subroutine ctu_normal_states(lo, hi, &
                               vlo, vhi, &
                               q, qd_lo, qd_hi, &
                               flatn, f_lo, f_hi, &
                               qaux, qa_lo, qa_hi, &
                               srcQ, src_lo, src_hi, &
                               shk, sk_lo, sk_hi, &
                               Ip, Ip_lo, Ip_hi, &
                               Im, Im_lo, Im_hi, &
                               Ip_src, Ips_lo, Ips_hi, &
                               Im_src, Ims_lo, Ims_hi, &
                               Ip_gc, Ipg_lo, Ipg_hi, &
                               Im_gc, Img_lo, Img_hi, &
                               dq, dq_lo, dq_hi, &
                               sm, sm_lo, sm_hi, &
                               sp, sp_lo, sp_hi, &
                               qxm, qxm_lo, qxm_hi, &
                               qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM >= 2
                               qym, qym_lo, qym_hi, &
                               qyp, qyp_lo, qyp_hi, &
#endif
#if AMREX_SPACEDIM == 3
                               qzm, qzm_lo, qzm_hi, &
                               qzp, qzp_lo, qzp_hi, &
#endif
                               dx, dt, &
#if AMREX_SPACEDIM < 3
                               dloga, dloga_lo, dloga_hi, &
#endif
                               domlo, domhi)

    ! everything in this routine has the same lo:hi requirements (we fill one ghost cell)

    use meth_params_module, only : QVAR, NQ, NVAR, &
                                   QFS, QFX, QTEMP, QREINT, &
                                   QC, QGAMC, NQAUX, QGAME, QREINT, &
                                   NGDNV, GDU, GDV, GDW, GDPRES, &
                                   ppm_type, ppm_predict_gammae, &
                                   plm_iorder, use_pslope, ppm_temp_fix, &
                                   hybrid_riemann
    use trace_plm_module, only : trace_plm
#if AMREX_SPACEDIM >= 2
    use transverse_module
#endif
    use ppm_module, only : ca_ppm_reconstruct, ppm_int_profile, ppm_reconstruct_with_eos
    use slope_module, only : uslope, pslope
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use trace_ppm_rad_module, only : trace_ppm_rad
#else
    use trace_ppm_module, only : trace_ppm
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
    use advection_util_module, only : ca_shock
    use prob_params_module, only : dg

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: sk_lo(3), sk_hi(3)
    integer, intent(in) :: Ip_lo(3), Ip_hi(3)
    integer, intent(in) :: Im_lo(3), Im_hi(3)
    integer, intent(in) :: Ips_lo(3), Ips_hi(3)
    integer, intent(in) :: Ims_lo(3), Ims_hi(3)
    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)
    integer, intent(in) :: dq_lo(3), dq_hi(3)
    integer, intent(in) :: sm_lo(3), sm_hi(3)
    integer, intent(in) :: sp_lo(3), sp_hi(3)
    integer, intent(in) :: qxm_lo(3), qxm_hi(3)
    integer, intent(in) :: qxp_lo(3), qxp_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: qym_lo(3), qym_hi(3)
    integer, intent(in) :: qyp_lo(3), qyp_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: qzm_lo(3), qzm_hi(3)
    integer, intent(in) :: qzp_lo(3), qzp_hi(3)
#endif
#if AMREX_SPACEDIM < 3
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    real(rt), intent(in) :: dx(3), dt
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)

    real(rt), intent(inout) :: shk(sk_lo(1):sk_hi(1), sk_lo(2):sk_hi(2), sk_lo(3):sk_hi(3))
    real(rt), intent(inout) :: Ip(Ip_lo(1):Ip_hi(1),Ip_lo(2):Ip_hi(2),Ip_lo(3):Ip_hi(3),1:AMREX_SPACEDIM,1:3,NQ)
    real(rt), intent(inout) :: Im(Im_lo(1):Im_hi(1),Im_lo(2):Im_hi(2),Im_lo(3):Im_hi(3),1:AMREX_SPACEDIM,1:3,NQ)
    real(rt), intent(inout) :: Ip_src(Ips_lo(1):Ips_hi(1),Ips_lo(2):Ips_hi(2),Ips_lo(3):Ips_hi(3),1:AMREX_SPACEDIM,1:3,QVAR)
    real(rt), intent(inout) :: Im_src(Ims_lo(1):Ims_hi(1),Ims_lo(2):Ims_hi(2),Ims_lo(3):Ims_hi(3),1:AMREX_SPACEDIM,1:3,QVAR)
    real(rt), intent(inout) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:AMREX_SPACEDIM,1:3,1)
    real(rt), intent(inout) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:AMREX_SPACEDIM,1:3,1)
    real(rt), intent(inout) :: dq(dq_lo(1):dq_hi(1), dq_lo(2):dq_hi(2), dq_lo(3):dq_hi(3), NQ, AMREX_SPACEDIM)
    real(rt), intent(inout) :: sm(sm_lo(1):sm_hi(1), sm_lo(2):sm_hi(2), sm_lo(3):sm_hi(3))
    real(rt), intent(inout) :: sp(sp_lo(1):sp_hi(1), sp_lo(2):sp_hi(2), sp_lo(3):sp_hi(3))

    real(rt), intent(inout) :: qxm(qxm_lo(1):qxm_hi(1), qxm_lo(2):qxm_hi(2), qxm_lo(3):qxm_hi(3), NQ)
    real(rt), intent(inout) :: qxp(qxp_lo(1):qxp_hi(1), qxp_lo(2):qxp_hi(2), qxp_lo(3):qxp_hi(3), NQ)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: qym(qym_lo(1):qym_hi(1), qym_lo(2):qym_hi(2), qym_lo(3):qym_hi(3), NQ)
    real(rt), intent(inout) :: qyp(qyp_lo(1):qyp_hi(1), qyp_lo(2):qyp_hi(2), qyp_lo(3):qyp_hi(3), NQ)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: qzm(qzm_lo(1):qzm_hi(1), qzm_lo(2):qzm_hi(2), qzm_lo(3):qzm_hi(3), NQ)
    real(rt), intent(inout) :: qzp(qzp_lo(1):qzp_hi(1), qzp_lo(2):qzp_hi(2), qzp_lo(3):qzp_hi(3), NQ)
#endif
#if AMREX_SPACEDIM < 3
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt) :: hdt
    integer :: i, j, k, n

    logical :: source_nonzero(QVAR)
    logical :: reconstruct_state(NQ)

    logical :: compute_shock

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
    if (ppm_temp_fix < 3) then
       reconstruct_state(QTEMP) = .false.
    endif

    ! preprocess the sources -- we don't want to trace under a source
    ! that is empty.  Note, we need to do this check over the entire
    ! grid, to be sure, e.g., use vlo:vhi.  We cannot rely on lo:hi,
    ! since that may just be a single zone on the GPU.
    if (ppm_type > 0) then
       do n = 1, QVAR
          if (minval(srcQ(vlo(1)-2:vhi(1)+2,vlo(2)-2*dg(2):vhi(2)+2*dg(2),vlo(3)-2*dg(3):vhi(3)+2*dg(3),n)) == ZERO .and. &
              maxval(srcQ(vlo(1)-2:vhi(1)+2,vlo(2)-2*dg(2):vhi(2)+2*dg(2),vlo(3)-2*dg(3):vhi(3)+2*dg(3),n)) == ZERO) then
             source_nonzero(n) = .false.
          else
             source_nonzero(n) = .true.
          endif
       enddo

       ! Compute Ip and Im -- this does the parabolic reconstruction,
       ! limiting, and returns the integral of each profile under each
       ! wave to each interface
       do n = 1, NQ
          if (.not. reconstruct_state(n)) cycle

          call ca_ppm_reconstruct(lo, hi, 0, &
                                  q, qd_lo, qd_hi, NQ, n, n, &
                                  flatn, f_lo, f_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  1, 1, 1)

          call ppm_int_profile(lo, hi, &
                               q, qd_lo, qd_hi, NQ, n, &
                               q, qd_lo, qd_hi, &
                               qaux, qa_lo, qa_hi, &
                               sm, sm_lo, sm_hi, &
                               sp, sp_lo, sp_hi, &
                               Ip, Ip_lo, Ip_hi, &
                               Im, Im_lo, Im_hi, NQ, n, &
                               dx, dt)
       end do


       if (ppm_temp_fix /= 1) then
          call ca_ppm_reconstruct(lo, hi, 0, &
                                  qaux, qa_lo, qa_hi, NQAUX, QGAMC, QGAMC, &
                                  flatn, f_lo, f_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  1, 1, 1)

          call ppm_int_profile(lo, hi, &
                               qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               q, qd_lo, qd_hi, &
                               qaux, qa_lo, qa_hi, &
                               sm, sm_lo, sm_hi, &
                               sp, sp_lo, sp_hi, &
                               Ip_gc, Ipg_lo, Ipg_hi, &
                               Im_gc, Img_lo, Img_hi, 1, 1, &
                               dx, dt)
       else

          ! temperature-based PPM
          call ppm_reconstruct_with_eos(lo, hi, &
                                        Ip, Ip_lo, Ip_hi, &
                                        Im, Im_lo, Im_hi, &
                                        Ip_gc, Ipg_lo, Ipg_hi, &
                                        Im_gc, Img_lo, Img_hi)

       end if


       ! source terms
       do n = 1, QVAR
          if (source_nonzero(n)) then
             call ca_ppm_reconstruct(lo, hi, 0, &
                                     srcQ, src_lo, src_hi, QVAR, n, n, &
                                     flatn, f_lo, f_hi, &
                                     sm, sm_lo, sm_hi, &
                                     sp, sp_lo, sp_hi, &
                                     1, 1, 1)

             call ppm_int_profile(lo, hi, &
                                  srcQ, src_lo, src_hi, QVAR, n, &
                                  q, qd_lo, qd_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  sm, sm_lo, sm_hi, &
                                  sp, sp_lo, sp_hi, &
                                  Ip_src, Ips_lo, Ips_hi, &
                                  Im_src, Ims_lo, Ims_hi, QVAR, n, &
                                  dx, dt)
          else
             Ip_src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:,:,n) = ZERO
             Im_src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:,:,n) = ZERO
          endif

       enddo


       ! compute the interface states

#ifdef RADIATION
       call trace_ppm_rad(lo, hi, &
                          1, q, qd_lo, qd_hi, &
                          qaux, qa_lo, qa_hi, &
                          Ip, Ip_lo, Ip_hi, &
                          Im, Im_lo, Im_hi, &
                          Ip_src, Ips_lo, Ips_hi, &
                          Im_src, Ims_lo, Ims_hi, &
                          qxm, qxm_lo, qxm_hi, &
                          qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM <= 2
                          dloga, dloga_lo, dloga_hi, &
#endif
                          vlo, vhi, domlo, domhi, &
                          dx, dt)

#if AMREX_SPACEDIM >= 2
       call trace_ppm_rad(lo, hi, &
                          2, q, qd_lo, qd_hi, &
                          qaux, qa_lo, qa_hi, &
                          Ip, Ip_lo, Ip_hi, &
                          Im, Im_lo, Im_hi, &
                          Ip_src, Ips_lo, Ips_hi, &
                          Im_src, Ims_lo, Ims_hi, &
                          qym, qym_lo, qym_hi, &
                          qyp, qyp_lo, qyp_hi, &
#if AMREX_SPACEDIM == 2
                          dloga, dloga_lo, dloga_hi, &
#endif
                          vlo, vhi, domlo, domhi, &
                          dx, dt)
#endif

#if AMREX_SPACEDIM == 3
       call trace_ppm_rad(lo, hi, &
                          3, q, qd_lo, qd_hi, &
                          qaux, qa_lo, qa_hi, &
                          Ip, Ip_lo, Ip_hi, &
                          Im, Im_lo, Im_hi, &
                          Ip_src, Ips_lo, Ips_hi, &
                          Im_src, Ims_lo, Ims_hi, &
                          qzm, qzm_lo, qzm_hi, &
                          qzp, qzp_lo, qzp_hi, &
                          vlo, vhi, domlo, domhi, &
                          dx, dt)
#endif

#else
       call trace_ppm(lo, hi, &
                      1, q, qd_lo, qd_hi, &
                      qaux, qa_lo, qa_hi, &
                      Ip, Ip_lo, Ip_hi, &
                      Im, Im_lo, Im_hi, &
                      Ip_src, Ips_lo, Ips_hi, &
                      Im_src, Ims_lo, Ims_hi, &
                      Ip_gc, Ipg_lo, Ipg_hi, &
                      Im_gc, Img_lo, Img_hi, &
                      qxm, qxm_lo, qxm_hi, &
                      qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM <= 2
                      dloga, dloga_lo, dloga_hi, &
#endif
                      vlo, vhi, domlo, domhi, &
                      dx, dt)

#if AMREX_SPACEDIM >= 2
       call trace_ppm(lo, hi, &
                      2, q, qd_lo, qd_hi, &
                      qaux, qa_lo, qa_hi, &
                      Ip, Ip_lo, Ip_hi, &
                      Im, Im_lo, Im_hi, &
                      Ip_src, Ips_lo, Ips_hi, &
                      Im_src, Ims_lo, Ims_hi, &
                      Ip_gc, Ipg_lo, Ipg_hi, &
                      Im_gc, Img_lo, Img_hi, &
                      qym, qym_lo, qym_hi, &
                      qyp, qyp_lo, qyp_hi, &
#if AMREX_SPACEDIM == 2
                      dloga, dloga_lo, dloga_hi, &
#endif
                      vlo, vhi, domlo, domhi, &
                      dx, dt)
#endif

#if AMREX_SPACEDIM == 3
       call trace_ppm(lo, hi, &
                      3, q, qd_lo, qd_hi, &
                      qaux, qa_lo, qa_hi, &
                      Ip, Ip_lo, Ip_hi, &
                      Im, Im_lo, Im_hi, &
                      Ip_src, Ips_lo, Ips_hi, &
                      Im_src, Ims_lo, Ims_hi, &
                      Ip_gc, Ipg_lo, Ipg_hi, &
                      Im_gc, Img_lo, Img_hi, &
                      qzm, qzm_lo, qzm_hi, &
                      qzp, qzp_lo, qzp_hi, &
                      vlo, vhi, domlo, domhi, &
                      dx, dt)
#endif

#endif
    else
       ! PLM

#ifdef RADIATION
#ifndef AMREX_USE_CUDA
       call amrex_error("ppm_type <=0 is not supported in with radiation")
#endif
#endif
       ! Compute all slopes
       do n = 1, NQ
          if (.not. reconstruct_state(n)) cycle
          call uslope(lo, hi, &
                      q, qd_lo, qd_hi, n, &
                      flatn, f_lo, f_hi, &
                      dq, dq_lo, dq_hi)
       end do

       if (use_pslope == 1) then
          call pslope(lo, hi, &
                      q, qd_lo, qd_hi, &
                      flatn, f_lo, f_hi, &
                      dq, dq_lo, dq_hi, &
                      srcQ, src_lo, src_hi, &
                      dx)
       endif


       ! compute the interface states

       call trace_plm(lo, hi, &
                      1, q, qd_lo, qd_hi, &
                      qaux, qa_lo, qa_hi, &
                      dq, dq_lo, dq_hi, &
                      qxm, qxm_lo, qxm_hi, &
                      qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM < 3
                      dloga, dloga_lo, dloga_hi, &
#endif
                      SrcQ, src_lo, src_hi, &
                      vlo, vhi, domlo, domhi, &
                      dx, dt)

#if AMREX_SPACEDIM >= 2
       call trace_plm(lo, hi, &
                      2, q, qd_lo, qd_hi, &
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
#endif

#if AMREX_SPACEDIM == 3
       call trace_plm(lo, hi, &
                      3, q, qd_lo, qd_hi, &
                      qaux, qa_lo, qa_hi, &
                      dq, dq_lo, dq_hi, &
                      qzm, qzm_lo, qzm_hi, &
                      qzp, qzp_lo, qzp_hi, &
                      SrcQ, src_lo, src_hi, &
                      vlo, vhi, domlo, domhi, &
                      dx, dt)
#endif

    end if  ! ppm test

  end subroutine ctu_normal_states


  !> @brief Compute hyperbolic fluxes using unsplit second
  !! order Godunov integrator.
  !!
  !! @todo we can get rid of the the different temporary q Godunov
  !! state arrays
  !!
  !! @param[in] q            (const)  input state, primitives
  !! @param[in] qaux         (const)  auxiliary hydro data
  !! @param[in] flatn        (const)  flattening parameter
  !! @param[in] srcQ         (const)  primitive variable source
  !! @param[in] dx           (const)  grid spacing in X, Y, Z direction
  !! @param[in] dt           (const)  time stepsize
  !! @param[inout] flux1        (modify) flux in X direction on X edges
  !! @param[inout] flux2        (modify) flux in Y direction on Y edges
  !! @param[inout] flux3        (modify) flux in Z direction on Z edges
  !! @param[inout] q1           (modify) Godunov interface state in X
  !! @param[inout] q2           (modify) Godunov interface state in Y
  !! @param[inout] q3           (modify) Godunov interface state in Z
  !!
  subroutine ctu_compute_fluxes(lo, hi, &
                                q, qd_lo, qd_hi, &
                                qaux, qa_lo, qa_hi, &
                                dx, dt, &
                                shk, sk_lo, sk_hi, &
                                qxm, qxm_lo, qxm_hi, &
                                qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM >= 2
                                qym, qym_lo, qym_hi, &
                                qyp, qyp_lo, qyp_hi, &
#endif
#if AMREX_SPACEDIM == 3
                                qzm, qzm_lo, qzm_hi, &
                                qzp, qzp_lo, qzp_hi, &
#endif
                                uout, uout_lo, uout_hi, &
                                flux1, f1_lo, f1_hi, &
#if AMREX_SPACEDIM >= 2
                                flux2, f2_lo, f2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                                flux3, f3_lo, f3_hi, &
#endif
#ifdef RADIATION
                                rflux1, rf1_lo, rf1_hi, &
#if AMREX_SPACEDIM >= 2
                                rflux2, rf2_lo, rf2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                                rflux3, rf3_lo, rf3_hi, &
#endif
#endif
                                q1, q1_lo, q1_hi, &
#if AMREX_SPACEDIM >= 2
                                q2, q2_lo, q2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                                q3, q3_lo, q3_hi, &
#endif
#if AMREX_SPACEDIM <= 2
                                area1, area1_lo, area1_hi, &
#endif
#if AMREX_SPACEDIM == 2
                                area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM <= 2
                                vol, vol_lo, vol_hi, &
#endif
                                domlo, domhi)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QVAR, NQ, NVAR, &
                                   QFS, QFX, QTEMP, QREINT, &
                                   QC, QGAMC, NQAUX, QGAME, QREINT, &
                                   NGDNV, GDU, GDV, GDW, GDPRES, &
                                   hybrid_riemann
#if AMREX_SPACEDIM >= 2
    use transverse_module
#endif
    use riemann_module, only: cmpflx
    use riemann_util_module, only : ca_store_godunov_state
#ifdef RADIATION
    use rad_params_module, only : ngroups
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
    use prob_params_module, only : dg

    implicit none

    integer, intent(in) :: lo(3), hi(3)

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    integer, intent(in) :: sk_lo(3), sk_hi(3)
    integer, intent(in) :: qxm_lo(3), qxm_hi(3)
    integer, intent(in) :: qxp_lo(3), qxp_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: qym_lo(3), qym_hi(3)
    integer, intent(in) :: qyp_lo(3), qyp_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: qzm_lo(3), qzm_hi(3)
    integer, intent(in) :: qzp_lo(3), qzp_hi(3)
#endif

    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: f1_lo(3), f1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: f2_lo(3), f2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: f3_lo(3), f3_hi(3)
    integer, intent(in) :: q3_lo(3), q3_hi(3)
#endif
#if AMREX_SPACEDIM <= 2
    integer, intent(in) :: area1_lo(3), area1_hi(3)
#endif
#if AMREX_SPACEDIM == 2
    integer, intent(in) :: area2_lo(3), area2_hi(3)
#endif
#if AMREX_SPACEDIM <= 2
    integer, intent(in) :: vol_lo(3), vol_hi(3)
#endif

#ifdef RADIATION
    integer, intent(in) :: rf1_lo(3), rf1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: rf2_lo(3), rf2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: rf3_lo(3), rf3_hi(3)
#endif
#endif
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: dx(3), dt

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: shk(sk_lo(1):sk_hi(1), sk_lo(2):sk_hi(2), sk_lo(3):sk_hi(3))

    real(rt), intent(inout) :: qxm(qxm_lo(1):qxm_hi(1), qxm_lo(2):qxm_hi(2), qxm_lo(3):qxm_hi(3), NQ)
    real(rt), intent(inout) :: qxp(qxp_lo(1):qxp_hi(1), qxp_lo(2):qxp_hi(2), qxp_lo(3):qxp_hi(3), NQ)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: qym(qym_lo(1):qym_hi(1), qym_lo(2):qym_hi(2), qym_lo(3):qym_hi(3), NQ)
    real(rt), intent(inout) :: qyp(qyp_lo(1):qyp_hi(1), qyp_lo(2):qyp_hi(2), qyp_lo(3):qyp_hi(3), NQ)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: qzm(qzm_lo(1):qzm_hi(1), qzm_lo(2):qzm_hi(2), qzm_lo(3):qzm_hi(3), NQ)
    real(rt), intent(inout) :: qzp(qzp_lo(1):qzp_hi(1), qzp_lo(2):qzp_hi(2), qzp_lo(3):qzp_hi(3), NQ)
#endif
    real(rt), intent(inout) ::  uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(inout) :: flux1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NVAR)
    real(rt), intent(inout) ::    q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: flux2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3),NVAR)
    real(rt), intent(inout) ::    q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: flux3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3),NVAR)
    real(rt), intent(inout) ::    q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
#endif

#ifdef RADIATION
    real(rt), intent(inout) :: rflux1(rf1_lo(1):rf1_hi(1),rf1_lo(2):rf1_hi(2),rf1_lo(3):rf1_hi(3),0:ngroups-1)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: rflux2(rf2_lo(1):rf2_hi(1),rf2_lo(2):rf2_hi(2),rf2_lo(3):rf2_hi(3),0:ngroups-1)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: rflux3(rf3_lo(1):rf3_hi(1),rf3_lo(2):rf3_hi(2),rf3_lo(3):rf3_hi(3),0:ngroups-1)
#endif
#endif
#if AMREX_SPACEDIM <= 2
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
#endif
#if AMREX_SPACEDIM == 2
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
#endif

    real(rt) :: dxinv, dyinv, dzinv
    real(rt) :: dtdx, dtdy, dtdz, hdt
#if AMREX_SPACEDIM == 3
    real(rt) :: cdtdx, cdtdy, cdtdz
#endif
    real(rt) :: hdtdx, hdtdy, hdtdz

    integer :: i, j, k, iwave, n

    real(rt), pointer :: q_int(:,:,:,:)

#ifdef RADIATION
    real(rt), pointer :: lambda_int(:,:,:,:)
#endif


    ! Left and right state arrays (edge centered, cell centered)
    double precision, dimension(:,:,:,:), pointer :: &
         ql, qr, &
         qmxy, qpxy, qmxz, qpxz, qmyx, qpyx, &
         qmyz, qpyz, qmzx, qpzx, qmzy, qpzy, &
         qxl, qxr, qyl, qyr, qzl, qzr

    double precision, dimension(:,:,:,:), pointer:: &
         fx, fy, fz, fxy, fxz, fyx, fyz, fzx, fzy

#ifdef RADIATION
    double precision, dimension(:,:,:,:), pointer:: &
         rfx, rfy, rfz, rfxy, rfxz, rfyx, rfyz, rfzx, rfzy
#endif

    double precision, dimension(:,:,:,:), pointer:: &
         qgdnvx, qgdnvy, qgdnvz, &
         qgdnvxy, qgdnvxz, &
         qgdnvyx, qgdnvyz, &
         qgdnvzx, qgdnvzy

    ! these will be the temporary arrays we actually allocate space for
    double precision, dimension(:,:,:,:), pointer :: ftmp1, ftmp2, rftmp1, rftmp2
    double precision, dimension(:,:,:,:), pointer :: qgdnvtmp1, qgdnvtmp2

    integer :: fglo(3), fghi(3), glo(3), ghi(3)

    fglo = lo - dg(:)  ! face + one ghost
    fghi = hi + 2*dg(:)

    glo = lo - dg(:)  ! one ghost,  this can be used for face-based arrays too
    ghi = hi + dg(:)

    ! Local constants
    dxinv = ONE/dx(1)
    dtdx = dt*dxinv
    hdtdx = HALF*dtdx

#if AMREX_SPACEDIM >= 2
    dyinv = ONE/dx(2)
    dtdy = dt*dyinv
    hdtdy = HALF*dtdy
#endif

#if AMREX_SPACEDIM == 3
    dzinv = ONE/dx(3)
    dtdz = dt*dzinv
    hdtdz = HALF*dtdz
    cdtdx = dtdx*THIRD
    cdtdy = dtdy*THIRD
    cdtdz = dtdz*THIRD
#endif

    hdt = HALF*dt

    integer, intent(in) :: sk_lo(3), sk_hi(3)
    integer, intent(in) :: qxm_lo(3), qxm_hi(3)
    integer, intent(in) :: qxp_lo(3), qxp_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: qym_lo(3), qym_hi(3)
    integer, intent(in) :: qyp_lo(3), qyp_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: qzm_lo(3), qzm_hi(3)
    integer, intent(in) :: qzp_lo(3), qzp_hi(3)
    integer, intent(in) :: qmxy_lo(3), qmxy_hi(3)
    integer, intent(in) :: qpxy_lo(3), qpxy_hi(3)
    integer, intent(in) :: qmxz_lo(3), qmxz_hi(3)
    integer, intent(in) :: qpxz_lo(3), qpxz_hi(3)
    integer, intent(in) :: qmyx_lo(3), qmyx_hi(3)
    integer, intent(in) :: qpyx_lo(3), qpyx_hi(3)
    integer, intent(in) :: qmyz_lo(3), qmyz_hi(3)
    integer, intent(in) :: qpyz_lo(3), qpyz_hi(3)
    integer, intent(in) :: qmzx_lo(3), qmzx_hi(3)
    integer, intent(in) :: qpzx_lo(3), qpzx_hi(3)
    integer, intent(in) :: qmzy_lo(3), qmzy_hi(3)
    integer, intent(in) :: qpzy_lo(3), qpzy_hi(3)
#endif
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: ft1_lo(3), ft1_hi(3)
    integer, intent(in) :: ft2_lo(3), ft2_hi(3)
#ifdef RADIATION
    integer, intent(in) :: rft1_lo(3), rft1_hi(3)
    integer, intent(in) :: rft2_lo(3), rft2_hi(3)
#endif
    integer, intent(in) :: qg1_lo(3), qg1_hi(3)
    integer, intent(in) :: qg2_lo(3), qg2_hi(3)
    integer, intent(in) :: ql_lo(3), ql_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
#endif
    integer, intent(in) :: qi_lo(3), qi_hi(3)
#ifdef RADIATION
    integer, intent(in) :: li_lo(3), li_hi(3)
#endif
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: f1_lo(3), f1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: f2_lo(3), f2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: f3_lo(3), f3_hi(3)
    integer, intent(in) :: q3_lo(3), q3_hi(3)
#endif
#if AMREX_SPACEDIM <= 2
    integer, intent(in) :: area1_lo(3), area1_hi(3)
#endif
#if AMREX_SPACEDIM == 2
    integer, intent(in) :: area2_lo(3), area2_hi(3)
#endif
#if AMREX_SPACEDIM <= 2
    integer, intent(in) :: vol_lo(3), vol_hi(3)
#endif

#ifdef RADIATION
    integer, intent(in) :: rf1_lo(3), rf1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: rf2_lo(3), rf2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: rf3_lo(3), rf3_hi(3)
#endif
#endif
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: dx(3), dt

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: shk(sk_lo(1):sk_hi(1), sk_lo(2):sk_hi(2), sk_lo(3):sk_hi(3))

    real(rt), intent(inout) :: qxm(qxm_lo(1):qxm_hi(1), qxm_lo(2):qxm_hi(2), qxm_lo(3):qxm_hi(3), NQ)
    real(rt), intent(inout) :: qxp(qxp_lo(1):qxp_hi(1), qxp_lo(2):qxp_hi(2), qxp_lo(3):qxp_hi(3), NQ)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: qym(qym_lo(1):qym_hi(1), qym_lo(2):qym_hi(2), qym_lo(3):qym_hi(3), NQ)
    real(rt), intent(inout) :: qyp(qyp_lo(1):qyp_hi(1), qyp_lo(2):qyp_hi(2), qyp_lo(3):qyp_hi(3), NQ)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: qzm(qzm_lo(1):qzm_hi(1), qzm_lo(2):qzm_hi(2), qzm_lo(3):qzm_hi(3), NQ)
    real(rt), intent(inout) :: qzp(qzp_lo(1):qzp_hi(1), qzp_lo(2):qzp_hi(2), qzp_lo(3):qzp_hi(3), NQ)
    real(rt), intent(inout) :: qmxy(qmxy_lo(1):qmxy_hi(1), qmxy_lo(2):qmxy_hi(2), qmxy_lo(3):qmxy_hi(3), NQ)
    real(rt), intent(inout) :: qpxy(qpxy_lo(1):qpxy_hi(1), qpxy_lo(2):qpxy_hi(2), qpxy_lo(3):qpxy_hi(3), NQ)
    real(rt), intent(inout) :: qmxz(qmxz_lo(1):qmxz_hi(1), qmxz_lo(2):qmxz_hi(2), qmxz_lo(3):qmxz_hi(3), NQ)
    real(rt), intent(inout) :: qpxz(qpxz_lo(1):qpxz_hi(1), qpxz_lo(2):qpxz_hi(2), qpxz_lo(3):qpxz_hi(3), NQ)
    real(rt), intent(inout) :: qmyx(qmyx_lo(1):qmyx_hi(1), qmyx_lo(2):qmyx_hi(2), qmyx_lo(3):qmyx_hi(3), NQ)
    real(rt), intent(inout) :: qpyx(qpyx_lo(1):qpyx_hi(1), qpyx_lo(2):qpyx_hi(2), qpyx_lo(3):qpyx_hi(3), NQ)
    real(rt), intent(inout) :: qmyz(qmyz_lo(1):qmyz_hi(1), qmyz_lo(2):qmyz_hi(2), qmyz_lo(3):qmyz_hi(3), NQ)
    real(rt), intent(inout) :: qpyz(qpyz_lo(1):qpyz_hi(1), qpyz_lo(2):qpyz_hi(2), qpyz_lo(3):qpyz_hi(3), NQ)
    real(rt), intent(inout) :: qmzx(qmzx_lo(1):qmzx_hi(1), qmzx_lo(2):qmzx_hi(2), qmzx_lo(3):qmzx_hi(3), NQ)
    real(rt), intent(inout) :: qpzx(qpzx_lo(1):qpzx_hi(1), qpzx_lo(2):qpzx_hi(2), qpzx_lo(3):qpzx_hi(3), NQ)
    real(rt), intent(inout) :: qmzy(qmzy_lo(1):qmzy_hi(1), qmzy_lo(2):qmzy_hi(2), qmzy_lo(3):qmzy_hi(3), NQ)
    real(rt), intent(inout) :: qpzy(qpzy_lo(1):qpzy_hi(1), qpzy_lo(2):qpzy_hi(2), qpzy_lo(3):qpzy_hi(3), NQ)
#endif
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout), target :: ftmp1(ft1_lo(1):ft1_hi(1), ft1_lo(2):ft1_hi(2), ft1_lo(3):ft1_hi(3), NVAR)
    real(rt), intent(inout), target :: ftmp2(ft2_lo(1):ft2_hi(1), ft2_lo(2):ft2_hi(2), ft2_lo(3):ft2_hi(3), NVAR)
#ifdef RADIATION
    real(rt), intent(inout), target :: rftmp1(rft1_lo(1):rft1_hi(1), rft1_lo(2):rft1_hi(2), rft1_lo(3):rft1_hi(3), 0:ngroups-1)
    real(rt), intent(inout), target :: rftmp2(rft2_lo(1):rft2_hi(1), rft2_lo(2):rft2_hi(2), rft2_lo(3):rft2_hi(3), 0:ngroups-1)
#endif
    real(rt), intent(inout), target :: qgdnvtmp1(qg1_lo(1):qg2_hi(1), qg1_lo(2):qg2_hi(2), qg1_lo(3):qg2_hi(3), NGDNV)
    real(rt), intent(inout), target :: qgdnvtmp2(qg2_lo(1):qg2_hi(1), qg2_lo(2):qg2_hi(2), qg2_lo(3):qg2_hi(3), NGDNV)
    real(rt), intent(inout), target :: ql(ql_lo(1):ql_hi(1), ql_lo(2):ql_hi(2), ql_lo(3):ql_hi(3), NQ)
    real(rt), intent(inout), target :: qr(qr_lo(1):qr_hi(1), qr_lo(2):qr_hi(2), qr_lo(3):qr_hi(3), NQ)
#endif
    real(rt), intent(inout) :: q_int(qi_lo(1):qi_hi(1), qi_lo(2):qi_hi(2), qi_lo(3):qi_hi(3), NQ)
#ifdef RADIATION
    real(rt), intent(inout) :: lambda_int(li_lo(1):li_hi(1), li_lo(2):li_hi(2), li_lo(3):li_hi(3), 0:ngroups-1)
#endif
    real(rt), intent(inout) ::  uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(inout) :: flux1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NVAR)
    real(rt), intent(inout) ::    q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: flux2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3),NVAR)
    real(rt), intent(inout) ::    q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: flux3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3),NVAR)
    real(rt), intent(inout) ::    q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
#endif

#ifdef RADIATION
    real(rt), intent(inout) :: rflux1(rf1_lo(1):rf1_hi(1),rf1_lo(2):rf1_hi(2),rf1_lo(3):rf1_hi(3),0:ngroups-1)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: rflux2(rf2_lo(1):rf2_hi(1),rf2_lo(2):rf2_hi(2),rf2_lo(3):rf2_hi(3),0:ngroups-1)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: rflux3(rf3_lo(1):rf3_hi(1),rf3_lo(2):rf3_hi(2),rf3_lo(3):rf3_hi(3),0:ngroups-1)
#endif
#endif
#if AMREX_SPACEDIM <= 2
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
#endif
#if AMREX_SPACEDIM == 2
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
#endif

    real(rt) :: dxinv, dyinv, dzinv
    real(rt) :: dtdx, dtdy, dtdz, hdt
#if AMREX_SPACEDIM == 3
    real(rt) :: cdtdx, cdtdy, cdtdz
#endif
    real(rt) :: hdtdx, hdtdy, hdtdz

    integer :: i, j, k, iwave, n


    integer :: fglo(3), fghi(3), glo(3), ghi(3)

    ! Left and right state arrays (edge centered, cell centered)
    double precision, dimension(:,:,:,:), pointer :: &
         qxl, qxr, qyl, qyr, qzl, qzr

    double precision, dimension(:,:,:,:), pointer:: &
         qgdnvx, qgdnvy, qgdnvz, &
         qgdnvxy, qgdnvxz, &
         qgdnvyx, qgdnvyz, &
         qgdnvzx, qgdnvzy

    double precision, dimension(:,:,:,:), pointer:: &
         fx, fy, fz, fxy, fxz, fyx, fyz, fzx, fzy

#ifdef RADIATION
    double precision, dimension(:,:,:,:), pointer:: &
         rfx, rfy, rfz, rfxy, rfxz, rfyx, rfyz, rfzx, rfzy
#endif

    fglo = lo - dg(:)  ! face + one ghost
    fghi = hi + 2*dg(:)

    glo = lo - dg(:)  ! one ghost,  this can be used for face-based arrays too
    ghi = hi + dg(:)

    ! Local constants
    dxinv = ONE/dx(1)
    dtdx = dt*dxinv
    hdtdx = HALF*dtdx

#if AMREX_SPACEDIM >= 2
    dyinv = ONE/dx(2)
    dtdy = dt*dyinv
    hdtdy = HALF*dtdy
#endif

#if AMREX_SPACEDIM == 3
    dzinv = ONE/dx(3)
    dtdz = dt*dzinv
    hdtdz = HALF*dtdz
    cdtdx = dtdx*THIRD
    cdtdy = dtdy*THIRD
    cdtdz = dtdz*THIRD
#endif

    hdt = HALF*dt


    !-------------------------------------------------------------------------!
    ! Some notes on the work index (i.e., lo and hi arguments near the end    !
    !                               of the argument list).                    !
    ! * For cmpflx, we use face index in the flux direction and cell-centered !
    !   index for others.                                                     !
    ! * For trans*, we use cell-centered index of the valid region.           !
    !-------------------------------------------------------------------------!

#if AMREX_SPACEDIM == 1
    !==========================================================================
    ! 1-d code path
    !==========================================================================
    ! Solve Riemann problem, compute xflux from improved predicted states
    call cmpflx(lo, hi+dg(:), &
                qxm, qxm_lo, qxm_hi, &
                qxp, qxp_lo, qxp_hi, 1, 1, &
                flux1, f1_lo, f1_hi, &
                q_int, fglo, fghi, &
#ifdef RADIATION
                rflux1, rf1_lo, rf1_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                1, domlo, domhi)

  call ca_store_godunov_state(lo, hi+dg(:), &
                              q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                              lambda_int, li_lo, li_hi, &
#endif
                              q1, q1_lo, q1_hi)

#endif

#if AMREX_SPACEDIM == 2
    !==========================================================================
    ! 2-d code path
    !==========================================================================

    fx     =>     ftmp1
#ifdef RADIATION
    rfx    =>    rftmp1
#endif
    qgdnvx  =>  qgdnvtmp1

    ! Compute F^x
    ! Inputs: qxm, qxp                     : xface, +-1 at y
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fx, ugdnvx, pgdnvx, gegdnvx : xface, +-1 at y
    call cmpflx([lo(1), lo(2)-1, 0], [hi(1)+1, hi(2)+1, 0], &
                qxm, qxm_lo, qxm_hi, &
                qxp, qxp_lo, qxp_hi, 1, 1, &
                fx, ft1_lo, ft1_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rfx, rft1_lo, rft1_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                1, domlo, domhi)

    call ca_store_godunov_state([lo(1), lo(2)-1, 0], [hi(1)+1, hi(2)+1, 0], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnvx, qg1_lo, qg1_hi)

    fy     =>     ftmp2
#ifdef RADIATION
    rfy    =>    rftmp2
#endif

    ! Compute F^y
    ! Inputs: qym, qyp                     : yface, +-1 at x
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fy, ugdnvy, pgdnvy, gegdnvy : yface, +-1 at x
    call cmpflx([lo(1)-1, lo(2), 0], [hi(1)+1, hi(2)+1, 0], &
                qym, qym_lo, qym_hi, &
                qyp, qyp_lo, qyp_hi, 1, 1, &
                fy, ft2_lo, ft2_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rfy, rft2_lo, rft2_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                2, domlo, domhi)

    call ca_store_godunov_state([lo(1)-1, lo(2), 0], [hi(1)+1, hi(2)+1, 0], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                q2, q2_lo, q2_hi)

    ! add the transverse flux difference in y to the x states
    ! Inputs: qxm, qxp                     : xface, +-1 at y
    !         fy, ugdnvy, pgdnvy, gegdnvy  : yface, +-1 at x
    !         gamc                         : +-4
    ! Outputs: qm, qp                      : xface, +-0 at y
    call transy([lo(1)-1, lo(2), 0], [hi(1)+1, hi(2), 0], &
                qxm, qxm_lo, qxm_hi, &
                ql, ql_lo, ql_hi, &
                qxp, qxp_lo, qxp_hi, &
                qr, qr_lo, qr_hi, &
                qaux, qa_lo, qa_hi, &
                fy, ft2_lo, ft2_hi, &
#ifdef RADIATION
                rfy, rft2_lo, rft2_hi, &
#endif
                q2, q2_lo, q2_hi, &
                hdtdy, &
                lo, hi)

    ! Solve the final Riemann problem across the x-interfaces with the
    ! full unsplit states.  The resulting flux through the x-interfaces
    ! is flux1
    call cmpflx([lo(1), lo(2), 0], [hi(1)+1, hi(2), 0], &
                ql, ql_lo, ql_hi, &
                qr, qr_lo, qr_hi, 1, 1, &
                flux1, f1_lo, f1_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rflux1, rf1_lo, rf1_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                1, domlo, domhi)

    call ca_store_godunov_state([lo(1), lo(2), 0], [hi(1)+1, hi(2), 0], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                q1, q1_lo, q1_hi)

    ! add the transverse flux difference in x to the y states
    ! Inputs: qym, qyp                     : yface, +-1 at x
    !         fx, ugdnvx, pgdnvx, gegdnvx  : xface, +-1 at y
    !         gamc                         : +-4
    ! Outputs: qm, qp                      : yface, +-0 at x
    call transx([lo(1), lo(2)-1, 0], [hi(1), hi(2)+1, 0], &
                qym, qym_lo, qym_hi, &
                ql, ql_lo, ql_hi, &
                qyp, qyp_lo. qyp_hi, &
                qr, qr_lo, qr_hi, &
                qaux, qa_lo, qa_hi, &
                fx, ft1_lo, ft1_hi, &
#ifdef RADIATION
                rf1, rft1_lo, rft1_hi, &
#endif
                qgdnvx, qg1_lo, qg1_hi, &
                area1, area1_lo, area1_hi, &
                vol, vol_lo, vol_hi, &
                hdt, hdtdx, &
                lo, hi)

    ! Solve the final Riemann problem across the y-interfaces with the
    ! full unsplit states.  The resulting flux through the y-interfaces
    ! is flux2
    call cmpflx([lo(1), lo(2), 0], [hi(1), hi(2)+1, 0], &
                ql, ql_lo, ql_hi, &
                qr, qr_lo, qr_hi, 1, 1, &
                flux2, f2_lo, f2_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rflux2, rf2_lo, rf2_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                2, domlo, domhi)

    call ca_store_godunov_state([lo(1), lo(2), 0], [hi(1), hi(2)+1, 0], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                q2, q2_lo, q2_hi)

    nullify(fx, fy, qgdnvx)
#ifdef RADIATION
    nullify(rfx, rfy)
#endif

#endif

#if AMREX_SPACEDIM == 3
    !==========================================================================
    ! 3-d code path
    !==========================================================================



    fx     =>     ftmp1
#ifdef RADIATION
    rfx    =>    rftmp1
#endif
    qgdnvx  =>  qgdnvtmp1

    ! Compute F^x
    ! Inputs: qxm, qxp                     : xface, +-1 at y & z
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fx, ugdnvx, pgdnvx, gegdnvx : xface, +-1 at y & z
    call cmpflx([lo(1), lo(2)-dg(2), lo(3)-dg(3)], [hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                qxm, qxm_lo, qxm_hi, &
                qxp, qxp_lo, qxp_hi, 1, 1, &
                fx, ft1_lo, ft1_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rfx, rft1_lo, rft1_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                1, domlo, domhi)

    call ca_store_godunov_state([lo(1), lo(2)-dg(2), lo(3)-dg(3)], [hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnvx, qg1_lo, qg1_hi)

    ! add the transverse flux difference in x to the y and z states
    ! Inputs: qym, qyp                     : yface, +-1 at x & z
    !         qzm, qzp                     : zface, +-1 at x & y
    !         fx, ugdnvx, pgdnvx, gegdnvx  : xface, +-1 at y & z
    !         gamc                         : +-4
    ! Outputs: qmyx, qpyx                  : yface, +-0 at x, +-1 at z
    !          qmzx, qpzx                  : zface, +-0 at x, +-1 at y
    call transx([lo(1), lo(2)-1, lo(3)-1], [hi(1), hi(2)+1, hi(3)+1], &
                qym, qym_lo, qym_hi, &
                qmyx, qmyx_lo, qmyx_hi, &
                qyp, qyp_lo, qyp_hi, &
                qpyx, qpyx_lo, qpyx_hi, &
                qzm, qzm_lo, qzm_hi, &
                qmzx, qmzx_lo, qmzx_hi, &
                qzp, qzp_lo, qzp_hi, &
                qpzx, qpzx_lo, qpzx_hi, &
                qaux, qa_lo, qa_hi, &
                fx, ft1_lo, ft1_hi, &
#ifdef RADIATION
                rfx, rft1_lo, rft1_hi, &
#endif
                qgdnvx, qg1_lo, qg1_hi, &
                hdt, cdtdx, &
                lo, hi)

    nullify(fx, qgdnvx)
#ifdef RADIATION
    nullify(rfx)
#endif

    fy     =>     ftmp1
#ifdef RADIATION
    rfy    =>    rftmp1
#endif
    qgdnvy => qgdnvtmp1

    ! Compute F^y
    ! Inputs: qym, qyp                     : yface, +-1 at x & z
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fy, ugdnvy, pgdnvy, gegdnvy : yface, +-1 at x & z
    call cmpflx([lo(1)-1, lo(2), lo(3)-dg(3)], [hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                qym, qym_lo, qym_hi, &
                qyp, qyp_lo, qyp_hi, 1, 1, &
                fy, ft1_lo, ft1_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rfy, rft1_lo, rft1_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                2, domlo, domhi)

    call ca_store_godunov_state([lo(1)-1, lo(2), lo(3)-dg(3)], [hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnvy, qg1_lo, qg1_hi)

    ! add the transverse flux difference in y to the x and z states
    ! Inputs: qxm, qxp                     : xface, +-1 at y & z
    !         qzm, qzp                     : zface, +-1 at x & y
    !         fy, ugdnvy, pgdnvy, gegdnvy  : yface, +-1 at x & z
    !         gamc                         : +-4
    ! Outputs: qmxy, qpxy                  : xface, +-0 at y, +-1 at z
    !          qmzy, qpzy                  : zface, +-0 at y, +-1 at x
    call transy([lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1], &
                qxm, qxm_lo, qxm_hi, &
                qmxy, qmxy_lo, qmxy_hi, &
                qxp, qxp_lo, qxp_hi, &
                qpxy, qpxy_lo, qpxy_hi, &
                qzm, qzm_lo, qzm_hi, &
                qmzy, qmzy_lo, qmzy_hi, &
                qzp, qzp_lo, qzp_hi, &
                qpzy, qpzy_lo, qpzy_hi, &
                qaux, qa_lo, qa_hi, &
                fy, ft1_lo, ft1_hi, &
#ifdef RADIATION
                rfy, rft1_lo, rft1_hi, &
#endif
                qgdnvy, qg1_lo, qg1_hi, &
                cdtdy, &
                lo, hi)

    nullify(fy, qgdnvy)
#ifdef RADIATION
    nullify(rfy)
#endif

    fz      =>     ftmp1
#ifdef RADIATION
    rfz     =>    rftmp1
#endif
    qgdnvz  =>  qgdnvtmp1

    ! Compute F^z
    ! Inputs: qzm, qzp                     : zface, +-1 at x & y
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fz, ugdnvz, pgdnvz, gegdnvz : zface, +-1 at x & y
    call cmpflx([lo(1)-1, lo(2)-dg(2), lo(3)], [ hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                qzm, qzm_lo, qzm_hi, &
                qzp, qzp_lo, qzp_hi, 1, 1, &
                fz, ft1_lo, ft1_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rfz, rft1_lo, rft1_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                3, domlo, domhi)

    call ca_store_godunov_state([lo(1)-1, lo(2)-dg(2), lo(3)], [ hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnvz, qg1_lo, qg1_hi)

    ! add the transverse flux difference in z to the x and y states
    ! Inputs: qxm, qxp                     : xface, +-1 at y & z
    !         qym, qyp                     : yface, +-1 at x & z
    !         fz, ugdnvz, pgdnvz, gegdnvz  : zface, +-1 at x & y
    !         gamc                         : +-4
    ! Outputs: qmxz, qpxz                  : xface, +-0 at z, +-1 at y
    !          qmyz, qpyz                  : yface, +-0 at z, +-1 at x
    call transz([lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)], &
                qxm, qxm_lo, qxm_hi, &
                qmxz, qmxz_lo, qmxz_hi, &
                qxp, qxp_lo, qxp_hi, &
                qpxz, qpxz_lo, qpxz_hi, &
                qym, qym_lo, qym_hi, &
                qmyz, qmyz_lo, qmyz_hi, &
                qyp, qyp_lo, qyp_hi, &
                qpyz, qpyz_lo, qpyz_hi, &
                qaux, qa_lo, qa_hi, &
                fz, ft1_lo, ft1_hi, &
#ifdef RADIATION
                rfz, rft1_lo, rft1_hi, &
#endif
                qgdnvz, qg1_lo, qg1_hi, &
                cdtdz, &
                lo, hi)

    nullify(fz, qgdnvz)
#ifdef RADIATION
    nullify(rfz)
#endif

    ! We now have qx?, qy?, qz?
    !         and q?zx, q?yx, q?zy, q?xy, q?yz, q?xz

    !
    ! Use qx?, q?yz, q?zy to compute final x-flux
    !

    fyz      =>      ftmp1
#ifdef RADIATION
    rfyz     =>     rftmp1
#endif
    qgdnvyz  =>  qgdnvtmp1

    ! Compute F^{y|z}
    ! Inputs: qmyz, qpyz                       : yface, +-1 at x, +-0 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fyz, ugdnvyz, pgdnvyz, gegdnvyz : yface, +-1 at x, +-0 at z
    call cmpflx([lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+dg(2), hi(3)], &
                qmyz, qmyz_lo, qmyz_hi, &
                qpyz, qpyz_lo, qpyz_hi, 1, 1, &
                fyz, ft1_lo, ft1_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rfyz, rft1_lo, rft1_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                2, domlo, domhi)

    call ca_store_godunov_state([lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+dg(2), hi(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnvyz, qg1_lo, qg1_hi)

    fzy      =>      ftmp2
#ifdef RADIATION
    rfzy     =>     rftmp2
#endif
    qgdnvzy  =>  qgdnvtmp2

    ! Compute F^{z|y}
    ! Inputs: qmzy, qpzy                       : zface, +-1 at x, +-0 at y
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fzy, ugdnvzy, pgdnvzy, gegdnvzy : zface, +-1 at x, +-0 at y
    call cmpflx([lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+dg(3)], &
                qmzy, qmzy_lo, qmzy_hi, &
                qpzy, qpzy_lo, qpzy_hi, 1, 1, &
                fzy, ft2_lo, ft2_hi, &
                q_int, qi_lo, q1_hi, &
#ifdef RADIATION
                rfzy, rft2_lo, rft2_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                3, domlo, domhi)

    call ca_store_godunov_state([lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+dg(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnvzy, qg2_lo, qg2_hi)

    qxl => ql
    qxr => qr

    ! Compute the corrected x interface states
    ! Inputs: qxm, qxp                        : xface, +-1 at y & z
    !         fyz, ugdnvyz, pgdnvyz, gegdnvyz : yface, +-1 at x, +-0 at z
    !         fzy, ugdnvzy, pgdnvzy, gegdnvzy : zface, +-1 at x, +-0 at y
    !         gamc, grav, rot                 : +-4
    !         srcQ                            : +-1
    ! Outputs: qxl, qxr                       : xface, +-0 at y & z
    call transyz([lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)], &
                 qxm, qxm_lo, qxm_hi, &
                 qxl, ql_lo, ql_hi, &
                 qxp, qxp_lo, qxp_hi, &
                 qxr, qr_lo, qr_hi, &
                 qaux, qa_lo, qa_hi, &
                 fyz, ft1_lo, ft1_hi, &
#ifdef RADIATION
                 rfyz, rft1_lo, rft1_hi, &
#endif
                 fzy, ft2_lo, ft2_hi, &
#ifdef RADIATION
                 rfzy, rft2_lo, rft2_hi, &
#endif
                 qgdnvyz, qg1_lo, qg1_hi, &
                 qgdnvzy, qg2_lo, qg2_hi, &
                 hdt, hdtdy, hdtdz, &
                 lo, hi)

    nullify(fyz, qgdnvyz)
    nullify(fzy, qgdnvzy)
#ifdef RADIATION
    nullify(rfyz, rfzy)
#endif

    qgdnvx  =>  qgdnvtmp1

    ! now compute the final x fluxes, F^x
    ! Inputs: qxl, qxr                        : xface, +-0 at y & z
    !         gamc, csml, c                   : +-4
    !         shk                             : +-1
    ! Outputs: flux1, ugdnvx, pgdnvx, gegdnvx : xface, +-0 at y & z
    call cmpflx(lo, [hi(1)+1, hi(2), hi(3)], &
                qxl, ql_lo, ql_hi, &
                qxr, qr_lo, qr_hi, 1, 1, &
                flux1, f1_lo, f1_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rflux1, rf1_lo, rf1_hi, &
                lambda_int, fglo, fghi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                1, domlo, domhi)

    call ca_store_godunov_state(lo, [hi(1)+1, hi(2), hi(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                q1, q1_lo, q1_hi)

    nullify(qgdnvx)
    nullify(qxl, qxr)

    !
    ! Use qy?, q?zx, q?xz to compute final y-flux
    !

    fzx      =>      ftmp1
#ifdef RADIATION
    rfzx     =>     rftmp1
#endif
    qgdnvzx  =>  qgdnvtmp1

    ! Compute F^{z|x}
    ! Inputs: qmzx, qpzx                       : zface, +-0 at x, +-1 at y
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fzx, ugdnvzx, pgdnvzx, gegdnvzx : zface, +-0 at x, +-1 at y
    call cmpflx([lo(1), lo(2)-dg(2), lo(3)], [hi(1), hi(2)+dg(2), hi(3)+dg(3)], &
                qmzx, qmzx_lo, qmzx_hi, &
                qpzx, qpzx_lo, qpzx_hi, 1, 1, &
                fzx, ft1_lo, ft1_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rfzx, rft1_lo, rft1_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                3, domlo, domhi)

    call ca_store_godunov_state([lo(1), lo(2)-dg(2), lo(3)], [hi(1), hi(2)+dg(2), hi(3)+dg(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnvzx, qg1_lo, qg1_hi)

    fxz      =>      ftmp2
#ifdef RADIATION
    rfxz     =>     rftmp2
#endif
    qgdnvxz  =>  qgdnvtmp2

    ! Compute F^{x|z} at km
    ! Inputs: qmxz, qpxz                       : xface, +-1 at y, +-0 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fxz, ugdnvxz, pgdnvxz, gegdnvxz : xface, +-1 at y, +-0 at z
    call cmpflx([lo(1), lo(2)-dg(2), lo(3)], [hi(1)+1, hi(2)+dg(2), hi(3)], &
                qmxz, qmxz_lo, qmxz_hi, &
                qpxz, qpxz_lo, qpxz_hi, 1, 1, &
                fxz, ft2_lo, ft2_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rfxz, rft2_lo, rft2_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                1, domlo, domhi)

    call ca_store_godunov_state([lo(1), lo(2)-dg(2), lo(3)], [hi(1)+1, hi(2)+dg(2), hi(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnvxz, qg2_lo, qg2_hi)

    qyl => ql
    qyr => qr

    ! Compute the corrected y interface states
    ! Inputs: qym, qyp                        : yface, +-1 at x & z
    !         fxz, ugdnvxz, pgdnvxz, gegdnvxz : xface, +-1 at y, +-0 at z
    !         fzx, ugdnvzx, pgdnvzx, gegdnvzx : zface, +-0 at x, +-1 at y
    !         gamc, grav, rot                 : +-4
    !         srcQ                            : +-1
    ! Outputs: qyl, qyr                       : yface, +-0 at x & z
    call transxz([lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)], &
                 qym, qym_lo, qym_hi, &
                 qyl, ql_lo, ql_hi, &
                 qyp, qyp_lo, qyp_hi, &
                 qyr, qr_lo, qr_hi, &
                 qaux, qa_lo, qa_hi, &
                 fxz, ft2_lo, ft2_hi, &
#ifdef RADIATION
                 rfxz, rft2_lo, rft2_hi, &
#endif
                 fzx, ft1_lo, ft1_hi, &
#ifdef RADIATION
                 rfzx, rft1_lo, rft1_hi, &
#endif
                 qgdnvxz, qg2_lo, qg2_hi, &
                 qgdnvzx, qg1_lo, qg1_hi, &
                 hdt, hdtdx, hdtdz, &
                 lo, hi)

    nullify(fzx, qgdnvzx)
    nullify(fxz, qgdnvxz)
#ifdef RADIATION
    nullify(rfzx, rfxz)
#endif

    qgdnvy  =>  qgdnvtmp1

    ! now compute the final y fluxes F^y
    ! Inputs: qyl, qyr                        : yface, +-0 at x & y
    !         gamc, csml, c                   : +-4
    !         shk                             : +-1
    ! Outputs: flux2, ugdnvy, pgdnvy, gegdnvy : yface, +-0 at x & y
    call cmpflx([lo(1), lo(2), lo(3)], [hi(1), hi(2)+dg(2), hi(3)], &
                qyl, ql_lo, ql_hi, &
                qyr, qr_lo, qr_hi, 1, 1, &
                flux2, f2_lo, f2_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rflux2, rf2_lo, rf2_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                2, domlo, domhi)

    call ca_store_godunov_state([lo(1), lo(2), lo(3)], [hi(1), hi(2)+dg(2), hi(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                q2, q2_lo, q2_hi)

    nullify(qgdnvy)
    nullify(qyl,qyr)

    !
    ! Use qz?, q?xy, q?yx to compute final z-flux
    !

    fxy      =>      ftmp1
#ifdef RADIATION
    rfxy     =>     rftmp1
#endif
    qgdnvxy  =>   qgdnvtmp1


    ! Compute F^{x|y}
    ! Inputs: qmxy, qpxy                       : xface, +-0 at y, +-1 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fxy, ugdnvxy, pgdnvxy, gegdnvxy : xface, +-0 at y, +-1 at z
    call cmpflx([lo(1), lo(2), lo(3)-dg(3)], [hi(1)+1, hi(2), hi(3)+dg(3)], &
                qmxy, qmxy_lo, qmxy_hi, &
                qpxy, qpxy_lo, qpxy_hi, 1, 1, &
                fxy, ft1_lo, ft1_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rfxy, rft1_lo, rft1_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                1, domlo, domhi)

    call ca_store_godunov_state([lo(1), lo(2), lo(3)-dg(3)], [hi(1)+1, hi(2), hi(3)+dg(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnvxy, qg1_lo, qg1_hi)


    fyx      =>      ftmp2
#ifdef RADIATION
    rfyx     =>     rftmp2
#endif
    qgdnvyx  =>  qgdnvtmp2

    ! Compute F^{y|x}
    ! Inputs: qmyx, qpyx                       : yface, +-0 at x, +-1 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fyx, ugdnvyx, pgdnvyx, gegdnvyx : yface, +-0 at x, +-1 at z
    call cmpflx([lo(1), lo(2), lo(3)-dg(3)], [hi(1), hi(2)+dg(2), hi(3)+dg(3)], &
                qmyx, qmyx_lo, qmyx_hi, &
                qpyx, qpyx_lo, qpyx_hi, 1, 1, &
                fyx, ft2_lo, ft2_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rfyx, rft2_lo, rft2_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                2, domlo, domhi)

    call ca_store_godunov_state([lo(1), lo(2), lo(3)-dg(3)], [hi(1), hi(2)+dg(2), hi(3)+dg(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                qgdnvyx, qg2_lo, qg2_hi)

    qzl => ql
    qzr => qr

    ! Compute the corrected z interface states
    ! Inputs: qzm, qzp                        : zface, +-1 at x & y
    !         fxy, ugdnvxy, pgdnvxy, gegdnvxy : xface, +-0 at y, +-1 at z
    !         fyx, ugdnvyx, pgdnvyx, gegdnvyx : yface, +-0 at x, +-1 at z
    !         gamc, grav, rot                 : +-4
    !         srcQ                            : +-1
    ! Outputs: qzl, qzr                       : zface, +-0 at x & y
    call transxy([lo(1), lo(2), lo(3)-1], [hi(1), hi(2), hi(3)+1], &
                 qzm, qzm_lo, qzm_hi, &
                 qzl, ql_lo, ql_hi, &
                 qzp, qzp_lo, qzp_hi, &
                 qzr, qr_lo, qr_hi, &
                 qaux, qa_lo, qa_hi, &
                 fxy, ft1_lo, ft1_hi, &
#ifdef RADIATION
                 rfxy, rft1_lo, rft1_hi, &
#endif
                 fyx, ft2_lo, ft2_hi, &
#ifdef RADIATION
                 rfyx, rft2_lo, rft2_hi, &
#endif
                 qgdnvxy, qg1_lo, qg1_hi, &
                 qgdnvyx, qg2_lo, qg2_hi, &
                 hdt, hdtdx, hdtdy, &
                 lo, hi)

    nullify(fxy, qgdnvxy)
    nullify(fyx, qgdnvyx)
#ifdef RADIATION
    nullify(rfxy, rfyx)
#endif

    qgdnvz  =>  qgdnvtmp1

    ! Compute the final z fluxes F^z
    ! Inputs: qzl, qzr                        : zface, +-0 at x & y
    !         gamc, csml, c                   : +-4
    !         shk                             : +-1
    ! Outputs: flux3, ugdnvz, pgdnvz, gegdnvz : zface, +-0 at x & y
    call cmpflx([lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+dg(3)], &
                qzl, ql_lo, ql_hi, &
                qzr, qr_lo, qr_hi, 1, 1, &
                flux3, f3_lo, f3_hi, &
                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                rflux3, rf3_lo, rf3_hi, &
                lambda_int, li_lo, li_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                3, domlo, domhi)

    call ca_store_godunov_state([lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+dg(3)], &
                                q_int, qi_lo, qi_hi, &
#ifdef RADIATION
                                lambda_int, li_lo, li_hi, &
#endif
                                q3, q3_lo, q3_hi)

    nullify(qgdnvz)
    nullify(qzl,qzr)
#endif


  end subroutine ctu_compute_fluxes


  subroutine consup(uin, uin_lo, uin_hi, &
                    q, q_lo, q_hi, &
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
                    vol,vol_lo,vol_hi, &
                    div, div_lo, div_hi, &
                    lo, hi, dx, dt, &
                    mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                    eden_lost, xang_lost, yang_lost, zang_lost, &
                    verbose)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, NGDNV, NQ, &
                                   GDPRES, &
#ifdef RADIATION
                                   fspace_type, comoving, &
                                   GDU, GDV, GDW, GDLAMS, GDERADS, &
#endif
                                   track_grid_losses, limit_fluxes_on_small_dens
    use advection_util_module, only : limit_hydro_fluxes_on_small_dens, normalize_species_fluxes, calc_pdivu, apply_av
    use castro_util_module, only : position, linear_to_angular_momentum
    use prob_params_module, only : mom_flux_has_p, domlo_level, domhi_level, center, dg, coord_type
    use amrinfo_module, only : amr_level
#ifdef RADIATION
    use rad_params_module, only : ngroups, nugroup, dlognu
    use radhydro_nd_module, only : advect_in_fspace
    use fluxlimiter_module, only : Edd_factor
    use advection_util_module, only : apply_av_rad
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
    integer, intent(in) ::   div_lo(3),   div_hi(3)
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

    integer, intent(in) :: verbose

    real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: shk(sk_lo(1):sk_hi(1),sk_lo(2):sk_hi(2),sk_lo(3):sk_hi(3))
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),updt_lo(3):updt_hi(3),NVAR)

    real(rt), intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt), intent(in) ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)

#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    real(rt), intent(in) ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
#endif

#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    real(rt), intent(in) ::    qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
#endif

    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(in) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in) :: dx(3), dt

#ifdef RADIATION
    real(rt)          Erin(Erin_lo(1):Erin_hi(1),Erin_lo(2):Erin_hi(2),Erin_lo(3):Erin_hi(3),0:ngroups-1)
    real(rt)         Erout(Erout_lo(1):Erout_hi(1),Erout_lo(2):Erout_hi(2),Erout_lo(3):Erout_hi(3),0:ngroups-1)
    real(rt)         radflux1(radflux1_lo(1):radflux1_hi(1),radflux1_lo(2):radflux1_hi(2),radflux1_lo(3):radflux1_hi(3),0:ngroups-1)
#if AMREX_SPACEDIM >= 2
    real(rt)         radflux2(radflux2_lo(1):radflux2_hi(1),radflux2_lo(2):radflux2_hi(2),radflux2_lo(3):radflux2_hi(3),0:ngroups-1)
#endif
#if AMREX_SPACEDIM == 3
    real(rt)         radflux3(radflux3_lo(1):radflux3_hi(1),radflux3_lo(2):radflux3_hi(2),radflux3_lo(3):radflux3_hi(3),0:ngroups-1)
#endif

#endif

    real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

    real(rt)         :: div1, volinv
    integer          :: i, j, g, k, n
    integer          :: domlo(3), domhi(3)
    real(rt)         :: loc(3), ang_mom(3)

#ifdef RADIATION
    real(rt)        , dimension(0:ngroups-1) :: Erscale
    real(rt)        , dimension(0:ngroups-1) :: ustar, af
    real(rt)         :: Eddf, Eddfxm, Eddfxp, Eddfym, Eddfyp, Eddfzm, Eddfzp
    real(rt)         :: f1, f2, f1xm, f1xp, f1ym, f1yp, f1zm, f1zp
    real(rt)         :: Gf1E(3)
    real(rt)         :: ux, uy, uz, divu, lamc, Egdc
    real(rt)         :: dudx(3), dudy(3), dudz(3), nhat(3), GnDotu(3), nnColonDotGu
    real(rt)         :: dprdx, dprdy, dprdz, ek1, ek2, dek, dpdx
    real(rt)         :: urho_new
    real(rt)         :: umx_new1, umy_new1, umz_new1
    real(rt)         :: umx_new2, umy_new2, umz_new2
#endif
    real(rt)        , pointer:: pdivu(:,:,:)

#ifdef RADIATION
    if (ngroups .gt. 1) then
       if (fspace_type .eq. 1) then
          Erscale = dlognu
       else
          Erscale = nugroup*dlognu
       end if
    end if
#endif

    call bl_allocate(pdivu, lo, hi)

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
                    dx, pdivu, lo, hi)

    ! zero out shock and temp fluxes -- these are physically meaningless here
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             flux1(i,j,k,UTEMP) = ZERO
#ifdef SHOCK_VAR
             flux1(i,j,k,USHK) = ZERO
#endif

#if AMREX_SPACEDIM >= 2
             flux2(i,j,k,UTEMP) = ZERO
#ifdef SHOCK_VAR
             flux2(i,j,k,USHK) = ZERO
#endif
#endif

#if AMREX_SPACEDIM == 3
             flux3(i,j,k,UTEMP) = ZERO
#ifdef SHOCK_VAR
             flux3(i,j,k,USHK) = ZERO
#endif
#endif

          end do
       end do
    end do

    call apply_av(lo, [hi(1)+1, hi(2), hi(3)], 1, dx, &
                  div, div_lo, div_hi, &
                  uin, uin_lo, uin_hi, &
                  flux1, flux1_lo, flux1_hi)

#if AMREX_SPACEDIM >= 2
    call apply_av(lo, [hi(1), hi(2)+1, hi(3)], 2, dx, &
                  div, div_lo, div_hi, &
                  uin, uin_lo, uin_hi, &
                  flux2, flux2_lo, flux2_hi)
#endif
#if AMREX_SPACEDIM == 3
    call apply_av(lo, [hi(1), hi(2), hi(3)+1], 3, dx, &
                  div, div_lo, div_hi, &
                  uin, uin_lo, uin_hi, &
                  flux3, flux3_lo, flux3_hi)
#endif

#ifdef RADIATION
   call apply_av_rad(lo, [hi(1)+1, hi(2), hi(3)], 1, dx, &
                      div, lo, hi+dg, &
                      Erin, Erin_lo, Erin_hi, &
                      radflux1, radflux1_lo, radflux1_hi)

#if AMREX_SPACEDIM >= 2
    call apply_av_rad(lo, [hi(1), hi(2)+1, hi(3)], 2, dx, &
                      div, lo, hi+dg, &
                      Erin, Erin_lo, Erin_hi, &
                      radflux2, radflux2_lo, radflux2_hi)
#endif

#if AMREX_SPACEDIM == 3
    call apply_av_rad(lo, [hi(1), hi(2), hi(3)+1], 3, dx, &
                      div, lo, hi+dg, &
                      Erin, Erin_lo, Erin_hi, &
                      radflux3, radflux3_lo, radflux3_hi)
#endif
#endif

    if (limit_fluxes_on_small_dens == 1) then
       call limit_hydro_fluxes_on_small_dens(uin,uin_lo,uin_hi, &
                                             q,q_lo,q_hi, &
                                             vol,vol_lo,vol_hi, &
                                             flux1,flux1_lo,flux1_hi, &
                                             area1,area1_lo,area1_hi, &
#if AMREX_SPACEDIM >= 2
                                             flux2,flux2_lo,flux2_hi, &
                                             area2,area2_lo,area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                                             flux3,flux3_lo,flux3_hi, &
                                             area3,area3_lo,area3_hi, &
#endif
                                             lo,hi,dt,dx)

    endif

    call normalize_species_fluxes(flux1_lo, flux1_hi, flux1, flux1_lo,flux1_hi)
#if AMREX_SPACEDIM >= 2
    call normalize_species_fluxes(flux2_lo, flux2_hi, flux2, flux2_lo,flux2_hi)
#endif
#if AMREX_SPACEDIM == 3
    call normalize_species_fluxes(flux3_lo, flux3_hi, flux3, flux3_lo, flux3_hi)
#endif

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

    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1) + 1
                flux1(i,j,k,n) = dt * flux1(i,j,k,n) * area1(i,j,k)
#if AMREX_SPACEDIM == 1
                ! Correct the momentum flux with the grad p part.
                if (coord_type .eq. 0 .and. n == UMX) then
                   flux1(i,j,k,n) = flux1(i,j,k,n) + dt * area1(i,j,k) * qx(i,j,k,GDPRES)
                endif
#endif
             enddo
          enddo
       enddo
    enddo

#if AMREX_SPACEDIM >= 2
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)
                flux2(i,j,k,n) = dt * flux2(i,j,k,n) * area2(i,j,k)
             enddo
          enddo
       enddo
    enddo
#endif

#if AMREX_SPACEDIM == 3
    do n = 1, NVAR
       do k = lo(3), hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux3(i,j,k,n) = dt * flux3(i,j,k,n) * area3(i,j,k)
             enddo
          enddo
       enddo
    enddo
#endif

#ifdef RADIATION
    do g = 0, ngroups-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1) + 1
                radflux1(i,j,k,g) = dt * radflux1(i,j,k,g) * area1(i,j,k)
             enddo
          enddo
       enddo
    enddo

#if AMREX_SPACEDIM >= 2
    do g = 0, ngroups-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)
                radflux2(i,j,k,g) = dt * radflux2(i,j,k,g) * area2(i,j,k)
             enddo
          enddo
       enddo
    enddo
#endif

#if AMREX_SPACEDIM == 3
    do g = 0, ngroups-1
       do k = lo(3), hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                radflux3(i,j,k,g) = dt * radflux3(i,j,k,g) * area3(i,j,k)
             enddo
          enddo
       enddo
    enddo
#endif

#endif


    ! Add up some diagnostic quantities. Note that we are not dividing by the cell volume.

    if (track_grid_losses .eq. 1) then

       domlo = domlo_level(:,amr_level)
       domhi = domhi_level(:,amr_level)

#if AMREX_SPACEDIM == 3
       if (lo(3) .le. domlo(3) .and. hi(3) .ge. domlo(3)) then

          k = domlo(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccz=.false.)

                mass_lost = mass_lost - flux3(i,j,k,URHO)
                xmom_lost = xmom_lost - flux3(i,j,k,UMX)
                ymom_lost = ymom_lost - flux3(i,j,k,UMY)
                zmom_lost = zmom_lost - flux3(i,j,k,UMZ)
                eden_lost = eden_lost - flux3(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux3(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(3) .le. domhi(3) .and. hi(3) .ge. domhi(3)) then

          k = domhi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccz=.false.)

                mass_lost = mass_lost + flux3(i,j,k,URHO)
                xmom_lost = xmom_lost + flux3(i,j,k,UMX)
                ymom_lost = ymom_lost + flux3(i,j,k,UMY)
                zmom_lost = zmom_lost + flux3(i,j,k,UMZ)
                eden_lost = eden_lost + flux3(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux3(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif
#endif

#if AMREX_SPACEDIM >= 2
       if (lo(2) .le. domlo(2) .and. hi(2) .ge. domlo(2)) then

          j = domlo(2)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccy=.false.)

                mass_lost = mass_lost - flux2(i,j,k,URHO)
                xmom_lost = xmom_lost - flux2(i,j,k,UMX)
                ymom_lost = ymom_lost - flux2(i,j,k,UMY)
                zmom_lost = zmom_lost - flux2(i,j,k,UMZ)
                eden_lost = eden_lost - flux2(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux2(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(2) .le. domhi(2) .and. hi(2) .ge. domhi(2)) then

          j = domhi(2) + 1
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccy=.false.)

                mass_lost = mass_lost + flux2(i,j,k,URHO)
                xmom_lost = xmom_lost + flux2(i,j,k,UMX)
                ymom_lost = ymom_lost + flux2(i,j,k,UMY)
                zmom_lost = zmom_lost + flux2(i,j,k,UMZ)
                eden_lost = eden_lost + flux2(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux2(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif
#endif

       if (lo(1) .le. domlo(1) .and. hi(1) .ge. domlo(1)) then

          i = domlo(1)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                loc = position(i,j,k,ccx=.false.)

                mass_lost = mass_lost - flux1(i,j,k,URHO)
                xmom_lost = xmom_lost - flux1(i,j,k,UMX)
                ymom_lost = ymom_lost - flux1(i,j,k,UMY)
                zmom_lost = zmom_lost - flux1(i,j,k,UMZ)
                eden_lost = eden_lost - flux1(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux1(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(1) .le. domhi(1) .and. hi(1) .ge. domhi(1)) then

          i = domhi(1) + 1
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                loc = position(i,j,k,ccx=.false.)

                mass_lost = mass_lost + flux1(i,j,k,URHO)
                xmom_lost = xmom_lost + flux1(i,j,k,UMX)
                ymom_lost = ymom_lost + flux1(i,j,k,UMY)
                zmom_lost = zmom_lost + flux1(i,j,k,UMZ)
                eden_lost = eden_lost + flux1(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux1(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif

    endif

    call bl_deallocate(pdivu)

  end subroutine consup


  subroutine ca_ctu_update(lo, hi, is_finest_level, time, &
                           domlo, domhi, &
                           uin, uin_lo, uin_hi, &
                           uout, uout_lo, uout_hi, &
#ifdef RADIATION
                           Erin, Erin_lo, Erin_hi, &
                           Erout, Erout_lo, Erout_hi, &
#endif
                           q, q_lo, q_hi, &
                           qaux, qa_lo, qa_hi, &
                           srcQ, srQ_lo, srQ_hi, &
                           update, updt_lo, updt_hi, &
                           dx, dt, &
                           flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
                           flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                           flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
                           radflux1, radflux1_lo, radflux1_hi, &
#if AMREX_SPACEDIM >= 2
                           radflux2, radflux2_lo, radflux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                           radflux3, radflux3_lo, radflux3_hi, &
#endif
#endif
                           area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                           area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                           area3, area3_lo, area3_hi, &
#endif
#if AMREX_SPACEDIM <= 2
                           pradial, p_lo, p_hi, &
                           dloga, dloga_lo, dloga_hi, &
#endif
                           vol, vol_lo, vol_hi, &
                           verbose, &
#ifdef RADIATION
                           nstep_fsp, &
#endif
                           mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                           eden_lost, xang_lost, yang_lost, zang_lost) bind(C, name="ca_ctu_update")

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : NQ, QVAR, QPRES, NQAUX, NVAR, NHYP, NGDNV, UMX, GDPRES, &
#ifdef RADIATION
                                   QPTOT, &
#endif
                                   use_flattening, first_order_hydro
    use advection_util_module, only : divu
    use amrex_constants_module, only : ZERO, ONE
    use flatten_module, only: ca_uflatten
    use prob_params_module, only : mom_flux_has_p, dg, coord_type
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use flatten_module, only : rad_flatten
#endif

    implicit none

#ifdef RADIATION
    integer, intent(inout) :: nstep_fsp
#endif
    integer, intent(in) :: is_finest_level
    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
#ifdef RADIATION
    integer, intent(in) :: Erin_lo(3), Erin_hi(3)
    integer, intent(in) :: Erout_lo(3), Erout_hi(3)
#endif
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: srQ_lo(3), srQ_hi(3)
    integer, intent(in) :: updt_lo(3), updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
#endif
#ifdef RADIATION
    integer, intent(in) :: radflux1_lo(3), radflux1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: radflux2_lo(3), radflux2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: radflux3_lo(3), radflux3_hi(3)
#endif
#endif
    integer, intent(in) :: area1_lo(3), area1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: area2_lo(3), area2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: area3_lo(3), area3_hi(3)
#endif
    integer, intent(in) :: vol_lo(3), vol_hi(3)
#if AMREX_SPACEDIM <= 2
    integer, intent(in) :: p_lo(3), p_hi(3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif

    real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
#ifdef RADIATION
    real(rt)        , intent(in) :: Erin(Erin_lo(1):Erin_hi(1), Erin_lo(2):Erin_hi(2), Erin_lo(3):Erin_hi(3), 0:ngroups-1)
    real(rt)        , intent(inout) :: Erout(Erout_lo(1):Erout_hi(1), Erout_lo(2):Erout_hi(2), Erout_lo(3):Erout_hi(3), 0:ngroups-1)
#endif
    real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt)        , intent(in) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt)        , intent(in) :: srcQ(srQ_lo(1):srQ_hi(1), srQ_lo(2):srQ_hi(2), srQ_lo(3):srQ_hi(3), QVAR)
    real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
    real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
#if AMREX_SPACEDIM >= 2
    real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
#endif
#if AMREX_SPACEDIM == 3
    real(rt)        , intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)
#endif
#ifdef RADIATION
    real(rt)        , intent(inout) :: radflux1(radflux1_lo(1):radflux1_hi(1), radflux1_lo(2):radflux1_hi(2), &
         radflux1_lo(3):radflux1_hi(3), 0:ngroups-1)
#if AMREX_SPACEDIM >= 2
    real(rt)        , intent(inout) :: radflux2(radflux2_lo(1):radflux2_hi(1), radflux2_lo(2):radflux2_hi(2), &
         radflux2_lo(3):radflux2_hi(3), 0:ngroups-1)
#endif
#if AMREX_SPACEDIM == 3
    real(rt)        , intent(inout) :: radflux3(radflux3_lo(1):radflux3_hi(1), radflux3_lo(2):radflux3_hi(2), &
         radflux3_lo(3):radflux3_hi(3), 0:ngroups-1)
#endif
#endif
    real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
#if AMREX_SPACEDIM >= 2
    real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
#endif
#if AMREX_SPACEDIM == 3
    real(rt)        , intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
#endif
    real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))

#if AMREX_SPACEDIM < 3
    real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
    real(rt)        , intent(inout) :: pradial(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
#endif

    real(rt)        , intent(in) :: dx(3), dt, time

    real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

    ! Automatic arrays for workspace
    real(rt)        , pointer:: flatn(:,:,:)
#ifdef RADIATION
    real(rt)        , pointer:: flatg(:,:,:)
#endif
    real(rt)        , pointer:: div(:,:,:)

    ! Edge-centered primitive variables (Riemann state)
    real(rt)        , pointer:: q1(:,:,:,:)
    real(rt)        , pointer:: q2(:,:,:,:)
    real(rt)        , pointer:: q3(:,:,:,:)

    integer :: q1_lo(3), q1_hi(3), q2_lo(3), q2_hi(3), q3_lo(3), q3_hi(3)

    integer :: fglo(3), fghi(3), glo(3), ghi(3)

    ! for computing the normal states
    real(rt), pointer :: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
    real(rt), pointer :: Ip_src(:,:,:,:,:,:), Im_src(:,:,:,:,:,:)
    real(rt), pointer :: Ip_gc(:,:,:,:,:,:), Im_gc(:,:,:,:,:,:)

    real(rt), pointer :: shk(:,:,:)

    real(rt), pointer :: sm(:,:,:,:), sp(:,:,:,:)
    real(rt), pointer :: dq(:,:,:,:,:)

    ! Left and right state arrays (edge centered, cell centered)
    double precision, dimension(:,:,:,:), pointer :: &
         qxm, qym, qzm, qxp, qyp, qzp

    ! for the transverse corrections
    real(rt), pointer :: q_int(:,:,:,:)

#ifdef RADIATION
    real(rt), pointer :: lambda_int(:,:,:,:)
#endif


    ! Left and right state arrays (edge centered, cell centered)
    double precision, dimension(:,:,:,:), pointer :: &
         ql, qr, &
         qmxy, qpxy, qmxz, qpxz, qmyx, qpyx, &
         qmyz, qpyz, qmzx, qpzx, qmzy, qpzy


    ! these will be the temporary arrays we actually allocate space for
    double precision, dimension(:,:,:,:), pointer :: ftmp1, ftmp2, rftmp1, rftmp2
    double precision, dimension(:,:,:,:), pointer :: qgdnvtmp1, qgdnvtmp2

    call bl_allocate(   div, lo, hi+dg)

    q1_lo = flux1_lo - dg
    q1_hi = flux1_hi + dg
#if AMREX_SPACEDIM >= 2
    q2_lo = flux2_lo - dg
    q2_hi = flux2_hi + dg
#endif
#if AMREX_SPACEDIM == 3
    q3_lo = flux3_lo - dg
    q3_hi = flux3_hi + dg
#endif

    call bl_allocate(q1, q1_lo, q1_hi, NGDNV)
#if AMREX_SPACEDIM >= 2
    call bl_allocate(q2, q2_lo, q2_hi, NGDNV)
#endif
#if AMREX_SPACEDIM == 3
    call bl_allocate(q3, q3_lo, q3_hi, NGDNV)
#endif

    ! Compute flattening coefficient for slope calculations.
    call bl_allocate( flatn, q_lo, q_hi)
#ifdef RADIATION
    call bl_allocate( flatg, q_lo, q_hi)
#endif
    if (first_order_hydro == 1) then
       flatn = ZERO
    elseif (use_flattening == 1) then
       call ca_uflatten(lo-dg, hi+dg, &
            q, q_lo, q_hi, &
            flatn, q_lo, q_hi, QPRES)
#ifdef RADIATION
       call ca_uflatten(lo-dg, hi+dg, &
            q, q_lo, q_hi, &
            flatg, q_lo, q_hi, QPTOT)

       call rad_flatten(lo-dg, hi+dg, &
            q, q_lo, q_hi, &
            flatn, q_lo, q_hi, &
            flatg, q_lo, q_hi)
#endif
    else
       flatn = ONE
    endif
#ifdef RADIATION
    call bl_deallocate(flatg)
#endif

    fglo = lo - dg(:)  ! face + one ghost
    fghi = hi + 2*dg(:)

    glo = lo - dg(:)  ! one ghost,  this can be used for face-based arrays too
    ghi = hi + dg(:)

    call bl_allocate ( qxm, fglo, fghi, NQ)
    call bl_allocate ( qxp, fglo, fghi, NQ)

#if AMREX_SPACEDIM >= 2
    call bl_allocate ( qym, fglo, fghi, NQ)
    call bl_allocate ( qyp, fglo, fghi, NQ)
#endif

#if AMREX_SPACEDIM == 3
    call bl_allocate ( qzm, fglo, fghi, NQ)
    call bl_allocate ( qzp, fglo, fghi, NQ)
#endif

    ! x-index, y-index, z-index, dim, characteristics, variables
    call bl_allocate ( Ip, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,AMREX_SPACEDIM,1,3,1,NQ)
    call bl_allocate ( Im, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,AMREX_SPACEDIM,1,3,1,NQ)

    ! for source terms
    call bl_allocate ( Ip_src, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,AMREX_SPACEDIM,1,3,1,QVAR)
    call bl_allocate ( Im_src, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,AMREX_SPACEDIM,1,3,1,QVAR)

    ! for gamc -- needed for the reference state in eigenvectors
    call bl_allocate ( Ip_gc, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,AMREX_SPACEDIM,1,3,1,1)
    call bl_allocate ( Im_gc, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,AMREX_SPACEDIM,1,3,1,1)

    ! only needed for PLM -- maybe we can piggyback on the others?
    call bl_allocate ( dq, glo(1),ghi(1), glo(2),ghi(2), glo(3),ghi(3), 1, NQ, 1, AMREX_SPACEDIM)

    ! for the hybrid Riemann solver
    call bl_allocate(shk, glo, ghi)

    call bl_allocate(sm, glo, ghi, AMREX_SPACEDIM)
    call bl_allocate(sp, glo, ghi, AMREX_SPACEDIM)

    call ctu_normal_states(lo-dg, hi+dg, &
                           lo, hi, &
                           q, q_lo, q_hi, &
                           flatn, q_lo, q_hi, &
                           qaux, qa_lo, qa_hi, &
                           srcQ, srQ_lo, srQ_hi, &
                           shk, glo, ghi, &
                           Ip, glo, ghi, &
                           Im, glo, ghi, &
                           Ip_src, glo, ghi, &
                           Im_src, glo, ghi, &
                           Ip_gc, glo, ghi, &
                           Im_gc, glo, ghi, &
                           dq, glo, ghi, &
                           sm, glo, ghi, &
                           sp, glo, ghi, &
                           qxm, fglo, fghi, &
                           qxp, fglo, fghi, &
#if AMREX_SPACEDIM >= 2
                           qym, fglo, fghi, &
                           qyp, fglo, fghi, &
#endif
#if AMREX_SPACEDIM == 3
                           qzm, fglo, fghi, &
                           qzp, fglo, fghi, &
#endif
                           dx, dt, &
#if AMREX_SPACEDIM < 3
                           dloga, dloga_lo, dloga_hi, &
#endif
                           domlo, domhi)

    call bl_deallocate ( Ip)
    call bl_deallocate ( Im)

    call bl_deallocate ( Ip_src)
    call bl_deallocate ( Im_src)

    call bl_deallocate ( Ip_gc)
    call bl_deallocate ( Im_gc)

    call bl_deallocate ( dq)

    call bl_deallocate(sm)
    call bl_deallocate(sp)

    call bl_deallocate( flatn)

#if AMREX_SPACEDIM == 3
    call bl_allocate( qmxy, fglo, fghi, NQ)
    call bl_allocate( qpxy, fglo, fghi, NQ)

    call bl_allocate( qmxz, fglo, fghi, NQ)
    call bl_allocate( qpxz, fglo, fghi, NQ)

    call bl_allocate( qmyx, fglo, fghi, NQ)
    call bl_allocate( qpyx, fglo, fghi, NQ)

    call bl_allocate( qmyz, fglo, fghi, NQ)
    call bl_allocate( qpyz, fglo, fghi, NQ)

    call bl_allocate( qmzx, fglo, fghi, NQ)
    call bl_allocate( qpzx, fglo, fghi, NQ)

    call bl_allocate( qmzy, fglo, fghi, NQ)
    call bl_allocate( qpzy, fglo, fghi, NQ)
#endif

#if AMREX_SPACEDIM >= 2
    call bl_allocate( ftmp1, glo, ghi, NVAR)
    call bl_allocate( ftmp2, glo, ghi, NVAR)
#ifdef RADIATION
    call bl_allocate ( rftmp1, glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3), 0, ngroups-1)
    call bl_allocate ( rftmp2, glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3), 0, ngroups-1)
#endif
    call bl_allocate ( qgdnvtmp1, fglo, fghi, NGDNV)
    call bl_allocate ( qgdnvtmp2, fglo, fghi, NGDNV)

    call bl_allocate( ql, fglo, fghi, NQ)
    call bl_allocate( qr, fglo, fghi, NQ)
#endif

    call bl_allocate(q_int, fglo, fghi, NQ)
#ifdef RADIATION
    call bl_allocate(lambda_int, fglo(1), fghi(1), fglo(2), fghi(2), fglo(3), fghi(3), 0, ngroups-1)
#endif

    ! Compute hyperbolic fluxes using unsplit Godunov
    call ctu_compute_fluxes(lo, hi, &
                            q, q_lo, q_hi, &
                            qaux, qa_lo, qa_hi, &
                            dx, dt, &
                            shk, glo, ghi, &
                            qxm, fglo, fghi, &
                            qxp, fglo, fghi, &
#if AMREX_SPACEDIM >= 2
                            qym, fglo, fghi, &
                            qyp, fglo, fghi, &
#endif
#if AMREX_SPACEDIM == 3
                            qzm, fglo, fghi, &
                            qzp, fglo, fghi, &
                            qmxy, fglo, fghi, &
                            qpxy, fglo, fghi, &
                            qmxz, fglo, fghi, &
                            qpxz, fglo, fghi, &
                            qmyx, fglo, fghi, &
                            qpyx, fglo, fghi, &
                            qmyz, fglo, fghi, &
                            qpyz, fglo, fghi, &
                            qmzx, fglo, fghi, &
                            qpzx, fglo, fghi, &
                            qmzy, fglo, fghi, &
                            qpzy, fglo, fghi, &
#endif
#if AMREX_SPACEDIM >= 2
                            ftmp1, glo, ghi, &
                            ftmp2, glo, ghi, &
#ifdef RADIATION
                            rftmp1, glo, ghi, &
                            rftmp2, glo, ghi, &
#endif
                            qgdnvtmp1, fglo, fghi, &
                            qgdnvtmp2, fglo, fghi, &
                            ql, fglo, fghi, &
                            qr, fglo, fghi, &
#endif
                            q_int, fglo, fghi, &
#ifdef RADIATION
                            lambda_int, fglo, fghi, &
#endif
                            uout, uout_lo, uout_hi, &
                            flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
                            flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                            flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
                            radflux1, radflux1_lo, radflux1_hi, &
#if AMREX_SPACEDIM >= 2
                            radflux2, radflux2_lo, radflux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                            radflux3, radflux3_lo, radflux3_hi, &
#endif
#endif
                            q1, q1_lo, q1_hi, &
#if AMREX_SPACEDIM >= 2
                            q2, q2_lo, q2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                            q3, q3_lo, q3_hi, &
#endif
#if AMREX_SPACEDIM < 3
                            area1, area1_lo, area1_hi, &
#endif
#if AMREX_SPACEDIM == 2
                            area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM < 3
                            vol, vol_lo, vol_hi, &
#endif
                            domlo, domhi)

    call bl_deallocate(qxm)
    call bl_deallocate(qxp)
#if AMREX_SPACEDIM >= 2
    call bl_deallocate(qym)
    call bl_deallocate(qyp)
#endif
#if AMREX_SPACEDIM == 3
    call bl_deallocate(qzm)
    call bl_deallocate(qzp)
#endif

    call bl_deallocate(q_int)
#ifdef RADIATION
    call bl_deallocate(lambda_int)
#endif

#if AMREX_SPACEDIM >= 2
    call bl_deallocate ( ql)
    call bl_deallocate ( qr)
#endif

#if AMREX_SPACEDIM == 3
    call bl_deallocate ( qmxy)
    call bl_deallocate ( qpxy)

    call bl_deallocate ( qmxz)
    call bl_deallocate ( qpxz)

    call bl_deallocate ( qmzx)
    call bl_deallocate ( qpzx)

    call bl_deallocate ( qmzy)
    call bl_deallocate ( qpzy)

    call bl_deallocate ( qmyx)
    call bl_deallocate ( qpyx)

    call bl_deallocate ( qmyz)
    call bl_deallocate ( qpyz)
#endif

#if AMREX_SPACEDIM >= 2
    call bl_deallocate ( ftmp1)
    call bl_deallocate ( ftmp2)

#ifdef RADIATION
    call bl_deallocate (rftmp1)
    call bl_deallocate (rftmp2)
#endif

    call bl_deallocate(qgdnvtmp1)
    call bl_deallocate(qgdnvtmp2)
#endif

    ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
    call divu(lo, hi+dg, q, q_lo, q_hi, dx, div, lo, hi+dg)

    ! Conservative update
    ! TODO: store the shock variable in uout

    call consup(uin, uin_lo, uin_hi, &
                q, q_lo, q_hi, &
                shk, glo, ghi, &
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
                q1, q1_lo, q1_hi, &
#if AMREX_SPACEDIM >= 2
                q2, q2_lo, q2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                q3, q3_lo, q3_hi, &
#endif
                area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                area3, area3_lo, area3_hi, &
#endif
                vol, vol_lo, vol_hi, &
                div, lo, hi+dg, &
                lo, hi, dx, dt, &
                mass_lost,xmom_lost,ymom_lost,zmom_lost, &
                eden_lost,xang_lost,yang_lost,zang_lost, &
                verbose)


    call bl_deallocate(shk)



#if AMREX_SPACEDIM == 1
    if (coord_type > 0) then
       pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = q1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
    end if
#endif
#if AMREX_SPACEDIM == 2
    if (.not. mom_flux_has_p(1)%comp(UMX)) then
       pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = q1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
    end if
#endif

    call bl_deallocate(   div)

    call bl_deallocate(    q1)
#if AMREX_SPACEDIM >= 2
    call bl_deallocate(    q2)
#endif
#if AMREX_SPACEDIM == 3
    call bl_deallocate(    q3)
#endif

  end subroutine ca_ctu_update

end module ctu_module
 
