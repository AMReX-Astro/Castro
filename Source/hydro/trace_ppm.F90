module trace_ppm_module
  !
  ! These routines do the characteristic tracing under the parabolic
  ! profiles in each zone to the edge / half-time.

  use prob_params_module, only : dg
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module, only : ZERO, HALF, ONE

  implicit none

contains

  subroutine trace_ppm(lo, hi, &
                       idir, &
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
                       qm_core, qcm_lo, qcm_hi, &
                       qp_core, qcp_lo, qcp_hi, &
                       qm_pass, qpm_lo, qpm_hi, &
                       qp_pass, qpp_lo, qpp_hi, &
#if (AMREX_SPACEDIM < 3)
                       dloga, dloga_lo, dloga_hi, &
#endif
                       vlo, vhi, domlo, domhi, &
                       dx, dt)
    ! here, lo and hi are the range we loop over -- this can include ghost cells
    ! vlo and vhi are the bounds of the valid box (no ghost cells)

    use network, only : nspec, naux
    use meth_params_module, only : NQC, NQAUX, NQC_SRC, NQP, ppm_predict_gammae, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                                   NQP_SRC, &
#endif
#ifdef RADIATION
                                   NQR, &
#endif
                                   ppm_temp_fix, QU, QV, QW, npassive, qpass_map
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qc_lo(3), qc_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)

    integer, intent(in) :: qa_lo(3), qa_hi(3)

    integer, intent(in) :: Icp_lo(3), Icp_hi(3)
    integer, intent(in) :: Icm_lo(3), Icm_hi(3)
    integer, intent(in) :: Ipp_lo(3), Ipp_hi(3)
    integer, intent(in) :: Ipm_lo(3), Ipm_hi(3)

    integer, intent(in) :: Icsp_lo(3), Icsp_hi(3)
    integer, intent(in) :: Icsm_lo(3), Icsm_hi(3)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    integer, intent(in) :: Ipsp_lo(3), Ipsp_hi(3)
    integer, intent(in) :: Ipsm_lo(3), Ipsm_hi(3)
#endif

    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)

    integer, intent(in) :: qcm_lo(3), qcm_hi(3)
    integer, intent(in) :: qcp_lo(3), qcp_hi(3)
    integer, intent(in) :: qpm_lo(3), qpm_hi(3)
    integer, intent(in) :: qpp_lo(3), qpp_hi(3)

#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q_core(qc_lo(1):qc_hi(1),qc_lo(2):qc_hi(2),qc_lo(3):qc_hi(3),NQC)
    real(rt), intent(in) :: q_pass(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQP)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip_core(Icp_lo(1):Icp_hi(1),Icp_lo(2):Icp_hi(2),Icp_lo(3):Icp_hi(3),1:3,NQC)
    real(rt), intent(in) :: Im_core(Icm_lo(1):Icm_hi(1),Icm_lo(2):Icm_hi(2),Icm_lo(3):Icm_hi(3),1:3,NQC)
    real(rt), intent(in) :: Ip_pass(Ipp_lo(1):Ipp_hi(1),Ipp_lo(2):Ipp_hi(2),Ipp_lo(3):Ipp_hi(3),1:3,NQP)
    real(rt), intent(in) :: Im_pass(Ipm_lo(1):Ipm_hi(1),Ipm_lo(2):Ipm_hi(2),Ipm_lo(3):Ipm_hi(3),1:3,NQP)

    real(rt), intent(in) :: Ip_core_src(Icsp_lo(1):Icsp_hi(1),Icsp_lo(2):Icsp_hi(2),Icsp_lo(3):Icsp_hi(3),1:3,NQC_SRC)
    real(rt), intent(in) :: Im_core_src(Icsm_lo(1):Icsm_hi(1),Icsm_lo(2):Icsm_hi(2),Icsm_lo(3):Icsm_hi(3),1:3,NQC_SRC)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    real(rt), intent(in) :: Ip_pass_src(Ipsp_lo(1):Ipsp_hi(1),Ipsp_lo(2):Ipsp_hi(2),Ipsp_lo(3):Ipsp_hi(3),1:3,NQP_SRC)
    real(rt), intent(in) :: Im_pass_src(Ipsm_lo(1):Ipsm_hi(1),Ipsm_lo(2):Ipsm_hi(2),Ipsm_lo(3):Ipsm_hi(3),1:3,NQP_SRC)
#endif

    real(rt), intent(in) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:3,1)
    real(rt), intent(in) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:3,1)

    real(rt), intent(inout) :: qm_core(qcm_lo(1):qcm_hi(1),qcm_lo(2):qcm_hi(2),qcm_lo(3):qcm_hi(3),NQC)
    real(rt), intent(inout) :: qp_core(qcp_lo(1):qcp_hi(1),qcp_lo(2):qcp_hi(2),qcp_lo(3):qcp_hi(3),NQC)
    real(rt), intent(inout) :: qm_pass(qpm_lo(1):qpm_hi(1),qpm_lo(2):qpm_hi(2),qpm_lo(3):qpm_hi(3),NQP)
    real(rt), intent(inout) :: qp_pass(qpp_lo(1):qpp_hi(1),qpp_lo(2):qpp_hi(2),qpp_lo(3):qpp_hi(3),NQP)

#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)


    real(rt) :: un
    integer :: ipassive, n, i, j, k

    !$gpu

    ! the passive stuff is the same regardless of the tracing
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Plus state on face i
                if ((idir == 1 .and. i >= vlo(1)) .or. &
                    (idir == 2 .and. j >= vlo(2)) .or. &
                    (idir == 3 .and. k >= vlo(3))) then

                   un = q_core(i,j,k,QU-1+idir)

                   ! We have
                   !
                   ! q_l = q_ref - Proj{(q_ref - I)}
                   !
                   ! and Proj{} represents the characteristic projection.
                   ! But for these, there is only 1-wave that matters, the u
                   ! wave, so no projection is needed.  Since we are not
                   ! projecting, the reference state doesn't matter

                   if (un > ZERO) then
                      qp_pass(i,j,k,n) = q_pass(i,j,k,n)
                   else
                      qp_pass(i,j,k,n) = Im_pass(i,j,k,2,n)
                   end if
#ifdef PRIM_SPECIES_HAVE_SOURCES
                   qp_pass(i,j,k,n) = qp_pass(i,j,k,n) + HALF*dt*Im_pass_src(i,j,k,2,n)
#endif

                end if

                ! Minus state on face i+1
                if (idir == 1 .and. i <= vhi(1)) then
                   un = q_core(i,j,k,QU-1+idir)
                   if (un > ZERO) then
                      qm_pass(i+1,j,k,n) = Ip_pass(i,j,k,2,n)
                   else
                      qm_pass(i+1,j,k,n) = q_pass(i,j,k,n)
                   end if
#ifdef PRIM_SPECIES_HAVE_SOURCES
                   qm_pass(i+1,j,k,n) = qm_pass(i+1,j,k,n) + HALF*dt*Ip_pass_src(i,j,k,2,n)
#endif

                else if (idir == 2 .and. j <= vhi(2)) then
                   un = q_core(i,j,k,QU-1+idir)
                   if (un > ZERO) then
                      qm_pass(i,j+1,k,n) = Ip_pass(i,j,k,2,n)
                   else
                      qm_pass(i,j+1,k,n) = q_pass(i,j,k,n)
                   end if
#ifdef PRIM_SPECIES_HAVE_SOURCES
                   qm_pass(i,j+1,k,n) = qm_pass(i,j+1,k,n) + HALF*dt*Ip_pass_src(i,j,k,2,n)
#endif

                else if (idir == 3 .and. k <= vhi(3)) then
                   un = q_core(i,j,k,QU-1+idir)
                   if (un > ZERO) then
                      qm_pass(i,j,k+1,n) = Ip_pass(i,j,k,2,n)
                   else
                      qm_pass(i,j,k+1,n) = q_pass(i,j,k,n)
                   end if
#ifdef PRIM_SPECIES_HAVE_SOURCES
                   qm_pass(i,j,k+1,n) = qm_pass(i,j,k+1,n) + HALF*dt*Ip_pass_src(i,j,k,2,n)
#endif
                end if

             end do

          end do
       end do

    end do

    if (ppm_temp_fix < 3) then
       if (ppm_predict_gammae == 0) then
          call trace_ppm_rhoe(lo, hi, &
                              idir, &
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
                              qm_core, qcm_lo, qcm_hi, &
                              qp_core, qcp_lo, qcp_hi, &
#if (AMREX_SPACEDIM < 3)
                              dloga, dloga_lo, dloga_hi, &
#endif
                              vlo, vhi, domlo, domhi, &
                              dx, dt)
       else
          call trace_ppm_gammae(lo, hi, &
                                idir, &
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
                                qm_core, qcm_lo, qcm_hi, &
                                qp_core, qcp_lo, qcp_hi, &
#if (AMREX_SPACEDIM < 3)
                                dloga, dloga_lo, dloga_hi, &
#endif
                                vlo, vhi, domlo, domhi, &
                                dx, dt)
       end if
    else
       call trace_ppm_temp(lo, hi, &
                           idir, &
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
                           qm_core, qcm_lo, qcm_hi, &
                           qp_core, qcp_lo, qcp_hi, &
                           qm_pass, qpm_lo, qpm_hi, &
                           qp_pass, qpp_lo, qpp_hi, &
#if (AMREX_SPACEDIM < 3)
                           dloga, dloga_lo, dloga_hi, &
#endif
                           vlo, vhi, domlo, domhi, &
                           dx, dt)
    end if

  end subroutine trace_ppm



  subroutine trace_ppm_rhoe(lo, hi, &
                            idir, &
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
                            qm_core, qcm_lo, qcm_hi, &
                            qp_core, qcp_lo, qcp_hi, &
#if (AMREX_SPACEDIM < 3)
                            dloga, dloga_lo, dloga_hi, &
#endif
                            vlo, vhi, domlo, domhi, &
                            dx, dt)

    use network, only : nspec, naux

    use meth_params_module, only : NQC, NQP, NQAUX, NQC_SRC, QRHO, QU, QV, QW, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                                   NQP_SRC, &
#endif
                                   QREINT, QPRES, QGAME, QC, QGAMC, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   ppm_reference_eigenvectors
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qc_lo(3), qc_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)

    integer, intent(in) :: qa_lo(3), qa_hi(3)

    integer, intent(in) :: Icp_lo(3), Icp_hi(3)
    integer, intent(in) :: Icm_lo(3), Icm_hi(3)
    integer, intent(in) :: Ipp_lo(3), Ipp_hi(3)
    integer, intent(in) :: Ipm_lo(3), Ipm_hi(3)

    integer, intent(in) :: Icsp_lo(3), Icsp_hi(3)
    integer, intent(in) :: Icsm_lo(3), Icsm_hi(3)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    integer, intent(in) :: Ipsp_lo(3), Ipsp_hi(3)
    integer, intent(in) :: Ipsm_lo(3), Ipsm_hi(3)
#endif

    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)

    integer, intent(in) :: qcm_lo(3), qcm_hi(3)
    integer, intent(in) :: qcp_lo(3), qcp_hi(3)

#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q_core(qc_lo(1):qc_hi(1),qc_lo(2):qc_hi(2),qc_lo(3):qc_hi(3),NQC)
    real(rt), intent(in) :: q_pass(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQP)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip_core(Icp_lo(1):Icp_hi(1),Icp_lo(2):Icp_hi(2),Icp_lo(3):Icp_hi(3),1:3,NQC)
    real(rt), intent(in) :: Im_core(Icm_lo(1):Icm_hi(1),Icm_lo(2):Icm_hi(2),Icm_lo(3):Icm_hi(3),1:3,NQC)
    real(rt), intent(in) :: Ip_pass(Ipp_lo(1):Ipp_hi(1),Ipp_lo(2):Ipp_hi(2),Ipp_lo(3):Ipp_hi(3),1:3,NQP)
    real(rt), intent(in) :: Im_pass(Ipm_lo(1):Ipm_hi(1),Ipm_lo(2):Ipm_hi(2),Ipm_lo(3):Ipm_hi(3),1:3,NQP)

    real(rt), intent(in) :: Ip_core_src(Icsp_lo(1):Icsp_hi(1),Icsp_lo(2):Icsp_hi(2),Icsp_lo(3):Icsp_hi(3),1:3,NQC_SRC)
    real(rt), intent(in) :: Im_core_src(Icsm_lo(1):Icsm_hi(1),Icsm_lo(2):Icsm_hi(2),Icsm_lo(3):Icsm_hi(3),1:3,NQC_SRC)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    real(rt), intent(in) :: Ip_pass_src(Ipsp_lo(1):Ipsp_hi(1),Ipsp_lo(2):Ipsp_hi(2),Ipsp_lo(3):Ipsp_hi(3),1:3,NQP_SRC)
    real(rt), intent(in) :: Im_pass_src(Ipsm_lo(1):Ipsm_hi(1),Ipsm_lo(2):Ipsm_hi(2),Ipsm_lo(3):Ipsm_hi(3),1:3,NQP_SRC)
#endif

    real(rt), intent(in) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:3,1)
    real(rt), intent(in) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:3,1)

    real(rt), intent(inout) :: qm_core(qcm_lo(1):qcm_hi(1),qcm_lo(2):qcm_hi(2),qcm_lo(3):qcm_hi(3),NQC)
    real(rt), intent(inout) :: qp_core(qcp_lo(1):qcp_hi(1),qcp_lo(2):qcp_hi(2),qcp_lo(3):qcp_hi(3),NQC)

#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j, k

    real(rt) :: hdt

    integer :: QUN, QUT, QUTT

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! gam_g : the gas Gamma_1
    ! game : gas gamma_e
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    real(rt) :: cc, csq
    real(rt) :: rho, un, ut, utt, p, rhoe_g, h_g
    real(rt) :: gam_g

    real(rt) :: drho, dptot, drhoe_g
    real(rt) :: dup, dptotp
    real(rt) :: dum, dptotm

    real(rt) :: rho_ref, un_ref, p_ref, rhoe_g_ref, h_g_ref

    real(rt) :: cc_ref, csq_ref, gam_g_ref
    real(rt) :: cc_ev, csq_ev, rho_ev, h_g_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp

#ifndef AMREX_USE_CUDA
    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    !$gpu

    hdt = HALF * dt

    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    !
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    !
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).
    !
    ! The indices are: Ip(i, j, k, dim, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.

    if (idir == 1) then
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! Trace to left and right edges using upwind PPM
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rho = q_core(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2

             un = q_core(i,j,k,QUN)
             ut = q_core(i,j,k,QUT)
             utt = q_core(i,j,k,QUTT)

             p = q_core(i,j,k,QPRES)
             rhoe_g = q_core(i,j,k,QREINT)
             h_g = ( (p + rhoe_g)/rho)/csq

             gam_g = qaux(i,j,k,QGAMC)


             !-------------------------------------------------------------------
             ! plus state on face i
             !-------------------------------------------------------------------

             if ((idir == 1 .and. i >= vlo(1)) .or. &
                 (idir == 2 .and. j >= vlo(2)) .or. &
                 (idir == 3 .and. k >= vlo(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im_core(i,j,k,1,QRHO)
                un_ref    = Im_core(i,j,k,1,QUN)

                p_ref    = Im_core(i,j,k,1,QPRES)
                rhoe_g_ref = Im_core(i,j,k,1,QREINT)

                gam_g_ref  = Im_gc(i,j,k,1,1)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! *m are the jumps carried by un-c
                ! *p are the jumps carried by un+c

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dum = un_ref - Im_core(i,j,k,1,QUN) - hdt*Im_core_src(i,j,k,1,QUN)
                dptotm = p_ref - Im_core(i,j,k,1,QPRES) - hdt*Im_core_src(i,j,k,1,QPRES)

                drho = rho_ref - Im_core(i,j,k,2,QRHO) - hdt*Im_core_src(i,j,k,2,QRHO)
                dptot = p_ref - Im_core(i,j,k,2,QPRES) - hdt*Im_core_src(i,j,k,2,QPRES)
                drhoe_g = rhoe_g_ref - Im_core(i,j,k,2,QREINT) - hdt*Im_core_src(i,j,k,2,QREINT)

                dup = un_ref - Im_core(i,j,k,3,QUN) - hdt*Im_core_src(i,j,k,3,QUN)
                dptotp = p_ref - Im_core(i,j,k,3,QPRES) - hdt*Im_core_src(i,j,k,3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   h_g_ev = h_g
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   h_g_ev = h_g_ref
                endif

                ! (rho, u, p, (rho e) eigensystem

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is
                ! simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dum)*(rho_ev/cc_ev)
                alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dup)*(rho_ev/cc_ev)
                alpha0r = drho - dptot/csq_ev
                alpha0e_g = drhoe_g - dptot*h_g_ev  ! note h_g has a 1/c**2 in it

                if (un-cc > ZERO) then
                   alpham = ZERO
                else
                   alpham = -alpham
                end if

                if (un+cc > ZERO) then
                   alphap = ZERO
                else
                   alphap = -alphap
                end if

                if (un > ZERO) then
                   alpha0r = ZERO
                else
                   alpha0r = -alpha0r
                end if

                if (un > ZERO) then
                   alpha0e_g = ZERO
                else
                   alpha0e_g = -alpha0e_g
                end if

                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                qp_core(i,j,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                qp_core(i,j,k,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                qp_core(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                qp_core(i,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)


                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp_core(i,j,k,QUT) = Im_core(i,j,k,2,QUT) + hdt*Im_core_src(i,j,k,2,QUT)
                qp_core(i,j,k,QUTT) = Im_core(i,j,k,2,QUTT) + hdt*Im_core_src(i,j,k,2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= vhi(1)) .or. &
                 (idir == 2 .and. j <= vhi(2)) .or. &
                 (idir == 3 .and. k <= vhi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip_core(i,j,k,3,QRHO)
                un_ref    = Ip_core(i,j,k,3,QUN)

                p_ref    = Ip_core(i,j,k,3,QPRES)
                rhoe_g_ref = Ip_core(i,j,k,3,QREINT)

                gam_g_ref  = Ip_gc(i,j,k,3,1)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip_core(i,j,k,1,QUN) - hdt*Ip_core_src(i,j,k,1,QUN)
                dptotm  = p_ref - Ip_core(i,j,k,1,QPRES) - hdt*Ip_core_src(i,j,k,1,QPRES)

                drho = rho_ref - Ip_core(i,j,k,2,QRHO) - hdt*Ip_core_src(i,j,k,2,QRHO)
                dptot = p_ref - Ip_core(i,j,k,2,QPRES) - hdt*Ip_core_src(i,j,k,2,QPRES)
                drhoe_g = rhoe_g_ref - Ip_core(i,j,k,2,QREINT) - hdt*Ip_core_src(i,j,k,2,QREINT)

                dup = un_ref - Ip_core(i,j,k,3,QUN) - hdt*Ip_core_src(i,j,k,3,QUN)
                dptotp = p_ref - Ip_core(i,j,k,3,QPRES) - hdt*Ip_core_src(i,j,k,3,QPRES)

                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   h_g_ev = h_g
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   h_g_ev = h_g_ref
                endif

                ! (rho, u, p, (rho e)) eigensystem

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is
                ! simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dum)*(rho_ev/cc_ev)
                alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dup)*(rho_ev/cc_ev)
                alpha0r = drho - dptot/csq_ev
                alpha0e_g = drhoe_g - dptot*h_g_ev  ! h_g has a 1/c**2 in it

                if (un-cc > ZERO) then
                   alpham = -alpham
                else
                   alpham = ZERO
                end if

                if (un+cc > ZERO) then
                   alphap = -alphap
                else
                   alphap = ZERO
                end if

                if (un > ZERO) then
                   alpha0r = -alpha0r
                else
                   alpha0r = ZERO
                end if

                if (un > ZERO) then
                   alpha0e_g = -alpha0e_g
                else
                   alpha0e_g = ZERO
                end if

                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                if (idir == 1) then
                   qm_core(i+1,j,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm_core(i+1,j,k,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                   qm_core(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                   qm_core(i+1,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)

                   ! transverse velocities
                   qm_core(i+1,j,k,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i+1,j,k,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                else if (idir == 2) then
                   qm_core(i,j+1,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm_core(i,j+1,k,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                   qm_core(i,j+1,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                   qm_core(i,j+1,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)

                   ! transverse velocities
                   qm_core(i,j+1,k,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i,j+1,k,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                else if (idir == 3) then
                   qm_core(i,j,k+1,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm_core(i,j,k+1,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                   qm_core(i,j,k+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                   qm_core(i,j,k+1,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)

                   ! transverse velocities
                   qm_core(i,j,k+1,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i,j,k+1,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)
                endif

             end if

             !-------------------------------------------------------------------
             ! geometry source terms
             !-------------------------------------------------------------------

#if (AMREX_SPACEDIM < 3)
             ! these only apply for x states (idir = 1)
             if (idir == 1 .and. dloga(i,j,k) /= 0) then
                courn = dt/dx(1)*(cc+abs(un))
                eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k)))
                dlogatmp = min(eta, ONE)*dloga(i,j,k)
                sourcr = -HALF*dt*rho*dlogatmp*un
                sourcp = sourcr*csq
                source = sourcp*h_g

                if (i <= vhi(1)) then

                   qm_core(i+1,j,k,QRHO) = qm_core(i+1,j,k,QRHO) + sourcr
                   qm_core(i+1,j,k,QRHO) = max(qm_core(i+1,j,k,QRHO), small_dens)
                   qm_core(i+1,j,k,QPRES) = qm_core(i+1,j,k,QPRES) + sourcp
                   qm_core(i+1,j,k,QREINT) = qm_core(i+1,j,k,QREINT) + source
                end if

                if (i >= vlo(1)) then

                   qp_core(i,j,k,QRHO) = qp_core(i,j,k,QRHO) + sourcr
                   qp_core(i,j,k,QRHO) = max(qp_core(i,j,k,QRHO), small_dens)
                   qp_core(i,j,k,QPRES) = qp_core(i,j,k,QPRES) + sourcp
                   qp_core(i,j,k,QREINT) = qp_core(i,j,k,QREINT) + source
                end if

             end if
#endif

          end do
       end do
    end do


  end subroutine trace_ppm_rhoe

  subroutine trace_ppm_gammae(lo, hi, &
                              idir, &
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
                              qm_core, qcm_lo, qcm_hi, &
                              qp_core, qcp_lo, qcp_hi, &
#if (AMREX_SPACEDIM < 3)
                              dloga, dloga_lo, dloga_hi, &
#endif
                              vlo, vhi, domlo, domhi, &
                              dx, dt)

    use network, only : nspec, naux
    use meth_params_module, only : NQC, NQP, NQAUX, NQC_SRC, QRHO, QU, QV, QW, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                                   NQP_SRC, &
#endif
                                   QREINT, QPRES, QGAME, QC, QGAMC, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   ppm_reference_eigenvectors
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qc_lo(3), qc_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)

    integer, intent(in) :: qa_lo(3), qa_hi(3)

    integer, intent(in) :: Icp_lo(3), Icp_hi(3)
    integer, intent(in) :: Icm_lo(3), Icm_hi(3)
    integer, intent(in) :: Ipp_lo(3), Ipp_hi(3)
    integer, intent(in) :: Ipm_lo(3), Ipm_hi(3)

    integer, intent(in) :: Icsp_lo(3), Icsp_hi(3)
    integer, intent(in) :: Icsm_lo(3), Icsm_hi(3)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    integer, intent(in) :: Ipsp_lo(3), Ipsp_hi(3)
    integer, intent(in) :: Ipsm_lo(3), Ipsm_hi(3)
#endif

    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)

    integer, intent(in) :: qcm_lo(3), qcm_hi(3)
    integer, intent(in) :: qcp_lo(3), qcp_hi(3)

#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q_core(qc_lo(1):qc_hi(1),qc_lo(2):qc_hi(2),qc_lo(3):qc_hi(3),NQC)
    real(rt), intent(in) :: q_pass(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQP)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip_core(Icp_lo(1):Icp_hi(1),Icp_lo(2):Icp_hi(2),Icp_lo(3):Icp_hi(3),1:3,NQC)
    real(rt), intent(in) :: Im_core(Icm_lo(1):Icm_hi(1),Icm_lo(2):Icm_hi(2),Icm_lo(3):Icm_hi(3),1:3,NQC)
    real(rt), intent(in) :: Ip_pass(Ipp_lo(1):Ipp_hi(1),Ipp_lo(2):Ipp_hi(2),Ipp_lo(3):Ipp_hi(3),1:3,NQP)
    real(rt), intent(in) :: Im_pass(Ipm_lo(1):Ipm_hi(1),Ipm_lo(2):Ipm_hi(2),Ipm_lo(3):Ipm_hi(3),1:3,NQP)

    real(rt), intent(in) :: Ip_core_src(Icsp_lo(1):Icsp_hi(1),Icsp_lo(2):Icsp_hi(2),Icsp_lo(3):Icsp_hi(3),1:3,NQC_SRC)
    real(rt), intent(in) :: Im_core_src(Icsm_lo(1):Icsm_hi(1),Icsm_lo(2):Icsm_hi(2),Icsm_lo(3):Icsm_hi(3),1:3,NQC_SRC)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    real(rt), intent(in) :: Ip_pass_src(Ipsp_lo(1):Ipsp_hi(1),Ipsp_lo(2):Ipsp_hi(2),Ipsp_lo(3):Ipsp_hi(3),1:3,NQP_SRC)
    real(rt), intent(in) :: Im_pass_src(Ipsm_lo(1):Ipsm_hi(1),Ipsm_lo(2):Ipsm_hi(2),Ipsm_lo(3):Ipsm_hi(3),1:3,NQP_SRC)
#endif

    real(rt), intent(in) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:3,1)
    real(rt), intent(in) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:3,1)

    real(rt), intent(inout) :: qm_core(qcm_lo(1):qcm_hi(1),qcm_lo(2):qcm_hi(2),qcm_lo(3):qcm_hi(3),NQC)
    real(rt), intent(inout) :: qp_core(qcp_lo(1):qcp_hi(1),qcp_lo(2):qcp_hi(2),qcp_lo(3):qcp_hi(3),NQC)

#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j, k

    real(rt) :: hdt

    integer :: QUN, QUT, QUTT

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! gam_g : the gas Gamma_1
    ! game : gas gamma_e
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    real(rt) :: cc, csq, Clag
    real(rt) :: rho, un, ut, utt, p, rhoe_g, h_g
    real(rt) :: gam_g, game

    real(rt) :: dptot
    real(rt) :: dge, dtau, dtaum, dtaup
    real(rt) :: dup, dptotp
    real(rt) :: dum, dptotm

    real(rt) :: rho_ref, un_ref, p_ref, rhoe_g_ref
    real(rt) :: tau_ref

    real(rt) :: Clag_ref, gam_g_ref, game_ref, gfactor
    real(rt) :: Clag_ev, tau_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp
    real(rt) :: tau_s

#ifndef AMREX_USE_CUDA
    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    !$gpu

    hdt = HALF * dt

    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    !
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    !
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).
    !
    ! The indices are: Ip(i, j, k, dim, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.

    if (idir == 1) then
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! Trace to left and right edges using upwind PPM
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             gfactor = ONE ! to help compiler resolve ANTI dependence

             rho = q_core(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2
             Clag = rho*cc

             un = q_core(i,j,k,QUN)
             ut = q_core(i,j,k,QUT)
             utt = q_core(i,j,k,QUTT)

             p = q_core(i,j,k,QPRES)
             rhoe_g = q_core(i,j,k,QREINT)

             gam_g = qaux(i,j,k,QGAMC)
             game = q_core(i,j,k,QGAME)


             !-------------------------------------------------------------------
             ! plus state on face i
             !-------------------------------------------------------------------

             if ((idir == 1 .and. i >= vlo(1)) .or. &
                 (idir == 2 .and. j >= vlo(2)) .or. &
                 (idir == 3 .and. k >= vlo(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im_core(i,j,k,1,QRHO)
                un_ref    = Im_core(i,j,k,1,QUN)

                p_ref    = Im_core(i,j,k,1,QPRES)
                rhoe_g_ref = Im_core(i,j,k,1,QREINT)

                tau_ref  = ONE/Im_core(i,j,k,1,QRHO)

                gam_g_ref  = Im_gc(i,j,k,1,1)
                game_ref = Im_core(i,j,k,1,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                Clag_ref = sqrt(gam_g_ref*p_ref*rho_ref)

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dum = un_ref - Im_core(i,j,k,1,QUN) - hdt*Im_core_src(i,j,k,1,QUN)
                dptotm = p_ref - Im_core(i,j,k,1,QPRES) - hdt*Im_core_src(i,j,k,1,QPRES)

                dptot = p_ref - Im_core(i,j,k,2,QPRES) - hdt*Im_core_src(i,j,k,2,QPRES)

                ! we are treating tau as 1/rho, but we could have reconstructed
                ! it separately
                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Im_core(i,j,k,1,QRHO) + hdt*Im_core_src(i,j,k,1,QRHO)/Im_core(i,j,k,1,QRHO)**2
                dtau  = tau_ref - ONE/Im_core(i,j,k,2,QRHO) + hdt*Im_core_src(i,j,k,2,QRHO)/Im_core(i,j,k,2,QRHO)**2
                dtaup = tau_ref - ONE/Im_core(i,j,k,3,QRHO) + hdt*Im_core_src(i,j,k,3,QRHO)/Im_core(i,j,k,3,QRHO)**2

                dup = un_ref - Im_core(i,j,k,3,QUN) - hdt*Im_core_src(i,j,k,3,QUN)
                dptotp = p_ref - Im_core(i,j,k,3,QPRES) - hdt*Im_core_src(i,j,k,3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   Clag_ev = Clag
                   tau_ev  = ONE/rho
                else
                   Clag_ev = Clag_ref
                   tau_ev  = tau_ref
                endif

                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system

                alpham = HALF*( dum - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                alphap = HALF*(-dup - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                dge   = game_ref - Im_core(i,j,k,2,QGAME)
                gfactor = (game - ONE)*(game - gam_g)
                alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                if (un-cc > ZERO) then
                   alpham = ZERO
                else
                   alpham = -alpham
                end if

                if (un+cc > ZERO) then
                   alphap = ZERO
                else
                   alphap = -alphap
                end if

                if (un > ZERO) then
                   alpha0r = ZERO
                else
                   alpha0r = -alpha0r
                end if

                if (un > ZERO) then
                   alpha0e_g = ZERO
                else
                   alpha0e_g = -alpha0e_g
                end if

                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                tau_s = tau_ref + alphap + alpham + alpha0r
                qp_core(i,j,k,QRHO) = max(small_dens, ONE/tau_s)

                qp_core(i,j,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                qp_core(i,j,k,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

                qp_core(i,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                qp_core(i,j,k,QREINT) = qp_core(i,j,k,QPRES )/(qp_core(i,j,k,QGAME) - ONE)


                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp_core(i,j,k,QUT) = Im_core(i,j,k,2,QUT) + hdt*Im_core_src(i,j,k,2,QUT)
                qp_core(i,j,k,QUTT) = Im_core(i,j,k,2,QUTT) + hdt*Im_core_src(i,j,k,2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= vhi(1)) .or. &
                 (idir == 2 .and. j <= vhi(2)) .or. &
                 (idir == 3 .and. k <= vhi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip_core(i,j,k,3,QRHO)
                un_ref    = Ip_core(i,j,k,3,QUN)

                p_ref    = Ip_core(i,j,k,3,QPRES)
                rhoe_g_ref = Ip_core(i,j,k,3,QREINT)

                tau_ref  = ONE/Ip_core(i,j,k,3,QRHO)

                gam_g_ref  = Ip_gc(i,j,k,3,1)
                game_ref = Ip_core(i,j,k,3,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                Clag_ref = sqrt(gam_g_ref*p_ref*rho_ref)

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip_core(i,j,k,1,QUN) - hdt*Ip_core_src(i,j,k,1,QUN)
                dptotm  = p_ref - Ip_core(i,j,k,1,QPRES) - hdt*Ip_core_src(i,j,k,1,QPRES)

                dptot = p_ref - Ip_core(i,j,k,2,QPRES) - hdt*Ip_core_src(i,j,k,2,QPRES)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Ip_core(i,j,k,1,QRHO) + hdt*Ip_core_src(i,j,k,1,QRHO)/Ip_core(i,j,k,1,QRHO)**2
                dtau = tau_ref - ONE/Ip_core(i,j,k,2,QRHO) + hdt*Ip_core_src(i,j,k,2,QRHO)/Ip_core(i,j,k,2,QRHO)**2
                dtaup = tau_ref - ONE/Ip_core(i,j,k,3,QRHO) + hdt*Ip_core_src(i,j,k,3,QRHO)/Ip_core(i,j,k,3,QRHO)**2

                dup = un_ref - Ip_core(i,j,k,3,QUN) - hdt*Ip_core_src(i,j,k,3,QUN)
                dptotp = p_ref - Ip_core(i,j,k,3,QPRES) - hdt*Ip_core_src(i,j,k,3,QPRES)

                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   Clag_ev = Clag
                   tau_ev  = ONE/rho
                else
                   Clag_ev = Clag_ref
                   tau_ev  = tau_ref
                endif

                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                alpham = HALF*( dum - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                alphap = HALF*(-dup - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                dge = game_ref - Ip_core(i,j,k,2,QGAME)
                gfactor = (game - ONE)*(game - gam_g)
                alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                if (un-cc > ZERO) then
                   alpham = -alpham
                else
                   alpham = ZERO
                end if

                if (un+cc > ZERO) then
                   alphap = -alphap
                else
                   alphap = ZERO
                end if

                if (un > ZERO) then
                   alpha0r = -alpha0r
                else
                   alpha0r = ZERO
                end if

                if (un > ZERO) then
                   alpha0e_g = -alpha0e_g
                else
                   alpha0e_g = ZERO
                end if


                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                if (idir == 1) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm_core(i+1,j,k,QRHO) = max(small_dens, ONE/tau_s)

                   qm_core(i+1,j,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm_core(i+1,j,k,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

                   qm_core(i+1,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                   qm_core(i+1,j,k,QREINT) = qm_core(i+1,j,k,QPRES )/(qm_core(i+1,j,k,QGAME) - ONE)

                   ! transverse velocities
                   qm_core(i+1,j,k,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i+1,j,k,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                else if (idir == 2) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm_core(i,j+1,k,QRHO) = max(small_dens, ONE/tau_s)

                   qm_core(i,j+1,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm_core(i,j+1,k,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

                   qm_core(i,j+1,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                   qm_core(i,j+1,k,QREINT) = qm_core(i,j+1,k,QPRES )/(qm_core(i,j+1,k,QGAME) - ONE)

                   ! transverse velocities
                   qm_core(i,j+1,k,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i,j+1,k,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                else if (idir == 3) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm_core(i,j,k+1,QRHO) = max(small_dens, ONE/tau_s)

                   qm_core(i,j,k+1,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm_core(i,j,k+1,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

                   qm_core(i,j,k+1,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                   qm_core(i,j,k+1,QREINT) = qm_core(i,j,k+1,QPRES )/(qm_core(i,j,k+1,QGAME) - ONE)

                   ! transverse velocities
                   qm_core(i,j,k+1,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i,j,k+1,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                end if
             end if

             !-------------------------------------------------------------------
             ! geometry source terms
             !-------------------------------------------------------------------

#if (AMREX_SPACEDIM < 3)
             ! these only apply for x states (dim = 1)
             if (idir == 1 .and. dloga(i,j,k) /= 0) then
                h_g = ( (p + rhoe_g)/rho)/csq
                courn = dt/dx(1)*(cc+abs(un))
                eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k)))
                dlogatmp = min(eta, ONE)*dloga(i,j,k)
                sourcr = -HALF*dt*rho*dlogatmp*un
                sourcp = sourcr*csq
                source = sourcp*h_g

                if (i <= vhi(1)) then
                   qm_core(i+1,j,k,QRHO) = qm_core(i+1,j,k,QRHO) + sourcr
                   qm_core(i+1,j,k,QRHO) = max(qm_core(i+1,j,k,QRHO), small_dens)
                   qm_core(i+1,j,k,QPRES) = qm_core(i+1,j,k,QPRES) + sourcp
                   qm_core(i+1,j,k,QREINT) = qm_core(i+1,j,k,QREINT) + source
                end if

                if (i >= vlo(1)) then
                   qp_core(i,j,k,QRHO) = qp_core(i,j,k,QRHO) + sourcr
                   qp_core(i,j,k,QRHO) = max(qp_core(i,j,k,QRHO), small_dens)
                   qp_core(i,j,k,QPRES) = qp_core(i,j,k,QPRES) + sourcp
                   qp_core(i,j,k,QREINT) = qp_core(i,j,k,QREINT) + source
                end if

             endif
#endif

          end do
       end do
    end do

  end subroutine trace_ppm_gammae


  subroutine trace_ppm_temp(lo, hi, &
                            idir, &
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
                            qm_core, qcm_lo, qcm_hi, &
                            qp_core, qcp_lo, qcp_hi, &
                            qm_pass, qpm_lo, qpm_hi, &
                            qp_pass, qpp_lo, qpp_hi, &
#if (AMREX_SPACEDIM < 3)
                            dloga, dloga_lo, dloga_hi, &
#endif
                            vlo, vhi, domlo, domhi, &
                            dx, dt)

    use network, only : nspec, naux
    use meth_params_module, only : NQC, NQP, NQAUX, NQC_SRC, QRHO, QU, QV, QW, &
#ifdef PRIM_SPECIES_HAVE_SOURCES
                                   NQP_SRC, &
#endif
                                   QREINT, QPRES, QTEMP, QGAME, QC, QGAMC, QFS, QFX, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   ppm_reference_eigenvectors
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qc_lo(3), qc_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)

    integer, intent(in) :: qa_lo(3), qa_hi(3)

    integer, intent(in) :: Icp_lo(3), Icp_hi(3)
    integer, intent(in) :: Icm_lo(3), Icm_hi(3)
    integer, intent(in) :: Ipp_lo(3), Ipp_hi(3)
    integer, intent(in) :: Ipm_lo(3), Ipm_hi(3)

    integer, intent(in) :: Icsp_lo(3), Icsp_hi(3)
    integer, intent(in) :: Icsm_lo(3), Icsm_hi(3)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    integer, intent(in) :: Ipsp_lo(3), Ipsp_hi(3)
    integer, intent(in) :: Ipsm_lo(3), Ipsm_hi(3)
#endif

    integer, intent(in) :: Ipg_lo(3), Ipg_hi(3)
    integer, intent(in) :: Img_lo(3), Img_hi(3)

    integer, intent(in) :: qcm_lo(3), qcm_hi(3)
    integer, intent(in) :: qcp_lo(3), qcp_hi(3)
    integer, intent(in) :: qpm_lo(3), qpm_hi(3)
    integer, intent(in) :: qpp_lo(3), qpp_hi(3)

#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q_core(qc_lo(1):qc_hi(1),qc_lo(2):qc_hi(2),qc_lo(3):qc_hi(3),NQC)
    real(rt), intent(in) :: q_pass(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQP)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip_core(Icp_lo(1):Icp_hi(1),Icp_lo(2):Icp_hi(2),Icp_lo(3):Icp_hi(3),1:3,NQC)
    real(rt), intent(in) :: Im_core(Icm_lo(1):Icm_hi(1),Icm_lo(2):Icm_hi(2),Icm_lo(3):Icm_hi(3),1:3,NQC)
    real(rt), intent(in) :: Ip_pass(Ipp_lo(1):Ipp_hi(1),Ipp_lo(2):Ipp_hi(2),Ipp_lo(3):Ipp_hi(3),1:3,NQP)
    real(rt), intent(in) :: Im_pass(Ipm_lo(1):Ipm_hi(1),Ipm_lo(2):Ipm_hi(2),Ipm_lo(3):Ipm_hi(3),1:3,NQP)

    real(rt), intent(in) :: Ip_core_src(Icsp_lo(1):Icsp_hi(1),Icsp_lo(2):Icsp_hi(2),Icsp_lo(3):Icsp_hi(3),1:3,NQC_SRC)
    real(rt), intent(in) :: Im_core_src(Icsm_lo(1):Icsm_hi(1),Icsm_lo(2):Icsm_hi(2),Icsm_lo(3):Icsm_hi(3),1:3,NQC_SRC)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    real(rt), intent(in) :: Ip_pass_src(Ipsp_lo(1):Ipsp_hi(1),Ipsp_lo(2):Ipsp_hi(2),Ipsp_lo(3):Ipsp_hi(3),1:3,NQP_SRC)
    real(rt), intent(in) :: Im_pass_src(Ipsm_lo(1):Ipsm_hi(1),Ipsm_lo(2):Ipsm_hi(2),Ipsm_lo(3):Ipsm_hi(3),1:3,NQP_SRC)
#endif

    real(rt), intent(in) :: Ip_gc(Ipg_lo(1):Ipg_hi(1),Ipg_lo(2):Ipg_hi(2),Ipg_lo(3):Ipg_hi(3),1:3,1)
    real(rt), intent(in) :: Im_gc(Img_lo(1):Img_hi(1),Img_lo(2):Img_hi(2),Img_lo(3):Img_hi(3),1:3,1)

    real(rt), intent(inout) :: qm_core(qcm_lo(1):qcm_hi(1),qcm_lo(2):qcm_hi(2),qcm_lo(3):qcm_hi(3),NQC)
    real(rt), intent(inout) :: qp_core(qcp_lo(1):qcp_hi(1),qcp_lo(2):qcp_hi(2),qcp_lo(3):qcp_hi(3),NQC)
    real(rt), intent(in) :: qm_pass(qpm_lo(1):qpm_hi(1),qpm_lo(2):qpm_hi(2),qpm_lo(3):qpm_hi(3),NQP)
    real(rt), intent(in) :: qp_pass(qpp_lo(1):qpp_hi(1),qpp_lo(2):qpp_hi(2),qpp_lo(3):qpp_hi(3),NQP)

#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j, k

    type(eos_t) :: eos_state

    real(rt) :: hdt

    integer :: QUN, QUT, QUTT

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! gam_g : the gas Gamma_1
    ! game : gas gamma_e
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    real(rt) :: cc, csq, Clag
    real(rt) :: rho, un, ut, utt, p, rhoe_g, h_g, temp
    real(rt) :: gam_g, game

    real(rt) :: drho, dptot
    real(rt) :: dtau, dtaum, dtaup
    real(rt) :: dup, dptotp
    real(rt) :: dum, dptotm
    real(rt) :: dT0, dTp, dTm
    real(rt) :: p_r, p_T

    real(rt) :: rho_ref, un_ref, p_ref, temp_ref
    real(rt) :: tau_ref

    real(rt) :: cc_ref, csq_ref, Clag_ref, gam_g_ref, game_ref, gfactor
    real(rt) :: cc_ev, csq_ev, Clag_ev, rho_ev, tau_ev, temp_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp
    real(rt) :: tau_s

#ifndef AMREX_USE_CUDA
    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    !$gpu

    hdt = HALF * dt


    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    !
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    !
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).
    !
    ! The indices are: Ip(i, j, k, dim, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.

    if (idir == 1) then
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! Trace to left and right edges using upwind PPM
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             gfactor = ONE ! to help compiler resolve ANTI dependence

             rho = q_core(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2
             Clag = rho*cc

             un = q_core(i,j,k,QUN)
             ut = q_core(i,j,k,QUT)
             utt = q_core(i,j,k,QUTT)

             p = q_core(i,j,k,QPRES)
             rhoe_g = q_core(i,j,k,QREINT)
             temp = q_core(i,j,k,QTEMP)

             gam_g = qaux(i,j,k,QGAMC)
             game = q_core(i,j,k,QGAME)


             !-------------------------------------------------------------------
             ! plus state on face i
             !-------------------------------------------------------------------

             if ((idir == 1 .and. i >= vlo(1)) .or. &
                 (idir == 2 .and. j >= vlo(2)) .or. &
                 (idir == 3 .and. k >= vlo(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im_core(i,j,k,1,QRHO)
                un_ref    = Im_core(i,j,k,1,QUN)

                p_ref    = Im_core(i,j,k,1,QPRES)
                temp_ref = Im_core(i,j,k,1,QTEMP)

                tau_ref  = ONE/Im_core(i,j,k,1,QRHO)

                gam_g_ref  = Im_gc(i,j,k,1,1)
                game_ref = Im_core(i,j,k,1,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dum = un_ref - Im_core(i,j,k,1,QUN) - hdt*Im_core_src(i,j,k,1,QUN)
                dptotm = p_ref - Im_core(i,j,k,1,QPRES) - hdt*Im_core_src(i,j,k,1,QPRES)

                drho = rho_ref - Im_core(i,j,k,2,QRHO) - hdt*Im_core_src(i,j,k,2,QRHO)
                dptot = p_ref - Im_core(i,j,k,2,QPRES) - hdt*Im_core_src(i,j,k,2,QPRES)

                ! TODO: need to figure sources for this out...
                dTm = temp_ref - Im_core(i,j,k,1,QTEMP)
                dT0 = temp_ref - Im_core(i,j,k,2,QTEMP)
                dTp = temp_ref - Im_core(i,j,k,3,QTEMP)

                ! we are treating tau as 1/rho, but we could have reconstructed
                ! it separately
                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Im_core(i,j,k,1,QRHO) + hdt*Im_core_src(i,j,k,1,QRHO)/Im_core(i,j,k,1,QRHO)**2
                dtau  = tau_ref - ONE/Im_core(i,j,k,2,QRHO) + hdt*Im_core_src(i,j,k,2,QRHO)/Im_core(i,j,k,2,QRHO)**2
                dtaup = tau_ref - ONE/Im_core(i,j,k,3,QRHO) + hdt*Im_core_src(i,j,k,3,QRHO)/Im_core(i,j,k,3,QRHO)**2

                dup = un_ref - Im_core(i,j,k,3,QUN) - hdt*Im_core_src(i,j,k,3,QUN)
                dptotp = p_ref - Im_core(i,j,k,3,QPRES) - hdt*Im_core_src(i,j,k,3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   Clag_ev = Clag
                   tau_ev  = ONE/rho
                   temp_ev = temp
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   Clag_ev = Clag_ref
                   tau_ev  = tau_ref
                   temp_ev = temp_ref
                endif

                ! (tau, u T) eigensystem

                ! eos to get some thermodynamics
                eos_state%T = temp_ev
                eos_state%rho = rho_ev
                eos_state%xn(:) = q_pass(i,j,k,QFS:QFS-1+nspec)
                eos_state%aux(:) = q_pass(i,j,k,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                p_r = eos_state%dpdr
                p_T = eos_state%dpdT

                alpham = HALF*(rho_ev**2*p_r*dtaum/Clag_ev + dum - p_T*dTm/Clag_ev)/Clag_ev
                alphap = HALF*(rho_ev**2*p_r*dtaup/Clag_ev - dup - p_T*dTp/Clag_ev)/Clag_ev
                alpha0r = dtau + (-rho_ev**2*p_r*dtau + p_T*dT0)/Clag_ev**2

                ! not used, but needed to prevent bad invalid ops
                alpha0e_g = ZERO

                if (un-cc > ZERO) then
                   alpham = ZERO
                else
                   alpham = -alpham
                end if

                if (un+cc > ZERO) then
                   alphap = ZERO
                else
                   alphap = -alphap
                end if

                if (un > ZERO) then
                   alpha0r = ZERO
                else
                   alpha0r = -alpha0r
                end if

                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                tau_s = tau_ref + alphap + alpham + alpha0r
                qp_core(i,j,k,QRHO) = max(small_dens, ONE/tau_s)

                qp_core(i,j,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                qp_core(i,j,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                     rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                ! we defer getting the pressure until later, once we do the species
                qp_core(i,j,k,QPRES) = small_pres ! just to make it defined

                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp_core(i,j,k,QUT) = Im_core(i,j,k,2,QUT) + hdt*Im_core_src(i,j,k,2,QUT)
                qp_core(i,j,k,QUTT) = Im_core(i,j,k,2,QUTT) + hdt*Im_core_src(i,j,k,2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= vhi(1)) .or. &
                 (idir == 2 .and. j <= vhi(2)) .or. &
                 (idir == 3 .and. k <= vhi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip_core(i,j,k,3,QRHO)
                un_ref    = Ip_core(i,j,k,3,QUN)

                p_ref    = Ip_core(i,j,k,3,QPRES)
                temp_ref = Ip_core(i,j,k,3,QTEMP)

                tau_ref  = ONE/Ip_core(i,j,k,3,QRHO)

                gam_g_ref  = Ip_gc(i,j,k,3,1)
                game_ref = Ip_core(i,j,k,3,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip_core(i,j,k,1,QUN) - hdt*Ip_core_src(i,j,k,1,QUN)
                dptotm  = p_ref - Ip_core(i,j,k,1,QPRES) - hdt*Ip_core_src(i,j,k,1,QPRES)

                drho = rho_ref - Ip_core(i,j,k,2,QRHO) - hdt*Ip_core_src(i,j,k,2,QRHO)
                dptot = p_ref - Ip_core(i,j,k,2,QPRES) - hdt*Ip_core_src(i,j,k,2,QPRES)

                dTm = temp_ref - Ip_core(i,j,k,1,QTEMP)
                dT0 = temp_ref - Ip_core(i,j,k,2,QTEMP)
                dTp = temp_ref - Ip_core(i,j,k,3,QTEMP)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Ip_core(i,j,k,1,QRHO) + hdt*Ip_core_src(i,j,k,1,QRHO)/Ip_core(i,j,k,1,QRHO)**2
                dtau = tau_ref - ONE/Ip_core(i,j,k,2,QRHO) + hdt*Ip_core_src(i,j,k,2,QRHO)/Ip_core(i,j,k,2,QRHO)**2
                dtaup = tau_ref - ONE/Ip_core(i,j,k,3,QRHO) + hdt*Ip_core_src(i,j,k,3,QRHO)/Ip_core(i,j,k,3,QRHO)**2

                dup = un_ref - Ip_core(i,j,k,3,QUN) - hdt*Ip_core_src(i,j,k,3,QUN)
                dptotp = p_ref - Ip_core(i,j,k,3,QPRES) - hdt*Ip_core_src(i,j,k,3,QPRES)

                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   Clag_ev = Clag
                   tau_ev  = ONE/rho
                   temp_ev = temp
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   Clag_ev = Clag_ref
                   tau_ev  = tau_ref
                   temp_ev = temp_ref
                endif

                ! (tau, u T) eigensystem

                ! eos to get some thermodynamics
                eos_state%T = temp_ev
                eos_state%rho = rho_ev
                eos_state%xn(:) = q_pass(i,j,k,QFS:QFS-1+nspec)
                eos_state%aux(:) = q_pass(i,j,k,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                p_r = eos_state%dpdr
                p_T = eos_state%dpdT

                alpham = HALF*(rho_ev**2*p_r*dtaum/Clag_ev + dum - p_T*dTm/Clag_ev)/Clag_ev
                alphap = HALF*(rho_ev**2*p_r*dtaup/Clag_ev - dup - p_T*dTp/Clag_ev)/Clag_ev
                alpha0r = dtau + (-rho_ev**2*p_r*dtau + p_T*dT0)/Clag_ev**2

                ! not used, but needed to prevent bad invalid ops
                alpha0e_g = ZERO

                if (un-cc > ZERO) then
                   alpham = -alpham
                else
                   alpham = ZERO
                end if

                if (un+cc > ZERO) then
                   alphap = -alphap
                else
                   alphap = ZERO
                end if

                if (un > ZERO) then
                   alpha0r = -alpha0r
                else
                   alpha0r = ZERO
                end if


                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                if (idir == 1) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm_core(i+1,j,k,QRHO) = max(small_dens, ONE/tau_s)

                   qm_core(i+1,j,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm_core(i+1,j,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once
                   ! we do the species
                   qm_core(i+1,j,k,QPRES) = small_pres ! just to make it defined

                   ! transverse velocities
                   qm_core(i+1,j,k,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i+1,j,k,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                else if (idir == 2) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm_core(i,j+1,k,QRHO) = max(small_dens, ONE/tau_s)

                   qm_core(i,j+1,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm_core(i,j+1,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once
                   ! we do the species
                   qm_core(i,j+1,k,QPRES) = small_pres ! just to make it defined

                   ! transverse velocities
                   qm_core(i,j+1,k,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i,j+1,k,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                else if (idir == 3) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm_core(i,j,k+1,QRHO) = max(small_dens, ONE/tau_s)

                   qm_core(i,j,k+1,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm_core(i,j,k+1,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once
                   ! we do the species
                   qm_core(i,j,k+1,QPRES) = small_pres ! just to make it defined

                   ! transverse velocities
                   qm_core(i,j,k+1,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i,j,k+1,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                endif

             end if

             !-------------------------------------------------------------------
             ! geometry source terms
             !-------------------------------------------------------------------

#if (AMREX_SPACEDIM < 3)
             ! these only apply for x states (idir = 1)
             if (idir == 1 .and. dloga(i,j,k) /= 0) then
                h_g = ( (p + rhoe_g)/rho)/csq
                courn = dt/dx(1)*(cc+abs(un))
                eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k)))
                dlogatmp = min(eta, ONE)*dloga(i,j,k)
                sourcr = -HALF*dt*rho*dlogatmp*un
                sourcp = sourcr*csq
                source = sourcp*h_g

                if (i <= vhi(1)) then
                   qm_core(i+1,j,k,QRHO) = qm_core(i+1,j,k,QRHO) + sourcr
                   qm_core(i+1,j,k,QRHO) = max(qm_core(i+1,j,k,QRHO), small_dens)
                   qm_core(i+1,j,k,QPRES) = qm_core(i+1,j,k,QPRES) + sourcp
                   qm_core(i+1,j,k,QREINT) = qm_core(i+1,j,k,QREINT) + source
                end if

                if (i >= vlo(1)) then
                   qp_core(i,j,k,QRHO) = qp_core(i,j,k,QRHO) + sourcr
                   qp_core(i,j,k,QRHO) = max(qp_core(i,j,k,QRHO), small_dens)
                   qp_core(i,j,k,QPRES) = qp_core(i,j,k,QPRES) + sourcp
                   qp_core(i,j,k,QREINT) = qp_core(i,j,k,QREINT) + source
                end if

             endif
#endif

          end do
       end do
    end do

    ! we predicted T, now make p, (rho e) consistent
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if ((idir == 1 .and. i >= vlo(1)) .or. &
                 (idir == 2 .and. j >= vlo(2)) .or. &
                 (idir == 3 .and. k >= vlo(3))) then

                ! plus face
                eos_state%T     = qp_core(i,j,k,QTEMP)
                eos_state%rho   = qp_core(i,j,k,QRHO)
                eos_state%xn(:) = qp_pass(i,j,k,QFS:QFS-1+nspec)
                eos_state%aux(:) = qp_pass(i,j,k,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                qp_core(i,j,k,QPRES) = eos_state%p
                qp_core(i,j,k,QREINT) = qp_core(i,j,k,QRHO)*eos_state%e

                qp_core(i,j,k,QPRES) = max(qp_core(i,j,k,QPRES), small_pres)
             end if

             if (idir == 1 .and. i <= vhi(1)) then

                ! minus face
                eos_state%T     = qm_core(i+1,j,k,QTEMP)
                eos_state%rho   = qm_core(i+1,j,k,QRHO)
                eos_state%xn(:) = qm_pass(i+1,j,k,QFS:QFS-1+nspec)
                eos_state%aux(:) = qm_pass(i+1,j,k,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                qm_core(i+1,j,k,QPRES) = max(small_pres, eos_state%p)
                qm_core(i+1,j,k,QREINT) = qm_core(i+1,j,k,QRHO)*eos_state%e

             else if (idir == 2 .and. j <= vhi(2)) then

                ! minus face
                eos_state%T     = qm_core(i,j+1,k,QTEMP)
                eos_state%rho   = qm_core(i,j+1,k,QRHO)
                eos_state%xn(:) = qm_pass(i,j+1,k,QFS:QFS-1+nspec)
                eos_state%aux(:) = qm_pass(i,j+1,k,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                qm_core(i,j+1,k,QPRES) = max(small_pres, eos_state%p)
                qm_core(i,j+1,k,QREINT) = qm_core(i,j+1,k,QRHO)*eos_state%e

             else if (idir == 3 .and. k <= vhi(3)) then

                ! minus face
                eos_state%T     = qm_core(i,j,k+1,QTEMP)
                eos_state%rho   = qm_core(i,j,k+1,QRHO)
                eos_state%xn(:) = qm_pass(i,j,k+1,QFS:QFS-1+nspec)
                eos_state%aux(:) = qm_pass(i,j,k+1,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                qm_core(i,j,k+1,QPRES) = max(small_pres, eos_state%p)
                qm_core(i,j,k+1,QREINT) = qm_core(i,j,k+1,QRHO)*eos_state%e

             end if

          end do
       end do
    end do

  end subroutine trace_ppm_temp

end module trace_ppm_module
