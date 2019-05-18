! These routines do the characteristic tracing under the parabolic
! profiles in each zone to the edge / half-time.

module trace_ppm_rad_module

  use prob_params_module, only : dg
  use amrex_error_module, only : amrex_error
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public trace_ppm_rad

contains


  subroutine trace_ppm_rad(lo, hi, &
                           idir, &
                           q_core, qc_lo, qc_hi, &
                           q_pass, qp_lo, qp_hi, &
                           q_rad, qr_lo, qr_hi, &
                           qaux, qa_lo, qa_hi, &
                           Ip_core, Icp_lo, Icp_hi, &
                           Im_core, Icm_lo, Icm_hi, &
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
                           qm_core, qcm_lo, qcm_hi, &
                           qp_core, qcp_lo, qcp_hi, &
                           qm_pass, qpm_lo, qpm_hi, &
                           qp_pass, qpp_lo, qpp_hi, &
                           qm_rad, qrm_lo, qrm_hi, &
                           qp_rad, qrp_lo, qrp_hi, &
#if (AMREX_SPACEDIM < 3)
                           dloga, dloga_lo, dloga_hi, &
#endif
                           vlo, vhi, domlo, domhi, &
                           dx, dt)

    use network, only : nspec, naux
    use meth_params_module, only : NQC, NQP, NQR, NQAUX, NQC_SRC, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QGAME, QC, QCG, QGAMC, QGAMCG, QLAMS, &
                                   qrad, qradhi, qptot, qreitot, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   ppm_reference_eigenvectors, ppm_predict_gammae, &
                                   npassive, qpass_map
    use rad_params_module, only : ngroups
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qc_lo(3), qc_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: Icp_lo(3), Icp_hi(3)
    integer, intent(in) :: Icm_lo(3), Icm_hi(3)
    integer, intent(in) :: Ipp_lo(3), Ipp_hi(3)
    integer, intent(in) :: Ipm_lo(3), Ipm_hi(3)
    integer, intent(in) :: Irp_lo(3), Irp_hi(3)
    integer, intent(in) :: Irm_lo(3), Irm_hi(3)
    integer, intent(in) :: Icsp_lo(3), Icsp_hi(3)
    integer, intent(in) :: Icsm_lo(3), Icsm_hi(3)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    integer, intent(in) :: Ipsp_lo(3), Ipsp_hi(3)
    integer, intent(in) :: Ipsm_lo(3), Ipsm_hi(3)
#endif
    integer, intent(in) :: qcm_lo(3), qcm_hi(3)
    integer, intent(in) :: qcp_lo(3), qcp_hi(3)
    integer, intent(in) :: qpm_lo(3), qpm_hi(3)
    integer, intent(in) :: qpp_lo(3), qpp_hi(3)
    integer, intent(in) :: qrm_lo(3), qrm_hi(3)
    integer, intent(in) :: qrp_lo(3), qrp_hi(3)
#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q_core(qc_lo(1):qc_hi(1),qc_lo(2):qc_hi(2),qc_lo(3):qc_hi(3),NQC)
    real(rt), intent(in) :: q_pass(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQP)
    real(rt), intent(in) :: q_rad(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),NQR)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip_core(Icp_lo(1):Icp_hi(1),Icp_lo(2):Icp_hi(2),Icp_lo(3):Icp_hi(3),1:3,NQC)
    real(rt), intent(in) :: Im_core(Icm_lo(1):Icm_hi(1),Icm_lo(2):Icm_hi(2),Icm_lo(3):Icm_hi(3),1:3,NQC)
    real(rt), intent(in) :: Ip_pass(Ipp_lo(1):Ipp_hi(1),Ipp_lo(2):Ipp_hi(2),Ipp_lo(3):Ipp_hi(3),1:3,NQP)
    real(rt), intent(in) :: Im_pass(Ipm_lo(1):Ipm_hi(1),Ipm_lo(2):Ipm_hi(2),Ipm_lo(3):Ipm_hi(3),1:3,NQP)
    real(rt), intent(in) :: Ip_rad(Irp_lo(1):Irp_hi(1),Irp_lo(2):Irp_hi(2),Irp_lo(3):Irp_hi(3),1:3,NQR)
    real(rt), intent(in) :: Im_rad(Irm_lo(1):Irm_hi(1),Irm_lo(2):Irm_hi(2),Irm_lo(3):Irm_hi(3),1:3,NQR)

    real(rt), intent(in) :: Ip_core_src(Icsp_lo(1):Icsp_hi(1),Icsp_lo(2):Icsp_hi(2),Icsp_lo(3):Icsp_hi(3),1:3,NQC_SRC)
    real(rt), intent(in) :: Im_core_src(Icsm_lo(1):Icsm_hi(1),Icsm_lo(2):Icsm_hi(2),Icsm_lo(3):Icsm_hi(3),1:3,NQC_SRC)
#ifdef PRIM_SPECIES_HAVE_SOURCES
    real(rt), intent(in) :: Ip_pass_src(Ipsp_lo(1):Ipsp_hi(1),Ipsp_lo(2):Ipsp_hi(2),Ipsp_lo(3):Ipsp_hi(3),1:3,NQP_SRC)
    real(rt), intent(in) :: Im_pass_src(Ipsm_lo(1):Ipsm_hi(1),Ipsm_lo(2):Ipsm_hi(2),Ipsm_lo(3):Ipsm_hi(3),1:3,NQP_SRC)
#endif
    real(rt), intent(inout) :: qm_core(qcm_lo(1):qcm_hi(1),qcm_lo(2):qcm_hi(2),qcm_lo(3):qcm_hi(3),NQC)
    real(rt), intent(inout) :: qp_core(qcp_lo(1):qcp_hi(1),qcp_lo(2):qcp_hi(2),qcp_lo(3):qcp_hi(3),NQC)
    real(rt), intent(inout) :: qm_pass(qcm_lo(1):qpm_hi(1),qpm_lo(2):qpm_hi(2),qpm_lo(3):qpm_hi(3),NQP)
    real(rt), intent(inout) :: qp_pass(qcp_lo(1):qpp_hi(1),qpp_lo(2):qpp_hi(2),qpp_lo(3):qpp_hi(3),NQP)
    real(rt), intent(inout) :: qm_rad(qrm_lo(1):qrm_hi(1),qrm_lo(2):qrm_hi(2),qrm_lo(3):qrm_hi(3),NQR)
    real(rt), intent(inout) :: qp_rad(qrp_lo(1):qrp_hi(1),qrp_lo(2):qrp_hi(2),qrp_lo(3):qrp_hi(3),NQR)

#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    integer :: QUN, QUT, QUTT

    ! Local variables
    integer :: i, j, k, g
    integer :: n, ipassive

    real(rt) :: hdt

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

    real(rt) :: cc, csq, cgassq, Clag
    real(rt) :: rho, un, p, rhoe_g, h_g, tau
    real(rt) :: ptot, gam_g, game

    real(rt) :: drho, dptot, drhoe_g
    real(rt) :: de, dge, dtau
    real(rt) :: dup, dvp, dptotp
    real(rt) :: dum, dvm, dptotm

    real(rt) :: rho_ref, un_ref, p_ref, rhoe_g_ref, h_g_ref
    real(rt) :: tau_ref
    real(rt) :: ptot_ref

    real(rt) :: gam_ref, game_ref, gfactor

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp, sourcer(0:ngroups-1)
    real(rt) :: tau_s, e_s

    real(rt), dimension(0:ngroups-1) :: er, der, alphar, qrtmp,hr
    real(rt), dimension(0:ngroups-1) :: lam0, lamp, lamm

    real(rt), dimension(0:ngroups-1) :: er_ref

    real(rt) :: er_foo

    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_rad_nd.f90 :: tracexy_ppm_rad")
    end if

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
    ! The indices are: Ip(i, j, k, wave, var)
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

             do g=0, ngroups-1
                lam0(g) = qaux(i,j,k,QLAMS+g)
                lamp(g) = qaux(i,j,k,QLAMS+g)
                lamm(g) = qaux(i,j,k,QLAMS+g)
             end do

             rho = q_core(i,j,k,QRHO)
             tau = ONE/rho

             ! cgassq is the gas soundspeed **2
             ! cc is the total soundspeed **2 (gas + radiation)
             cgassq = qaux(i,j,k,QCG)**2
             cc = qaux(i,j,k,QC)
             csq = cc**2
             Clag = rho*cc

             un = q_core(i,j,k,QUN)

             p = q_core(i,j,k,QPRES)
             rhoe_g = q_core(i,j,k,QREINT)
             h_g = ( (p+rhoe_g)/rho)/csq

             gam_g = qaux(i,j,k,QGAMCG)
             game = q_core(i,j,k,QGAME)

             ptot = q_rad(i,j,k,qptot)

             er(:) = q_rad(i,j,k,qrad:qradhi)
             hr(:) = (lam0+ONE)*er/rho

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

                !gam_g_ref  = Im_gc(i,j,k,1,1)
                game_ref = Im_core(i,j,k,1,QGAME)

                ptot_ref = Im_rad(i,j,k,1,QPTOT)

                er_ref(:) = Im_rad(i,j,k,1,QRAD:QRADHI)


                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)

                ! we also add the sources here so they participate in the tracing
                dum    = un_ref    - Im_core(i,j,k,1,QUN) - hdt*Im_core_src(i,j,k,1,QUN)
                dptotm = ptot_ref - Im_rad(i,j,k,1,qptot) - hdt*Im_core_src(i,j,k,1,QPRES)

                drho    = rho_ref    - Im_core(i,j,k,2,QRHO) - hdt*Im_core_src(i,j,k,2,QRHO)
                dptot   = ptot_ref   - Im_rad(i,j,k,2,qptot) - hdt*Im_core_src(i,j,k,2,QPRES)
                drhoe_g = rhoe_g_ref - Im_core(i,j,k,2,QREINT) - hdt*Im_core_src(i,j,k,2,QREINT)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtau  = tau_ref  - ONE/Im_core(i,j,k,2,QRHO) + hdt*Im_core_src(i,j,k,2,QRHO)/Im_core(i,j,k,2,QRHO)**2
                der(:)  = er_ref(:)  - Im_rad(i,j,k,2,qrad:qradhi)

                dup    = un_ref    - Im_core(i,j,k,3,QUN) - hdt*Im_core_src(i,j,k,3,QUN)
                dptotp = ptot_ref - Im_rad(i,j,k,3,qptot) - hdt*Im_core_src(i,j,k,3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors -- NOT YET IMPLEMENTED

                if (ppm_predict_gammae == 0) then

                   ! (rho, u, p, (rho e) eigensystem

                   ! These are analogous to the beta's from the original PPM
                   ! paper (except we work with rho instead of tau).  This is
                   ! simply (l . dq), where dq = qref - I(q)

                   alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
                   alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
                   alpha0r = drho - dptot/csq
                   alpha0e_g = drhoe_g - dptot*h_g
                else

                   ! (tau, u, p, game) eigensystem

                   ! This is the way things were done in the original PPM
                   ! paper -- here we work with tau in the characteristic
                   ! system

                   alpham = HALF*( dum - dptotm*(ONE/Clag))*(ONE/Clag)
                   alphap = HALF*(-dup - dptotp*(ONE/Clag))*(ONE/Clag)
                   alpha0r = dtau + dptot*(ONE/Clag)**2

                   dge   = game_ref - Im_core(i,j,k,2,QGAME)
                   gfactor = (game - ONE)*(game - gam_g)
                   alpha0e_g = gfactor*dptot/(tau*Clag**2) + dge

                endif    ! which tracing method

                alphar(:) = der(:) - dptot/csq*hr

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
                   alphar(:) = ZERO
                else
                   alphar(:) = -alphar(:)
                end if

                if (un > ZERO) then
                   alpha0e_g = ZERO
                else
                   alpha0e_g = -alpha0e_g
                end if


                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                if (ppm_predict_gammae == 0) then
                   qp_core(i,j,k,QRHO) = rho_ref + alphap + alpham + alpha0r
                   qp_core(i,j,k,QUN) = un_ref + (alphap - alpham)*cc/rho
                   qp_core(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                   qp_core(i,j,k,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

                   qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                   qp_rad(i,j,k,qrad:qradhi) = qrtmp

                   qp_rad(i,j,k,qptot) = ptot_ref + (alphap + alpham)*csq
                   qp_rad(i,j,k,qreitot) = qp_core(i,j,k,QREINT) + sum(qrtmp)

                else
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qp_core(i,j,k,QRHO  ) = ONE/tau_s

                   qp_core(i,j,k,QUN    ) = un_ref + (alpham - alphap)*Clag
                   qp_core(i,j,k,QPRES ) = p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamp(:)*alphar(:))

                   qp_core(i,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                   qp_core(i,j,k,QREINT) = qp_core(i,j,k,QPRES )/(qp_core(i,j,k,QGAME) - ONE)

                   qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                   qp_rad(i,j,k,qrad:qradhi) = qrtmp

                   qp_rad(i,j,k,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                   qp_rad(i,j,k,qreitot) = qp_core(i,j,k,QREINT) + sum(qrtmp)

                endif

                ! Enforce small_*
                qp_core(i,j,k,QRHO) = max(qp_core(i,j,k,QRHO), small_dens)
                qp_core(i,j,k,QPRES) = max(qp_core(i,j,k,QPRES), small_pres)

                do g = 0, ngroups-1
                   if (qp_rad(i,j,k,qrad+g) < ZERO) then
                      er_foo = - qp_rad(i,j,k,qrad+g)
                      qp_rad(i,j,k,qrad+g) = ZERO
                      qp_rad(i,j,k,qptot) = qp_rad(i,j,k,qptot) + lamp(g) * er_foo
                      qp_rad(i,j,k,qreitot) = qp_rad(i,j,k,qreitot) + er_foo
                   end if
                end do


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

                !gam_g_ref  = Ip_gc(i,j,k,3,1)
                game_ref = Ip_core(i,j,k,3,QGAME)


                ptot_ref = Ip_rad(i,j,k,3,QPTOT)

                er_ref(:) = Ip_rad(i,j,k,3,QRAD:QRADHI)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum    = un_ref    - Ip_core(i,j,k,1,QUN) - hdt*Ip_core_src(i,j,k,1,QUN)
                dptotm = ptot_ref - Ip_rad(i,j,k,1,qptot) - hdt*Ip_core_src(i,j,k,1,QPRES)

                drho    = rho_ref    - Ip_core(i,j,k,2,QRHO) - hdt*Ip_core_src(i,j,k,2,QRHO)
                dptot   = ptot_ref   - Ip_rad(i,j,k,2,qptot) - hdt*Ip_core_src(i,j,k,2,QPRES)
                drhoe_g = rhoe_g_ref - Ip_core(i,j,k,2,QREINT) - hdt*Ip_core_src(i,j,k,2,QREINT)
                dtau  = tau_ref  - ONE/Ip_core(i,j,k,2,QRHO) + hdt*Ip_core_src(i,j,k,2,QRHO)/Ip_core(i,j,k,2,QRHO)**2
                der(:)  = er_ref(:)  - Ip_rad(i,j,k,2,qrad:qradhi)

                dup    = un_ref    - Ip_core(i,j,k,3,QUN) - hdt*Ip_core_src(i,j,k,3,QUN)
                dptotp = ptot_ref - Ip_rad(i,j,k,3,qptot) - hdt*Ip_core_src(i,j,k,3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors -- NOT YET IMPLEMENTED

                if (ppm_predict_gammae == 0) then

                   ! (rho, u, p, (rho e)) eigensystem

                   ! These are analogous to the beta's from the original PPM
                   ! paper (except we work with rho instead of tau).  This is
                   ! simply (l . dq), where dq = qref - I(q)

                   alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
                   alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
                   alpha0r = drho - dptot/csq
                   alpha0e_g = drhoe_g - dptot*h_g

                else

                   ! (tau, u, p, game) eigensystem

                   ! This is the way things were done in the original PPM
                   ! paper -- here we work with tau in the characteristic
                   ! system
                   alpham = HALF*( dum - dptotm*(ONE/Clag))*(ONE/Clag)
                   alphap = HALF*(-dup - dptotp*(ONE/Clag))*(ONE/Clag)
                   alpha0r = dtau + dptot*(ONE/Clag)**2

                   dge = game_ref - Ip_core(i,j,k,2,QGAME)
                   gfactor = (game - ONE)*(game - gam_g)
                   alpha0e_g = gfactor*dptot/(tau*Clag**2) + dge

                endif

                alphar(:) = der(:) - dptot/csq*hr

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
                   alphar(:) = -alphar(:)
                else
                   alphar(:) = ZERO
                end if

                if (un > ZERO) then
                   alpha0e_g = -alpha0e_g
                else
                   alpha0e_g = ZERO
                end if

                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                if (ppm_predict_gammae == 0) then

                   if (idir == 1) then
                      qm_core(i+1,j,k,QRHO) = max(small_dens, rho_ref + alphap + alpham + alpha0r)
                      qm_core(i+1,j,k,QUN) = un_ref + (alphap - alpham)*cc/rho
                      qm_core(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                      qm_core(i+1,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:)))

                      qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                      qm_rad(i+1,j,k,qrad:qradhi) = qrtmp

                      qm_rad(i+1,j,k,qptot) = ptot_ref + (alphap + alpham)*csq
                      qm_rad(i+1,j,k,qreitot) = qm_core(i+1,j,k,QREINT) + sum(qrtmp)

                   else if (idir == 2) then
                      qm_core(i,j+1,k,QRHO) = max(small_dens, rho_ref + alphap + alpham + alpha0r)
                      qm_core(i,j+1,k,QUN) = un_ref + (alphap - alpham)*cc/rho
                      qm_core(i,j+1,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                      qm_core(i,j+1,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:)))

                      qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                      qm_rad(i,j+1,k,qrad:qradhi) = qrtmp

                      qm_rad(i,j+1,k,qptot) = ptot_ref + (alphap + alpham)*csq
                      qm_rad(i,j+1,k,qreitot) = qm_core(i,j+1,k,QREINT) + sum(qrtmp)

                   else if (idir == 3) then
                      qm_core(i,j,k+1,QRHO) = max(small_dens, rho_ref + alphap + alpham + alpha0r)
                      qm_core(i,j,k+1,QUN) = un_ref + (alphap - alpham)*cc/rho
                      qm_core(i,j,k+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                      qm_core(i,j,k+1,QPRES) = max(small_pres, p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:)))

                      qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                      qm_rad(i,j,k+1,qrad:qradhi) = qrtmp

                      qm_rad(i,j,k+1,qptot) = ptot_ref + (alphap + alpham)*csq
                      qm_rad(i,j,k+1,qreitot) = qm_core(i,j,k+1,QREINT) + sum(qrtmp)

                   end if

                else

                   if (idir == 1) then
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qm_core(i+1,j,k,QRHO  ) = max(small_dens, ONE/tau_s)

                      qm_core(i+1,j,k,QUN    ) = un_ref + (alpham - alphap)*Clag
                      qm_core(i+1,j,k,QPRES ) = max(small_pres, p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamm(:)*alphar(:)))

                      qm_core(i+1,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                      qm_core(i+1,j,k,QREINT) = qm_core(i+1,j,k,QPRES )/(qm_core(i+1,j,k,QGAME) - ONE)

                      qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                      qm_rad(i+1,j,k,qrad:qradhi) = qrtmp

                      qm_rad(i+1,j,k,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                      qm_rad(i+1,j,k,qreitot) = qm_core(i+1,j,k,QREINT) + sum(qrtmp)

                   else if (idir == 2) then
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qm_core(i,j+1,k,QRHO  ) = max(small_dens, ONE/tau_s)

                      qm_core(i,j+1,k,QUN    ) = un_ref + (alpham - alphap)*Clag
                      qm_core(i,j+1,k,QPRES ) = max(small_pres, p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamm(:)*alphar(:)))

                      qm_core(i,j+1,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                      qm_core(i,j+1,k,QREINT) = qm_core(i,j+1,k,QPRES )/(qm_core(i,j+1,k,QGAME) - ONE)

                      qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                      qm_rad(i,j+1,k,qrad:qradhi) = qrtmp

                      qm_rad(i,j+1,k,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                      qm_rad(i,j+1,k,qreitot) = qm_core(i,j+1,k,QREINT) + sum(qrtmp)

                   else if (idir == 3) then
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qm_core(i,j,k+1,QRHO  ) = max(small_dens, ONE/tau_s)

                      qm_core(i,j,k+1,QUN    ) = un_ref + (alpham - alphap)*Clag
                      qm_core(i,j,k+1,QPRES ) = max(small_pres, p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamm(:)*alphar(:)))

                      qm_core(i,j,k+1,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                      qm_core(i,j,k+1,QREINT) = qm_core(i,j,k+1,QPRES )/(qm_core(i,j,k+1,QGAME) - ONE)

                      qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                      qm_rad(i,j,k+1,qrad:qradhi) = qrtmp

                      qm_rad(i,j,k+1,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                      qm_rad(i,j,k+1,qreitot) = qm_core(i,j,k+1,QREINT) + sum(qrtmp)

                   end if

                endif

                if (idir == 1) then
                   do g=0,ngroups-1
                      if (qm_rad(i+1,j,k,qrad+g) < ZERO) then
                         er_foo = - qm_rad(i+1,j,k,qrad+g)
                         qm_rad(i+1,j,k,qrad+g) = ZERO
                         qm_rad(i+1,j,k,qptot) = qm_rad(i+1,j,k,qptot) + lamm(g) * er_foo
                         qm_rad(i+1,j,k,qreitot) = qm_rad(i+1,j,k,qreitot) + er_foo
                      end if
                   end do

                   ! transverse velocities
                   qm_core(i+1,j,k,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i+1,j,k,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                else if (idir == 2) then
                   do g=0,ngroups-1
                      if (qm_rad(i,j+1,k,qrad+g) < ZERO) then
                         er_foo = - qm_rad(i,j+1,k,qrad+g)
                         qm_rad(i,j+1,k,qrad+g) = ZERO
                         qm_rad(i,j+1,k,qptot) = qm_rad(i,j+1,k,qptot) + lamm(g) * er_foo
                         qm_rad(i,j+1,k,qreitot) = qm_rad(i,j+1,k,qreitot) + er_foo
                      end if
                   end do

                   ! transverse velocities
                   qm_core(i,j+1,k,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i,j+1,k,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                else if (idir == 3) then
                   do g=0,ngroups-1
                      if (qm_rad(i,j,k+1,qrad+g) < ZERO) then
                         er_foo = - qm_rad(i,j,k+1,qrad+g)
                         qm_rad(i,j,k+1,qrad+g) = ZERO
                         qm_rad(i,j,k+1,qptot) = qm_rad(i,j,k+1,qptot) + lamm(g) * er_foo
                         qm_rad(i,j,k+1,qreitot) = qm_rad(i,j,k+1,qreitot) + er_foo
                      end if
                   end do

                   ! transverse velocities
                   qm_core(i,j,k+1,QUT) = Ip_core(i,j,k,2,QUT) + hdt*Ip_core_src(i,j,k,2,QUT)
                   qm_core(i,j,k+1,QUTT) = Ip_core(i,j,k,2,QUTT) + hdt*Ip_core_src(i,j,k,2,QUTT)

                end if

             end if


             !-------------------------------------------------------------------
             ! geometry source terms
             !-------------------------------------------------------------------

#if (AMREX_SPACEDIM < 3)
             if (idir == 1 .and. dloga(i,j,k) /= 0) then
                courn = dt/dx(1)*(cc+abs(un))
                eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k)))
                dlogatmp = min(eta, ONE)*dloga(i,j,k)
                sourcr = -HALF*dt*rho*dlogatmp*un
                sourcp = sourcr*cgassq
                source = sourcp*h_g
                sourcer(:) = -HALF*dt*dlogatmp*un*(lam0(:)+ONE)*er(:)

                if (i <= vhi(1)) then
                   qm_core(i+1,j,k,QRHO  ) = qm_core(i+1,j,k,QRHO  ) + sourcr
                   qm_core(i+1,j,k,QRHO  ) = max(qm_core(i+1,j,k,QRHO), small_dens)
                   qm_core(i+1,j,k,QPRES ) = qm_core(i+1,j,k,QPRES ) + sourcp
                   qm_core(i+1,j,k,QREINT) = qm_core(i+1,j,k,QREINT) + source
                   qm_rad(i+1,j,k,qrad:qradhi) = qm_rad(i+1,j,k,qrad:qradhi) + sourcer(:)
                   ! qm(i+1,j,k,qptot ) = sum(lamm(:)*qm(i+1,j,k,qrad:qradhi)) + qm(i+1,j,k,QPRES)
                   qm_rad(i+1,j,k,qptot) = qm_rad(i+1,j,k,qptot) + sum(lamm(:)*sourcer(:)) + sourcp
                   qm_rad(i+1,j,k,qreitot) = sum(qm_rad(i+1,j,k,qrad:qradhi))  + qm_core(i+1,j,k,QREINT)
                end if

                if (i >= vlo(1)) then
                   qp_core(i,j,k,QRHO  ) = qp_core(i,j,k,QRHO  ) + sourcr
                   qp_core(i,j,k,QRHO  ) = max(qp_core(i,j,k,QRHO), small_dens)
                   qp_core(i,j,k,QPRES ) = qp_core(i,j,k,QPRES ) + sourcp
                   qp_core(i,j,k,QREINT) = qp_core(i,j,k,QREINT) + source
                   qp_rad(i,j,k,qrad:qradhi) = qp_rad(i,j,k,qrad:qradhi) + sourcer(:)
                   ! qp(i  ,qptot ) = sum(lamp(:)*qp(i,qrad:qradhi)) + qp(i,QPRES)
                   qp_rad(i,j,k,qptot) = qp_rad(i,j,k,qptot) + sum(lamp(:)*sourcer(:)) + sourcp
                   qp_rad(i,j,k,qreitot) = sum(qp_rad(i,j,k,qrad:qradhi))  + qp_core(i,j,k,QREINT)
                end if
             endif
#endif

          end do
       end do
    end do


    !-------------------------------------------------------------------------
    ! passively advected quantities
    !-------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
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

  end subroutine trace_ppm_rad

end module trace_ppm_rad_module
