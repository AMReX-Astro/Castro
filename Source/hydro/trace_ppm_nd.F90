! These routines do the characteristic tracing under the parabolic
! profiles in each zone to the edge / half-time.

module trace_ppm_module

  use prob_params_module, only : dg
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public tracexy_ppm, tracez_ppm

contains

  subroutine tracexy_ppm(q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                         qxm, qxp, qym, qyp, qs_lo, qs_hi, &
#if (AMREX_SPACEDIM < 3)
                         dloga, dloga_lo, dloga_hi, &
#endif
                         lo, hi, domlo, domhi, &
                         dx, dt)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QTEMP, QGAME, QC, QGAMC, QFX, QFS, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   ppm_reference_eigenvectors, ppm_predict_gammae, &
                                   npassive, qpass_map, ppm_temp_fix, &
                                   fix_mass_flux
    use amrex_constants_module, only : ZERO, HALF, ONE
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qs_lo(3),qs_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: I_lo(3), I_hi(3)
#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,NQ)
    real(rt), intent(in) :: Im(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,NQ)

    real(rt), intent(in) :: Ip_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,QVAR)
    real(rt), intent(in) :: Im_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,QVAR)

    real(rt), intent(in) :: Ip_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,1)
    real(rt), intent(in) :: Im_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,1)

    real(rt), intent(inout) :: qxm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qxp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qym(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qyp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j, k
    integer :: n, ipassive

    type(eos_t) :: eos_state

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
    real(rt) :: rho, u, v, w, p, rhoe_g, h_g, temp
    real(rt) :: gam_g, game

    real(rt) :: drho, dptot, drhoe_g
    real(rt) :: de, dge, dtau, dtaum, dtaup
    real(rt) :: dup, dvp, dptotp
    real(rt) :: dum, dvm, dptotm
    real(rt) :: dT0, dTp, dTm
    real(rt) :: p_r, p_T

    real(rt) :: rho_ref, u_ref, v_ref, p_ref, rhoe_g_ref, h_g_ref, temp_ref
    real(rt) :: tau_ref

    real(rt) :: cc_ref, csq_ref, Clag_ref, gam_g_ref, game_ref, gfactor
    real(rt) :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, h_g_ev, tau_ev, temp_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp
    real(rt) :: tau_s, e_s

    logical :: fix_mass_flux_lo, fix_mass_flux_hi

#ifndef AMREX_USE_CUDA    
    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    hdt = HALF * dt


    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (lo(1) == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (hi(1) == domhi(1))


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


    !-------------------------------------------------------------------------
    ! x-direction
    !-------------------------------------------------------------------------

    ! Trace to left and right edges using upwind PPM
    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-1, hi(1)+1

             gfactor = ONE ! to help compiler resolve ANTI dependence

             rho = q(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2
             Clag = rho*cc

             u = q(i,j,k,QU)
             v = q(i,j,k,QV)
             w = q(i,j,k,QW)

             p = q(i,j,k,QPRES)
             rhoe_g = q(i,j,k,QREINT)
             h_g = ( (p + rhoe_g)/rho)/csq
             temp = q(i,j,k,QTEMP)

             gam_g = qaux(i,j,k,QGAMC)
             game = q(i,j,k,QGAME)


             !-------------------------------------------------------------------
             ! plus state on face i
             !-------------------------------------------------------------------

             if (i >= lo(1)) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(i,j,k,1,1,QRHO)
                u_ref    = Im(i,j,k,1,1,QU)

                p_ref    = Im(i,j,k,1,1,QPRES)
                rhoe_g_ref = Im(i,j,k,1,1,QREINT)
                temp_ref = Im(i,j,k,1,1,QTEMP)

                tau_ref  = ONE/Im(i,j,k,1,1,QRHO)

                gam_g_ref  = Im_gc(i,j,k,1,1,1)
                game_ref = Im(i,j,k,1,1,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dum = u_ref - Im(i,j,k,1,1,QU) - hdt*Im_src(i,j,k,1,1,QU)
                dptotm = p_ref - Im(i,j,k,1,1,QPRES) - hdt*Im_src(i,j,k,1,1,QPRES)

                drho = rho_ref - Im(i,j,k,1,2,QRHO) - hdt*Im_src(i,j,k,1,2,QRHO)
                dptot = p_ref - Im(i,j,k,1,2,QPRES) - hdt*Im_src(i,j,k,1,2,QPRES)
                drhoe_g = rhoe_g_ref - Im(i,j,k,1,2,QREINT) - hdt*Im_src(i,j,k,1,2,QREINT)

                ! TODO: need to figure sources for this out...
                dTm = temp_ref - Im(i,j,k,1,1,QTEMP)
                dT0 = temp_ref - Im(i,j,k,1,2,QTEMP)
                dTp = temp_ref - Im(i,j,k,1,3,QTEMP)

                ! we are treating tau as 1/rho, but we could have reconstructed
                ! it separately
                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Im(i,j,k,1,1,QRHO) + hdt*Im_src(i,j,k,1,1,QRHO)/Im(i,j,k,1,1,QRHO)**2
                dtau  = tau_ref - ONE/Im(i,j,k,1,2,QRHO) + hdt*Im_src(i,j,k,1,2,QRHO)/Im(i,j,k,1,2,QRHO)**2
                dtaup = tau_ref - ONE/Im(i,j,k,1,3,QRHO) + hdt*Im_src(i,j,k,1,3,QRHO)/Im(i,j,k,1,3,QRHO)**2

                dup = u_ref - Im(i,j,k,1,3,QU) - hdt*Im_src(i,j,k,1,3,QU)
                dptotp = p_ref - Im(i,j,k,1,3,QPRES) - hdt*Im_src(i,j,k,1,3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   Clag_ev = Clag
                   h_g_ev = h_g
                   p_ev    = p
                   tau_ev  = ONE/rho
                   temp_ev = temp
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   Clag_ev = Clag_ref
                   h_g_ev = h_g_ref
                   p_ev    = p_ref
                   tau_ev  = tau_ref
                   temp_ev = temp_ref
                endif

                if (ppm_temp_fix < 3) then

                   ! we are generally working with p in our eigensystem

                   if (ppm_predict_gammae == 0) then

                      ! (rho, u, p, (rho e) eigensystem

                      ! These are analogous to the beta's from the original PPM
                      ! paper (except we work with rho instead of tau).  This is
                      ! simply (l . dq), where dq = qref - I(q)

                      alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dum)*(rho_ev/cc_ev)
                      alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dup)*(rho_ev/cc_ev)
                      alpha0r = drho - dptot/csq_ev
                      alpha0e_g = drhoe_g - dptot*h_g_ev  ! note h_g has a 1/c**2 in it

                   else

                      ! (tau, u, p, game) eigensystem

                      ! This is the way things were done in the original PPM
                      ! paper -- here we work with tau in the characteristic
                      ! system

                      alpham = HALF*( dum - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alphap = HALF*(-dup - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                      dge   = game_ref - Im(i,j,k,1,2,QGAME)
                      gfactor = (game - ONE)*(game - gam_g)
                      alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                   endif

                else

                   ! (tau, u T) eigensystem

                   ! eos to get some thermodynamics
                   eos_state%T = temp_ev
                   eos_state%rho = rho_ev
                   eos_state%xn(:) = q(i,j,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = q(i,j,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   p_r = eos_state%dpdr
                   p_T = eos_state%dpdT

                   alpham = HALF*(rho_ev**2*p_r*dtaum/Clag_ev + dum - p_T*dTm/Clag_ev)/Clag_ev
                   alphap = HALF*(rho_ev**2*p_r*dtaup/Clag_ev - dup - p_T*dTp/Clag_ev)/Clag_ev
                   alpha0r = dtau + (-rho_ev**2*p_r*dtau + p_T*dT0)/Clag_ev**2

                   ! not used, but needed to prevent bad invalid ops
                   alpha0e_g = ZERO

                endif    ! which tracing method

                if (u-cc > ZERO) then
                   alpham = ZERO
                else if (u-cc < ZERO) then
                   alpham = -alpham
                else
                   alpham = -HALF*alpham
                endif

                if (u+cc > ZERO) then
                   alphap = ZERO
                else if (u+cc < ZERO) then
                   alphap = -alphap
                else
                   alphap = -HALF*alphap
                endif

                if (u > ZERO) then
                   alpha0r = ZERO
                   alpha0e_g = ZERO
                else if (u < ZERO) then
                   alpha0r = -alpha0r
                   alpha0e_g = -alpha0e_g
                else
                   alpha0r = -HALF*alpha0r
                   alpha0e_g = -HALF*alpha0e_g
                endif

                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                if (ppm_temp_fix < 3) then

                   if (ppm_predict_gammae == 0) then
                      qxp(i,j,k,QRHO  ) =  rho_ref +  alphap + alpham + alpha0r
                      qxp(i,j,k,QU    ) =    u_ref + (alphap - alpham)*cc_ev/rho_ev
                      qxp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                      qxp(i,j,k,QPRES ) =    p_ref + (alphap + alpham)*csq_ev

                   else
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qxp(i,j,k,QRHO  ) = ONE/tau_s

                      qxp(i,j,k,QU    ) = u_ref + (alpham - alphap)*Clag_ev
                      qxp(i,j,k,QPRES ) = p_ref - (alphap + alpham)*Clag_ev**2

                      qxp(i,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                      qxp(i,j,k,QREINT) = qxp(i,j,k,QPRES )/(qxp(i,j,k,QGAME) - ONE)
                   endif

                else
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qxp(i,j,k,QRHO ) = ONE/tau_s

                   qxp(i,j,k,QU   ) = u_ref + (alpham - alphap)*Clag_ev
                   qxp(i,j,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once we do the species
                   qxp(i,j,k,QPRES) = small_pres ! just to make it defined
                endif

                ! Enforce small_*
                qxp(i,j,k,QRHO ) = max(qxp(i,j,k,QRHO ), small_dens)
                qxp(i,j,k,QPRES) = max(qxp(i,j,k,QPRES), small_pres)

                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qxp(i,j,k,QV) = Im(i,j,k,1,2,QV) + hdt*Im_src(i,j,k,1,2,QV)
                qxp(i,j,k,QW) = Im(i,j,k,1,2,QW) + hdt*Im_src(i,j,k,1,2,QW)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if (i <= hi(1)) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,k,1,3,QRHO)
                u_ref    = Ip(i,j,k,1,3,QU)

                p_ref    = Ip(i,j,k,1,3,QPRES)
                rhoe_g_ref = Ip(i,j,k,1,3,QREINT)
                temp_ref = Ip(i,j,k,1,3,QTEMP)

                tau_ref  = ONE/Ip(i,j,k,1,3,QRHO)

                gam_g_ref  = Ip_gc(i,j,k,1,3,1)
                game_ref = Ip(i,j,k,1,3,QGAME)

                rho_ref = max(rho_ref,small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = u_ref - Ip(i,j,k,1,1,QU) - hdt*Ip_src(i,j,k,1,1,QU)
                dptotm  = p_ref - Ip(i,j,k,1,1,QPRES) - hdt*Ip_src(i,j,k,1,1,QPRES)

                drho = rho_ref - Ip(i,j,k,1,2,QRHO) - hdt*Ip_src(i,j,k,1,2,QRHO)
                dptot = p_ref - Ip(i,j,k,1,2,QPRES) - hdt*Ip_src(i,j,k,1,2,QPRES)
                drhoe_g = rhoe_g_ref - Ip(i,j,k,1,2,QREINT) - hdt*Ip_src(i,j,k,1,2,QREINT)

                dTm = temp_ref - Ip(i,j,k,1,1,QTEMP)
                dT0 = temp_ref - Ip(i,j,k,1,2,QTEMP)
                dTp = temp_ref - Ip(i,j,k,1,3,QTEMP)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Ip(i,j,k,1,1,QRHO) + hdt*Ip_src(i,j,k,1,1,QRHO)/Ip(i,j,k,1,1,QRHO)**2
                dtau = tau_ref - ONE/Ip(i,j,k,1,2,QRHO) + hdt*Ip_src(i,j,k,1,2,QRHO)/Ip(i,j,k,1,2,QRHO)**2
                dtaup = tau_ref - ONE/Ip(i,j,k,1,3,QRHO) + hdt*Ip_src(i,j,k,1,3,QRHO)/Ip(i,j,k,1,3,QRHO)**2

                dup = u_ref - Ip(i,j,k,1,3,QU) - hdt*Ip_src(i,j,k,1,3,QU)
                dptotp = p_ref - Ip(i,j,k,1,3,QPRES) - hdt*Ip_src(i,j,k,1,3,QPRES)

                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   Clag_ev = Clag
                   h_g_ev = h_g
                   p_ev    = p
                   tau_ev  = ONE/rho
                   temp_ev = temp
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   Clag_ev = Clag_ref
                   h_g_ev = h_g_ref
                   p_ev    = p_ref
                   tau_ev  = tau_ref
                   temp_ev = temp_ref
                endif

                if (ppm_temp_fix < 3) then

                   if (ppm_predict_gammae == 0) then

                      ! (rho, u, p, (rho e)) eigensystem

                      ! These are analogous to the beta's from the original PPM
                      ! paper (except we work with rho instead of tau).  This is
                      ! simply (l . dq), where dq = qref - I(q)

                      alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dum)*(rho_ev/cc_ev)
                      alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dup)*(rho_ev/cc_ev)
                      alpha0r = drho - dptot/csq_ev
                      alpha0e_g = drhoe_g - dptot*h_g_ev  ! h_g has a 1/c**2 in it

                   else

                      ! (tau, u, p, game) eigensystem

                      ! This is the way things were done in the original PPM
                      ! paper -- here we work with tau in the characteristic
                      ! system
                      alpham = HALF*( dum - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alphap = HALF*(-dup - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                      dge = game_ref - Ip(i,j,k,1,2,QGAME)
                      gfactor = (game - ONE)*(game - gam_g)
                      alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                   endif

                else

                   ! (tau, u T) eigensystem

                   ! eos to get some thermodynamics
                   eos_state%T = temp_ev
                   eos_state%rho = rho_ev
                   eos_state%xn(:) = q(i,j,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = q(i,j,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   p_r = eos_state%dpdr
                   p_T = eos_state%dpdT

                   alpham = HALF*(rho_ev**2*p_r*dtaum/Clag_ev + dum - p_T*dTm/Clag_ev)/Clag_ev
                   alphap = HALF*(rho_ev**2*p_r*dtaup/Clag_ev - dup - p_T*dTp/Clag_ev)/Clag_ev
                   alpha0r = dtau + (-rho_ev**2*p_r*dtau + p_T*dT0)/Clag_ev**2

                   ! not used, but needed to prevent bad invalid ops
                   alpha0e_g = ZERO

                endif

                if (u-cc > ZERO) then
                   alpham = -alpham
                else if (u-cc < ZERO) then
                   alpham = ZERO
                else
                   alpham = -HALF*alpham
                endif

                if (u+cc > ZERO) then
                   alphap = -alphap
                else if (u+cc < ZERO) then
                   alphap = ZERO
                else
                   alphap = -HALF*alphap
                endif

                if (u > ZERO) then
                   alpha0r = -alpha0r
                   alpha0e_g = -alpha0e_g
                else if (u < ZERO) then
                   alpha0r = ZERO
                   alpha0e_g = ZERO
                else
                   alpha0r = -HALF*alpha0r
                   alpha0e_g = -HALF*alpha0e_g
                endif


                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                if (ppm_temp_fix < 3) then

                   if (ppm_predict_gammae == 0) then
                      qxm(i+1,j,k,QRHO  ) =  rho_ref +  alphap + alpham + alpha0r
                      qxm(i+1,j,k,QU    ) =    u_ref + (alphap - alpham)*cc_ev/rho_ev
                      qxm(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                      qxm(i+1,j,k,QPRES ) =    p_ref + (alphap + alpham)*csq_ev
                   else
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qxm(i+1,j,k,QRHO  ) = ONE/tau_s

                      qxm(i+1,j,k,QU    ) = u_ref + (alpham - alphap)*Clag_ev
                      qxm(i+1,j,k,QPRES ) = p_ref - (alphap + alpham)*Clag_ev**2

                      qxm(i+1,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                      qxm(i+1,j,k,QREINT) = qxm(i+1,j,k,QPRES )/(qxm(i+1,j,k,QGAME) - ONE)
                   endif

                else
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qxm(i+1,j,k,QRHO ) = ONE/tau_s

                   qxm(i+1,j,k,QU   ) = u_ref + (alpham - alphap)*Clag_ev
                   qxm(i+1,j,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once we do the species
                   qxm(i+1,j,k,QPRES) = small_pres ! just to make it defined
                endif

                ! Enforce small_*
                qxm(i+1,j,k,QRHO ) = max(qxm(i+1,j,k,QRHO ),small_dens)
                qxm(i+1,j,k,QPRES) = max(qxm(i+1,j,k,QPRES),small_pres)

                ! transverse velocities
                qxm(i+1,j,k,QV    ) = Ip(i,j,k,1,2,QV) + hdt*Ip_src(i,j,k,1,2,QV)
                qxm(i+1,j,k,QW    ) = Ip(i,j,k,1,2,QW) + hdt*Ip_src(i,j,k,1,2,QW)

             end if

             !-------------------------------------------------------------------
             ! geometry source terms
             !-------------------------------------------------------------------

#if (AMREX_SPACEDIM < 3)
             if (dloga(i,j,k) /= 0) then
                courn = dt/dx(1)*(cc+abs(u))
                eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k)))
                dlogatmp = min(eta, ONE)*dloga(i,j,k)
                sourcr = -HALF*dt*rho*dlogatmp*u
                sourcp = sourcr*csq
                source = sourcp*h_g

                if (i <= hi(1)) then
                   qxm(i+1,j,k,QRHO) = qxm(i+1,j,k,QRHO) + sourcr
                   qxm(i+1,j,k,QRHO) = max(qxm(i+1,j,k,QRHO), small_dens)
                   qxm(i+1,j,k,QPRES) = qxm(i+1,j,k,QPRES) + sourcp
                   qxm(i+1,j,k,QREINT) = qxm(i+1,j,k,QREINT) + source
                end if

                if (i >= lo(1)) then
                   qxp(i,j,k,QRHO) = qxp(i,j,k,QRHO) + sourcr
                   qxp(i,j,k,QRHO) = max(qxp(i,j,k,QRHO), small_dens)
                   qxp(i,j,k,QPRES) = qxp(i,j,k,QPRES) + sourcp
                   qxp(i,j,k,QREINT) = qxp(i,j,k,QREINT) + source
                end if

             endif
#endif

#if (AMREX_SPACEDIM == 1)
             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_lo) then
                qxm(lo(1),j,k,QRHO  ) = q(domlo(1)-1,j,k,QRHO)
                qxm(lo(1),j,k,QU    ) = q(domlo(1)-1,j,k,QU  )
                qxm(lo(1),j,k,QPRES ) = q(domlo(1)-1,j,k,QPRES)
                qxm(lo(1),j,k,QREINT) = q(domlo(1)-1,j,k,QREINT)
             end if

             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_hi) then
                qxp(hi(1)+1,j,k,QRHO  ) = q(domhi(1)+1,j,k,QRHO)
                qxp(hi(1)+1,j,k,QU    ) = q(domhi(1)+1,j,k,QU  )
                qxp(hi(1)+1,j,k,QPRES ) = q(domhi(1)+1,j,k,QPRES)
                qxp(hi(1)+1,j,k,QREINT) = q(domhi(1)+1,j,k,QREINT)
             end if
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

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  But we already dealt with all 3 velocity
       ! components above, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2)-dg(2), hi(2)+dg(2)

             ! Plus state on face i
             do i = lo(1), hi(1)+1
                u = q(i,j,k,QU)

                ! We have
                !
                ! q_l = q_ref - Proj{(q_ref - I)}
                !
                ! and Proj{} represents the characteristic projection.
                ! But for these, there is only 1-wave that matters, the u
                ! wave, so no projection is needed.  Since we are not
                ! projecting, the reference state doesn't matter

                if (u > ZERO) then
                   qxp(i,j,k,n) = q(i,j,k,n)
                else if (u < ZERO) then
                   qxp(i,j,k,n) = Im(i,j,k,1,2,n)
                else
                   qxp(i,j,k,n) = q(i,j,k,n) + HALF*(Im(i,j,k,1,2,n) - q(i,j,k,n))
                endif
             enddo

             ! Minus state on face i+1
             do i = lo(1)-1, hi(1)
                u = q(i,j,k,QU)

                if (u > ZERO) then
                   qxm(i+1,j,k,n) = Ip(i,j,k,1,2,n)
                else if (u < ZERO) then
                   qxm(i+1,j,k,n) = q(i,j,k,n)
                else
                   qxm(i+1,j,k,n) = q(i,j,k,n) + HALF*(Ip(i,j,k,1,2,n) - q(i,j,k,n))
                endif
             end do

#if AMREX_SPACEDIM == 1
             if (fix_mass_flux_hi) qxp(hi(1)+1,j,k,n) = q(hi(1)+1,j,k,n)
             if (fix_mass_flux_lo) qxm(lo(1),j,k,n) = q(lo(1)-1,j,k,n)
#endif

          end do
       end do
    end do

    if (ppm_temp_fix == 3) then

       ! we predicted T, now make p, (rho e) consistent
       do k = lo(3)-dg(3), hi(3)-dg(3)
          do j = lo(2)-dg(2), hi(2)+dg(2)
             do i = lo(1)-1, hi(1)+1

                if (i >= lo(1)) then
                   ! plus face
                   eos_state%T     = qxp(i,j,k,QTEMP)
                   eos_state%rho   = qxp(i,j,k,QRHO)
                   eos_state%xn(:) = qxp(i,j,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = qxp(i,j,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   qxp(i,j,k,QPRES) = eos_state%p
                   qxp(i,j,k,QREINT) = qxp(i,j,k,QRHO)*eos_state%e

                   qxp(i,j,k,QPRES) = max(qxp(i,j,k,QPRES), small_pres)
                endif

                if (i <= hi(1)) then
                   ! minus face
                   eos_state%T     = qxm(i+1,j,k,QTEMP)
                   eos_state%rho   = qxm(i+1,j,k,QRHO)
                   eos_state%xn(:) = qxm(i+1,j,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = qxm(i+1,j,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   qxm(i+1,j,k,QPRES) = eos_state%p
                   qxm(i+1,j,k,QREINT) = qxm(i+1,j,k,QRHO)*eos_state%e

                   qxm(i+1,j,k,QPRES) = max(qxm(i+1,j,k,QPRES), small_pres)
                endif

             end do
          end do
       end do
    endif


#if (AMREX_SPACEDIM >= 2)
    !-------------------------------------------------------------------------
    ! y-direction
    !-------------------------------------------------------------------------

    ! Trace to bottom and top edges using upwind PPM

    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             gfactor = ONE ! to help compiler resolve ANTI dependence

             rho = q(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2
             Clag = rho*cc

             u = q(i,j,k,QU)
             v = q(i,j,k,QV)
             w = q(i,j,k,QW)

             p = q(i,j,k,QPRES)
             rhoe_g = q(i,j,k,QREINT)
             h_g = ( (p + rhoe_g)/rho )/csq
             temp = q(i,j,k,QTEMP)

             gam_g = qaux(i,j,k,QGAMC)
             game = q(i,j,k,QGAME)

             !-------------------------------------------------------------------
             ! plus state on face j
             !-------------------------------------------------------------------

             if (j >= lo(2)) then

                ! Set the reference state
                ! This will be the fastest moving state to the left
                rho_ref  = Im(i,j,k,2,1,QRHO)
                v_ref    = Im(i,j,k,2,1,QV)

                p_ref    = Im(i,j,k,2,1,QPRES)
                rhoe_g_ref = Im(i,j,k,2,1,QREINT)
                temp_ref = Im(i,j,k,2,1,QTEMP)

                tau_ref  = ONE/Im(i,j,k,2,1,QRHO)

                gam_g_ref  = Im_gc(i,j,k,2,1,1)
                game_ref = Im(i,j,k,2,1,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! *m are the jumps carried by v-c
                ! *p are the jumps carried by v+c

                dvm = v_ref - Im(i,j,k,2,1,QV) - hdt*Im_src(i,j,k,2,1,QV)
                dptotm = p_ref - Im(i,j,k,2,1,QPRES) - hdt*Im_src(i,j,k,2,1,QPRES)

                drho = rho_ref - Im(i,j,k,2,2,QRHO) - hdt*Im_src(i,j,k,2,2,QRHO)
                dptot = p_ref - Im(i,j,k,2,2,QPRES) - hdt*Im_src(i,j,k,2,2,QPRES)
                drhoe_g = rhoe_g_ref - Im(i,j,k,2,2,QREINT) - hdt*Im_src(i,j,k,2,2,QREINT)

                dTm = temp_ref - Im(i,j,k,2,1,QTEMP)
                dT0 = temp_ref - Im(i,j,k,2,2,QTEMP)
                dTp = temp_ref - Im(i,j,k,2,3,QTEMP)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Im(i,j,k,2,1,QRHO) + hdt*Im_src(i,j,k,2,1,QRHO)/Im(i,j,k,2,1,QRHO)**2
                dtau = tau_ref - ONE/Im(i,j,k,2,2,QRHO) + hdt*Im_src(i,j,k,2,2,QRHO)/Im(i,j,k,2,2,QRHO)**2
                dtaup = tau_ref - ONE/Im(i,j,k,2,3,QRHO) + hdt*Im_src(i,j,k,2,2,QRHO)/Im(i,j,k,2,3,QRHO)**2

                dvp = v_ref - Im(i,j,k,2,3,QV) - hdt*Im_src(i,j,k,2,3,QV)
                dptotp = p_ref - Im(i,j,k,2,3,QPRES) - hdt*Im_src(i,j,k,2,3,QPRES)

                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   Clag_ev = Clag
                   h_g_ev = h_g
                   p_ev    = p
                   tau_ev  = ONE/rho
                   temp_ev = temp
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   Clag_ev = Clag_ref
                   h_g_ev = h_g_ref
                   p_ev    = p_ref
                   tau_ev  = tau_ref
                   temp_ev = temp_ref
                endif

                if (ppm_temp_fix < 3) then

                   if (ppm_predict_gammae == 0) then

                      ! (rho, u, p, (rho e) eigensystem

                      ! These are analogous to the beta's from the original PPM
                      ! paper (except we work with rho instead of tau).  This
                      ! is simply (l . dq), where dq = qref - I(q)

                      alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dvm)*(rho_ev/cc_ev)
                      alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dvp)*(rho_ev/cc_ev)
                      alpha0r = drho - dptot/csq_ev
                      alpha0e_g = drhoe_g - dptot*h_g_ev

                   else

                      ! (tau, u, p, game) eigensystem

                      ! This is the way things were done in the original PPM
                      ! paper -- here we work with tau in the characteristic
                      ! system

                      alpham = HALF*( dvm - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alphap = HALF*(-dvp - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                      dge = game_ref - Im(i,j,k,2,2,QGAME)
                      gfactor = (game - ONE)*(game - gam_g)
                      alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                   end if
                else

                   ! (tau, u T) eigensystem

                   ! eos to get some thermodynamics
                   eos_state%T = temp_ev
                   eos_state%rho = rho_ev
                   eos_state%xn(:) = q(i,j,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = q(i,j,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   p_r = eos_state%dpdr
                   p_T = eos_state%dpdT

                   alpham = HALF*(rho_ev**2*p_r*dtaum/Clag_ev + dvm - p_T*dTm/Clag_ev)/Clag_ev
                   alphap = HALF*(rho_ev**2*p_r*dtaup/Clag_ev - dvp - p_T*dTp/Clag_ev)/Clag_ev
                   alpha0r = dtau + (-rho_ev**2*p_r*dtau + p_T*dT0)/Clag_ev**2

                   ! not used, but needed to prevent bad invalid ops
                   alpha0e_g = ZERO

                endif    ! which tracing method

                if (v-cc > ZERO) then
                   alpham = ZERO
                else if (v-cc < ZERO) then
                   alpham = -alpham
                else
                   alpham = -HALF*alpham
                endif

                if (v+cc > ZERO) then
                   alphap = ZERO
                else if (v+cc < ZERO) then
                   alphap = -alphap
                else
                   alphap = -HALF*alphap
                endif

                if (v > ZERO) then
                   alpha0r = ZERO
                   alpha0e_g = ZERO
                else if (v < ZERO) then
                   alpha0r = -alpha0r
                   alpha0e_g = -alpha0e_g
                else
                   alpha0r = -HALF*alpha0r
                   alpha0e_g = -HALF*alpha0e_g
                endif

                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                if (ppm_temp_fix < 3) then

                   if (ppm_predict_gammae == 0) then
                      qyp(i,j,k,QRHO  ) = rho_ref + alphap + alpham + alpha0r
                      qyp(i,j,k,QV    ) = v_ref + (alphap - alpham)*cc_ev/rho_ev
                      qyp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                      qyp(i,j,k,QPRES ) = p_ref + (alphap + alpham)*csq_ev

                   else
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qyp(i,j,k,QRHO  ) = ONE/tau_s

                      qyp(i,j,k,QV    ) = v_ref + (alpham - alphap)*Clag_ev
                      qyp(i,j,k,QPRES ) = p_ref - (alphap + alpham)*Clag_ev**2

                      qyp(i,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                      qyp(i,j,k,QREINT) = qyp(i,j,k,QPRES )/(qyp(i,j,k,QGAME) - ONE)
                   endif

                else
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qyp(i,j,k,QRHO ) = ONE/tau_s

                   qyp(i,j,k,QV   ) = v_ref + (alpham - alphap)*Clag_ev
                   qyp(i,j,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once we do the species
                   qyp(i,j,k,QPRES) = small_pres ! just to make it defined

                endif

                ! Enforce small_*
                qyp(i,j,k,QRHO ) = max(qyp(i,j,k,QRHO ), small_dens)
                qyp(i,j,k,QPRES) = max(qyp(i,j,k,QPRES), small_pres)

                ! transverse velocities
                qyp(i,j,k,QU    ) = Im(i,j,k,2,2,QU) + hdt*Im_src(i,j,k,2,2,QU)
                qyp(i,j,k,QW    ) = Im(i,j,k,2,2,QW) + hdt*Im_src(i,j,k,2,2,QW)

             end if

             !-------------------------------------------------------------------
             ! minus state on face j+1
             !-------------------------------------------------------------------

             if (j <= hi(2)) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,k,2,3,QRHO)
                v_ref    = Ip(i,j,k,2,3,QV)

                p_ref    = Ip(i,j,k,2,3,QPRES)
                rhoe_g_ref = Ip(i,j,k,2,3,QREINT)
                temp_ref = Ip(i,j,k,2,3,QTEMP)

                tau_ref  = ONE/Ip(i,j,k,2,3,QRHO)

                gam_g_ref  = Ip_gc(i,j,k,2,3,1)
                game_ref = Ip(i,j,k,2,3,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! *m are the jumps carried by v-c
                ! *p are the jumps carried by v+c

                dvm = v_ref - Ip(i,j,k,2,1,QV) - hdt*Ip_src(i,j,k,2,1,QV)
                dptotm = p_ref - Ip(i,j,k,2,1,QPRES) - hdt*Ip_src(i,j,k,2,1,QPRES)

                drho = rho_ref - Ip(i,j,k,2,2,QRHO) - hdt*Ip_src(i,j,k,2,2,QRHO)
                dptot = p_ref - Ip(i,j,k,2,2,QPRES) - hdt*Ip_src(i,j,k,2,2,QPRES)
                drhoe_g = rhoe_g_ref - Ip(i,j,k,2,2,QREINT) - hdt*Ip_src(i,j,k,2,2,QREINT)

                dTm = temp_ref - Ip(i,j,k,2,1,QTEMP)
                dT0 = temp_ref - Ip(i,j,k,2,2,QTEMP)
                dTp = temp_ref - Ip(i,j,k,2,3,QTEMP)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Ip(i,j,k,2,1,QRHO) + hdt*Ip_src(i,j,k,2,2,QRHO)/Ip(i,j,k,2,2,QRHO)**2
                dtau = tau_ref - ONE/Ip(i,j,k,2,2,QRHO) + hdt*Ip_src(i,j,k,2,2,QRHO)/Ip(i,j,k,2,2,QRHO)**2
                dtaup = tau_ref - ONE/Ip(i,j,k,2,3,QRHO) + hdt*Ip_src(i,j,k,2,2,QRHO)/Ip(i,j,k,2,2,QRHO)**2

                dvp = v_ref - Ip(i,j,k,2,3,QV) - hdt*Ip_src(i,j,k,2,3,QV)
                dptotp = p_ref - Ip(i,j,k,2,3,QPRES) - hdt*Ip_src(i,j,k,2,3,QPRES)

                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   Clag_ev = Clag
                   h_g_ev = h_g
                   p_ev    = p
                   tau_ev  = ONE/rho
                   temp_ev = temp
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   Clag_ev = Clag_ref
                   h_g_ev = h_g_ref
                   p_ev    = p_ref
                   tau_ev  = tau_ref
                   temp_ev = temp_ref
                endif

                if (ppm_temp_fix < 3) then

                   if (ppm_predict_gammae == 0) then

                      ! (rho, u, p, (rho e) eigensystem

                      ! These are analogous to the beta's from the original PPM
                      ! paper.  This is simply (l . dq), where dq = qref - I(q)

                      alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dvm)*(rho_ev/cc_ev)
                      alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dvp)*(rho_ev/cc_ev)
                      alpha0r = drho - dptot/csq_ev
                      alpha0e_g = drhoe_g - dptot*h_g_ev

                   else

                      ! (tau, u, p, game) eigensystem

                      ! This is the way things were done in the original PPM
                      ! paper -- here we work with tau in the characteristic
                      ! system

                      alpham = HALF*( dvm - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alphap = HALF*(-dvp - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                      dge = game_ref - Ip(i,j,k,2,2,QGAME)
                      gfactor = (game - ONE)*(game - gam_g)
                      alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                   end if
                else

                   ! (tau, u T) eigensystem

                   ! eos to get some thermodynamics
                   eos_state%T = temp_ev
                   eos_state%rho = rho_ev
                   eos_state%xn(:) = q(i,j,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = q(i,j,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   p_r = eos_state%dpdr
                   p_T = eos_state%dpdT

                   alpham = HALF*(rho_ev**2*p_r*dtaum/Clag_ev + dvm - p_T*dTm/Clag_ev)/Clag_ev
                   alphap = HALF*(rho_ev**2*p_r*dtaup/Clag_ev - dvp - p_T*dTp/Clag_ev)/Clag_ev
                   alpha0r = dtau + (-rho_ev**2*p_r*dtau + p_T*dT0)/Clag_ev**2

                   ! not used, but needed to prevent bad invalid ops
                   alpha0e_g = ZERO

                endif

                if (v-cc > ZERO) then
                   alpham = -alpham
                else if (v-cc < ZERO) then
                   alpham = ZERO
                else
                   alpham = -HALF*alpham
                endif

                if (v+cc > ZERO) then
                   alphap = -alphap
                else if (v+cc < ZERO) then
                   alphap = ZERO
                else
                   alphap = -HALF*alphap
                endif

                if (v > ZERO) then
                   alpha0r = -alpha0r
                   alpha0e_g = -alpha0e_g
                else if (v < ZERO) then
                   alpha0r = ZERO
                   alpha0e_g = ZERO
                else
                   alpha0r = -HALF*alpha0r
                   alpha0e_g = -HALF*alpha0e_g
                endif

                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above has the minus already
                if (ppm_temp_fix < 3) then

                   if (ppm_predict_gammae == 0) then
                      qym(i,j+1,k,QRHO  ) = rho_ref + alphap + alpham + alpha0r
                      qym(i,j+1,k,QV    ) = v_ref + (alphap - alpham)*cc_ev/rho_ev
                      qym(i,j+1,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                      qym(i,j+1,k,QPRES ) = p_ref + (alphap + alpham)*csq_ev

                   else
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qym(i,j+1,k,QRHO  ) = ONE/tau_s

                      qym(i,j+1,k,QV    ) = v_ref + (alpham - alphap)*Clag_ev
                      qym(i,j+1,k,QPRES ) = p_ref - (alphap + alpham)*Clag_ev**2

                      qym(i,j+1,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                      qym(i,j+1,k,QREINT) = qym(i,j+1,k,QPRES )/(qym(i,j+1,k,QGAME) - ONE)

                   endif
                else
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qym(i,j+1,k,QRHO ) = ONE/tau_s

                   qym(i,j+1,k,QV   ) = v_ref + (alpham - alphap)*Clag_ev
                   qym(i,j+1,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once we do the species
                   qym(i,j+1,k,QPRES) = small_pres ! just to make it defined

                endif

                ! Enforce small_*
                qym(i,j+1,k,QRHO ) = max(qym(i,j+1,k,QRHO ), small_dens)
                qym(i,j+1,k,QPRES) = max(qym(i,j+1,k,QPRES), small_pres)

                ! transverse velocities
                qym(i,j+1,k,QU    ) = Ip(i,j,k,2,2,QU) + hdt*Ip_src(i,j,k,2,2,QU)
                qym(i,j+1,k,QW    ) = Ip(i,j,k,2,2,QW) + hdt*Ip_src(i,j,k,2,2,QW)
             end if

          end do
       end do
    end do

    !-------------------------------------------------------------------------
    ! Passively advected quantities
    !-------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  But we already dealt with all 3 velocity
       ! components above, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       ! Plus state on face j
       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2), hi(2)+1
             do i = lo(1)-1, hi(1)+1
                v = q(i,j,k,QV)

                if (v > ZERO) then
                   qyp(i,j,k,n) = q(i,j,k,n)
                else if (v < ZERO) then
                   qyp(i,j,k,n) = Im(i,j,k,2,2,n)
                else
                   qyp(i,j,k,n) = q(i,j,k,n) + HALF*(Im(i,j,k,2,2,n) - q(i,j,k,n))
                endif
             end do
          end do
       end do

       ! Minus state on face j+1
       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2)-1, hi(2)
             do i = lo(1)-1, hi(1)+1
                v = q(i,j,k,QV)

                if (v > ZERO) then
                   qym(i,j+1,k,n) = Ip(i,j,k,2,2,n)
                else if (v < ZERO) then
                   qym(i,j+1,k,n) = q(i,j,k,n)
                else
                   qym(i,j+1,k,n) = q(i,j,k,n) + HALF*(Ip(i,j,k,2,2,n) - q(i,j,k,n))
                endif
             end do
          end do
       end do
    end do

    if (ppm_temp_fix == 3) then
       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

                if (j >= lo(2)) then
                   ! plus face
                   eos_state%T     = qyp(i,j,k,QTEMP)
                   eos_state%rho   = qyp(i,j,k,QRHO)
                   eos_state%xn(:) = qyp(i,j,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = qyp(i,j,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   qyp(i,j,k,QPRES) = eos_state%p
                   qyp(i,j,k,QREINT) = qyp(i,j,k,QRHO)*eos_state%e

                   qyp(i,j,k,QPRES) = max(qyp(i,j,k,QPRES), small_pres)
                endif

                if (j <= hi(2)) then
                   ! minus face
                   eos_state%T     = qym(i,j+1,k,QTEMP)
                   eos_state%rho   = qym(i,j+1,k,QRHO)
                   eos_state%xn(:) = qym(i,j+1,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = qym(i,j+1,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   qym(i,j+1,k,QPRES) = eos_state%p
                   qym(i,j+1,k,QREINT) = qym(i,j+1,k,QRHO)*eos_state%e

                   qym(i,j+1,k,QPRES) = max(qym(i,j+1,k,QPRES), small_pres)
                endif

             end do
          end do
       end do
    endif
#endif

  end subroutine tracexy_ppm



  subroutine tracez_ppm(q, qd_lo, qd_hi, &
                        qaux, qa_lo, qa_hi, &
                        Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                        qzm, qzp, qs_lo, qs_hi, &
                        lo, hi, domlo, domhi, &
                        dt)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QTEMP, QGAME, QC, QGAMC, QFS, QFX, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   ppm_reference_eigenvectors, ppm_predict_gammae, &
                                   ppm_temp_fix, &
                                   npassive, qpass_map
    use amrex_constants_module, only : ZERO, HALF, ONE
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qs_lo(3),qs_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: I_lo(3), I_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,NQ)
    real(rt), intent(in) :: Im(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,NQ)

    real(rt), intent(in) :: Ip_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)
    real(rt), intent(in) :: Im_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)

    real(rt), intent(in) :: Ip_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,1)
    real(rt), intent(in) :: Im_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,1)

    real(rt), intent(inout) :: qzm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qzp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)

    real(rt), intent(in) :: dt

    !     Local variables
    integer :: i, j, k
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
    real(rt) :: rho, u, v, w, p, rhoe_g, h_g, temp
    real(rt) :: gam_g, game

    real(rt) :: drho, dptot, drhoe_g
    real(rt) :: de, dge, dtau, dtaum, dtaup
    real(rt) :: dwp, dptotp
    real(rt) :: dwm, dptotm
    real(rt) :: dT0, dTp, dTm
    real(rt) :: p_r, p_T

    real(rt) :: rho_ref, w_ref, p_ref, rhoe_g_ref, h_g_ref, temp_ref
    real(rt) :: tau_ref

    real(rt) :: cc_ref, csq_ref, Clag_ref, gam_g_ref, game_ref, gfactor
    real(rt) :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, h_g_ev, tau_ev, temp_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g

    real(rt) :: tau_s, e_s

    type(eos_t) :: eos_state

    hdt = HALF * dt

#ifndef AMREX_USE_CUDA
    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_3d.f90 :: tracez_ppm")
    end if
#endif

    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! Trace to left and right edges using upwind PPM
    !
    ! Note: in contrast to the above code for x and y, here the loop
    ! is over interfaces, not over cell-centers.

    !-------------------------------------------------------------------------
    ! z-direction
    !-------------------------------------------------------------------------

    ! Trace to left and right edges using upwind PPM
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             gfactor = ONE ! to help compiler resolve ANTI dependence

             rho  = q(i,j,k,QRHO)

             cc   = qaux(i,j,k,QC)
             csq  = cc**2
             Clag = rho*cc

             u    = q(i,j,k,QU)
             v    = q(i,j,k,QV)
             w    = q(i,j,k,QW)

             p    = q(i,j,k,QPRES)
             rhoe_g = q(i,j,k,QREINT)
             h_g = ( (p+rhoe_g)/rho )/csq
             temp = q(i,j,k,QTEMP)

             gam_g = qaux(i,j,k,QGAMC)
             game = q(i,j,k,QGAME)


             !-------------------------------------------------------------------
             ! plus state on face k
             !-------------------------------------------------------------------

             if (k >= lo(3)) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward user
                rho_ref  = Im(i,j,k,3,1,QRHO)
                w_ref    = Im(i,j,k,3,1,QW)

                p_ref    = Im(i,j,k,3,1,QPRES)
                rhoe_g_ref = Im(i,j,k,3,1,QREINT)
                temp_ref = Im(i,j,k,3,1,QTEMP)

                tau_ref  = ONE/Im(i,j,k,3,1,QRHO)

                gam_g_ref  = Im_gc(i,j,k,3,1,1)
                game_ref = Im(i,j,k,3,1,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! Note: for the transverse velocities, the jump is carried
                !       only by the w wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dwm = w_ref - Im(i,j,k,3,1,QW) - hdt*Im_src(i,j,k,3,1,QW)
                dptotm = p_ref - Im(i,j,k,3,1,QPRES) - hdt*Im_src(i,j,k,3,1,QPRES)

                drho = rho_ref - Im(i,j,k,3,2,QRHO) - hdt*Im_src(i,j,k,3,2,QRHO)
                dptot = p_ref - Im(i,j,k,3,2,QPRES) - hdt*Im_src(i,j,k,3,2,QPRES)
                drhoe_g = rhoe_g_ref - Im(i,j,k,3,2,QREINT) - hdt*Im_src(i,j,k,3,2,QREINT)

                ! TODO: need to figure sources for this out...
                dTm = temp_ref - Im(i,j,k,3,1,QTEMP)
                dT0 = temp_ref - Im(i,j,k,3,2,QTEMP)
                dTp = temp_ref - Im(i,j,k,3,3,QTEMP)

                ! we are treating tau as 1/rho, but we could have reconstructed
                ! it separately
                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Im(i,j,k,3,1,QRHO) + hdt*Im_src(i,j,k,3,1,QRHO)/Im(i,j,k,3,1,QRHO)**2
                dtau = tau_ref - ONE/Im(i,j,k,3,2,QRHO) + hdt*Im_src(i,j,k,3,2,QRHO)/Im(i,j,k,3,2,QRHO)**2
                dtaup = tau_ref - ONE/Im(i,j,k,3,3,QRHO) + hdt*Im_src(i,j,k,3,3,QRHO)/Im(i,j,k,3,3,QRHO)**2

                dwp = w_ref - Im(i,j,k,3,3,QW) - hdt*Im_src(i,j,k,3,3,QW)
                dptotp = p_ref - Im(i,j,k,3,3,QPRES) - hdt*Im_src(i,j,k,3,3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   Clag_ev = Clag
                   h_g_ev = h_g
                   p_ev    = p
                   tau_ev  = ONE/rho
                   temp_ev = temp
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   Clag_ev = Clag_ref
                   h_g_ev = h_g_ref
                   p_ev    = p_ref
                   tau_ev  = tau_ref
                   temp_ev = temp_ref
                endif

                if (ppm_temp_fix < 3) then

                   if (ppm_predict_gammae == 0) then

                      ! (rho, u, p, (rho e) eigensystem

                      ! These are analogous to the beta's from the original PPM
                      ! paper (except we work with rho instead of tau).  This is
                      ! simply (l . dq), where dq = qref - I(q)

                      alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dwm)*(rho_ev/cc_ev)
                      alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dwp)*(rho_ev/cc_ev)
                      alpha0r = drho - dptot/csq_ev
                      alpha0e_g = drhoe_g - dptot*h_g_ev

                   else

                      ! (tau, u, p, game) eigensystem

                      ! This is the way things were done in the original PPM
                      ! paper -- here we work with tau in the characteristic
                      ! system

                      alpham = HALF*( dwm - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alphap = HALF*(-dwp - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                      dge = game_ref - Im(i,j,k,3,2,QGAME)
                      gfactor = (game - ONE)*(game - gam_g)
                      alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                   endif

                else

                   ! (tau, u T) eigensystem

                   ! eos to get some thermodynamics
                   eos_state%T = temp_ev
                   eos_state%rho = rho_ev
                   eos_state%xn(:) = q(i,j,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = q(i,j,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   p_r = eos_state%dpdr
                   p_T = eos_state%dpdT

                   alpham = HALF*(rho_ev**2*p_r*dtaum/Clag_ev + dwm - p_T*dTm/Clag_ev)/Clag_ev
                   alphap = HALF*(rho_ev**2*p_r*dtaup/Clag_ev - dwp - p_T*dTp/Clag_ev)/Clag_ev
                   alpha0r = dtau + (-rho_ev**2*p_r*dtau + p_T*dT0)/Clag_ev**2

                   ! not used, but needed to prevent bad invalid ops
                   alpha0e_g = ZERO

                endif

                if (w-cc > ZERO) then
                   alpham = ZERO
                else if (w-cc < ZERO) then
                   alpham = -alpham
                else
                   alpham = -HALF*alpham
                endif

                if (w+cc > ZERO) then
                   alphap = ZERO
                else if (w+cc < ZERO) then
                   alphap = -alphap
                else
                   alphap = -HALF*alphap
                endif

                if (w > ZERO) then
                   alpha0r = ZERO
                   alpha0e_g = ZERO
                else if (w < ZERO) then
                   alpha0r = -alpha0r
                   alpha0e_g = -alpha0e_g
                else
                   alpha0r = -HALF*alpha0r
                   alpha0e_g = -HALF*alpha0e_g
                endif

                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                if (ppm_temp_fix < 3) then

                   if (ppm_predict_gammae == 0) then
                      qzp(i,j,k,QRHO  ) = rho_ref + alphap + alpham + alpha0r
                      qzp(i,j,k,QW    ) = w_ref + (alphap - alpham)*cc_ev/rho_ev
                      qzp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                      qzp(i,j,k,QPRES ) = p_ref + (alphap + alpham)*csq_ev

                   else
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qzp(i,j,k,QRHO  ) = ONE/tau_s

                      qzp(i,j,k,QW    ) = w_ref + (alpham - alphap)*Clag_ev
                      qzp(i,j,k,QPRES ) = p_ref - (alphap + alpham)*Clag_ev**2

                      qzp(i,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                      qzp(i,j,k,QREINT) = qzp(i,j,k,QPRES )/(qzp(i,j,k,QGAME) - ONE)

                   endif
                else
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qzp(i,j,k,QRHO ) = ONE/tau_s

                   qzp(i,j,k,QW   ) = w_ref + (alpham - alphap)*Clag_ev
                   qzp(i,j,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once we do the species
                   qzp(i,j,k,QPRES) = small_pres ! just to make it defined

                endif

                ! Enforce small_*
                qzp(i,j,k,QRHO ) = max(qzp(i,j,k,QRHO ), small_dens)
                qzp(i,j,k,QPRES) = max(qzp(i,j,k,QPRES), small_pres)

                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qzp(i,j,k,QU    ) = Im(i,j,k,3,2,QU) + hdt*Im_src(i,j,k,3,2,QU)
                qzp(i,j,k,QV    ) = Im(i,j,k,3,2,QV) + hdt*Im_src(i,j,k,3,2,QV)

             endif

             !-------------------------------------------------------------------
             ! minus state on face k + 1
             !-------------------------------------------------------------------
             if (k <= hi(3)) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,k,3,3,QRHO)
                w_ref    = Ip(i,j,k,3,3,QW)

                p_ref    = Ip(i,j,k,3,3,QPRES)
                rhoe_g_ref = Ip(i,j,k,3,3,QREINT)
                temp_ref = Ip(i,j,k,3,3,QTEMP)

                tau_ref  = ONE/Ip(i,j,k,3,3,QRHO)

                gam_g_ref  = Ip_gc(i,j,k,3,3,1)
                game_ref = Ip(i,j,k,3,3,QGAME)

                rho_ref = max(rho_ref,small_dens)
                p_ref = max(p_ref,small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! *m are the jumps carried by w-c
                ! *p are the jumps carried by w+c

                dwm = w_ref - Ip(i,j,k,3,1,QW) - hdt*Ip_src(i,j,k,3,1,QW)
                dptotm = p_ref - Ip(i,j,k,3,1,QPRES) - hdt*Ip_src(i,j,k,3,1,QPRES)

                drho = rho_ref - Ip(i,j,k,3,2,QRHO) - hdt*Ip_src(i,j,k,3,2,QRHO)
                dptot = p_ref - Ip(i,j,k,3,2,QPRES) - hdt*Ip_src(i,j,k,3,2,QPRES)
                drhoe_g = rhoe_g_ref - Ip(i,j,k,3,2,QREINT) - hdt*Ip_src(i,j,k,3,2,QREINT)

                dTm = temp_ref - Ip(i,j,k,3,1,QTEMP)
                dT0 = temp_ref - Ip(i,j,k,3,2,QTEMP)
                dTp = temp_ref - Ip(i,j,k,3,3,QTEMP)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Ip(i,j,k,3,1,QRHO) + hdt*Ip_src(i,j,k,3,1,QRHO)/Ip(i,j,k,3,1,QRHO)**2
                dtau = tau_ref - ONE/Ip(i,j,k,3,2,QRHO) + hdt*Ip_src(i,j,k,3,2,QRHO)/Ip(i,j,k,3,2,QRHO)**2
                dtaup = tau_ref - ONE/Ip(i,j,k,3,3,QRHO) + hdt*Ip_src(i,j,k,3,3,QRHO)/Ip(i,j,k,3,3,QRHO)**2

                dwp = w_ref - Ip(i,j,k,3,3,QW) - hdt*Ip_src(i,j,k,3,3,QW)
                dptotp = p_ref - Ip(i,j,k,3,3,QPRES) - hdt*Ip_src(i,j,k,3,3,QPRES)

                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   rho_ev  = rho
                   cc_ev   = cc
                   csq_ev  = csq
                   Clag_ev = Clag
                   h_g_ev = h_g
                   p_ev    = p
                   tau_ev  = ONE/rho
                   temp_ev = temp
                else
                   rho_ev  = rho_ref
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
                   Clag_ev = Clag_ref
                   h_g_ev = h_g_ref
                   p_ev    = p_ref
                   tau_ev  = tau_ref
                   temp_ev = temp_ref
                endif

                if (ppm_temp_fix < 3) then

                   if (ppm_predict_gammae == 0) then

                      ! (rho, u, p, (rho e)) eigensystem

                      ! These are analogous to the beta's from the original PPM
                      ! paper (except we work with rho instead of tau).  This is
                      ! simply (l . dq), where dq = qref - I(q)

                      alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dwm)*(rho_ev/cc_ev)
                      alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dwp)*(rho_ev/cc_ev)
                      alpha0r = drho - dptot/csq_ev
                      alpha0e_g = drhoe_g - dptot*h_g_ev

                   else

                      ! (tau, u, p, game) eigensystem

                      ! This is the way things were done in the original PPM
                      ! paper -- here we work with tau in the characteristic
                      ! system

                      alpham = HALF*( dwm - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alphap = HALF*(-dwp - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                      alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                      dge = game_ref - Ip(i,j,k,3,2,QGAME)
                      gfactor = (game - ONE)*(game - gam_g)
                      alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                   endif
                else
                   ! (tau, u T) eigensystem

                   ! eos to get some thermodynamics
                   eos_state%T = temp_ev
                   eos_state%rho = rho_ev
                   eos_state%xn(:) = q(i,j,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = q(i,j,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   p_r = eos_state%dpdr
                   p_T = eos_state%dpdT

                   alpham = HALF*(rho_ev**2*p_r*dtaum/Clag_ev + dwm - p_T*dTm/Clag_ev)/Clag_ev
                   alphap = HALF*(rho_ev**2*p_r*dtaup/Clag_ev - dwp - p_T*dTp/Clag_ev)/Clag_ev
                   alpha0r = dtau + (-rho_ev**2*p_r*dtau + p_T*dT0)/Clag_ev**2

                   ! not used, but needed to prevent bad invalid ops
                   alpha0e_g = ZERO

                endif

                if (w-cc > ZERO) then
                   alpham = -alpham
                else if (w-cc < ZERO) then
                   alpham = ZERO
                else
                   alpham = -HALF*alpham
                endif

                if (w+cc > ZERO) then
                   alphap = -alphap
                else if (w+cc < ZERO) then
                   alphap = ZERO
                else
                   alphap = -HALF*alphap
                endif

                if (w > ZERO) then
                   alpha0r = -alpha0r
                   alpha0e_g = -alpha0e_g
                else if (w < ZERO) then
                   alpha0r = ZERO
                   alpha0e_g = ZERO
                else
                   alpha0r = -HALF*alpha0r
                   alpha0e_g = -HALF*alpha0e_g
                endif

                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                if (ppm_temp_fix < 3) then
                   if (ppm_predict_gammae == 0) then
                      qzm(i,j,k+1,QRHO  ) = rho_ref + alphap + alpham + alpha0r
                      qzm(i,j,k+1,QW    ) = w_ref + (alphap - alpham)*cc_ev/rho_ev
                      qzm(i,j,k+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                      qzm(i,j,k+1,QPRES ) = p_ref + (alphap + alpham)*csq_ev

                   else
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qzm(i,j,k+1,QRHO  ) = ONE/tau_s

                      qzm(i,j,k+1,QW    ) = w_ref + (alpham - alphap)*Clag_ev
                      qzm(i,j,k+1,QPRES ) = p_ref - (alphap + alpham)*Clag_ev**2

                      qzm(i,j,k+1,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                      qzm(i,j,k+1,QREINT) = qzm(i,j,k+1,QPRES )/(qzm(i,j,k+1,QGAME) - ONE)

                   endif
                else
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qzm(i,j,k+1,QRHO ) = ONE/tau_s

                   qzm(i,j,k+1,QW   ) = w_ref + (alpham - alphap)*Clag_ev
                   qzm(i,j,k+1,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once we do the species
                   qzm(i,j,k+1,QPRES) = small_pres ! just to make it defined
                endif

                ! Enforce small_*
                qzm(i,j,k+1,QRHO ) = max(qzm(i,j,k+1,QRHO ),small_dens)
                qzm(i,j,k+1,QPRES) = max(qzm(i,j,k+1,QPRES),small_pres)

                ! Transverse velocity
                qzm(i,j,k+1,QU    ) = Ip(i,j,k,3,2,QU) + hdt*Ip_src(i,j,k,3,2,QU)
                qzm(i,j,k+1,QV    ) = Ip(i,j,k,3,2,QV) + hdt*Ip_src(i,j,k,3,2,QV)

             end if

          end do
       end do
    end do

    !-------------------------------------------------------------------------
    ! passively advected quantities
    !-------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  But we already dealt with all 3 velocity
       ! components above, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       ! Plus state on face k
       do k = lo(3), hi(3)+1
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

                w = q(i,j,k,QW)

                if (w > ZERO) then
                   qzp(i,j,k,n) = q(i,j,k,n)
                else if (w < ZERO) then
                   qzp(i,j,k,n) = Im(i,j,k,3,2,n)
                else
                   qzp(i,j,k,n) = q(i,j,k,n) + HALF*(Im(i,j,k,3,2,n) - q(i,j,k,n))
                endif
             end do
          end do
       end do

       ! Minus state on face k+1
       do k = lo(3)-1, hi(3)
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

                w = q(i,j,k,QW)

                if (w > ZERO) then
                   qzm(i,j,k+1,n) = Ip(i,j,k,3,2,n)
                else if (w < ZERO) then
                   qzm(i,j,k+1,n) = q(i,j,k,n)
                else
                   qzm(i,j,k+1,n) = q(i,j,k,n) + HALF*(Ip(i,j,k,3,2,n) - q(i,j,k,n))
                endif

             end do
          end do
       end do
    end do

    if (ppm_temp_fix == 3) then
       ! we predicted T, now make p, (rho e) consistent

       do k = lo(3)-1, hi(3)+1
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1

                if (k >= lo(3)) then
                   ! plus face
                   eos_state%T     = qzp(i,j,k,QTEMP)
                   eos_state%rho   = qzp(i,j,k,QRHO)
                   eos_state%xn(:) = qzp(i,j,k,QFS:QFS-1+nspec)
                   eos_state%aux(:) = qzp(i,j,k,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   qzp(i,j,k,QPRES) = eos_state%p
                   qzp(i,j,k,QREINT) = qzp(i,j,k,QRHO)*eos_state%e
                   qzp(i,j,k,QPRES) = max(qzp(i,j,k,QPRES), small_pres)
                endif

                if (k <= hi(3)) then
                   ! minus face
                   eos_state%T     = qzm(i,j,k+1,QTEMP)
                   eos_state%rho   = qzm(i,j,k+1,QRHO)
                   eos_state%xn(:) = qzm(i,j,k+1,QFS:QFS-1+nspec)
                   eos_state%aux(:) = qzm(i,j,k+1,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   qzm(i,j,k+1,QPRES) = eos_state%p
                   qzm(i,j,k+1,QREINT) = qzm(i,j,k+1,QRHO)*eos_state%e

                   qzm(i,j,k+1,QPRES) = max(qzm(i,j,k+1,QPRES), small_pres)
                end if

             end do
          end do
       end do
    endif

  end subroutine tracez_ppm

end module trace_ppm_module
