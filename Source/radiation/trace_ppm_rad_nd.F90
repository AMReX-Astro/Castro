! These routines do the characteristic tracing under the parabolic
! profiles in each zone to the edge / half-time.

module trace_ppm_rad_module

  use prob_params_module, only : dg
  use castro_error_module, only : castro_error
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public trace_ppm_rad

contains


  subroutine trace_ppm_rad(lo, hi, &
                           idir, &
                           q, qd_lo, qd_hi, &
                           qaux, qa_lo, qa_hi, &
                           srcQ, src_lo, src_hi, &
                           flatn, f_lo, f_hi, &
                           qm, qm_lo, qm_hi, &
                           qp, qp_lo, qp_hi, &
#if (AMREX_SPACEDIM < 3)
                           dloga, dloga_lo, dloga_hi, &
#endif
                           vlo, vhi, domlo, domhi, &
                           dx, dt)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, NQSRC, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QGAME, QC, QCG, QGAMC, QGAMCG, QLAMS, QTEMP, &
                                   qrad, qradhi, qptot, qreitot, &
                                   small_dens, small_pres, &
                                   ppm_type, ppm_temp_fix, &
                                   ppm_reference_eigenvectors, ppm_predict_gammae

    use rad_params_module, only : ngroups
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow
    use trace_ppm_module, only : trace_ppm_species
    use amrex_fort_module, only : rt => amrex_real
    use ppm_module, only : ppm_reconstruct, ppm_int_profile, ppm_reconstruct_with_eos

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NQSRC)
    real(rt), intent(in) ::  flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)

#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    integer :: QUN, QUT, QUTT

    ! Local variables
    integer :: i, j, k, g
    integer :: n, ipassive

    real(rt) :: hdt, dtdx

    real(rt) :: sm, sp

    real(rt) :: s(-2:2)
    real(rt) :: Ip(1:3,NQ)
    real(rt) :: Im(1:3,NQ)

    real(rt) :: Ip_src(1:3,NQSRC)
    real(rt) :: Im_src(1:3,NQSRC)


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

    logical :: source_nonzero(NQSRC)
    logical :: reconstruct_state(NQ)

    !$gpu

    hdt = HALF * dt
    dtdx = dt / dx(idir)

    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call castro_error("Error:: trace_ppm_rad_nd.f90 :: tracexy_ppm_rad")
    end if


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

             rho = q(i,j,k,QRHO)
             tau = ONE/rho

             ! cgassq is the gas soundspeed **2
             ! cc is the total soundspeed **2 (gas + radiation)
             cgassq = qaux(i,j,k,QCG)**2
             cc = qaux(i,j,k,QC)
             csq = cc**2
             Clag = rho*cc

             un = q(i,j,k,QUN)

             p = q(i,j,k,QPRES)
             rhoe_g = q(i,j,k,QREINT)
             h_g = ( (p+rhoe_g)/rho)/csq

             gam_g = qaux(i,j,k,QGAMCG)
             game = q(i,j,k,QGAME)

             ptot = q(i,j,k,qptot)

             er(:) = q(i,j,k,qrad:qradhi)
             hr(:) = (lam0+ONE)*er/rho


             ! do the parabolic reconstruction and compute the
             ! integrals under the characteristic waves
             do n = 1, NQ
                if (.not. reconstruct_state(n)) cycle

                if (idir == 1) then
                   s(:) = q(i-2:i+2,j,k,n)
                else if (idir == 2) then
                   s(:) = q(i,j-2:j+2,k,n)
                else
                   s(:) = q(i,j,k-2:k+2,n)
                end if

                call ppm_reconstruct(s, flatn(i,j,k), sm, sp)

                call ppm_int_profile(sm, sp, s(0), un, cc, dtdx, Ip(:,n), Im(:,n))
             end do

             ! source terms
             do n = 1, NQSRC
                if (source_nonzero(n)) then

                   if (idir == 1) then
                      s(:) = srcQ(i-2:i+2,j,k,n)
                   else if (idir == 2) then
                      s(:) = srcQ(i,j-2:j+2,k,n)
                   else
                      s(:) = srcQ(i,j,k-2:k+2,n)
                   end if

                   call ppm_reconstruct(s, flatn(i,j,k), sm, sp)

                   call ppm_int_profile(sm, sp, s(0), un, cc, dtdx, Ip_src(:,n), Im_src(:,n))
                else
                   Ip_src(:,n) = ZERO
                   Im_src(:,n) = ZERO
                end if

             end do


             ! do the passives separately
             call trace_ppm_species(i, j, k, &
                                    idir, &
                                    q, qd_lo, qd_hi, &
                                    Ip, Im, &
                                    Ip_src, Im_src, &
                                    qm, qm_lo, qm_hi, &
                                    qp, qp_lo, qp_hi, &
                                    vlo, vhi, domlo, domhi, &
                                    dx, dt)


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
                rho_ref  = Im(1,QRHO)
                un_ref    = Im(1,QUN)

                p_ref    = Im(1,QPRES)
                rhoe_g_ref = Im(1,QREINT)

                tau_ref  = ONE/Im(1,QRHO)

                !gam_g_ref  = Im_gc(i,j,k,1,1)
                game_ref = Im(1,QGAME)

                ptot_ref = Im(1,QPTOT)

                er_ref(:) = Im(1,QRAD:QRADHI)


                rho_ref = max(rho_ref,small_dens)
                p_ref = max(p_ref,small_pres)

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)

                ! we also add the sources here so they participate in the tracing
                dum    = un_ref    - Im(1,QUN) - hdt*Im_src(1,QUN)
                dptotm = ptot_ref - Im(1,qptot) - hdt*Im_src(1,QPRES)

                drho    = rho_ref    - Im(2,QRHO) - hdt*Im_src(2,QRHO)
                dptot   = ptot_ref   - Im(2,qptot) - hdt*Im_src(2,QPRES)
                drhoe_g = rhoe_g_ref - Im(2,QREINT) - hdt*Im_src(2,QREINT)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtau  = tau_ref  - ONE/Im(2,QRHO) + hdt*Im_src(2,QRHO)/Im(2,QRHO)**2
                der(:)  = er_ref(:)  - Im(2,qrad:qradhi)

                dup    = un_ref    - Im(3,QUN) - hdt*Im_src(3,QUN)
                dptotp = ptot_ref - Im(3,qptot) - hdt*Im_src(3,QPRES)


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

                   dge   = game_ref - Im(2,QGAME)
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
                   qp(i,j,k,QRHO) = rho_ref + alphap + alpham + alpha0r
                   qp(i,j,k,QUN) = un_ref + (alphap - alpham)*cc/rho
                   qp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                   qp(i,j,k,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

                   qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                   qp(i,j,k,qrad:qradhi) = qrtmp

                   qp(i,j,k,qptot) = ptot_ref + (alphap + alpham)*csq
                   qp(i,j,k,qreitot) = qp(i,j,k,QREINT) + sum(qrtmp)

                else
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qp(i,j,k,QRHO  ) = ONE/tau_s

                   qp(i,j,k,QUN    ) = un_ref + (alpham - alphap)*Clag
                   qp(i,j,k,QPRES ) = p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamp(:)*alphar(:))

                   qp(i,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                   qp(i,j,k,QREINT) = qp(i,j,k,QPRES )/(qp(i,j,k,QGAME) - ONE)

                   qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                   qp(i,j,k,qrad:qradhi) = qrtmp

                   qp(i,j,k,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                   qp(i,j,k,qreitot) = qp(i,j,k,QREINT) + sum(qrtmp)

                endif

                ! Enforce small_*
                qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                qp(i,j,k,QPRES) = max(qp(i,j,k,QPRES),small_pres)

                do g = 0, ngroups-1
                   if (qp(i,j,k,qrad+g) < ZERO) then
                      er_foo = - qp(i,j,k,qrad+g)
                      qp(i,j,k,qrad+g) = ZERO
                      qp(i,j,k,qptot) = qp(i,j,k,qptot) + lamp(g) * er_foo
                      qp(i,j,k,qreitot) = qp(i,j,k,qreitot) + er_foo
                   end if
                end do


                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp(i,j,k,QUT) = Im(2,QUT) + hdt*Im_src(2,QUT)
                qp(i,j,k,QUTT) = Im(2,QUTT) + hdt*Im_src(2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= vhi(1)) .or. &
                 (idir == 2 .and. j <= vhi(2)) .or. &
                 (idir == 3 .and. k <= vhi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(3,QRHO)
                un_ref    = Ip(3,QUN)

                p_ref    = Ip(3,QPRES)
                rhoe_g_ref = Ip(3,QREINT)

                tau_ref  = ONE/Ip(3,QRHO)

                !gam_g_ref  = Ip_gc(i,j,k,3,1)
                game_ref = Ip(3,QGAME)


                ptot_ref = Ip(3,QPTOT)

                er_ref(:) = Ip(3,QRAD:QRADHI)

                rho_ref = max(rho_ref,small_dens)
                p_ref = max(p_ref,small_pres)

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum    = un_ref    - Ip(1,QUN) - hdt*Ip_src(1,QUN)
                dptotm = ptot_ref - Ip(1,qptot) - hdt*Ip_src(1,QPRES)

                drho    = rho_ref    - Ip(2,QRHO) - hdt*Ip_src(2,QRHO)
                dptot   = ptot_ref   - Ip(2,qptot) - hdt*Ip_src(2,QPRES)
                drhoe_g = rhoe_g_ref - Ip(2,QREINT) - hdt*Ip_src(2,QREINT)
                dtau  = tau_ref  - ONE/Ip(2,QRHO) + hdt*Ip_src(2,QRHO)/Ip(2,QRHO)**2
                der(:)  = er_ref(:)  - Ip(2,qrad:qradhi)

                dup    = un_ref    - Ip(3,QUN) - hdt*Ip_src(3,QUN)
                dptotp = ptot_ref - Ip(3,qptot) - hdt*Ip_src(3,QPRES)


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

                   dge = game_ref - Ip(2,QGAME)
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
                      qm(i+1,j,k,QRHO) = max(small_dens, rho_ref + alphap + alpham + alpha0r)
                      qm(i+1,j,k,QUN) = un_ref + (alphap - alpham)*cc/rho
                      qm(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                      qm(i+1,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:)))

                      qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                      qm(i+1,j,k,qrad:qradhi) = qrtmp

                      qm(i+1,j,k,qptot) = ptot_ref + (alphap + alpham)*csq
                      qm(i+1,j,k,qreitot) = qm(i+1,j,k,QREINT) + sum(qrtmp)

                   else if (idir == 2) then
                      qm(i,j+1,k,QRHO) = max(small_dens, rho_ref + alphap + alpham + alpha0r)
                      qm(i,j+1,k,QUN) = un_ref + (alphap - alpham)*cc/rho
                      qm(i,j+1,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                      qm(i,j+1,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:)))

                      qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                      qm(i,j+1,k,qrad:qradhi) = qrtmp

                      qm(i,j+1,k,qptot) = ptot_ref + (alphap + alpham)*csq
                      qm(i,j+1,k,qreitot) = qm(i,j+1,k,QREINT) + sum(qrtmp)

                   else if (idir == 3) then
                      qm(i,j,k+1,QRHO) = max(small_dens, rho_ref + alphap + alpham + alpha0r)
                      qm(i,j,k+1,QUN) = un_ref + (alphap - alpham)*cc/rho
                      qm(i,j,k+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                      qm(i,j,k+1,QPRES) = max(small_pres, p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:)))

                      qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                      qm(i,j,k+1,qrad:qradhi) = qrtmp

                      qm(i,j,k+1,qptot) = ptot_ref + (alphap + alpham)*csq
                      qm(i,j,k+1,qreitot) = qm(i,j,k+1,QREINT) + sum(qrtmp)

                   end if

                else

                   if (idir == 1) then
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qm(i+1,j,k,QRHO  ) = max(small_dens, ONE/tau_s)

                      qm(i+1,j,k,QUN    ) = un_ref + (alpham - alphap)*Clag
                      qm(i+1,j,k,QPRES ) = max(small_pres, p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamm(:)*alphar(:)))

                      qm(i+1,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                      qm(i+1,j,k,QREINT) = qm(i+1,j,k,QPRES )/(qm(i+1,j,k,QGAME) - ONE)

                      qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                      qm(i+1,j,k,qrad:qradhi) = qrtmp

                      qm(i+1,j,k,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                      qm(i+1,j,k,qreitot) = qm(i+1,j,k,QREINT) + sum(qrtmp)

                   else if (idir == 2) then
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qm(i,j+1,k,QRHO  ) = max(small_dens, ONE/tau_s)

                      qm(i,j+1,k,QUN    ) = un_ref + (alpham - alphap)*Clag
                      qm(i,j+1,k,QPRES ) = max(small_pres, p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamm(:)*alphar(:)))

                      qm(i,j+1,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                      qm(i,j+1,k,QREINT) = qm(i,j+1,k,QPRES )/(qm(i,j+1,k,QGAME) - ONE)

                      qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                      qm(i,j+1,k,qrad:qradhi) = qrtmp

                      qm(i,j+1,k,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                      qm(i,j+1,k,qreitot) = qm(i,j+1,k,QREINT) + sum(qrtmp)

                   else if (idir == 3) then
                      tau_s = tau_ref + alphap + alpham + alpha0r
                      qm(i,j,k+1,QRHO  ) = max(small_dens, ONE/tau_s)

                      qm(i,j,k+1,QUN    ) = un_ref + (alpham - alphap)*Clag
                      qm(i,j,k+1,QPRES ) = max(small_pres, p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamm(:)*alphar(:)))

                      qm(i,j,k+1,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                      qm(i,j,k+1,QREINT) = qm(i,j,k+1,QPRES )/(qm(i,j,k+1,QGAME) - ONE)

                      qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                      qm(i,j,k+1,qrad:qradhi) = qrtmp

                      qm(i,j,k+1,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                      qm(i,j,k+1,qreitot) = qm(i,j,k+1,QREINT) + sum(qrtmp)

                   end if

                endif

                if (idir == 1) then
                   do g=0,ngroups-1
                      if (qm(i+1,j,k,qrad+g) < ZERO) then
                         er_foo = - qm(i+1,j,k,qrad+g)
                         qm(i+1,j,k,qrad+g) = ZERO
                         qm(i+1,j,k,qptot) = qm(i+1,j,k,qptot) + lamm(g) * er_foo
                         qm(i+1,j,k,qreitot) = qm(i+1,j,k,qreitot) + er_foo
                      end if
                   end do

                   ! transverse velocities
                   qm(i+1,j,k,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i+1,j,k,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                else if (idir == 2) then
                   do g=0,ngroups-1
                      if (qm(i,j+1,k,qrad+g) < ZERO) then
                         er_foo = - qm(i,j+1,k,qrad+g)
                         qm(i,j+1,k,qrad+g) = ZERO
                         qm(i,j+1,k,qptot) = qm(i,j+1,k,qptot) + lamm(g) * er_foo
                         qm(i,j+1,k,qreitot) = qm(i,j+1,k,qreitot) + er_foo
                      end if
                   end do

                   ! transverse velocities
                   qm(i,j+1,k,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i,j+1,k,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                else if (idir == 3) then
                   do g=0,ngroups-1
                      if (qm(i,j,k+1,qrad+g) < ZERO) then
                         er_foo = - qm(i,j,k+1,qrad+g)
                         qm(i,j,k+1,qrad+g) = ZERO
                         qm(i,j,k+1,qptot) = qm(i,j,k+1,qptot) + lamm(g) * er_foo
                         qm(i,j,k+1,qreitot) = qm(i,j,k+1,qreitot) + er_foo
                      end if
                   end do

                   ! transverse velocities
                   qm(i,j,k+1,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i,j,k+1,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

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
                   qm(i+1,j,k,QRHO  ) = qm(i+1,j,k,QRHO  ) + sourcr
                   qm(i+1,j,k,QRHO  ) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES ) = qm(i+1,j,k,QPRES ) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                   qm(i+1,j,k,qrad:qradhi) = qm(i+1,j,k,qrad:qradhi) + sourcer(:)
                   ! qm(i+1,j,k,qptot ) = sum(lamm(:)*qm(i+1,j,k,qrad:qradhi)) + qm(i+1,j,k,QPRES)
                   qm(i+1,j,k,qptot) = qm(i+1,j,k,qptot) + sum(lamm(:)*sourcer(:)) + sourcp
                   qm(i+1,j,k,qreitot) = sum(qm(i+1,j,k,qrad:qradhi))  + qm(i+1,j,k,QREINT)
                end if

                if (i >= vlo(1)) then
                   qp(i,j,k,QRHO  ) = qp(i,j,k,QRHO  ) + sourcr
                   qp(i,j,k,QRHO  ) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES ) = qp(i,j,k,QPRES ) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                   qp(i,j,k,qrad:qradhi) = qp(i,j,k,qrad:qradhi) + sourcer(:)
                   ! qp(i  ,qptot ) = sum(lamp(:)*qp(i,qrad:qradhi)) + qp(i,QPRES)
                   qp(i,j,k,qptot) = qp(i,j,k,qptot) + sum(lamp(:)*sourcer(:)) + sourcp
                   qp(i,j,k,qreitot) = sum(qp(i,j,k,qrad:qradhi))  + qp(i,j,k,QREINT)
                end if
             endif
#endif

          end do
       end do
    end do

  end subroutine trace_ppm_rad

end module trace_ppm_rad_module
