! These routines do the characteristic tracing under the parabolic
! profiles in each zone to the edge / half-time.

module trace_ppm_module

  use prob_params_module, only : dg
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public trace_ppm, trace_ppm_gammae, trace_ppm_temp

contains

  subroutine trace_ppm(idir, q, qd_lo, qd_hi, &
                       qaux, qa_lo, qa_hi, &
                       Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                       qm, qp, qs_lo, qs_hi, &
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

    integer, intent(in) :: idir
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

    real(rt), intent(inout) :: qm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j, k
    integer :: n, ipassive

    type(eos_t) :: eos_state

    real(rt) :: hdt

    integer :: ix, iy, iz, QUN, QUT, QUTT

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

    real(rt) :: cc, csq, cgassq
    real(rt) :: rho, un, ut, utt, p, rhoe_g, h_g
    real(rt) :: gam_g

    real(rt) :: drho, dptot, drhoe_g
    real(rt) :: dup, dvp, dptotp
    real(rt) :: dum, dvm, dptotm

    real(rt) :: rho_ref, un_ref, p_ref, rhoe_g_ref, h_g_ref

    real(rt) :: cc_ref, csq_ref, gam_g_ref
    real(rt) :: cc_ev, csq_ev, rho_ev, h_g_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp

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

    if (idir == 1) then
       ix = 1
       iy = 0
       iz = 0
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       ix = 0
       iy = 1
       iz = 0
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       ix = 0
       iy = 0
       iz = 1
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! Trace to left and right edges using upwind PPM
    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-1, hi(1)+1

             rho = q(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2

             un = q(i,j,k,QUN)
             ut = q(i,j,k,QUT)
             utt = q(i,j,k,QUTT)

             p = q(i,j,k,QPRES)
             rhoe_g = q(i,j,k,QREINT)
             h_g = ( (p + rhoe_g)/rho)/csq

             gam_g = qaux(i,j,k,QGAMC)


             !-------------------------------------------------------------------
             ! plus state on face i
             !-------------------------------------------------------------------

             if ((idir == 1 .and. i >= lo(1)) .or. &
                 (idir == 2 .and. j >= lo(2)) .or. &
                 (idir == 3 .and. k >= lo(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(i,j,k,idir,1,QRHO)
                un_ref    = Im(i,j,k,idir,1,QUN)

                p_ref    = Im(i,j,k,idir,1,QPRES)
                rhoe_g_ref = Im(i,j,k,idir,1,QREINT)

                gam_g_ref  = Im_gc(i,j,k,idir,1,1)

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
                dum = un_ref - Im(i,j,k,idir,1,QUN) - hdt*Im_src(i,j,k,idir,1,QUN)
                dptotm = p_ref - Im(i,j,k,idir,1,QPRES) - hdt*Im_src(i,j,k,idir,1,QPRES)

                drho = rho_ref - Im(i,j,k,idir,2,QRHO) - hdt*Im_src(i,j,k,idir,2,QRHO)
                dptot = p_ref - Im(i,j,k,idir,2,QPRES) - hdt*Im_src(i,j,k,idir,2,QPRES)
                drhoe_g = rhoe_g_ref - Im(i,j,k,idir,2,QREINT) - hdt*Im_src(i,j,k,idir,2,QREINT)

                dup = un_ref - Im(i,j,k,idir,3,QUN) - hdt*Im_src(i,j,k,idir,3,QUN)
                dptotp = p_ref - Im(i,j,k,idir,3,QPRES) - hdt*Im_src(i,j,k,idir,3,QPRES)


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

                alpham = merge(ZERO, -alpham, un-cc > ZERO)
                alphap = merge(ZERO, -alphap, un+cc > ZERO)
                alpha0r = merge(ZERO, -alpha0r, un > ZERO)
                alpha0e_g = merge(ZERO, -alpha0e_g, un > ZERO)

                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                qp(i,j,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                qp(i,j,k,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                qp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                qp(i,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)


                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp(i,j,k,QUT) = Im(i,j,k,idir,2,QUT) + hdt*Im_src(i,j,k,idir,2,QUT)
                qp(i,j,k,QUTT) = Im(i,j,k,idir,2,QUTT) + hdt*Im_src(i,j,k,idir,2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= hi(1)) .or. &
                 (idir == 2 .and. j <= hi(2)) .or. &
                 (idir == 3 .and. k <= hi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,k,idir,3,QRHO)
                un_ref    = Ip(i,j,k,idir,3,QUN)

                p_ref    = Ip(i,j,k,idir,3,QPRES)
                rhoe_g_ref = Ip(i,j,k,idir,3,QREINT)

                gam_g_ref  = Ip_gc(i,j,k,idir,3,1)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip(i,j,k,idir,1,QUN) - hdt*Ip_src(i,j,k,idir,1,QUN)
                dptotm  = p_ref - Ip(i,j,k,idir,1,QPRES) - hdt*Ip_src(i,j,k,idir,1,QPRES)

                drho = rho_ref - Ip(i,j,k,idir,2,QRHO) - hdt*Ip_src(i,j,k,idir,2,QRHO)
                dptot = p_ref - Ip(i,j,k,idir,2,QPRES) - hdt*Ip_src(i,j,k,idir,2,QPRES)
                drhoe_g = rhoe_g_ref - Ip(i,j,k,idir,2,QREINT) - hdt*Ip_src(i,j,k,idir,2,QREINT)

                dup = un_ref - Ip(i,j,k,idir,3,QUN) - hdt*Ip_src(i,j,k,idir,3,QUN)
                dptotp = p_ref - Ip(i,j,k,idir,3,QPRES) - hdt*Ip_src(i,j,k,idir,3,QPRES)

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

                alpham = merge(-alpham, ZERO, un-cc > ZERO)
                alphap = merge(-alphap, ZERO, un+cc > ZERO)
                alpha0r = merge(-alpha0r, ZERO, un > ZERO)
                alpha0e_g = merge(-alpha0e_g, ZERO, un > ZERO)

                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                qm(i+ix,j+iy,k+iz,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                qm(i+ix,j+iy,k+iz,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                qm(i+ix,j+iy,k+iz,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                qm(i+ix,j+iy,k+iz,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)

                ! transverse velocities
                qm(i+ix,j+iy,k+iz,QUT) = Ip(i,j,k,idir,2,QUT) + hdt*Ip_src(i,j,k,idir,2,QUTT)
                qm(i+ix,j+iy,k+iz,QUTT) = Ip(i,j,k,idir,2,QUTT) + hdt*Ip_src(i,j,k,idir,2,QUTT)

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

                if (i <= hi(1)) then

                   qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr
                   qm(i+1,j,k,QRHO) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                end if

                if (i >= lo(1)) then

                   qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr
                   qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                end if

             end if
#endif

#if (AMREX_SPACEDIM == 1)
             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_lo) then
                qm(lo(1),j,k,QRHO  ) = q(domlo(1)-1,j,k,QRHO)
                qm(lo(1),j,k,QUN   ) = q(domlo(1)-1,j,k,QUN )
                qm(lo(1),j,k,QPRES ) = q(domlo(1)-1,j,k,QPRES)
                qm(lo(1),j,k,QREINT) = q(domlo(1)-1,j,k,QREINT)
             end if

             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_hi) then
                qp(hi(1)+1,j,k,QRHO  ) = q(domhi(1)+1,j,k,QRHO)
                qp(hi(1)+1,j,k,QUN   ) = q(domhi(1)+1,j,k,QUN  )
                qp(hi(1)+1,j,k,QPRES ) = q(domhi(1)+1,j,k,QPRES)
                qp(hi(1)+1,j,k,QREINT) = q(domhi(1)+1,j,k,QREINT)
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
             do i = lo(1)-1, hi(1)+1

                ! Plus state on face i
                if ((idir == 1 .and. i >= lo(1)) .or. &
                    (idir == 2 .and. j >= lo(2)) .or. &
                    (idir == 3 .and. k >= lo(3))) then

                   un = q(i,j,k,QUN)

                   ! We have
                   !
                   ! q_l = q_ref - Proj{(q_ref - I)}
                   !
                   ! and Proj{} represents the characteristic projection.
                   ! But for these, there is only 1-wave that matters, the u
                   ! wave, so no projection is needed.  Since we are not
                   ! projecting, the reference state doesn't matter

                   qp(i,j,k,n) = merge(q(i,j,k,n), Im(i,j,k,idir,2,n), un > ZERO)
                endif

                ! Minus state on face i+1
                if (idir == 1 .and. i <= hi(1)) then

                   un = q(i,j,k,QUN)
                   qm(i+1,j,k,n) = merge(Ip(i,j,k,idir,2,n), q(i,j,k,n), un > ZERO)

                else if (idir == 2 .and. j <= hi(2)) then

                   un = q(i,j,k,QUN)
                   qm(i,j+1,k,n) = merge(Ip(i,j,k,idir,2,n), q(i,j,k,n), un > ZERO)

                else if (idir == 3 .and. k <= hi(3)) then

                   un = q(i,j,k,QUN)
                   qm(i,j,k+1,n) = merge(Ip(i,j,k,idir,2,n), q(i,j,k,n), un > ZERO)

                endif

             end do

#if AMREX_SPACEDIM == 1
             if (fix_mass_flux_hi) qp(hi(1)+1,j,k,n) = q(hi(1)+1,j,k,n)
             if (fix_mass_flux_lo) qm(lo(1),j,k,n) = q(lo(1)-1,j,k,n)
#endif
          end do
       end do

    end do

  end subroutine trace_ppm


  subroutine trace_ppm_gammae(idir, q, qd_lo, qd_hi, &
                              qaux, qa_lo, qa_hi, &
                              Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                              qm, qp, qs_lo, qs_hi, &
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

    integer, intent(in) :: idir
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

    real(rt), intent(inout) :: qm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j, k
    integer :: n, ipassive

    type(eos_t) :: eos_state

    real(rt) :: hdt

    integer :: ix, iy, iz, QUN, QUT, QUTT

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
    real(rt) :: rho, un, ut, utt, p, rhoe_g, h_g
    real(rt) :: gam_g, game

    real(rt) :: dptot
    real(rt) :: dge, dtau, dtaum, dtaup
    real(rt) :: dup, dvp, dptotp
    real(rt) :: dum, dvm, dptotm

    real(rt) :: rho_ref, un_ref, p_ref, rhoe_g_ref
    real(rt) :: tau_ref

    real(rt) :: cc_ref, csq_ref, Clag_ref, gam_g_ref, game_ref, gfactor
    real(rt) :: cc_ev, csq_ev, Clag_ev, tau_ev

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

    if (idir == 1) then
       ix = 1
       iy = 0
       iz = 0
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       ix = 0
       iy = 1
       iz = 0
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       ix = 0
       iy = 0
       iz = 1
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! Trace to left and right edges using upwind PPM
    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-1, hi(1)+1

             gfactor = ONE ! to help compiler resolve ANTI dependence

             rho = q(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2
             Clag = rho*cc

             un = q(i,j,k,QUN)
             ut = q(i,j,k,QUT)
             utt = q(i,j,k,QUTT)

             p = q(i,j,k,QPRES)
             rhoe_g = q(i,j,k,QREINT)

             gam_g = qaux(i,j,k,QGAMC)
             game = q(i,j,k,QGAME)


             !-------------------------------------------------------------------
             ! plus state on face i
             !-------------------------------------------------------------------

             if ((idir == 1 .and. i >= lo(1)) .or. &
                 (idir == 2 .and. j >= lo(2)) .or. &
                 (idir == 3 .and. k >= lo(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(i,j,k,idir,1,QRHO)
                un_ref    = Im(i,j,k,idir,1,QUN)

                p_ref    = Im(i,j,k,idir,1,QPRES)
                rhoe_g_ref = Im(i,j,k,idir,1,QREINT)

                tau_ref  = ONE/Im(i,j,k,idir,1,QRHO)

                gam_g_ref  = Im_gc(i,j,k,idir,1,1)
                game_ref = Im(i,j,k,idir,1,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dum = un_ref - Im(i,j,k,idir,1,QUN) - hdt*Im_src(i,j,k,idir,1,QUN)
                dptotm = p_ref - Im(i,j,k,idir,1,QPRES) - hdt*Im_src(i,j,k,idir,1,QPRES)

                dptot = p_ref - Im(i,j,k,idir,2,QPRES) - hdt*Im_src(i,j,k,idir,2,QPRES)

                ! we are treating tau as 1/rho, but we could have reconstructed
                ! it separately
                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Im(i,j,k,idir,1,QRHO) + hdt*Im_src(i,j,k,idir,1,QRHO)/Im(i,j,k,idir,1,QRHO)**2
                dtau  = tau_ref - ONE/Im(i,j,k,idir,2,QRHO) + hdt*Im_src(i,j,k,idir,2,QRHO)/Im(i,j,k,idir,2,QRHO)**2
                dtaup = tau_ref - ONE/Im(i,j,k,idir,3,QRHO) + hdt*Im_src(i,j,k,idir,3,QRHO)/Im(i,j,k,idir,3,QRHO)**2

                dup = un_ref - Im(i,j,k,idir,3,QUN) - hdt*Im_src(i,j,k,idir,3,QUN)
                dptotp = p_ref - Im(i,j,k,idir,3,QPRES) - hdt*Im_src(i,j,k,idir,3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   cc_ev   = cc
                   csq_ev  = csq
                   Clag_ev = Clag
                   tau_ev  = ONE/rho
                else
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
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

                dge   = game_ref - Im(i,j,k,idir,2,QGAME)
                gfactor = (game - ONE)*(game - gam_g)
                alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                alpham = merge(ZERO, -alpham, un-cc > ZERO)
                alphap = merge(ZERO, -alphap, un+cc > ZERO)
                alpha0r = merge(ZERO, -alpha0r, un > ZERO)
                alpha0e_g = merge(ZERO, -alpha0e_g, un > ZERO)

                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                tau_s = tau_ref + alphap + alpham + alpha0r
                qp(i,j,k,QRHO) = max(small_dens, ONE/tau_s)

                qp(i,j,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                qp(i,j,k,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

                qp(i,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                qp(i,j,k,QREINT) = qp(i,j,k,QPRES )/(qp(i,j,k,QGAME) - ONE)


                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp(i,j,k,QUT) = Im(i,j,k,idir,2,QUT) + hdt*Im_src(i,j,k,idir,2,QUT)
                qp(i,j,k,QUTT) = Im(i,j,k,idir,2,QUTT) + hdt*Im_src(i,j,k,idir,2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= hi(1)) .or. &
                 (idir == 2 .and. j <= hi(2)) .or. &
                 (idir == 3 .and. k <= hi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,k,idir,3,QRHO)
                un_ref    = Ip(i,j,k,idir,3,QUN)

                p_ref    = Ip(i,j,k,idir,3,QPRES)
                rhoe_g_ref = Ip(i,j,k,idir,3,QREINT)

                tau_ref  = ONE/Ip(i,j,k,idir,3,QRHO)

                gam_g_ref  = Ip_gc(i,j,k,idir,3,1)
                game_ref = Ip(i,j,k,idir,3,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip(i,j,k,idir,1,QUN) - hdt*Ip_src(i,j,k,idir,1,QUN)
                dptotm  = p_ref - Ip(i,j,k,idir,1,QPRES) - hdt*Ip_src(i,j,k,idir,1,QPRES)

                dptot = p_ref - Ip(i,j,k,idir,2,QPRES) - hdt*Ip_src(i,j,k,idir,2,QPRES)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Ip(i,j,k,idir,1,QRHO) + hdt*Ip_src(i,j,k,idir,1,QRHO)/Ip(i,j,k,idir,1,QRHO)**2
                dtau = tau_ref - ONE/Ip(i,j,k,idir,2,QRHO) + hdt*Ip_src(i,j,k,idir,2,QRHO)/Ip(i,j,k,idir,2,QRHO)**2
                dtaup = tau_ref - ONE/Ip(i,j,k,idir,3,QRHO) + hdt*Ip_src(i,j,k,idir,3,QRHO)/Ip(i,j,k,idir,3,QRHO)**2

                dup = un_ref - Ip(i,j,k,idir,3,QUN) - hdt*Ip_src(i,j,k,idir,3,QUN)
                dptotp = p_ref - Ip(i,j,k,idir,3,QPRES) - hdt*Ip_src(i,j,k,idir,3,QPRES)

                ! Optionally use the reference state in evaluating the
                ! eigenvectors
                if (ppm_reference_eigenvectors == 0) then
                   cc_ev   = cc
                   csq_ev  = csq
                   Clag_ev = Clag
                   tau_ev  = ONE/rho
                else
                   cc_ev   = cc_ref
                   csq_ev  = csq_ref
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

                dge = game_ref - Ip(i,j,k,idir,2,QGAME)
                gfactor = (game - ONE)*(game - gam_g)
                alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

                alpham = merge(-alpham, ZERO, un-cc > ZERO)
                alphap = merge(-alphap, ZERO, un+cc > ZERO)
                alpha0r = merge(-alpha0r, ZERO, un > ZERO)
                alpha0e_g = merge(-alpha0e_g, ZERO, un > ZERO)


                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                tau_s = tau_ref + alphap + alpham + alpha0r
                qm(i+ix,j+iy,k+iz,QRHO) = max(small_dens, ONE/tau_s)

                qm(i+ix,j+iy,k+iz,QUN) = un_ref + (alpham - alphap)*Clag_ev
                qm(i+ix,j+iy,k+iz,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

                qm(i+ix,j+iy,k+iz,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                qm(i+ix,j+iy,k+iz,QREINT) = qm(i+ix,j+iy,k+iz,QPRES )/(qm(i+ix,j+iy,k+iz,QGAME) - ONE)

                ! transverse velocities
                qm(i+ix,j+iy,k+iz,QUT) = Ip(i,j,k,idir,2,QUT) + hdt*Ip_src(i,j,k,idir,2,QUTT)
                qm(i+ix,j+iy,k+iz,QUTT) = Ip(i,j,k,idir,2,QUTT) + hdt*Ip_src(i,j,k,idir,2,QUTT)

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

                if (i <= hi(1)) then
                   qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr
                   qm(i+1,j,k,QRHO) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                end if

                if (i >= lo(1)) then
                   qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr
                   qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                end if

             endif
#endif

#if (AMREX_SPACEDIM == 1)
             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_lo) then
                qm(lo(1),j,k,QRHO  ) = q(domlo(1)-1,j,k,QRHO)
                qm(lo(1),j,k,QUN   ) = q(domlo(1)-1,j,k,QUN )
                qm(lo(1),j,k,QPRES ) = q(domlo(1)-1,j,k,QPRES)
                qm(lo(1),j,k,QREINT) = q(domlo(1)-1,j,k,QREINT)
             end if

             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_hi) then
                qp(hi(1)+1,j,k,QRHO  ) = q(domhi(1)+1,j,k,QRHO)
                qp(hi(1)+1,j,k,QUN   ) = q(domhi(1)+1,j,k,QUN  )
                qp(hi(1)+1,j,k,QPRES ) = q(domhi(1)+1,j,k,QPRES)
                qp(hi(1)+1,j,k,QREINT) = q(domhi(1)+1,j,k,QREINT)
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
             do i = lo(1)-1, hi(1)+1

                ! Plus state on face i
                if ((idir == 1 .and. i >= lo(1)) .or. &
                    (idir == 2 .and. j >= lo(2)) .or. &
                    (idir == 3 .and. k >= lo(3))) then

                   un = q(i,j,k,QUN)

                   ! We have
                   !
                   ! q_l = q_ref - Proj{(q_ref - I)}
                   !
                   ! and Proj{} represents the characteristic projection.
                   ! But for these, there is only 1-wave that matters, the u
                   ! wave, so no projection is needed.  Since we are not
                   ! projecting, the reference state doesn't matter

                   if (un > ZERO) then
                      qp(i,j,k,n) = q(i,j,k,n)
                   else if (un < ZERO) then
                      qp(i,j,k,n) = Im(i,j,k,idir,2,n)
                   else
                      qp(i,j,k,n) = q(i,j,k,n) + HALF*(Im(i,j,k,idir,2,n) - q(i,j,k,n))
                   endif
                endif

                ! Minus state on face i+1
                if ((idir == 1 .and. i <= hi(1)) .or. &
                    (idir == 2 .and. j <= hi(2)) .or. &
                    (idir == 3 .and. k <= hi(3))) then

                   un = q(i,j,k,QUN)

                   if (un > ZERO) then
                      qm(i+ix,j+iy,k+iz,n) = Ip(i,j,k,idir,2,n)
                   else if (un < ZERO) then
                      qm(i+ix,j+iy,k+iz,n) = q(i,j,k,n)
                   else
                      qm(i+ix,j+iy,k+iz,n) = q(i,j,k,n) + HALF*(Ip(i,j,k,idir,2,n) - q(i,j,k,n))
                   endif
                endif

             end do

#if AMREX_SPACEDIM == 1
             if (fix_mass_flux_hi) qp(hi(1)+1,j,k,n) = q(hi(1)+1,j,k,n)
             if (fix_mass_flux_lo) qm(lo(1),j,k,n) = q(lo(1)-1,j,k,n)
#endif
          end do
       end do

    end do

  end subroutine trace_ppm_gammae


  subroutine trace_ppm_temp(idir, q, qd_lo, qd_hi, &
                       qaux, qa_lo, qa_hi, &
                       Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                       qm, qp, qs_lo, qs_hi, &
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

    integer, intent(in) :: idir
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

    real(rt), intent(inout) :: qm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j, k
    integer :: n, ipassive

    type(eos_t) :: eos_state

    real(rt) :: hdt

    integer :: ix, iy, iz, QUN, QUT, QUTT

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
    real(rt) :: rho, un, ut, utt, p, rhoe_g, h_g, temp
    real(rt) :: gam_g, game

    real(rt) :: drho, dptot, drhoe_g
    real(rt) :: de, dge, dtau, dtaum, dtaup
    real(rt) :: dup, dvp, dptotp
    real(rt) :: dum, dvm, dptotm
    real(rt) :: dT0, dTp, dTm
    real(rt) :: p_r, p_T

    real(rt) :: rho_ref, un_ref, p_ref, rhoe_g_ref, temp_ref
    real(rt) :: tau_ref

    real(rt) :: cc_ref, csq_ref, Clag_ref, gam_g_ref, game_ref, gfactor
    real(rt) :: cc_ev, csq_ev, Clag_ev, rho_ev, tau_ev, temp_ev

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

    if (idir == 1) then
       ix = 1
       iy = 0
       iz = 0
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       ix = 0
       iy = 1
       iz = 0
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       ix = 0
       iy = 0
       iz = 1
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! Trace to left and right edges using upwind PPM
    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-1, hi(1)+1

             gfactor = ONE ! to help compiler resolve ANTI dependence

             rho = q(i,j,k,QRHO)

             cc = qaux(i,j,k,QC)
             csq = cc**2
             Clag = rho*cc

             un = q(i,j,k,QUN)
             ut = q(i,j,k,QUT)
             utt = q(i,j,k,QUTT)

             p = q(i,j,k,QPRES)
             rhoe_g = q(i,j,k,QREINT)
             temp = q(i,j,k,QTEMP)

             gam_g = qaux(i,j,k,QGAMC)
             game = q(i,j,k,QGAME)


             !-------------------------------------------------------------------
             ! plus state on face i
             !-------------------------------------------------------------------

             if ((idir == 1 .and. i >= lo(1)) .or. &
                 (idir == 2 .and. j >= lo(2)) .or. &
                 (idir == 3 .and. k >= lo(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(i,j,k,idir,1,QRHO)
                un_ref    = Im(i,j,k,idir,1,QUN)

                p_ref    = Im(i,j,k,idir,1,QPRES)
                rhoe_g_ref = Im(i,j,k,idir,1,QREINT)
                temp_ref = Im(i,j,k,idir,1,QTEMP)

                tau_ref  = ONE/Im(i,j,k,idir,1,QRHO)

                gam_g_ref  = Im_gc(i,j,k,idir,1,1)
                game_ref = Im(i,j,k,idir,1,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dum = un_ref - Im(i,j,k,idir,1,QUN) - hdt*Im_src(i,j,k,idir,1,QUN)
                dptotm = p_ref - Im(i,j,k,idir,1,QPRES) - hdt*Im_src(i,j,k,idir,1,QPRES)

                drho = rho_ref - Im(i,j,k,idir,2,QRHO) - hdt*Im_src(i,j,k,idir,2,QRHO)
                dptot = p_ref - Im(i,j,k,idir,2,QPRES) - hdt*Im_src(i,j,k,idir,2,QPRES)
                drhoe_g = rhoe_g_ref - Im(i,j,k,idir,2,QREINT) - hdt*Im_src(i,j,k,idir,2,QREINT)

                ! TODO: need to figure sources for this out...
                dTm = temp_ref - Im(i,j,k,idir,1,QTEMP)
                dT0 = temp_ref - Im(i,j,k,idir,2,QTEMP)
                dTp = temp_ref - Im(i,j,k,idir,3,QTEMP)

                ! we are treating tau as 1/rho, but we could have reconstructed
                ! it separately
                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Im(i,j,k,idir,1,QRHO) + hdt*Im_src(i,j,k,idir,1,QRHO)/Im(i,j,k,idir,1,QRHO)**2
                dtau  = tau_ref - ONE/Im(i,j,k,idir,2,QRHO) + hdt*Im_src(i,j,k,idir,2,QRHO)/Im(i,j,k,idir,2,QRHO)**2
                dtaup = tau_ref - ONE/Im(i,j,k,idir,3,QRHO) + hdt*Im_src(i,j,k,idir,3,QRHO)/Im(i,j,k,idir,3,QRHO)**2

                dup = un_ref - Im(i,j,k,idir,3,QUN) - hdt*Im_src(i,j,k,idir,3,QUN)
                dptotp = p_ref - Im(i,j,k,idir,3,QPRES) - hdt*Im_src(i,j,k,idir,3,QPRES)


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

                alpham = merge(ZERO, -alpham, un-cc > ZERO)
                alphap = merge(ZERO, -alphap, un+cc > ZERO)
                alpha0r = merge(ZERO, -alpha0r, un > ZERO)

                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                tau_s = tau_ref + alphap + alpham + alpha0r
                qp(i,j,k,QRHO) = max(small_dens, ONE/tau_s)

                qp(i,j,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                qp(i,j,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                     rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                ! we defer getting the pressure until later, once we do the species
                qp(i,j,k,QPRES) = small_pres ! just to make it defined

                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp(i,j,k,QUT) = Im(i,j,k,idir,2,QUT) + hdt*Im_src(i,j,k,idir,2,QUT)
                qp(i,j,k,QUTT) = Im(i,j,k,idir,2,QUTT) + hdt*Im_src(i,j,k,idir,2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= hi(1)) .or. &
                 (idir == 2 .and. j <= hi(2)) .or. &
                 (idir == 3 .and. k <= hi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,k,idir,3,QRHO)
                un_ref    = Ip(i,j,k,idir,3,QUN)

                p_ref    = Ip(i,j,k,idir,3,QPRES)
                rhoe_g_ref = Ip(i,j,k,idir,3,QREINT)
                temp_ref = Ip(i,j,k,idir,3,QTEMP)

                tau_ref  = ONE/Ip(i,j,k,idir,3,QRHO)

                gam_g_ref  = Ip_gc(i,j,k,idir,3,1)
                game_ref = Ip(i,j,k,idir,3,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip(i,j,k,idir,1,QUN) - hdt*Ip_src(i,j,k,idir,1,QUN)
                dptotm  = p_ref - Ip(i,j,k,idir,1,QPRES) - hdt*Ip_src(i,j,k,idir,1,QPRES)

                drho = rho_ref - Ip(i,j,k,idir,2,QRHO) - hdt*Ip_src(i,j,k,idir,2,QRHO)
                dptot = p_ref - Ip(i,j,k,idir,2,QPRES) - hdt*Ip_src(i,j,k,idir,2,QPRES)
                drhoe_g = rhoe_g_ref - Ip(i,j,k,idir,2,QREINT) - hdt*Ip_src(i,j,k,idir,2,QREINT)

                dTm = temp_ref - Ip(i,j,k,idir,1,QTEMP)
                dT0 = temp_ref - Ip(i,j,k,idir,2,QTEMP)
                dTp = temp_ref - Ip(i,j,k,idir,3,QTEMP)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Ip(i,j,k,idir,1,QRHO) + hdt*Ip_src(i,j,k,idir,1,QRHO)/Ip(i,j,k,idir,1,QRHO)**2
                dtau = tau_ref - ONE/Ip(i,j,k,idir,2,QRHO) + hdt*Ip_src(i,j,k,idir,2,QRHO)/Ip(i,j,k,idir,2,QRHO)**2
                dtaup = tau_ref - ONE/Ip(i,j,k,idir,3,QRHO) + hdt*Ip_src(i,j,k,idir,3,QRHO)/Ip(i,j,k,idir,3,QRHO)**2

                dup = un_ref - Ip(i,j,k,idir,3,QUN) - hdt*Ip_src(i,j,k,idir,3,QUN)
                dptotp = p_ref - Ip(i,j,k,idir,3,QPRES) - hdt*Ip_src(i,j,k,idir,3,QPRES)

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

                alpham = merge(-alpham, ZERO, un-cc > ZERO)
                alphap = merge(-alphap, ZERO, un+cc > ZERO)
                alpha0r = merge(-alpha0r, ZERO, un > ZERO)


                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already
                tau_s = tau_ref + alphap + alpham + alpha0r
                qm(i+ix,j+iy,k+iz,QRHO) = max(small_dens, ONE/tau_s)

                qm(i+ix,j+iy,k+iz,QUN) = un_ref + (alpham - alphap)*Clag_ev
                qm(i+ix,j+iy,k+iz,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                     rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                ! we defer getting the pressure until later, once we do the species
                qm(i+ix,j+iy,k+iz,QPRES) = small_pres ! just to make it defined

                ! transverse velocities
                qm(i+ix,j+iy,k+iz,QUT) = Ip(i,j,k,idir,2,QUT) + hdt*Ip_src(i,j,k,idir,2,QUTT)
                qm(i+ix,j+iy,k+iz,QUTT) = Ip(i,j,k,idir,2,QUTT) + hdt*Ip_src(i,j,k,idir,2,QUTT)

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

                if (i <= hi(1)) then
                   qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr
                   qm(i+1,j,k,QRHO) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                end if

                if (i >= lo(1)) then
                   qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr
                   qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                end if

             endif
#endif

#if (AMREX_SPACEDIM == 1)
             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_lo) then
                qm(lo(1),j,k,QRHO  ) = q(domlo(1)-1,j,k,QRHO)
                qm(lo(1),j,k,QUN   ) = q(domlo(1)-1,j,k,QUN )
                qm(lo(1),j,k,QPRES ) = q(domlo(1)-1,j,k,QPRES)
                qm(lo(1),j,k,QREINT) = q(domlo(1)-1,j,k,QREINT)
             end if

             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_hi) then
                qp(hi(1)+1,j,k,QRHO  ) = q(domhi(1)+1,j,k,QRHO)
                qp(hi(1)+1,j,k,QUN   ) = q(domhi(1)+1,j,k,QUN  )
                qp(hi(1)+1,j,k,QPRES ) = q(domhi(1)+1,j,k,QPRES)
                qp(hi(1)+1,j,k,QREINT) = q(domhi(1)+1,j,k,QREINT)
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
             do i = lo(1)-1, hi(1)+1

                ! Plus state on face i
                if ((idir == 1 .and. i >= lo(1)) .or. &
                    (idir == 2 .and. j >= lo(2)) .or. &
                    (idir == 3 .and. k >= lo(3))) then

                   un = q(i,j,k,QUN)

                   ! We have
                   !
                   ! q_l = q_ref - Proj{(q_ref - I)}
                   !
                   ! and Proj{} represents the characteristic projection.
                   ! But for these, there is only 1-wave that matters, the u
                   ! wave, so no projection is needed.  Since we are not
                   ! projecting, the reference state doesn't matter

                   if (un > ZERO) then
                      qp(i,j,k,n) = q(i,j,k,n)
                   else if (un < ZERO) then
                      qp(i,j,k,n) = Im(i,j,k,idir,2,n)
                   else
                      qp(i,j,k,n) = q(i,j,k,n) + HALF*(Im(i,j,k,idir,2,n) - q(i,j,k,n))
                   endif
                endif

                ! Minus state on face i+1
                if ((idir == 1 .and. i <= hi(1)) .or. &
                    (idir == 2 .and. j <= hi(2)) .or. &
                    (idir == 3 .and. k <= hi(3))) then

                   un = q(i,j,k,QUN)

                   if (un > ZERO) then
                      qm(i+ix,j+iy,k+iz,n) = Ip(i,j,k,idir,2,n)
                   else if (un < ZERO) then
                      qm(i+ix,j+iy,k+iz,n) = q(i,j,k,n)
                   else
                      qm(i+ix,j+iy,k+iz,n) = q(i,j,k,n) + HALF*(Ip(i,j,k,idir,2,n) - q(i,j,k,n))
                   endif
                endif

             end do

#if AMREX_SPACEDIM == 1
             if (fix_mass_flux_hi) qp(hi(1)+1,j,k,n) = q(hi(1)+1,j,k,n)
             if (fix_mass_flux_lo) qm(lo(1),j,k,n) = q(lo(1)-1,j,k,n)
#endif
          end do
       end do

    end do

    ! we predicted T, now make p, (rho e) consistent
    do k = lo(3)-dg(3), hi(3)-dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-1, hi(1)+1

             if ((idir == 1 .and. i >= lo(1)) .or. &
                 (idir == 2 .and. j >= lo(2)) .or. &
                 (idir == 3 .and. k >= lo(3))) then

                ! plus face
                eos_state%T     = qp(i,j,k,QTEMP)
                eos_state%rho   = qp(i,j,k,QRHO)
                eos_state%xn(:) = qp(i,j,k,QFS:QFS-1+nspec)
                eos_state%aux(:) = qp(i,j,k,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                qp(i,j,k,QPRES) = eos_state%p
                qp(i,j,k,QREINT) = qp(i,j,k,QRHO)*eos_state%e

                qp(i,j,k,QPRES) = max(qp(i,j,k,QPRES), small_pres)
             endif

             if ((idir == 1 .and. i <= hi(1)) .or. &
                 (idir == 2 .and. j <= hi(2)) .or. &
                 (idir == 3 .and. k <= hi(3))) then

                ! minus face
                eos_state%T     = qm(i+ix,j+iy,k+iz,QTEMP)
                eos_state%rho   = qm(i+ix,j+iy,k+iz,QRHO)
                eos_state%xn(:) = qm(i+ix,j+iy,k+iz,QFS:QFS-1+nspec)
                eos_state%aux(:) = qm(i+ix,j+iy,k+iz,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                qm(i+ix,j+iy,k+iz,QPRES) = eos_state%p
                qm(i+ix,j+iy,k+iz,QREINT) = qm(i+ix,j+iy,k+iz,QRHO)*eos_state%e

                qm(i+ix,j+iy,k+iz,QPRES) = max(qm(i+ix,j+iy,k+iz,QPRES), small_pres)
             endif

          end do
       end do
    end do

  end subroutine trace_ppm_temp

end module trace_ppm_module
