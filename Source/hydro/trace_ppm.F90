module trace_ppm_module
  !
  ! These routines do the characteristic tracing under the parabolic
  ! profiles in each zone to the edge / half-time.

  use prob_params_module, only : dg
  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module, only : ZERO, HALF, ONE

  implicit none

contains

  subroutine trace_ppm(lo, hi, &
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
    ! here, lo and hi are the range we loop over -- this can include ghost cells
    ! vlo and vhi are the bounds of the valid box (no ghost cells)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, NQSRC, ppm_predict_gammae, &
                                   ppm_temp_fix, QU, QV, QW, QGAME, QREINT, QTEMP, npassive, qpass_map
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

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
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    integer :: n, i, j, k

    logical :: source_nonzero(NQSRC)
    logical :: reconstruct_state(NQ)

    !$gpu

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


    if (ppm_temp_fix < 3) then
       if (ppm_predict_gammae == 0) then
          call trace_ppm_rhoe(lo, hi, &
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
                              reconstruct_state, source_nonzero, &
                              dx, dt)
       else
          call trace_ppm_gammae(lo, hi, &
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
                                reconstruct_state, source_nonzero, &
                                dx, dt)
       end if
    else
       call trace_ppm_temp(lo, hi, &
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
                           reconstruct_state, source_nonzero, &
                           dx, dt)
    end if

  end subroutine trace_ppm


  subroutine trace_ppm_species(i, j, k, &
                               idir, &
                               q, qd_lo, qd_hi, &
                               Ip, Im, &
                               Ip_src, Im_src, &
                               qm, qm_lo, qm_hi, &
                               qp, qp_lo, qp_hi, &
                               vlo, vhi, domlo, domhi, &
                               dx, dt)
    ! here, lo and hi are the range we loop over -- this can include ghost cells
    ! vlo and vhi are the bounds of the valid box (no ghost cells)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQSRC, ppm_predict_gammae, &
                                   ppm_temp_fix, QU, QV, QW, npassive, qpass_map
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)

    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    real(rt), intent(in) :: Ip(1:3,NQ)
    real(rt), intent(in) :: Im(1:3,NQ)

    real(rt), intent(in) :: Ip_src(1:3,NQSRC)
    real(rt), intent(in) :: Im_src(1:3,NQSRC)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)
    real(rt), intent(in) :: dt, dx(3)

    integer, intent(in) :: i, j, k

    integer :: ipassive, n

    !$gpu

    ! the passive stuff is the same regardless of the tracing
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! Plus state on face i
       if ((idir == 1 .and. i >= vlo(1)) .or. &
           (idir == 2 .and. j >= vlo(2)) .or. &
           (idir == 3 .and. k >= vlo(3))) then

          ! We have
          !
          ! q_l = q_ref - Proj{(q_ref - I)}
          !
          ! and Proj{} represents the characteristic projection.
          ! But for these, there is only 1-wave that matters, the u
          ! wave, so no projection is needed.  Since we are not
          ! projecting, the reference state doesn't matter

          qp(i,j,k,n) = Im(2,n)
          if (n <= NQSRC) qp(i,j,k,n) = qp(i,j,k,n) + HALF*dt*Im_src(2,n)

       end if

       ! Minus state on face i+1
       if (idir == 1 .and. i <= vhi(1)) then
          qm(i+1,j,k,n) = Ip(2,n)
          if (n <= NQSRC) qm(i+1,j,k,n) = qm(i+1,j,k,n) + HALF*dt*Ip_src(2,n)

       else if (idir == 2 .and. j <= vhi(2)) then
          qm(i,j+1,k,n) = Ip(2,n)
          if (n <= NQSRC) qm(i,j+1,k,n) = qm(i,j+1,k,n) + HALF*dt*Ip_src(2,n)

       else if (idir == 3 .and. k <= vhi(3)) then
          qm(i,j,k+1,n) = Ip(2,n)
          if (n <= NQSRC) qm(i,j,k+1,n) = qm(i,j,k+1,n) + HALF*dt*Ip_src(2,n)
       end if

    end do
  end subroutine trace_ppm_species


  subroutine trace_ppm_rhoe(lo, hi, &
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
                            reconstruct_state, source_nonzero, &
                            dx, dt)

    use network, only : nspec, naux

    use meth_params_module, only : NQ, NQAUX, NQSRC, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QGAME, QC, QGAMC, &
                                   small_dens, small_pres, &
                                   ppm_type, ppm_temp_fix, &
                                   ppm_reference_eigenvectors
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow
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
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    logical, intent(in) :: source_nonzero(NQSRC)
    logical, intent(in) :: reconstruct_state(NQ)

    ! Local variables
    integer :: i, j, k, n

    real(rt) :: hdt, dtdx

    real(rt) :: sm, sp

    real(rt) :: s(-2:2)
    real(rt) :: Ip(1:3,NQ), Im(1:3,NQ)
    real(rt) :: Ip_src(1:3,NQSRC), Im_src(1:3,NQSRC)
    real(rt) :: Ip_gc(1:3,1), Im_gc(1:3,1)

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
       call castro_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    !$gpu

    hdt = HALF * dt
    dtdx = dt / dx(idir)

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


             if (idir == 1) then
                s(:) = qaux(i-2:i+2,j,k,QGAMC)
             else if (idir == 2) then
                s(:) = qaux(i,j-2:j+2,k,QGAMC)
             else
                s(:) = qaux(i,j,k-2:k+2,QGAMC)
             end if

             call ppm_reconstruct(s, flatn(i,j,k), sm, sp)

             call ppm_int_profile(sm, sp, s(0), un, cc, dtdx, Ip_gc, Im_gc)

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

                gam_g_ref  = Im_gc(1,1)

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
                dum = un_ref - Im(1,QUN) - hdt*Im_src(1,QUN)
                dptotm = p_ref - Im(1,QPRES) - hdt*Im_src(1,QPRES)

                drho = rho_ref - Im(2,QRHO) - hdt*Im_src(2,QRHO)
                dptot = p_ref - Im(2,QPRES) - hdt*Im_src(2,QPRES)
                drhoe_g = rhoe_g_ref - Im(2,QREINT) - hdt*Im_src(2,QREINT)

                dup = un_ref - Im(3,QUN) - hdt*Im_src(3,QUN)
                dptotp = p_ref - Im(3,QPRES) - hdt*Im_src(3,QPRES)


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

                gam_g_ref  = Ip_gc(3,1)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip(1,QUN) - hdt*Ip_src(1,QUN)
                dptotm  = p_ref - Ip(1,QPRES) - hdt*Ip_src(1,QPRES)

                drho = rho_ref - Ip(2,QRHO) - hdt*Ip_src(2,QRHO)
                dptot = p_ref - Ip(2,QPRES) - hdt*Ip_src(2,QPRES)
                drhoe_g = rhoe_g_ref - Ip(2,QREINT) - hdt*Ip_src(2,QREINT)

                dup = un_ref - Ip(3,QUN) - hdt*Ip_src(3,QUN)
                dptotp = p_ref - Ip(3,QPRES) - hdt*Ip_src(3,QPRES)

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
                   qm(i+1,j,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm(i+1,j,k,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                   qm(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                   qm(i+1,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)

                   ! transverse velocities
                   qm(i+1,j,k,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i+1,j,k,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                else if (idir == 2) then
                   qm(i,j+1,k,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm(i,j+1,k,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                   qm(i,j+1,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                   qm(i,j+1,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)

                   ! transverse velocities
                   qm(i,j+1,k,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i,j+1,k,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                else if (idir == 3) then
                   qm(i,j,k+1,QRHO) = max(small_dens, rho_ref +  alphap + alpham + alpha0r)
                   qm(i,j,k+1,QUN) = un_ref + (alphap - alpham)*cc_ev/rho_ev
                   qm(i,j,k+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
                   qm(i,j,k+1,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ev)

                   ! transverse velocities
                   qm(i,j,k+1,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i,j,k+1,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)
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

                   qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr
                   qm(i+1,j,k,QRHO) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                end if

                if (i >= vlo(1)) then

                   qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr
                   qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                end if

             end if
#endif

          end do
       end do
    end do


  end subroutine trace_ppm_rhoe

  subroutine trace_ppm_gammae(lo, hi, &
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
                              reconstruct_state, source_nonzero, &
                              dx, dt)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, NQSRC, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QGAME, QC, QGAMC, &
                                   small_dens, small_pres, &
                                   ppm_type, ppm_temp_fix, &
                                   ppm_reference_eigenvectors
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow
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
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    logical, intent(in) :: source_nonzero(NQSRC)
    logical, intent(in) :: reconstruct_state(NQ)

    ! Local variables
    integer :: i, j, k, n

    real(rt) :: hdt, dtdx

    real(rt) :: sm, sp

    real(rt) :: s(-2:2)
    real(rt) :: Ip(1:3,NQ), Im(1:3,NQ)
    real(rt) :: Ip_src(1:3,NQSRC), Im_src(1:3,NQSRC)
    real(rt) :: Ip_gc(1:3,1), Im_gc(1:3,1)

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
       call castro_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    !$gpu

    hdt = HALF * dt
    dtdx = dt / dx(idir)

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


             if (ppm_temp_fix /= 1) then

                if (idir == 1) then
                   s(:) = qaux(i-2:i+2,j,k,QGAMC)
                else if (idir == 2) then
                   s(:) = qaux(i,j-2:j+2,k,QGAMC)
                else
                   s(:) = qaux(i,j,k-2:k+2,QGAMC)
                end if

                call ppm_reconstruct(s, flatn(i,j,k), sm, sp)

                call ppm_int_profile(sm, sp, s(0), un, cc, dtdx, Ip_gc, Im_gc)
             else

                ! temperature-based PPM
                call ppm_reconstruct_with_eos(Ip, Im, Ip_gc, Im_gc)

             end if


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

                gam_g_ref  = Im_gc(1,1)
                game_ref = Im(1,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                Clag_ref = sqrt(gam_g_ref*p_ref*rho_ref)

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dum = un_ref - Im(1,QUN) - hdt*Im_src(1,QUN)
                dptotm = p_ref - Im(1,QPRES) - hdt*Im_src(1,QPRES)

                dptot = p_ref - Im(2,QPRES) - hdt*Im_src(2,QPRES)

                ! we are treating tau as 1/rho, but we could have reconstructed
                ! it separately
                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Im(1,QRHO) + hdt*Im_src(1,QRHO)/Im(1,QRHO)**2
                dtau  = tau_ref - ONE/Im(2,QRHO) + hdt*Im_src(2,QRHO)/Im(2,QRHO)**2
                dtaup = tau_ref - ONE/Im(3,QRHO) + hdt*Im_src(3,QRHO)/Im(3,QRHO)**2

                dup = un_ref - Im(3,QUN) - hdt*Im_src(3,QUN)
                dptotp = p_ref - Im(3,QPRES) - hdt*Im_src(3,QPRES)


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

                dge   = game_ref - Im(2,QGAME)
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

                gam_g_ref  = Ip_gc(3,1)
                game_ref = Ip(3,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                Clag_ref = sqrt(gam_g_ref*p_ref*rho_ref)

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip(1,QUN) - hdt*Ip_src(1,QUN)
                dptotm  = p_ref - Ip(1,QPRES) - hdt*Ip_src(1,QPRES)

                dptot = p_ref - Ip(2,QPRES) - hdt*Ip_src(2,QPRES)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Ip(1,QRHO) + hdt*Ip_src(1,QRHO)/Ip(1,QRHO)**2
                dtau = tau_ref - ONE/Ip(2,QRHO) + hdt*Ip_src(2,QRHO)/Ip(2,QRHO)**2
                dtaup = tau_ref - ONE/Ip(3,QRHO) + hdt*Ip_src(3,QRHO)/Ip(3,QRHO)**2

                dup = un_ref - Ip(3,QUN) - hdt*Ip_src(3,QUN)
                dptotp = p_ref - Ip(3,QPRES) - hdt*Ip_src(3,QPRES)

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

                dge = game_ref - Ip(2,QGAME)
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
                   qm(i+1,j,k,QRHO) = max(small_dens, ONE/tau_s)

                   qm(i+1,j,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm(i+1,j,k,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

                   qm(i+1,j,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QPRES )/(qm(i+1,j,k,QGAME) - ONE)

                   ! transverse velocities
                   qm(i+1,j,k,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i+1,j,k,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                else if (idir == 2) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm(i,j+1,k,QRHO) = max(small_dens, ONE/tau_s)

                   qm(i,j+1,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm(i,j+1,k,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

                   qm(i,j+1,k,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                   qm(i,j+1,k,QREINT) = qm(i,j+1,k,QPRES )/(qm(i,j+1,k,QGAME) - ONE)

                   ! transverse velocities
                   qm(i,j+1,k,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i,j+1,k,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                else if (idir == 3) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm(i,j,k+1,QRHO) = max(small_dens, ONE/tau_s)

                   qm(i,j,k+1,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm(i,j,k+1,QPRES) = max(small_pres, p_ref - (alphap + alpham)*Clag_ev**2)

                   qm(i,j,k+1,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
                   qm(i,j,k+1,QREINT) = qm(i,j,k+1,QPRES )/(qm(i,j,k+1,QGAME) - ONE)

                   ! transverse velocities
                   qm(i,j,k+1,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i,j,k+1,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

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
                   qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr
                   qm(i+1,j,k,QRHO) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                end if

                if (i >= vlo(1)) then
                   qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr
                   qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                end if

             endif
#endif

          end do
       end do
    end do

  end subroutine trace_ppm_gammae


  subroutine trace_ppm_temp(lo, hi, &
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
                            reconstruct_state, source_nonzero, &
                            dx, dt)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, NQSRC, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QTEMP, QGAME, QC, QGAMC, QFS, QFX, &
                                   small_dens, small_pres, &
                                   ppm_type, ppm_temp_fix, &
                                   ppm_reference_eigenvectors
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
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
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    logical, intent(in) :: source_nonzero(NQSRC)
    logical, intent(in) :: reconstruct_state(NQ)

    ! Local variables
    type(eos_t) :: eos_state

    real(rt) :: hdt, dtdx

    real(rt) :: sm, sp

    real(rt) :: s(-2:2)
    real(rt) :: Ip(1:3,NQ), Im(1:3,NQ)
    real(rt) :: Ip_src(1:3,NQSRC), Im_src(1:3,NQSRC)
    real(rt) :: Ip_gc(1:3,1), Im_gc(1:3,1)

    integer :: i, j, k, n

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
       call castro_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    !$gpu

    hdt = HALF * dt
    dtdx = dt / dx(idir)

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


             if (ppm_temp_fix /= 1) then
                
                if (idir == 1) then
                   s(:) = qaux(i-2:i+2,j,k,QGAMC)
                else if (idir == 2) then
                   s(:) = qaux(i,j-2:j+2,k,QGAMC)
                else
                   s(:) = qaux(i,j,k-2:k+2,QGAMC)
                end if

                call ppm_reconstruct(s, flatn(i,j,k), sm, sp)

                call ppm_int_profile(sm, sp, s(0), un, cc, dtdx, Ip_gc, Im_gc)
             else

                ! temperature-based PPM
                call ppm_reconstruct_with_eos(Ip, Im, Ip_gc, Im_gc)

             end if


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
                temp_ref = Im(1,QTEMP)

                tau_ref  = ONE/Im(1,QRHO)

                gam_g_ref  = Im_gc(1,1)
                game_ref = Im(1,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref

                ! Note: for the transverse velocities, the jump is carried
                !       only by the u wave (the contact)


                ! we also add the sources here so they participate in the tracing
                dum = un_ref - Im(1,QUN) - hdt*Im_src(1,QUN)
                dptotm = p_ref - Im(1,QPRES) - hdt*Im_src(1,QPRES)

                drho = rho_ref - Im(2,QRHO) - hdt*Im_src(2,QRHO)
                dptot = p_ref - Im(2,QPRES) - hdt*Im_src(2,QPRES)

                ! TODO: need to figure sources for this out...
                dTm = temp_ref - Im(1,QTEMP)
                dT0 = temp_ref - Im(2,QTEMP)
                dTp = temp_ref - Im(3,QTEMP)

                ! we are treating tau as 1/rho, but we could have reconstructed
                ! it separately
                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Im(1,QRHO) + hdt*Im_src(1,QRHO)/Im(1,QRHO)**2
                dtau  = tau_ref - ONE/Im(2,QRHO) + hdt*Im_src(2,QRHO)/Im(2,QRHO)**2
                dtaup = tau_ref - ONE/Im(3,QRHO) + hdt*Im_src(3,QRHO)/Im(3,QRHO)**2

                dup = un_ref - Im(3,QUN) - hdt*Im_src(3,QUN)
                dptotp = p_ref - Im(3,QPRES) - hdt*Im_src(3,QPRES)


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
                temp_ref = Ip(3,QTEMP)

                tau_ref  = ONE/Ip(3,QRHO)

                gam_g_ref  = Ip_gc(3,1)
                game_ref = Ip(3,QGAME)

                rho_ref = max(rho_ref, small_dens)
                p_ref = max(p_ref, small_pres)

                ! For tracing (optionally)
                csq_ref = gam_g_ref*p_ref/rho_ref
                cc_ref = sqrt(csq_ref)
                Clag_ref = rho_ref*cc_ref

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum = un_ref - Ip(1,QUN) - hdt*Ip_src(1,QUN)
                dptotm  = p_ref - Ip(1,QPRES) - hdt*Ip_src(1,QPRES)

                drho = rho_ref - Ip(2,QRHO) - hdt*Ip_src(2,QRHO)
                dptot = p_ref - Ip(2,QPRES) - hdt*Ip_src(2,QPRES)

                dTm = temp_ref - Ip(1,QTEMP)
                dT0 = temp_ref - Ip(2,QTEMP)
                dTp = temp_ref - Ip(3,QTEMP)

                ! since d(rho)/dt = S_rho, d(tau**{-1})/dt = S_rho, so d(tau)/dt = -S_rho*tau**2
                dtaum = tau_ref - ONE/Ip(1,QRHO) + hdt*Ip_src(1,QRHO)/Ip(1,QRHO)**2
                dtau = tau_ref - ONE/Ip(2,QRHO) + hdt*Ip_src(2,QRHO)/Ip(2,QRHO)**2
                dtaup = tau_ref - ONE/Ip(3,QRHO) + hdt*Ip_src(3,QRHO)/Ip(3,QRHO)**2

                dup = un_ref - Ip(3,QUN) - hdt*Ip_src(3,QUN)
                dptotp = p_ref - Ip(3,QPRES) - hdt*Ip_src(3,QPRES)

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
                   qm(i+1,j,k,QRHO) = max(small_dens, ONE/tau_s)

                   qm(i+1,j,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm(i+1,j,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once
                   ! we do the species
                   qm(i+1,j,k,QPRES) = small_pres ! just to make it defined

                   ! transverse velocities
                   qm(i+1,j,k,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i+1,j,k,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                else if (idir == 2) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm(i,j+1,k,QRHO) = max(small_dens, ONE/tau_s)

                   qm(i,j+1,k,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm(i,j+1,k,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once
                   ! we do the species
                   qm(i,j+1,k,QPRES) = small_pres ! just to make it defined

                   ! transverse velocities
                   qm(i,j+1,k,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i,j+1,k,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                else if (idir == 3) then
                   tau_s = tau_ref + alphap + alpham + alpha0r
                   qm(i,j,k+1,QRHO) = max(small_dens, ONE/tau_s)

                   qm(i,j,k+1,QUN) = un_ref + (alpham - alphap)*Clag_ev
                   qm(i,j,k+1,QTEMP) = temp_ref + (-Clag_ev**2 - rho_ev**2*p_r)*alpham/p_T + &
                        rho_ev**2*p_r*alpha0r/p_T - (-Clag_ev**2 - rho_ev**2*p_r)*alphap/p_T

                   ! we defer getting the pressure until later, once
                   ! we do the species
                   qm(i,j,k+1,QPRES) = small_pres ! just to make it defined

                   ! transverse velocities
                   qm(i,j,k+1,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i,j,k+1,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

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
                   qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr
                   qm(i+1,j,k,QRHO) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                end if

                if (i >= vlo(1)) then
                   qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr
                   qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
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
                eos_state%T     = qp(i,j,k,QTEMP)
                eos_state%rho   = qp(i,j,k,QRHO)
                eos_state%xn(:) = qp(i,j,k,QFS:QFS-1+nspec)
                eos_state%aux(:) = qp(i,j,k,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                qp(i,j,k,QPRES) = eos_state%p
                qp(i,j,k,QREINT) = qp(i,j,k,QRHO)*eos_state%e

                qp(i,j,k,QPRES) = max(qp(i,j,k,QPRES), small_pres)
             end if

             if (idir == 1 .and. i <= vhi(1)) then

                ! minus face
                eos_state%T     = qm(i+1,j,k,QTEMP)
                eos_state%rho   = qm(i+1,j,k,QRHO)
                eos_state%xn(:) = qm(i+1,j,k,QFS:QFS-1+nspec)
                eos_state%aux(:) = qm(i+1,j,k,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                qm(i+1,j,k,QPRES) = max(small_pres, eos_state%p)
                qm(i+1,j,k,QREINT) = qm(i+1,j,k,QRHO)*eos_state%e

             else if (idir == 2 .and. j <= vhi(2)) then

                ! minus face
                eos_state%T     = qm(i,j+1,k,QTEMP)
                eos_state%rho   = qm(i,j+1,k,QRHO)
                eos_state%xn(:) = qm(i,j+1,k,QFS:QFS-1+nspec)
                eos_state%aux(:) = qm(i,j+1,k,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                qm(i,j+1,k,QPRES) = max(small_pres, eos_state%p)
                qm(i,j+1,k,QREINT) = qm(i,j+1,k,QRHO)*eos_state%e

             else if (idir == 3 .and. k <= vhi(3)) then

                ! minus face
                eos_state%T     = qm(i,j,k+1,QTEMP)
                eos_state%rho   = qm(i,j,k+1,QRHO)
                eos_state%xn(:) = qm(i,j,k+1,QFS:QFS-1+nspec)
                eos_state%aux(:) = qm(i,j,k+1,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                qm(i,j,k+1,QPRES) = max(small_pres, eos_state%p)
                qm(i,j,k+1,QREINT) = qm(i,j,k+1,QRHO)*eos_state%e

             end if

          end do
       end do
    end do

  end subroutine trace_ppm_temp

end module trace_ppm_module
