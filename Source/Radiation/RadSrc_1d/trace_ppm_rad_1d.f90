! These routines do the characteristic tracing under the parabolic
! profiles in each zone to the edge / half-time.

module trace_ppm_rad_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public trace_ppm_rad

contains

  subroutine trace_ppm_rad(q, qaux, flatn, qd_l1, qd_h1, &
                           dloga, dloga_l1, dloga_h1, &
                           srcQ, src_l1, src_h1,&
                           qxm, qxp, qpd_l1, qpd_h1, &
                           ilo, ihi, domlo, domhi, dx, dt)

    use bl_constants_module
    use meth_params_module, only : QVAR, QRHO, QU, QREINT, QPRES, &
         NQ, NQAUX, qrad, qradhi, qptot, qreitot, &
         QGAME, QGAMC, QGAMCG, QC, QCG, QLAMS, &
         small_dens, small_pres, fix_mass_flux, &
         ppm_type, ppm_trace_sources, ppm_temp_fix, &
         ppm_predict_gammae, ppm_reference_eigenvectors, &
         npassive, qpass_map
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow
    use rad_params_module, only : ngroups
    use ppm_module, only : ppm_reconstruct, ppm_int_profile

    implicit none

    integer ilo,ihi
    integer domlo(1),domhi(1)
    integer    qd_l1,   qd_h1
    integer dloga_l1,dloga_h1
    integer   qpd_l1,  qpd_h1
    integer   src_l1,  src_h1
    real(rt)         dx, dt
    real(rt)             q( qd_l1: qd_h1, NQ)
    real(rt)         :: qaux(qd_l1:qd_h1, NQAUX)
    real(rt)         flatn(qd_l1:qd_h1)
    real(rt)          srcQ(src_l1:src_h1,QVAR)
    real(rt)         dloga(dloga_l1:dloga_h1)

    real(rt)          qxm( qpd_l1: qpd_h1, NQ)
    real(rt)          qxp( qpd_l1: qpd_h1, NQ)

    ! Local variables
    integer :: i, g
    integer :: n, ipassive

    real(rt)         :: hdt,dtdx

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

    real(rt)         :: cc, csq, cgassq, Clag
    real(rt)         :: rho, u, p, rhoe_g, h_g, tau
    real(rt)         :: ptot, gam_g, game

    real(rt)         :: drho, dptot, drhoe_g
    real(rt)         :: dge, dtau
    real(rt)         :: dup, dptotp
    real(rt)         :: dum, dptotm

    real(rt)         :: rho_ref, u_ref, p_ref, rhoe_g_ref, h_g_ref
    real(rt)         :: ptot_ref
    real(rt)         :: tau_ref

    real(rt)         :: gam_g_ref, game_ref, gfactor

    real(rt)         :: alpham, alphap, alpha0r, alpha0e_g
    real(rt)         :: sourcr, sourcp, source, courn, eta, dlogatmp

    real(rt)         :: tau_s

    logical :: fix_mass_flux_lo, fix_mass_flux_hi

    real(rt)        , dimension(0:ngroups-1) :: er, der, alphar, sourcer, qrtmp, hr
    real(rt)        , dimension(0:ngroups-1) :: lam0, lamp, lamm

    real(rt)        , dimension(0:ngroups-1) :: er_ref
    real(rt)         :: er_foo

    real(rt)        , allocatable :: Ip(:,:,:)
    real(rt)        , allocatable :: Im(:,:,:)

    real(rt)        , allocatable :: Ip_src(:,:,:)
    real(rt)        , allocatable :: Im_src(:,:,:)

    ! temporary interface values of the parabola
    real(rt), allocatable :: sxm(:), sxp(:)


    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (ilo == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (ihi == domhi(1))

    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
       call bl_error("Error:: RadHydro_1d.f90 :: trace_ppm_rad")
    end if

    if (ppm_temp_fix > 0) then
       call bl_error("ERROR: ppm_temp_fix > 0 not implemented with radiation")
    endif

    if (ppm_reference_eigenvectors == 1) then
       call bl_error("ERROR: ppm_reference_eigenvectors not implemented with radiation")
    endif

    hdt = HALF * dt
    dtdx = dt/dx

    allocate(Ip(ilo-1:ihi+1,3, NQ))
    allocate(Im(ilo-1:ihi+1,3, NQ))

    if (ppm_trace_sources == 1) then
       allocate(Ip_src(ilo-1:ihi+1,3,QVAR))
       allocate(Im_src(ilo-1:ihi+1,3,QVAR))
    endif

    allocate(sxm(qd_l1:qd_h1))
    allocate(sxp(qd_l1:qd_h1))


    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    !
    ! We first fill the Im and Ip arrays -- these are the averages of
    ! the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    !
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left state at that interface).
    !
    ! The indices are: Ip(i, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.


    ! Compute Ip and Im -- this does the parabolic reconstruction,
    ! limiting, and returns the integral of each profile under each
    ! wave to each interface
    do n = 1, NQ
       call ppm_reconstruct(q(:,n), qd_l1, qd_h1, &
                            flatn, &
                            sxm, sxp, &
                            ilo, ihi, dx)

       call ppm_int_profile(q(:,n), qd_l1, qd_h1, &
                            q(:,QU), qaux(:,QC), &
                            sxm, sxp, &
                            Ip(:,:,n), Im(:,:,n), &
                            ilo, ihi, dx, dt)
    end do

    if (ppm_trace_sources == 1) then
       do n = 1, QVAR
          call ppm_reconstruct(srcQ(:,n), src_l1, src_h1, &
                               flatn, &
                               sxm, sxp, &
                               ilo, ihi, dx)

          call ppm_int_profile(srcQ(:,n), src_l1, src_h1, &
                               q(:,QU), qaux(:,QC), &
                               sxm, sxp, &
                               Ip_src(:,:,n), Im_src(:,:,n), &
                               ilo, ihi, dx, dt)
       enddo
    endif

    deallocate(sxm, sxp)

    !-------------------------------------------------------------------------
    ! x-direction
    !-------------------------------------------------------------------------

    ! Trace to left and right edges using upwind PPM
    do i = ilo-1, ihi+1

       do g=0, ngroups-1
          lam0(g) = qaux(i,QLAMS+g)
          lamp(g) = qaux(i,QLAMS+g)
          lamm(g) = qaux(i,QLAMS+g)
       end do

       ! cgassq is the gas soundspeed **2
       ! cc is the total soundspeed **2 (gas + radiation)
       cgassq = qaux(i,QCG)**2
       cc = qaux(i,QC)
       csq = cc**2

       rho = q(i,QRHO)
       tau = ONE/rho  ! should not be needed once reference ev is implemented

       u = q(i,QU)
       p = q(i,QPRES)
       rhoe_g = q(i,QREINT)
       h_g = ( (p + rhoe_g)/rho )/csq

       ptot = q(i,qptot)

       er(:) = q(i,qrad:qradhi)
       hr(:) = (lam0+1.e0_rt)*er/rho

       Clag = rho*cc

       gam_g = qaux(i,QGAMCG)
       game = q(i,QGAME)

       !----------------------------------------------------------------------
       ! plus state on face i
       !----------------------------------------------------------------------

       ! set the reference state
       ! this will be the fastest moving state to the left --
       ! this is the method that Miller & Colella and Colella &
       ! Woodward use
       rho_ref = Im(i,1,QRHO)
       u_ref = Im(i,1,QU)

       p_ref = Im(i,1,QPRES)
       rhoe_g_ref = Im(i,1,QREINT)

       tau_ref = ONE/Im(i,1,QRHO)

       game_ref = Im(i,1,QGAME)

       ptot_ref = Im(i,1,qptot)

       er_ref(:) = Im(i,1,qrad:qradhi)

       rho_ref = max(rho_ref,small_dens)
       p_ref = max(p_ref,small_pres)

       ! *m are the jumps carried by u-c
       ! *p are the jumps carried by u+c

       ! note, for radiation, we use ptot here.  This is
       !    ptot = pgas + lambda Er
       ! (where lambda is the flux limiter)
       ! This is what appears in the betas (l . dq) directly

       dum    = u_ref    - Im(i,1,QU)
       dptotm = ptot_ref - Im(i,1,qptot)

       drho    = rho_ref    - Im(i,2,QRHO)
       drhoe_g = rhoe_g_ref - Im(i,2,QREINT)
       dptot   = ptot_ref   - Im(i,2,qptot)
       dtau    = tau_ref    - ONE/Im(i,2,QRHO)
       der(:)  = er_ref(:)  - Im(i,2,qrad:qradhi)

       dup    = u_ref    - Im(i,3,QU)
       dptotp = ptot_ref - Im(i,3,qptot)

       ! if we are doing source term tracing, then we add the force to
       ! the velocity here, otherwise we will deal with this later
       if (ppm_trace_sources == 1) then
          dum = dum - hdt*Im_src(i,1,QU)
          dup = dup - hdt*Im_src(i,3,QU)
       endif

       ! optionally use the reference state in evaluating the
       ! eigenvectors

       if (ppm_predict_gammae == 0) then

          ! (rho, u, p, (rho e)) eigensystem

          ! these are analogous to the beta's from the original
          ! PPM paper (except we work with rho instead of tau).
          ! This is simply (l . dq), where dq = qref - I(q)

          alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
          alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
          alpha0r = drho - dptot/csq
          alpha0e_g = drhoe_g - dptot*h_g  ! note h_g has a 1/c**2 in it

       else

          ! (tau, u, p, game) eigensystem

          ! this is the way things were done in the original PPM
          ! paper -- here we work with tau in the characteristic
          ! system.

          alpham = HALF*( dum - dptotm/Clag)/Clag
          alphap = HALF*(-dup - dptotp/Clag)/Clag
          alpha0r = dtau + dptot/Clag**2

          dge = game_ref - Im(i,2,QGAME)
          gfactor = (game - ONE)*(game - gam_g)
          alpha0e_g = gfactor*dptot/(tau*Clag**2) + dge

       endif ! which tracing method

       alphar(:) = der(:) - dptot/csq*hr

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
          alphar(:) = ZERO
       else if (u < ZERO) then
          alpha0r = -alpha0r
          alpha0e_g = -alpha0e_g
          alphar(:) = -alphar(:)
       else
          alpha0r = -HALF*alpha0r
          alpha0e_g = -HALF*alpha0e_g
          alphar(:) = -HALF*alphar(:)
       endif

       ! the final interface states are just
       ! q_s = q_ref - sum (l . dq) r
       if (i >= ilo) then

          if (ppm_predict_gammae == 0) then
             qxp(i,QRHO)   = rho_ref + alphap + alpham + alpha0r
             qxp(i,QU)     = u_ref + (alphap - alpham)*cc/rho
             qxp(i,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
             qxp(i,QPRES)  = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

             qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
             qxp(i,qrad:qradhi) = qrtmp

             qxp(i,qptot) = ptot_ref + (alphap + alpham)*csq
             qxp(i,qreitot) = qxp(i,QREINT) + sum(qrtmp)

          else

             tau_s = tau_ref + alphap + alpham + alpha0r
             qxp(i,QRHO)   = ONE/tau_s

             qxp(i,QU)     = u_ref + (alpham - alphap)*Clag
             qxp(i,QPRES)  = p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamp(:)*alphar(:))

             qxp(i,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
             qxp(i,QREINT) = qxp(i,QPRES )/(qxp(i,QGAME) - ONE)

             qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
             qxp(i,qrad:qradhi) = qrtmp

             qxp(i,qptot) = ptot_ref - (alphap + alpham)*Clag**2
             qxp(i,qreitot) = qxp(i,QREINT) + sum(qrtmp)

          endif

          ! enforce small_*
          qxp(i,QRHO) = max(small_dens,qxp(i,QRHO))
          qxp(i,QPRES) = max(small_pres,qxp(i,QPRES))

          ! add source terms
          qxp(i  ,QRHO  )  = qxp(i,QRHO   ) + hdt*srcQ(i,QRHO)
          qxp(i  ,QRHO  )  = max(small_dens,qxp(i,QRHO))
          qxp(i  ,QREINT)  = qxp(i,QREINT ) + hdt*srcQ(i,QREINT)
          qxp(i  ,QPRES )  = qxp(i,QPRES  ) + hdt*srcQ(i,QPRES)
          qxp(i  ,qptot )  = qxp(i,qptot  ) + hdt*srcQ(i,QPRES)
          qxp(i  ,qreitot) = qxp(i,qreitot) + hdt*srcQ(i,QREINT)

          ! add traced source terms as needed
          if (ppm_trace_sources == 0) then
             qxp(i  ,QU) = qxp(i,QU) + hdt*srcQ(i,QU)
          endif

          do g=0, ngroups-1
             if (qxp(i,qrad+g) < ZERO) then
                er_foo = - qxp(i,qrad+g)
                qxp(i,qrad+g) = ZERO
                qxp(i,qptot) = qxp(i,qptot) + lamp(g) * er_foo
                qxp(i,qreitot) = qxp(i,qreitot) + er_foo
             end if
          end do

          if ( qxp(i,QPRES) < ZERO ) then
             qxp(i,QPRES) = p
          end if

       end if


       !----------------------------------------------------------------------
       ! minus state on face i+1
       !----------------------------------------------------------------------

       ! set the reference state
       ! this will be the fastest moving state to the right
       rho_ref  = Ip(i,3,QRHO)
       u_ref    = Ip(i,3,QU)

       p_ref      = Ip(i,3,QPRES)
       rhoe_g_ref = Ip(i,3,QREINT)

       tau_ref = ONE/Ip(i,3,QRHO)

       game_ref = Ip(i,3,QGAME)

       ptot_ref = Ip(i,3,qptot)

       er_ref(:) = Ip(i,3,qrad:qradhi)

       rho_ref = max(rho_ref,small_dens)
       p_ref = max(p_ref,small_pres)

       !  *m are the jumps carried by u-c
       !  *p are the jumps carried by u+c

       dum    = u_ref    - Ip(i,1,QU)
       dptotm = ptot_ref - Ip(i,1,qptot)

       drho    = rho_ref    - Ip(i,2,QRHO)
       dptot   = ptot_ref   - Ip(i,2,qptot)
       drhoe_g = rhoe_g_ref - Ip(i,2,QREINT)
       dtau    = tau_ref    - ONE/Ip(i,2,QRHO)
       der(:)  = er_ref(:)  - Ip(i,2,qrad:qradhi)

       dup    = u_ref    - Ip(i,3,QU)
       dptotp = ptot_ref - Ip(i,3,qptot)

       ! if we are doing source term tracing, then we add the force to
       ! the velocity here, otherwise we will deal with this in the
       ! trans_X routines
       if (ppm_trace_sources == 1) then
          dum = dum - hdt*Ip_src(i,1,QU)
          dup = dup - hdt*Ip_src(i,3,QU)
       endif

       ! optionally use the reference state in evaluating the
       ! eigenvectors

       if (ppm_predict_gammae == 0) then

          ! (rho, u, p, (rho e)) eigensystem

          ! these are analogous to the beta's from the original
          ! PPM paper (except we work with rho instead of tau).
          ! This is simply (l . dq), where dq = qref - I(q)
          alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
          alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
          alpha0r = drho - dptot/csq
          alpha0e_g = drhoe_g - dptot*h_g

       else

          ! (tau, u, p, game) eigensystem

          ! this is the way things were done in the original PPM
          ! paper -- here we work with tau in the characteristic
          ! system.

          alpham = HALF*( dum - dptotm/Clag)/Clag
          alphap = HALF*(-dup - dptotp/Clag)/Clag
          alpha0r = dtau + dptot/Clag**2

          dge = game_ref - Ip(i,2,QGAME)
          gfactor = (game - ONE)*(game - gam_g)
          alpha0e_g = gfactor*dptot/(tau*Clag**2) + dge

       endif

       alphar(:) = der(:)- dptot/csq*hr

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
          alphar(:) = -alphar(:)
       else if (u < ZERO) then
          alpha0r = ZERO
          alpha0e_g = ZERO
          alphar(:) = ZERO
       else
          alpha0r = -HALF*alpha0r
          alpha0e_g = -HALF*alpha0e_g
          alphar(:) = -HALF*alphar(:)
       endif

       ! the final interface states are just
       ! q_s = q_ref - sum (l . dq) r
       if (i <= ihi) then

          if (ppm_predict_gammae == 0) then
             qxm(i+1,QRHO)   = rho_ref + alphap + alpham + alpha0r
             qxm(i+1,QU)     = u_ref + (alphap - alpham)*cc/rho
             qxm(i+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
             qxm(i+1,QPRES)  = p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))

             qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
             qxm(i+1,qrad:qradhi) = qrtmp

             qxm(i+1,qptot) = ptot_ref + (alphap + alpham)*csq
             qxm(i+1,qreitot) = qxm(i+1,QREINT) + sum(qrtmp)

          else

             tau_s = tau_ref + alphap + alpham + alpha0r
             qxm(i+1,QRHO)   = ONE/tau_s

             qxm(i+1,QU)     = u_ref + (alpham - alphap)*Clag
             qxm(i+1,QPRES)  = p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamm(:)*alphar(:))

             qxm(i+1,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
             qxm(i+1,QREINT) = qxm(i+1,QPRES )/(qxm(i+1,QGAME) - ONE)

             qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
             qxm(i+1,qrad:qradhi) = qrtmp

             qxm(i+1,qptot) = ptot_ref - (alphap + alpham)*Clag**2
             qxm(i+1,qreitot) = qxm(i+1,QREINT) + sum(qrtmp)

          endif

          ! enforce small_*
          qxm(i+1,QRHO) = max(qxm(i+1,QRHO),small_dens)
          qxm(i+1,QPRES) = max(qxm(i+1,QPRES),small_pres)

          ! add source terms
          qxm(i+1,QRHO   ) = qxm(i+1,QRHO   ) + hdt*srcQ(i,QRHO)
          qxm(i+1,QRHO   ) = max(small_dens, qxm(i+1,QRHO))
          qxm(i+1,QREINT ) = qxm(i+1,QREINT ) + hdt*srcQ(i,QREINT)
          qxm(i+1,QPRES  ) = qxm(i+1,QPRES  ) + hdt*srcQ(i,QPRES)
          qxm(i+1,qptot  ) = qxm(i+1,qptot  ) + hdt*srcQ(i,QPRES)
          qxm(i+1,qreitot) = qxm(i+1,qreitot) + hdt*srcQ(i,QREINT)

          ! add traced source terms as needed
          if (ppm_trace_sources == 0) then
             qxm(i+1,QU) = qxm(i+1,QU) + hdt*srcQ(i,QU)
          endif

          do g=0, ngroups-1
             if (qxm(i+1,qrad+g) < ZERO) then
                er_foo = - qxm(i+1,qrad+g)
                qxm(i+1,qrad+g) = ZERO
                qxm(i+1,qptot) = qxm(i+1,qptot) + lamm(g) * er_foo
                qxm(i+1,qreitot) = qxm(i+1,qreitot) + er_foo
             end if
          end do

          if ( qxm(i+1,QPRES) < ZERO ) then
             qxm(i+1,QPRES) = p
          end if

       end if


       !----------------------------------------------------------------------
       ! geometry source terms
       !----------------------------------------------------------------------

       if (dloga(i) /= 0) then
          courn = dtdx*(cc+abs(u))
          eta = (ONE-courn)/(cc*dt*abs(dloga(i)))
          dlogatmp = min(eta,ONE)*dloga(i)
          sourcr = -HALF*dt*rho*dlogatmp*u
          sourcp = sourcr*cgassq
          source = sourcp*h_g
          sourcer(:) = -HALF*dt*dlogatmp*u*(lam0(:)+ONE)*er(:)

          if (i <= ihi) then
             qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + sourcr
             qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
             qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + sourcp
             qxm(i+1,QREINT) = qxm(i+1,QREINT) + source
             qxm(i+1,qrad:qradhi) = qxm(i+1,qrad:qradhi) + sourcer(:)
             !           qxm(i+1,qptot ) = sum(lamm(:)*qxm(i+1,qrad:qradhi)) + qxm(i+1,QPRES)
             qxm(i+1,qptot) = qxm(i+1,qptot) + sum(lamm(:)*sourcer(:)) + sourcp
             qxm(i+1,qreitot) = sum(qxm(i+1,qrad:qradhi))  + qxm(i+1,QREINT)
          end if
          if (i >= ilo) then
             qxp(i  ,QRHO  ) = qxp(i  ,QRHO  ) + sourcr
             qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
             qxp(i  ,QPRES ) = qxp(i  ,QPRES ) + sourcp
             qxp(i  ,QREINT) = qxp(i  ,QREINT) + source
             qxp(i  ,qrad:qradhi) = qxp(i  ,qrad:qradhi) + sourcer(:)
             !           qxp(i  ,qptot ) = sum(lamp(:)*qxp(i,qrad:qradhi)) + qxp(i,QPRES)
             qxp(i,qptot) = qxp(i,qptot) + sum(lamp(:)*sourcer(:)) + sourcp
             qxp(i  ,qreitot) = sum(qxp(i,qrad:qradhi))  + qxp(i,QREINT)
          end if
       endif

    end do

    ! Enforce constant mass flux rate if specified
    if (fix_mass_flux_lo) then
       qxm(ilo,QRHO   ) = q(domlo(1)-1,QRHO)
       qxm(ilo,QU     ) = q(domlo(1)-1,QU  )
       qxm(ilo,QPRES  ) = q(domlo(1)-1,QPRES)
       qxm(ilo,QREINT ) = q(domlo(1)-1,QREINT)
       qxm(ilo,qrad:qradhi) = q(domlo(1)-1,qrad:qradhi)
       qxm(ilo,qptot  ) = q(domlo(1)-1,qptot)
       qxm(ilo,qreitot) = q(domlo(1)-1,qreitot)
    end if

    ! Enforce constant mass flux rate if specified
    if (fix_mass_flux_hi) then
       qxp(ihi+1,QRHO   ) = q(domhi(1)+1,QRHO)
       qxp(ihi+1,QU     ) = q(domhi(1)+1,QU  )
       qxp(ihi+1,QPRES  ) = q(domhi(1)+1,QPRES)
       qxp(ihi+1,QREINT ) = q(domhi(1)+1,QREINT)
       qxp(ihi+1,qrad:qradhi) = q(domhi(1)+1,qrad:qradhi)
       qxp(ihi+1,qptot  ) = q(domhi(1)+1,qptot)
       qxp(ihi+1,qreitot) = q(domhi(1)+1,qreitot)
    end if


    !-------------------------------------------------------------------------
    ! Now do the passively advected quantities
    !-------------------------------------------------------------------------

    ! We do all passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! plus state on face i
       do i = ilo, ihi+1
          u = q(i,QU)

          ! We want to do
          !
          ! q_l = q_ref - Proj{(q_ref - I)}
          !
          ! and Proj{} represents the characteristic projection.
          ! But for these, there is only 1-wave that matters, the u
          ! wave, so no projection is needed.  Since we are not
          ! projecting, the reference state doesn't matter

          if (u > ZERO) then
             qxp(i,n) = q(i,n)    ! we might want to change this to
                                  ! the limit of the parabola
          else if (u < ZERO) then
             qxp(i,n) = Im(i,2,n)
          else
             qxp(i,n) = q(i,n) + HALF*(Im(i,2,n) - q(i,n))
          endif
       enddo

       ! minus state on face i+1
       do i = ilo-1, ihi
          u = q(i,QU)

          if (u > ZERO) then
             qxm(i+1,n) = Ip(i,2,n)
          else if (u < ZERO) then
             qxm(i+1,n) = q(i,n)
          else
             qxm(i+1,n) = q(i,n) + HALF*(Ip(i,2,n) - q(i,n))
          endif
       enddo

       if (fix_mass_flux_hi) qxp(ihi+1,n) = q(ihi+1,n)
       if (fix_mass_flux_lo) qxm(ilo,n) = q(ilo-1,n)

    enddo

    deallocate(Ip,Im)
    if (ppm_trace_sources == 1) then
       deallocate(Ip_src,Im_src)
    endif

  end subroutine trace_ppm_rad

end module trace_ppm_rad_module
