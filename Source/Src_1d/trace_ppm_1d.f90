! These routines do the characteristic tracing under the parabolic
! profiles in each zone to the edge / half-time.

module trace_ppm_module

  use amrex_fort_module, only : rt => c_real
  implicit none

  private

  public trace_ppm

contains

  subroutine trace_ppm(q,c,flatn,gamc,qd_l1,qd_h1, &
                       dloga,dloga_l1,dloga_h1, &
                       srcQ,src_l1,src_h1,&
                       qxm,qxp,qpd_l1,qpd_h1, &
                       ilo,ihi,domlo,domhi,dx,dt)

    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use bl_constants_module
    use meth_params_module, only : QVAR, QRHO, QU, QREINT, QPRES, &
         QTEMP, QFS, QFX, QGAME, &
         small_dens, small_pres, fix_mass_flux, &
         ppm_type, ppm_trace_sources, ppm_temp_fix, &
         ppm_predict_gammae, ppm_reference_eigenvectors, &
         npassive, qpass_map
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow
    use ppm_module, only : ppm

    use amrex_fort_module, only : rt => c_real
    implicit none

    integer ilo,ihi
    integer domlo(1),domhi(1)
    integer    qd_l1,   qd_h1
    integer dloga_l1,dloga_h1
    integer   qpd_l1,  qpd_h1
    integer   src_l1,  src_h1
    real(rt)         dx, dt

    real(rt)             q( qd_l1: qd_h1,QVAR)
    real(rt)          srcQ(src_l1:src_h1,QVAR)
    real(rt)          gamc(qd_l1:qd_h1)
    real(rt)         flatn(qd_l1:qd_h1)
    real(rt)             c(qd_l1:qd_h1)
    real(rt)         dloga(dloga_l1:dloga_h1)

    real(rt)          qxm( qpd_l1: qpd_h1,QVAR)
    real(rt)          qxp( qpd_l1: qpd_h1,QVAR)

    ! Local variables
    integer :: i, iwave
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

    real(rt)         :: cc, csq, Clag
    real(rt)         :: rho, u, p, rhoe_g, h_g
    real(rt)         :: gam_g, game

    real(rt)         :: drho, dptot, drhoe_g
    real(rt)         :: dge, dtau
    real(rt)         :: dup, dptotp
    real(rt)         :: dum, dptotm

    real(rt)         :: rho_ref, u_ref, p_ref, rhoe_g_ref, h_g_ref
    real(rt)         :: tau_ref

    real(rt)         :: cc_ref, csq_ref, Clag_ref, gam_g_ref, game_ref, gfactor
    real(rt)         :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, h_g_ev, tau_ev

    real(rt)         :: alpham, alphap, alpha0r, alpha0e_g
    real(rt)         :: sourcr, sourcp, source, courn, eta, dlogatmp

    real(rt)         :: tau_s

    logical :: fix_mass_flux_lo, fix_mass_flux_hi

    real(rt)        , allocatable :: Ip(:,:,:)
    real(rt)        , allocatable :: Im(:,:,:)

    real(rt)        , allocatable :: Ip_src(:,:,:)
    real(rt)        , allocatable :: Im_src(:,:,:)

    real(rt)        , allocatable :: Ip_gc(:,:,:)
    real(rt)        , allocatable :: Im_gc(:,:,:)

    type (eos_t) :: eos_state

    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (ilo == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (ihi == domhi(1))

    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
       call bl_error("Error:: ppm_1d.f90 :: trace_ppm")
    end if

    hdt = HALF * dt
    dtdx = dt/dx

    allocate(Ip(ilo-1:ihi+1,3,QVAR))
    allocate(Im(ilo-1:ihi+1,3,QVAR))

    allocate(Ip_gc(ilo-1:ihi+1,3,1))
    allocate(Im_gc(ilo-1:ihi+1,3,1))

    if (ppm_trace_sources == 1) then
       allocate(Ip_src(ilo-1:ihi+1,3,QVAR))
       allocate(Im_src(ilo-1:ihi+1,3,QVAR))
    endif

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
    do n = 1, QVAR
       call ppm(q(:,n),qd_l1,qd_h1, &
                q(:,QU),c, &
                flatn, &
                Ip(:,:,n),Im(:,:,n), &
                ilo,ihi,dx,dt)
    enddo

    ! temperature-based PPM -- if desired, take the Ip(T)/Im(T)
    ! constructed above and use the EOS to overwrite Ip(p)/Im(p)
    if (ppm_temp_fix == 1) then
       do i = ilo-1, ihi+1
          do iwave = 1, 3
             eos_state%rho   = Ip(i,iwave,QRHO)
             eos_state%T     = Ip(i,iwave,QTEMP)
             eos_state%xn(:) = Ip(i,iwave,QFS:QFS-1+nspec)
             eos_state%aux   = Ip(i,iwave,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

             Ip(i,iwave,QPRES) = eos_state%p
             Ip(i,iwave,QREINT) = Ip(i,iwave,QRHO)*eos_state%e
             Ip_gc(i,iwave,1) = eos_state%gam1

             eos_state%rho   = Im(i,iwave,QRHO)
             eos_state%T     = Im(i,iwave,QTEMP)
             eos_state%xn(:) = Im(i,iwave,QFS:QFS-1+nspec)
             eos_state%aux   = Im(i,iwave,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

             Im(i,iwave,QPRES) = eos_state%p
             Im(i,iwave,QREINT) = Im(i,iwave,QRHO)*eos_state%e
             Im_gc(i,iwave,1) = eos_state%gam1
          enddo
       enddo

    endif


    ! get an edge-based gam1 here if we didn't get it from the EOS
    ! call above (for ppm_temp_fix = 1)
    if (ppm_temp_fix /= 1) then
       call ppm(gamc(:),qd_l1,qd_h1, &
                 q(:,QU),c, &
                 flatn, &
                 Ip_gc(:,:,1),Im_gc(:,:,1), &
                 ilo,ihi,dx,dt)
    endif

    if (ppm_trace_sources == 1) then
       do n=1,QVAR
          call ppm(srcQ(:,n),src_l1,src_h1, &
                   q(:,QU),c, &
                   flatn, &
                   Ip_src(:,:,n),Im_src(:,:,n), &
                   ilo,ihi,dx,dt)
       enddo
    endif

    !-------------------------------------------------------------------------
    ! x-direction
    !-------------------------------------------------------------------------

    ! Trace to left and right edges using upwind PPM
    do i = ilo-1, ihi+1

       cc = c(i)
       csq = cc**2

       rho = q(i,QRHO)
       u = q(i,QU)

       p = q(i,QPRES)
       rhoe_g = q(i,QREINT)
       h_g = ( (p + rhoe_g)/rho )/csq

       Clag = rho*cc

       gam_g = gamc(i)
       game = q(i,QGAME)

       !----------------------------------------------------------------------
       ! plus state on face i
       !----------------------------------------------------------------------

       ! set the reference state
       ! this will be the fastest moving state to the left --
       ! this is the method that Miller & Colella and Colella &
       ! Woodward use
       rho_ref  = Im(i,1,QRHO)
       u_ref    = Im(i,1,QU)

       p_ref    = Im(i,1,QPRES)
       rhoe_g_ref = Im(i,1,QREINT)

       tau_ref = ONE/Im(i,1,QRHO)

       gam_g_ref = Im_gc(i,1,1)
       game_ref = Im(i,1,QGAME)

       rho_ref = max(rho_ref,small_dens)
       p_ref = max(p_ref,small_pres)

       ! for tracing (optionally)
       cc_ref = sqrt(gam_g_ref*p_ref/rho_ref)
       csq_ref = cc_ref**2
       Clag_ref = rho_ref*cc_ref
       h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref )/csq_ref

       ! *m are the jumps carried by u-c
       ! *p are the jumps carried by u+c

       dum    = u_ref    - Im(i,1,QU)
       dptotm    = p_ref    - Im(i,1,QPRES)

       drho  = rho_ref  - Im(i,2,QRHO)
       drhoe_g = rhoe_g_ref - Im(i,2,QREINT)
       dptot    = p_ref    - Im(i,2,QPRES)
       dtau = tau_ref - ONE/Im(i,2,QRHO)

       dup    = u_ref    - Im(i,3,QU)
       dptotp    = p_ref    - Im(i,3,QPRES)

       ! if we are doing source term tracing, then we add the force to
       ! the velocity here, otherwise we will deal with this later
       if (ppm_trace_sources == 1) then
          dum = dum - hdt*Im_src(i,1,QU)
          dup = dup - hdt*Im_src(i,3,QU)
       endif

       ! optionally use the reference state in evaluating the
       ! eigenvectors
       if (ppm_reference_eigenvectors == 0) then
          rho_ev  = rho
          cc_ev   = cc
          csq_ev  = csq
          Clag_ev = Clag
          h_g_ev = h_g
          p_ev    = p
          tau_ev  = ONE/rho
       else
          rho_ev  = rho_ref
          cc_ev   = cc_ref
          csq_ev  = csq_ref
          Clag_ev = Clag_ref
          h_g_ev = h_g_ref
          p_ev    = p_ref
          tau_ev  = tau_ref
       endif


       if (ppm_predict_gammae == 0) then

          ! (rho, u, p, (rho e)) eigensystem

          ! these are analogous to the beta's from the original
          ! PPM paper (except we work with rho instead of tau).
          ! This is simply (l . dq), where dq = qref - I(q)

          alpham = HALF*(dptotm/(rho_ev*cc_ev) - dum)*rho_ev/cc_ev
          alphap = HALF*(dptotp/(rho_ev*cc_ev) + dup)*rho_ev/cc_ev
          alpha0r = drho - dptot/csq_ev
          alpha0e_g = drhoe_g - dptot*h_g_ev  ! note h_g has a 1/c**2 in it

       else

          ! (tau, u, p, game) eigensystem

          ! this is the way things were done in the original PPM
          ! paper -- here we work with tau in the characteristic
          ! system.

          alpham = HALF*( dum - dptotm/Clag_ev)/Clag_ev
          alphap = HALF*(-dup - dptotp/Clag_ev)/Clag_ev
          alpha0r = dtau + dptot/Clag_ev**2

          dge = game_ref - Im(i,2,QGAME)
          gfactor = (game - ONE)*(game - gam_g)  ! why not ref?
          alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

       endif ! which tracing method

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

       ! the final interface states are just
       ! q_s = q_ref - sum (l . dq) r
       if (i >= ilo) then

          if (ppm_predict_gammae == 0) then
             qxp(i,QRHO)   = rho_ref  + alphap + alpham + alpha0r
             qxp(i,QU)     = u_ref + (alphap - alpham)*cc_ev/rho_ev
             qxp(i,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
             qxp(i,QPRES)  = p_ref + (alphap + alpham)*csq_ev
          else
             tau_s = tau_ref + alphap + alpham + alpha0r
             qxp(i,QRHO)   = ONE/tau_s

             qxp(i,QU)     = u_ref + (alpham - alphap)*Clag_ev
             qxp(i,QPRES)  = p_ref - (alphap + alpham)*Clag_ev**2

             qxp(i,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
             qxp(i,QREINT) = qxp(i,QPRES )/(qxp(i,QGAME) - ONE)
          endif


          ! enforce small_*
          qxp(i,QRHO) = max(small_dens,qxp(i,QRHO))
          qxp(i,QPRES) = max(small_pres,qxp(i,QPRES))

          ! add source terms
          qxp(i  ,QRHO  ) = qxp(i,QRHO  ) + hdt*srcQ(i,QRHO)
          qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
          qxp(i  ,QREINT) = qxp(i,QREINT) + hdt*srcQ(i,QREINT)
          qxp(i  ,QPRES ) = qxp(i,QPRES ) + hdt*srcQ(i,QPRES)

          ! if we didn't add the gravity source via tracing above, then add it now
          if (ppm_trace_sources == 0) then
             qxp(i  ,QU) = qxp(i,QU) + hdt*srcQ(i,QU)
          endif

       endif


       !----------------------------------------------------------------------
       ! minus state on face i+1
       !----------------------------------------------------------------------

       ! set the reference state
       ! this will be the fastest moving state to the right
       rho_ref  = Ip(i,3,QRHO)
       u_ref    = Ip(i,3,QU)

       p_ref    = Ip(i,3,QPRES)
       rhoe_g_ref = Ip(i,3,QREINT)

       tau_ref = ONE/Ip(i,3,QRHO)

       gam_g_ref  = Ip_gc(i,3,1)
       game_ref = Ip(i,3,QGAME)

       rho_ref = max(rho_ref,small_dens)
       p_ref = max(p_ref,small_pres)

       ! for tracing (optionally)
       cc_ref = sqrt(gam_g_ref*p_ref/rho_ref)
       csq_ref = cc_ref**2
       Clag_ref = rho_ref*cc_ref
       h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref )/csq_ref

       ! *m are the jumps carried by u-c
       ! *p are the jumps carried by u+c

       dum    = u_ref    - Ip(i,1,QU)
       dptotm    = p_ref    - Ip(i,1,QPRES)

       drho  = rho_ref  - Ip(i,2,QRHO)
       dptot    = p_ref    - Ip(i,2,QPRES)
       drhoe_g = rhoe_g_ref - Ip(i,2,QREINT)
       dtau = tau_ref - ONE/Ip(i,2,QRHO)

       dup    = u_ref    - Ip(i,3,QU)
       dptotp    = p_ref    - Ip(i,3,QPRES)

       ! if we are doing source term tracing, then we add the force to
       ! the velocity here, otherwise we will deal with this in the
       ! trans_X routines
       if (ppm_trace_sources == 1) then
          dum = dum - hdt*Ip_src(i,1,QU)
          dup = dup - hdt*Ip_src(i,3,QU)
       endif

       ! optionally use the reference state in evaluating the
       ! eigenvectors
       if (ppm_reference_eigenvectors == 0) then
          rho_ev  = rho
          cc_ev   = cc
          csq_ev  = csq
          Clag_ev = Clag
          h_g_ev = h_g
          p_ev    = p
          tau_ev = ONE/rho
       else
          rho_ev  = rho_ref
          cc_ev   = cc_ref
          csq_ev  = csq_ref
          Clag_ev = Clag_ref
          h_g_ev = h_g_ref
          p_ev    = p_ref
          tau_ev  = tau_ref
       endif

       if (ppm_predict_gammae == 0) then

          ! (rho, u, p, (rho e)) eigensystem

          ! these are analogous to the beta's from the original
          ! PPM paper (except we work with rho instead of tau).
          ! This is simply (l . dq), where dq = qref - I(q)
          alpham = HALF*(dptotm/(rho_ev*cc_ev) - dum)*rho_ev/cc_ev
          alphap = HALF*(dptotp/(rho_ev*cc_ev) + dup)*rho_ev/cc_ev
          alpha0r = drho - dptot/csq_ev
          alpha0e_g = drhoe_g - dptot*h_g_ev  ! h_g has a 1/c**2 in it

       else

          ! (tau, u, p, game) eigensystem

          ! this is the way things were done in the original PPM
          ! paper -- here we work with tau in the characteristic
          ! system.

          alpham = HALF*( dum - dptotm/Clag_ev)/Clag_ev
          alphap = HALF*(-dup - dptotp/Clag_ev)/Clag_ev
          alpha0r = dtau + dptot/Clag_ev**2

          dge = game_ref - Ip(i,2,QGAME)
          gfactor = (game - ONE)*(game - gam_g)
          alpha0e_g = gfactor*dptot/(tau_ev*Clag_ev**2) + dge

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

       ! the final interface states are just
       ! q_s = q_ref - sum (l .dq) r
       if (i <= ihi) then

          if (ppm_predict_gammae == 0) then
             qxm(i+1,QRHO)   = rho_ref + alphap + alpham + alpha0r
             qxm(i+1,QU)     = u_ref + (alphap - alpham)*cc_ev/rho_ev
             qxm(i+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
             qxm(i+1,QPRES)  = p_ref + (alphap + alpham)*csq_ev

          else
             tau_s = tau_ref + (alphap + alpham + alpha0r)
             qxm(i+1,QRHO)   = ONE/tau_s

             qxm(i+1,QU)     = u_ref + (alpham - alphap)*Clag_ev
             qxm(i+1,QPRES)  = p_ref - (alphap + alpham)*Clag_ev**2

             qxm(i+1,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e_g
             qxm(i+1,QREINT) = qxm(i+1,QPRES )/(qxm(i+1,QGAME) - ONE)
          endif

          ! enforce small_*
          qxm(i+1,QRHO) = max(qxm(i+1,QRHO),small_dens)
          qxm(i+1,QPRES) = max(qxm(i+1,QPRES),small_pres)

          ! add source terms
          qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + hdt*srcQ(i,QRHO)
          qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
          qxm(i+1,QREINT) = qxm(i+1,QREINT) + hdt*srcQ(i,QREINT)
          qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + hdt*srcQ(i,QPRES)

          ! if we didn't add the gravity source via tracing above, then add it now
          if (ppm_trace_sources == 0) then
             qxm(i+1,QU) = qxm(i+1,QU) + hdt*srcQ(i,QU)
          endif

       end if


       !----------------------------------------------------------------------
       ! geometry source terms
       !----------------------------------------------------------------------

       if (dloga(i) /= 0) then
          courn = dtdx*(cc+abs(u))
          eta = (ONE-courn)/(cc*dt*abs(dloga(i)))
          dlogatmp = min(eta,ONE)*dloga(i)
          sourcr = -HALF*dt*rho*dlogatmp*u
          sourcp = sourcr*csq
          source = sourcp*h_g

          if (i <= ihi) then
             qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + sourcr
             qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
             qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + sourcp
             qxm(i+1,QREINT) = qxm(i+1,QREINT) + source
          end if
          if (i >= ilo) then
             qxp(i  ,QRHO  ) = qxp(i  ,QRHO  ) + sourcr
             qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
             qxp(i  ,QPRES ) = qxp(i  ,QPRES ) + sourcp
             qxp(i  ,QREINT) = qxp(i  ,QREINT) + source
          end if
       endif

    end do

    ! Enforce constant mass flux rate if specified
    if (fix_mass_flux_lo) then
       qxm(ilo,QRHO  ) = q(domlo(1)-1,QRHO)
       qxm(ilo,QU    ) = q(domlo(1)-1,QU  )
       qxm(ilo,QPRES ) = q(domlo(1)-1,QPRES)
       qxm(ilo,QREINT) = q(domlo(1)-1,QREINT)
    end if

    ! Enforce constant mass flux rate if specified
    if (fix_mass_flux_hi) then
       qxp(ihi+1,QRHO  ) = q(domhi(1)+1,QRHO)
       qxp(ihi+1,QU    ) = q(domhi(1)+1,QU  )
       qxp(ihi+1,QPRES ) = q(domhi(1)+1,QPRES)
       qxp(ihi+1,QREINT) = q(domhi(1)+1,QREINT)
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
             qxp(i,n) = q(i,n)
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
    deallocate(Ip_gc,Im_gc)
    if (ppm_trace_sources == 1) then
       deallocate(Ip_src,Im_src)
    endif

  end subroutine trace_ppm

end module trace_ppm_module
