module trace_ppm_module

  implicit none

  private

  public tracexy_ppm, tracez_ppm

contains

  subroutine tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         Ip,Im,Ip_g,Im_g,Ip_gc,Im_gc, &
                         qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         gamc,gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3, &
                         ilo1,ilo2,ihi1,ihi2,dt,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QPRES, QGAME, &
         small_dens, small_pres, &
         ppm_type, ppm_reference, ppm_trace_grav, &
         ppm_tau_in_tracing, ppm_reference_eigenvectors, &
         ppm_reference_edge_limit, ppm_flatten_before_integrals, &
         ppm_predict_gammae, &
         npassive, qpass_map
    use bl_constants_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer kc,k3d

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    double precision   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    double precision   Ip_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)
    double precision   Im_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)

    double precision   Ip_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,1)
    double precision   Im_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,1)

    double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)

    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2,gc_l3:gc_h3)

    double precision dt

    ! Local variables
    integer i, j
    integer n
    integer ipassive

    double precision cc, csq, Clag, rho, u, v, w, p, rhoe

    double precision drho, du, dv, dw, dp, drhoe, de, dge, dtau
    double precision dup, dvp, dpp
    double precision dum, dvm, dpm

    double precision :: rho_ref, u_ref, v_ref, p_ref, rhoe_ref, tau_ref
    double precision :: tau_s, e_s

    double precision :: cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref, game_ref, gfactor
    double precision :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev, tau_ev
    double precision :: gam, game

    double precision enth, alpham, alphap, alpha0r, alpha0e
    double precision apright, amright, azrright, azeright
    double precision apleft, amleft, azrleft, azeleft

    double precision xi, xi1
    double precision halfdt

    integer, parameter :: igx = 1
    integer, parameter :: igy = 2
    integer, parameter :: igz = 3

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call bl_error("Error:: trace_ppm_3d.f90 :: tracexy_ppm")
    end if

    halfdt = HALF * dt


    !==========================================================================
    ! PPM CODE
    !==========================================================================

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
    !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth) &
    !$OMP PRIVATE(rho_ref, u_ref, p_ref, rhoe_ref) &
    !$OMP PRIVATE(drho,dv,dw,dp,drhoe,de,dge,dum,dpm,dup,dpp,alpham,alphap,alpha0r) &
    !$OMP PRIVATE(alpha0e,amright,apright,azrright,azeright) &
    !$OMP PRIVATE(amleft,apleft,azrleft,azeleft,xi,xi1) &
    !$OMP PRIVATE(cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref, game_ref,gfactor) &
    !$OMP PRIVATE(cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev, tau_ev) &
    !$OMP PRIVATE(gam, game, tau_ref, dtau, tau_s, e_s)

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho = q(i,j,k3d,QRHO)

          cc = c(i,j,k3d)
          csq = cc**2
          Clag = rho*cc

          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)

          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)
          enth = ( (rhoe+p)/rho )/csq          

          game = q(i,j,k3d,QGAME)

          gam = gamc(i,j,k3d)


          !--------------------------------------------------------------------
          ! plus state on face i
          !--------------------------------------------------------------------

          if (i .ge. ilo1) then

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. u - cc >= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                u_ref    = u

                p_ref    = p
                rhoe_ref = rhoe

                tau_ref  = ONE/rho

                gam_ref  = gam
                
                game_ref = game

             else
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(i,j,kc,1,1,QRHO)
                u_ref    = Im(i,j,kc,1,1,QU)

                p_ref    = Im(i,j,kc,1,1,QPRES)
                rhoe_ref = Im(i,j,kc,1,1,QREINT)

                tau_ref  = ONE/Im(i,j,kc,1,1,QRHO)

                gam_ref  = Im_gc(i,j,kc,1,1,1)

                game_ref = Im(i,j,kc,1,1,QGAME)
             endif

             ! For tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)
   
             dum   = u_ref    - Im(i,j,kc,1,1,QU)
             dpm   = p_ref    - Im(i,j,kc,1,1,QPRES)
   
             drho  = rho_ref  - Im(i,j,kc,1,2,QRHO)
             dp    = p_ref    - Im(i,j,kc,1,2,QPRES)
             drhoe = rhoe_ref - Im(i,j,kc,1,2,QREINT)
             dtau  = tau_ref  - ONE/Im(i,j,kc,1,2,QRHO)

             dup   = u_ref    - Im(i,j,kc,1,3,QU)
             dpp   = p_ref    - Im(i,j,kc,1,3,QPRES)


             ! If we are doing gravity tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_grav .eq. 1) then
                dum = dum - halfdt*Im_g(i,j,kc,1,1,igx)
                dup = dup - halfdt*Im_g(i,j,kc,1,3,igx)
             endif
   
             ! Optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
                tau_ev  = 1.0d0/rho
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
                tau_ev  = tau_ref
             endif

             if (ppm_tau_in_tracing == 0) then

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is 
                ! simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dpm/(rho_ev*cc_ev) - dum)*rho_ev/cc_ev
                alphap = HALF*(dpp/(rho_ev*cc_ev) + dup)*rho_ev/cc_ev
                alpha0r = drho - dp/csq_ev
                alpha0e = drhoe - dp*enth_ev  ! note enth has a 1/c**2 in it

             else

                ! (tau, u, p, e) eigensystem
                ! or
                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_ref/rho_ref - Im(i,j,kc,1,2,QREINT)/Im(i,j,kc,1,2,QRHO))
                dge   = game_ref - Im(i,j,kc,1,2,QGAME)

                alpham = HALF*( dum - dpm/Clag_ev)/Clag_ev
                alphap = HALF*(-dup - dpp/Clag_ev)/Clag_ev
                alpha0r = dtau + dp/Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   alpha0e = de - dp*p_ev/Clag_ev**2
                else
                   gfactor = (game + 1.0d0)*(game - gam)
                   alpha0e = gfactor*dp/(tau_ev*Clag_ev**2) + dge
                endif

             endif    ! which tracing method 
                
             if (u-cc .gt. ZERO) then
                amright = ZERO
             else if (u-cc .lt. ZERO) then
                amright = -alpham
             else
                amright = -HALF*alpham
             endif

             if (u+cc .gt. ZERO) then
                apright = ZERO
             else if (u+cc .lt. ZERO) then
                apright = -alphap
             else
                apright = -HALF*alphap
             endif

             if (u .gt. ZERO) then
                azrright = ZERO
                azeright = ZERO
             else if (u .lt. ZERO) then
                azrright = -alpha0r
                azeright = -alpha0e
             else
                azrright = -HALF*alpha0r
                azeright = -HALF*alpha0e
             endif
                
             ! The final interface states are just
             ! q_s = q_ref - sum(l . dq) r
             ! note that the a{mpz}right as defined above have the minus already
             if (ppm_tau_in_tracing == 0) then
                qxp(i,j,kc,QRHO  ) =  rho_ref +  apright + amright + azrright
                qxp(i,j,kc,QU    ) =    u_ref + (apright - amright)*cc_ev/rho_ev
                qxp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enth_ev*csq_ev + azeright
                qxp(i,j,kc,QPRES ) =    p_ref + (apright + amright)*csq_ev
             else
                tau_s = tau_ref + apright + amright + azrright
                qxp(i,j,kc,QRHO  ) = ONE/tau_s
                qxp(i,j,kc,QU    ) = u_ref + (amright - apright)*Clag_ev

                qxp(i,j,kc,QPRES ) = p_ref + (-apright - amright)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_ref/rho_ref + (azeright - p_ev*amright - p_ev*apright)
                   qxp(i,j,kc,QREINT) = e_s/tau_s
                else
                   qxp(i,j,kc,QGAME) = game_ref + gfactor*(amright + apright)/tau_ev + azeright 
                   qxp(i,j,kc,QREINT) = qxp(i,j,kc,QPRES )/(qxp(i,j,kc,QGAME) - 1.0d0)
                endif
             endif    

             ! Enforce small_*
             qxp(i,j,kc,QRHO ) = max(qxp(i,j,kc,QRHO ),small_dens)
             qxp(i,j,kc,QPRES) = max(qxp(i,j,kc,QPRES),small_pres)

             ! Transverse velocities -- there's no projection here, so
             ! we don't need a reference state.  We only care about
             ! the state traced under the middle wave
             dv = Im(i,j,kc,1,2,QV)
             dw = Im(i,j,kc,1,2,QW)

             if (ppm_trace_grav .eq. 1) then
                dv  = dv  + halfdt*Im_g(i,j,kc,1,2,igy)
                dw  = dw  + halfdt*Im_g(i,j,kc,1,2,igz)
             endif

             ! Recall that I already takes the limit of the parabola
             ! in the event that the wave is not moving toward the
             ! interface
             if (u > ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qxp(i,j,kc,QV    ) = Im(i,j,kc,1,2,QV)
                   qxp(i,j,kc,QW    ) = Im(i,j,kc,1,2,QW)
                else
                   qxp(i,j,kc,QV    ) = v
                   qxp(i,j,kc,QW    ) = w
                endif
             else ! wave moving toward the interface 
                qxp(i,j,kc,QV    ) = dv
                qxp(i,j,kc,QW    ) = dw
             endif

             ! We may have already dealt with the flattening in the construction
             ! of the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE-flatn(i,j,k3d)
                xi = flatn(i,j,k3d)

                qxp(i,j,kc,QRHO  ) = xi1*rho  + xi*qxp(i,j,kc,QRHO  )
                qxp(i,j,kc,QU    ) = xi1*u    + xi*qxp(i,j,kc,QU    )
                qxp(i,j,kc,QV    ) = xi1*v    + xi*qxp(i,j,kc,QV    )
                qxp(i,j,kc,QW    ) = xi1*w    + xi*qxp(i,j,kc,QW    )
                qxp(i,j,kc,QREINT) = xi1*rhoe + xi*qxp(i,j,kc,QREINT)
                qxp(i,j,kc,QPRES ) = xi1*p    + xi*qxp(i,j,kc,QPRES )
             endif

          end if

          !--------------------------------------------------------------------
          ! minus state on face i + 1
          !--------------------------------------------------------------------
          if (i .le. ihi1) then

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. u + cc <= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc values
                rho_ref  = rho
                u_ref    = u

                p_ref    = p
                rhoe_ref = rhoe

                tau_ref  = ONE/rho

                gam_ref  = gam

                game_ref = game

             else
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,kc,1,3,QRHO)
                u_ref    = Ip(i,j,kc,1,3,QU)

                p_ref    = Ip(i,j,kc,1,3,QPRES)
                rhoe_ref = Ip(i,j,kc,1,3,QREINT)

                tau_ref  = ONE/Ip(i,j,kc,1,3,QRHO)

                gam_ref  = Ip_gc(i,j,kc,1,3,1)

                game_ref = Ip(i,j,kc,1,3,QGAME)
             endif

             ! For tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c
   
             dum   = u_ref    - Ip(i,j,kc,1,1,QU)
             dpm   = p_ref    - Ip(i,j,kc,1,1,QPRES)
   
             drho  = rho_ref  - Ip(i,j,kc,1,2,QRHO)
             dp    = p_ref    - Ip(i,j,kc,1,2,QPRES)
             drhoe = rhoe_ref - Ip(i,j,kc,1,2,QREINT)
             dtau  = tau_ref  - ONE/Ip(i,j,kc,1,2,QRHO)

             dup   = u_ref    - Ip(i,j,kc,1,3,QU)
             dpp   = p_ref    - Ip(i,j,kc,1,3,QPRES)

             ! If we are doing gravity tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_grav .eq. 1) then
                dum = dum - halfdt*Ip_g(i,j,kc,1,1,igx)
                dup = dup - halfdt*Ip_g(i,j,kc,1,3,igx)
             endif

             ! Optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
                tau_ev  = 1.0d0/rho
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
                tau_ev  = tau_ref
             endif

             if (ppm_tau_in_tracing == 0) then

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is 
                ! simply (l . dq), where dq = qref - I(q)
                
                alpham = HALF*(dpm/(rho_ev*cc_ev) - dum)*rho_ev/cc_ev
                alphap = HALF*(dpp/(rho_ev*cc_ev) + dup)*rho_ev/cc_ev
                alpha0r = drho - dp/csq_ev
                alpha0e = drhoe - dp*enth_ev  ! enth has a 1/c**2 in it

             else
                ! (tau, u, p, e) eigensystem
                ! or
                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_ref/rho_ref - Ip(i,j,kc,1,2,QREINT)/Ip(i,j,kc,1,2,QRHO))
                dge = game_ref - Ip(i,j,kc,1,2,QGAME)

                alpham = HALF*( dum - dpm/Clag_ev)/Clag_ev
                alphap = HALF*(-dup - dpp/Clag_ev)/Clag_ev
                alpha0r = dtau + dp/Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   alpha0e = de - dp*p_ev/Clag_ev**2
                else
                   gfactor = (game + 1.0d0)*(game - gam)
                   alpha0e = gfactor*dp/(tau_ev*Clag_ev**2) + dge
                endif

             end if
                
             if (u-cc .gt. ZERO) then
                amleft = -alpham
             else if (u-cc .lt. ZERO) then
                amleft = ZERO
             else
                amleft = -HALF*alpham
             endif
                
             if (u+cc .gt. ZERO) then
                apleft = -alphap
             else if (u+cc .lt. ZERO) then
                apleft = ZERO
             else
                apleft = -HALF*alphap
             endif
                
             if (u .gt. ZERO) then
                azrleft = -alpha0r
                azeleft = -alpha0e
             else if (u .lt. ZERO) then
                azrleft = ZERO
                azeleft = ZERO
             else
                azrleft = -HALF*alpha0r
                azeleft = -HALF*alpha0e
             endif

                
             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}left as defined above have the minus already
             if (ppm_tau_in_tracing == 0) then
                qxm(i+1,j,kc,QRHO  ) =  rho_ref +  apleft + amleft + azrleft
                qxm(i+1,j,kc,QU    ) =    u_ref + (apleft - amleft)*cc_ev/rho_ev
                qxm(i+1,j,kc,QREINT) = rhoe_ref + (apleft + amleft)*enth_ev*csq_ev + azeleft
                qxm(i+1,j,kc,QPRES ) =    p_ref + (apleft + amleft)*csq_ev
             else
                tau_s = tau_ref + (apleft + amleft + azrleft)
                qxm(i+1,j,kc,QRHO  ) = ONE/tau_s
                qxm(i+1,j,kc,QU    ) = u_ref + (amleft - apleft)*Clag_ev

                qxm(i+1,j,kc,QPRES ) = p_ref + (-apleft - amleft)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_ref/rho_ref + (azeleft - p_ev*amleft - p_ev*apleft)
                   qxm(i+1,j,kc,QREINT) = e_s/tau_s
                else
                   qxm(i+1,j,kc,QGAME) = game_ref + gfactor*(amleft + apleft)/tau_ev + azeleft
                   qxm(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QPRES )/(qxm(i+1,j,kc,QGAME) - 1.0d0)
                endif
             endif

             ! Enforce small_*
             qxm(i+1,j,kc,QRHO  ) = max(qxm(i+1,j,kc,QRHO ),small_dens)
             qxm(i+1,j,kc,QPRES)  = max(qxm(i+1,j,kc,QPRES),small_pres)

             ! transverse velocities
             dv    = Ip(i,j,kc,1,2,QV)
             dw    = Ip(i,j,kc,1,2,QW)

             if (ppm_trace_grav .eq. 1) then
                dv  = dv  + halfdt*Ip_g(i,j,kc,1,2,igy)
                dw  = dw  + halfdt*Ip_g(i,j,kc,1,2,igz)
             endif

             if (u < ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qxm(i+1,j,kc,QV    ) = Ip(i,j,kc,1,2,QV)
                   qxm(i+1,j,kc,QW    ) = Ip(i,j,kc,1,2,QW)
                else
                   qxm(i+1,j,kc,QV    ) = v
                   qxm(i+1,j,kc,QW    ) = w
                endif
             else ! wave moving toward interface
                qxm(i+1,j,kc,QV    ) = dv
                qxm(i+1,j,kc,QW    ) = dw
             endif

             ! We may have already dealt with flattening in the parabolas
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE - flatn(i,j,k3d)
                xi = flatn(i,j,k3d)

                qxm(i+1,j,kc,QRHO  ) = xi1*rho  + xi*qxm(i+1,j,kc,QRHO  )
                qxm(i+1,j,kc,QU    ) = xi1*u    + xi*qxm(i+1,j,kc,QU    )
                qxm(i+1,j,kc,QV    ) = xi1*v    + xi*qxm(i+1,j,kc,QV    )
                qxm(i+1,j,kc,QW    ) = xi1*w    + xi*qxm(i+1,j,kc,QW    )
                qxm(i+1,j,kc,QREINT) = xi1*rhoe + xi*qxm(i+1,j,kc,QREINT)
                qxm(i+1,j,kc,QPRES ) = xi1*p    + xi*qxm(i+1,j,kc,QPRES )
             endif

          end if

       end do
    end do
    !$OMP END PARALLEL DO

    
    !--------------------------------------------------------------------------
    ! passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    !$OMP parallel do private(ipassive,i,j,u,n,xi) IF(npassive .gt. 1)
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1

          ! Plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,QU)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = ONE
             endif

             ! the flattening here is a little confusing.  If 
             ! ppm_flatten_before_integrals = 0, then we are blending
             ! the cell centered state and the edge state here through
             ! the flattening procedure.  Otherwise, we've already
             ! took care of flattening.  What we want to do is: 
             !                                                                  
             ! q_l*  (1-xi)*q_i + xi*q_l                                        
             !                                                                  
             ! where                                                            
             !                                                                  
             ! q_l = q_ref - Proj{(q_ref - I)}                                  
             !                                                                  
             ! and Proj{} represents the characteristic projection.             
             ! But for these, there is only 1-wave that matters, the u          
             ! wave, so no projection is needed.  Since we are not              
             ! projecting, the reference state doesn't matter, so we            
             ! take it to be q_i, therefore, we reduce to                       
             !                                                                  
             ! q_l* = (1-xi)*q_i + xi*[q_i - (q_i - I)]                         
             !      = q_i + xi*(I - q_i)       

             if (u .gt. ZERO) then
                qxp(i,j,kc,n) = q(i,j,k3d,n)
             else if (u .lt. ZERO) then
                qxp(i,j,kc,n) = q(i,j,k3d,n) + xi*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else
                qxp(i,j,kc,n) = q(i,j,k3d,n) + HALF*xi*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo

          ! Minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = ONE
             endif

             if (u .gt. ZERO) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) + xi*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else if (u .lt. ZERO) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n)
             else
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) + HALF*xi*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo
       enddo
    enddo
    !$OMP end parallel do


    !--------------------------------------------------------------------------
    ! y-direction
    !--------------------------------------------------------------------------

    ! Trace to bottom and top edges using upwind PPM
    !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth) &
    !$OMP PRIVATE(rho_ref, v_ref, p_ref, rhoe_ref) &
    !$OMP PRIVATE(drho,du,dw,dp,drhoe,de,dge,dvm,dpm,dvp,dpp,alpham,alphap,alpha0r) &
    !$OMP PRIVATE(alpha0e,amright,apright,azrright,azeright,amleft) &
    !$OMP PRIVATE(apleft,azrleft,azeleft,xi,xi1) &
    !$OMP PRIVATE(cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref, game_ref,gfactor) &
    !$OMP PRIVATE(cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev, tau_ev) &
    !$OMP PRIVATE(gam, game, dtau, tau_ref, tau_s, e_s)

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho = q(i,j,k3d,QRHO)

          cc = c(i,j,k3d)
          csq = cc**2
          Clag = rho*cc

          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)

          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)
          enth = ( (rhoe+p)/rho )/csq

          gam = gamc(i,j,k3d)

          game = q(i,j,k3d,QGAME)

          !--------------------------------------------------------------------
          ! plus state on face j
          !--------------------------------------------------------------------

          if (j .ge. ilo2) then

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. v - cc >= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                v_ref    = v

                p_ref    = p
                rhoe_ref = rhoe

                tau_ref  = ONE/rho

                gam_ref  = gam

                game_ref = game

             else
                ! This will be the fastest moving state to the left
                rho_ref  = Im(i,j,kc,2,1,QRHO)
                v_ref    = Im(i,j,kc,2,1,QV)

                p_ref    = Im(i,j,kc,2,1,QPRES)
                rhoe_ref = Im(i,j,kc,2,1,QREINT)

                tau_ref  = ONE/Im(i,j,kc,2,1,QRHO)
                gam_ref  = Im_gc(i,j,kc,2,1,1)

                game_ref = Im(i,j,kc,2,1,QGAME)
             endif

             ! For tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             dvm   = v_ref    - Im(i,j,kc,2,1,QV)
             dpm   = p_ref    - Im(i,j,kc,2,1,QPRES)
   
             drho  = rho_ref  - Im(i,j,kc,2,2,QRHO)
             dp    = p_ref    - Im(i,j,kc,2,2,QPRES)
             drhoe = rhoe_ref - Im(i,j,kc,2,2,QREINT)
             dtau  = tau_ref  - ONE/Im(i,j,kc,2,2,QRHO)

             dvp   = v_ref    - Im(i,j,kc,2,3,QV)
             dpp   = p_ref    - Im(i,j,kc,2,3,QPRES)

             ! If we are doing gravity tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_grav .eq. 1) then
                dvm = dvm - halfdt*Im_g(i,j,kc,2,1,igy)
                dvp = dvp - halfdt*Im_g(i,j,kc,2,3,igy)
             endif

             ! Optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
                tau_ev  = 1.0d0/rho
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
                tau_ev  = tau_ref
             endif

             if (ppm_tau_in_tracing == 0) then

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This
                ! is simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dpm/(rho_ev*cc_ev) - dvm)*rho_ev/cc_ev
                alphap = HALF*(dpp/(rho_ev*cc_ev) + dvp)*rho_ev/cc_ev
                alpha0r = drho - dp/csq_ev
                alpha0e = drhoe - dp*enth_ev

             else

                ! (tau, u, p, e) eigensystem
                ! or
                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_ref/rho_ref - Im(i,j,kc,2,2,QREINT)/Im(i,j,kc,2,2,QRHO))
                dge = game_ref - Im(i,j,kc,2,2,QGAME)

                alpham = HALF*( dvm - dpm/Clag_ev)/Clag_ev
                alphap = HALF*(-dvp - dpp/Clag_ev)/Clag_ev
                alpha0r = dtau + dp/Clag_ev**2
                
                if (ppm_predict_gammae == 0) then
                   alpha0e = de - dp*p_ev/Clag_ev**2
                else
                   gfactor = (game + 1.0d0)*(game - gam)
                   alpha0e = gfactor*dp/(tau_ev*Clag_ev**2) + dge
                endif

             end if
   
             if (v-cc .gt. ZERO) then
                amright = ZERO
             else if (v-cc .lt. ZERO) then
                amright = -alpham
             else
                amright = -HALF*alpham
             endif
             
             if (v+cc .gt. ZERO) then
                apright = ZERO
             else if (v+cc .lt. ZERO) then
                apright = -alphap
             else
                apright = -HALF*alphap
             endif
             
             if (v .gt. ZERO) then
                azrright = ZERO
                azeright = ZERO
             else if (v .lt. ZERO) then
                azrright = -alpha0r
                azeright = -alpha0e
             else
                azrright = -HALF*alpha0r
                azeright = -HALF*alpha0e
             endif
                
             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}right as defined above have the minus already
             if (ppm_tau_in_tracing == 0) then
                qyp(i,j,kc,QRHO  ) = rho_ref + apright + amright + azrright
                qyp(i,j,kc,QV    ) = v_ref + (apright - amright)*cc_ev/rho_ev
                qyp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enth_ev*csq_ev + azeright
                qyp(i,j,kc,QPRES ) = p_ref + (apright + amright)*csq_ev
             else
                tau_s = tau_ref + apright + amright + azrright                
                qyp(i,j,kc,QRHO  ) = ONE/tau_s
                qyp(i,j,kc,QV    ) = v_ref + (amright - apright)*Clag_ev

                qyp(i,j,kc,QPRES ) = p_ref + (-apright - amright)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_ref/rho_ref + (azeright - p_ev*amright - p_ev*apright)
                   qyp(i,j,kc,QREINT) = e_s/tau_s
                else
                   qyp(i,j,kc,QGAME) = game_ref + gfactor*(amright + apright)/tau_ev + azeright
                   qyp(i,j,kc,QREINT) = qyp(i,j,kc,QPRES )/(qyp(i,j,kc,QGAME) - 1.0d0)
                endif
             endif

             ! Enforce small_*
             qyp(i,j,kc,QRHO ) = max(qyp(i,j,kc,QRHO ),small_dens)
             qyp(i,j,kc,QPRES) = max(qyp(i,j,kc,QPRES),small_pres)

             ! transverse velocities
             du    = Im(i,j,kc,2,2,QU)
             dw    = Im(i,j,kc,2,2,QW)

             if (ppm_trace_grav .eq. 1) then
                du  = du  + halfdt*Im_g(i,j,kc,2,2,igx)
                dw  = dw  + halfdt*Im_g(i,j,kc,2,2,igz)
             endif

             if (v > ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qyp(i,j,kc,QU    ) = Im(i,j,kc,2,2,QU)
                   qyp(i,j,kc,QW    ) = Im(i,j,kc,2,2,QW)
                else
                   qyp(i,j,kc,QU    ) = u
                   qyp(i,j,kc,QW    ) = w
                endif
             else ! wave moving toward the interface
                qyp(i,j,kc,QU    ) = du
                qyp(i,j,kc,QW    ) = dw
             endif

             ! We may have already dealt with flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE - flatn(i,j,k3d)
                xi = flatn(i,j,k3d)

                qyp(i,j,kc,QRHO  ) = xi1*rho  + xi*qyp(i,j,kc,QRHO  )
                qyp(i,j,kc,QV    ) = xi1*v    + xi*qyp(i,j,kc,QV    )
                qyp(i,j,kc,QU    ) = xi1*u    + xi*qyp(i,j,kc,QU    )
                qyp(i,j,kc,QW    ) = xi1*w    + xi*qyp(i,j,kc,QW    )
                qyp(i,j,kc,QREINT) = xi1*rhoe + xi*qyp(i,j,kc,QREINT)
                qyp(i,j,kc,QPRES ) = xi1*p    + xi*qyp(i,j,kc,QPRES )
             endif

          end if

          !--------------------------------------------------------------------
          ! minus state on face j+1
          !--------------------------------------------------------------------

          if (j .le. ihi2) then

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. v + cc <= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                v_ref    = v

                p_ref    = p
                rhoe_ref = rhoe

                tau_ref  = ONE/rho

                gam_ref  = gam

                game_ref = game

             else
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,kc,2,3,QRHO)
                v_ref    = Ip(i,j,kc,2,3,QV)

                p_ref    = Ip(i,j,kc,2,3,QPRES)
                rhoe_ref = Ip(i,j,kc,2,3,QREINT)

                tau_ref  = ONE/Ip(i,j,kc,2,3,QRHO)

                gam_ref  = Ip_gc(i,j,kc,2,3,1)

                game_ref = Ip(i,j,kc,2,3,QGAME)
             endif

             ! For tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             dvm   = v_ref    - Ip(i,j,kc,2,1,QV)
             dpm   = p_ref    - Ip(i,j,kc,2,1,QPRES)
   
             drho  = rho_ref  - Ip(i,j,kc,2,2,QRHO)
             dp    = p_ref    - Ip(i,j,kc,2,2,QPRES)
             drhoe = rhoe_ref - Ip(i,j,kc,2,2,QREINT)
             dtau  = tau_ref  - ONE/Ip(i,j,kc,2,2,QRHO)

             dvp   = v_ref    - Ip(i,j,kc,2,3,QV)
             dpp   = p_ref    - Ip(i,j,kc,2,3,QPRES)

             ! If we are doing gravity tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_grav .eq. 1) then
                dvm = dvm - halfdt*Ip_g(i,j,kc,2,1,igy)
                dvp = dvp - halfdt*Ip_g(i,j,kc,2,3,igy)
             endif

             ! Optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
                tau_ev  = 1.0d0/rho
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
                tau_ev  = tau_ref
             endif

             if (ppm_tau_in_tracing == 0) then

                ! These are analogous to the beta's from the original PPM
                ! paper.  This is simply (l . dq), where dq = qref - I(q)
 
                alpham = HALF*(dpm/(rho_ev*cc_ev) - dvm)*rho_ev/cc_ev
                alphap = HALF*(dpp/(rho_ev*cc_ev) + dvp)*rho_ev/cc_ev
                alpha0r = drho - dp/csq_ev
                alpha0e = drhoe - dp*enth_ev
                
             else
                ! (tau, u, p, e) eigensystem
                ! or
                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_ref/rho_ref - Ip(i,j,kc,2,2,QREINT)/Ip(i,j,kc,2,2,QRHO))
                dge = game_ref - Ip(i,j,kc,2,2,QGAME)

                alpham = HALF*( dvm - dpm/Clag_ev)/Clag_ev
                alphap = HALF*(-dvp - dpp/Clag_ev)/Clag_ev
                alpha0r = dtau + dp/Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   alpha0e = de - dp*p_ev/Clag_ev**2
                else
                   gfactor = (game + 1.0d0)*(game - gam)
                   alpha0e = gfactor*dp/(tau_ev*Clag_ev**2) + dge
                endif
                
             end if

             if (v-cc .gt. ZERO) then
                amleft = -alpham
             else if (v-cc .lt. ZERO) then
                amleft = ZERO
             else
                amleft = -HALF*alpham
             endif
             
             if (v+cc .gt. ZERO) then
                apleft = -alphap
             else if (v+cc .lt. ZERO) then
                apleft = ZERO
             else
                apleft = -HALF*alphap
             endif
             
             if (v .gt. ZERO) then
                azrleft = -alpha0r
                azeleft = -alpha0e
             else if (v .lt. ZERO) then
                azrleft = ZERO
                azeleft = ZERO
             else
                azrleft = -HALF*alpha0r
                azeleft = -HALF*alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}left as defined above has the minus already
             if (ppm_tau_in_tracing == 0) then
                qym(i,j+1,kc,QRHO  ) = rho_ref + apleft + amleft + azrleft
                qym(i,j+1,kc,QV    ) = v_ref + (apleft - amleft)*cc_ev/rho_ev
                qym(i,j+1,kc,QREINT) = rhoe_ref + (apleft + amleft)*enth_ev*csq_ev + azeleft
                qym(i,j+1,kc,QPRES ) = p_ref + (apleft + amleft)*csq_ev
             else
                tau_s = tau_ref + apleft + amleft + azrleft
                qym(i,j+1,kc,QRHO  ) = ONE/tau_s
                qym(i,j+1,kc,QV    ) = v_ref + (amleft - apleft)*Clag_ev

                qym(i,j+1,kc,QPRES ) = p_ref + (-apleft - amleft)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_ref/rho_ref + (azeleft - p_ev*amleft - p_ev*apleft)
                   qym(i,j+1,kc,QREINT) = e_s/tau_s
                else
                   qym(i,j+1,kc,QGAME) = game_ref + gfactor*(amleft + apleft)/tau_ev + azeleft 
                   qym(i,j+1,kc,QREINT) = qym(i,j+1,kc,QPRES )/(qym(i,j+1,kc,QGAME) - 1.0d0)
                endif

             endif

             ! Enforce small_*
             qym(i,j+1,kc,QRHO ) = max(qym(i,j+1,kc,QRHO ),small_dens)
             qym(i,j+1,kc,QPRES) = max(qym(i,j+1,kc,QPRES),small_pres)
             
             ! transverse velocities
             du    = Ip(i,j,kc,2,2,QU)
             dw    = Ip(i,j,kc,2,2,QW)

             if (ppm_trace_grav .eq. 1) then
                du  = du  + halfdt*Ip_g(i,j,kc,2,2,igx)
                dw  = dw  + halfdt*Ip_g(i,j,kc,2,2,igz)
             endif

             if (v < ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qym(i,j+1,kc,QU    ) = Ip(i,j,kc,2,2,QU)
                   qym(i,j+1,kc,QW    ) = Ip(i,j,kc,2,2,QW)
                else
                   qym(i,j+1,kc,QU    ) = u
                   qym(i,j+1,kc,QW    ) = w
                endif
             else ! wave is moving toward the interface
                qym(i,j+1,kc,QU    ) = du
                qym(i,j+1,kc,QW    ) = dw
             endif

             ! We may have already dealt with flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE - flatn(i,j,k3d)
                xi = flatn(i,j,k3d)

                qym(i,j+1,kc,QRHO  ) = xi1*rho  + xi*qym(i,j+1,kc,QRHO  )
                qym(i,j+1,kc,QV    ) = xi1*v    + xi*qym(i,j+1,kc,QV    )
                qym(i,j+1,kc,QU    ) = xi1*u    + xi*qym(i,j+1,kc,QU    )
                qym(i,j+1,kc,QW    ) = xi1*w    + xi*qym(i,j+1,kc,QW    )
                qym(i,j+1,kc,QREINT) = xi1*rhoe + xi*qym(i,j+1,kc,QREINT)
                qym(i,j+1,kc,QPRES ) = xi1*p    + xi*qym(i,j+1,kc,QPRES )
             endif

          end if

       end do
    end do
    !$OMP END PARALLEL DO

    !--------------------------------------------------------------------------
    ! Passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    !$OMP parallel do private(n,i,j,v,ipassive,xi) IF(npassive .gt. 1)
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do i = ilo1-1, ihi1+1

          ! Plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,k3d,QV)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = ONE
             endif

             if (v .gt. ZERO) then
                qyp(i,j,kc,n) = q(i,j,k3d,n)
             else if (v .lt. ZERO) then
                qyp(i,j,kc,n) = q(i,j,k3d,n) + xi*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else
                qyp(i,j,kc,n) = q(i,j,k3d,n) + HALF*xi*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo
          
          ! Minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,k3d,QV)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = ONE
             endif

             if (v .gt. ZERO) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n) + xi*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else if (v .lt. ZERO) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n)
             else
                qym(i,j+1,kc,n) = q(i,j,k3d,n) + HALF*xi*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo
          
       enddo
    enddo

  end subroutine tracexy_ppm

  subroutine tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        Ip,Im,Ip_g,Im_g,Ip_gc,Im_gc, &
                        qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        gamc,gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3, &
                        ilo1,ilo2,ihi1,ihi2,dt,km,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QPRES, QGAME, &
         small_dens, small_pres, &
         ppm_type, ppm_reference, ppm_trace_grav, &
         ppm_tau_in_tracing, ppm_reference_eigenvectors, &
         ppm_reference_edge_limit, ppm_flatten_before_integrals, &
         ppm_predict_gammae, &
         npassive, qpass_map
    use bl_constants_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer km,kc,k3d

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    double precision   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    double precision   Ip_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)
    double precision   Im_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,3)

    double precision   Ip_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,1)
    double precision   Im_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,1)

    double precision qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)

    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2,gc_l3:gc_h3)

    double precision dt

    !     Local variables
    integer i, j
    integer n
    integer ipassive

    double precision cc, csq, Clag, rho, u, v, w, p, rhoe

    double precision drho, du, dv, dp, drhoe, de, dge, dtau
    double precision dwp, dpp
    double precision dwm, dpm

    double precision :: rho_ref, w_ref, p_ref, rhoe_ref, tau_ref
    double precision :: tau_s, e_s

    double precision :: cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref, game_ref, gfactor
    double precision :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev, tau_ev
    double precision :: gam, game

    double precision enth, alpham, alphap, alpha0r, alpha0e
    double precision apright, amright, azrright, azeright
    double precision apleft, amleft, azrleft, azeleft

    double precision xi, xi1
    double precision halfdt

    integer, parameter :: igx = 1
    integer, parameter :: igy = 2
    integer, parameter :: igz = 3

    halfdt = HALF * dt

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call bl_error("Error:: trace_ppm_3d.f90 :: tracez_ppm")
    end if

    !==========================================================================
    ! PPM CODE
    !==========================================================================

    ! Trace to left and right edges using upwind PPM
    !
    ! Note: in contrast to the above code for x and y, here the loop
    ! is over interfaces, not over cell-centers.


    !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth) &
    !$OMP PRIVATE(rho_ref,w_ref,p_ref,rhoe_ref) &
    !$OMP PRIVATE(drho,du,dv,dp,drhoe,de,dge,dwm,dpm,dwp,dpp,alpham,alphap,alpha0r,alpha0e) &
    !$OMP PRIVATE(amright,apright,azrright,azeright,amleft,apleft)&
    !$OMP PRIVATE(azrleft,azeleft,xi,xi1) &
    !$OMP PRIVATE(cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref, game_ref,gfactor) &
    !$OMP PRIVATE(cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev, tau_ev) &
    !$OMP PRIVATE(gam, game, dtau, tau_ref, tau_s, e_s)

    !--------------------------------------------------------------------------
    ! construct qzp  -- plus state on face kc
    !--------------------------------------------------------------------------
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho  = q(i,j,k3d,QRHO)

          cc   = c(i,j,k3d)
          csq  = cc**2
          Clag = rho*cc

          u    = q(i,j,k3d,QU)
          v    = q(i,j,k3d,QV)
          w    = q(i,j,k3d,QW)

          p    = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)
          enth = ( (rhoe+p)/rho )/csq

          gam = gamc(i,j,k3d)

          game = q(i,j,k3d,QGAME)

          ! Set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. w - cc >= ZERO .and. &
                ppm_reference_edge_limit == 0) ) then
             ! original Castro way -- cc value
             rho_ref  = rho
             w_ref    = w

             p_ref    = p
             rhoe_ref = rhoe

             tau_ref  = ONE/rho

             gam_ref  = gam

             game_ref = game

          else
             ! This will be the fastest moving state to the left
             rho_ref  = Im(i,j,kc,3,1,QRHO)
             w_ref    = Im(i,j,kc,3,1,QW)

             p_ref    = Im(i,j,kc,3,1,QPRES)
             rhoe_ref = Im(i,j,kc,3,1,QREINT)

             tau_ref  = ONE/Im(i,j,kc,3,1,QRHO)
             gam_ref  = Im_gc(i,j,kc,3,1,1)

             game_ref = Im(i,j,kc,3,1,QGAME)
          endif

          ! For tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   = w_ref    - Im(i,j,kc,3,1,QW)
          dpm   = p_ref    - Im(i,j,kc,3,1,QPRES)

          drho  = rho_ref  - Im(i,j,kc,3,2,QRHO)
          dp    = p_ref    - Im(i,j,kc,3,2,QPRES)
          drhoe = rhoe_ref - Im(i,j,kc,3,2,QREINT)
          dtau  = tau_ref  - ONE/Im(i,j,kc,3,2,QRHO)

          dwp   = w_ref    - Im(i,j,kc,3,3,QW)
          dpp   = p_ref    - Im(i,j,kc,3,3,QPRES)

          ! If we are doing gravity tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_grav .eq. 1) then
             dwm = dwm - halfdt*Im_g(i,j,kc,3,1,igz)
             dwp = dwp - halfdt*Im_g(i,j,kc,3,3,igz)
          endif

          ! Optionally use the reference state in evaluating the
          ! eigenvectors
          if (ppm_reference_eigenvectors == 0) then
             rho_ev  = rho
             cc_ev   = cc
             csq_ev  = csq
             Clag_ev = Clag
             enth_ev = enth
             p_ev    = p
             tau_ev  = 1.0d0/rho
          else
             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             enth_ev = enth_ref
             p_ev    = p_ref
             tau_ev  = tau_ref
          endif
          
          if (ppm_tau_in_tracing == 0) then

             ! These are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)
             alpham = HALF*(dpm/(rho_ev*cc_ev) - dwm)*rho_ev/cc_ev
             alphap = HALF*(dpp/(rho_ev*cc_ev) + dwp)*rho_ev/cc_ev
             alpha0r = drho - dp/csq_ev
             alpha0e = drhoe - dp*enth_ev

          else 
             ! (tau, u, p, e) eigensystem
             ! or
             ! (tau, u, p, game) eigensystem

             ! This is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system
             de = (rhoe_ref/rho_ref - Im(i,j,kc,3,2,QREINT)/Im(i,j,kc,3,2,QRHO))
             dge = game_ref - Im(i,j,kc,3,2,QGAME)

             alpham = HALF*( dwm - dpm/Clag_ev)/Clag_ev
             alphap = HALF*(-dwp - dpp/Clag_ev)/Clag_ev
             alpha0r = dtau + dp/Clag_ev**2

             if (ppm_predict_gammae == 0) then
                alpha0e = de - dp*p_ev/Clag_ev**2
             else
                gfactor = (game + 1.0d0)*(game - gam)
                alpha0e = gfactor*dp/(tau_ev*Clag_ev**2) + dge
             endif

          endif

          if (w-cc .gt. ZERO) then
             amright = ZERO
          else if (w-cc .lt. ZERO) then
             amright = -alpham
          else
             amright = -HALF*alpham
          endif
          if (w+cc .gt. ZERO) then
             apright = ZERO
          else if (w+cc .lt. ZERO) then
             apright = -alphap
          else
             apright = -HALF*alphap
          endif
          if (w .gt. ZERO) then
             azrright = ZERO
             azeright = ZERO
          else if (w .lt. ZERO) then
             azrright = -alpha0r
             azeright = -alpha0e
          else
             azrright = -HALF*alpha0r
             azeright = -HALF*alpha0e
          endif

          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          ! note that the a{mpz}right as defined above have the minus already
          if (ppm_tau_in_tracing == 0) then
             qzp(i,j,kc,QRHO  ) = rho_ref + apright + amright + azrright
             qzp(i,j,kc,QW    ) = w_ref + (apright - amright)*cc_ev/rho_ev
             qzp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enth_ev*csq_ev + azeright
             qzp(i,j,kc,QPRES ) = p_ref + (apright + amright)*csq_ev
          else
             tau_s = tau_ref + apright + amright + azrright
             qzp(i,j,kc,QRHO  ) = ONE/tau_s
             qzp(i,j,kc,QW    ) = w_ref + (amright - apright)*Clag_ev

             qzp(i,j,kc,QPRES ) = p_ref + (-apright - amright)*Clag_ev**2

             if (ppm_predict_gammae == 0) then
                e_s = rhoe_ref/rho_ref + (azeright - p_ev*amright - p_ev*apright)           
                qzp(i,j,kc,QREINT) = e_s/tau_s
             else
                qzp(i,j,kc,QGAME) = game_ref + gfactor*(amright + apright)/tau_ev + azeright
                qzp(i,j,kc,QREINT) = qzp(i,j,kc,QPRES )/(qzp(i,j,kc,QGAME) - 1.0d0)
             endif

          endif

          ! Enforce small_*
          qzp(i,j,kc,QRHO ) = max(qzp(i,j,kc,QRHO ),small_dens)
          qzp(i,j,kc,QPRES) = max(qzp(i,j,kc,QPRES),small_pres)

          ! transverse velocities
          du    = Im(i,j,kc,3,2,QU)
          dv    = Im(i,j,kc,3,2,QV)

          if (ppm_trace_grav .eq. 1) then
             du  = du  + halfdt*Im_g(i,j,kc,3,2,igx)
             dv  = dv  + halfdt*Im_g(i,j,kc,3,2,igy)
          endif

          if (w > ZERO) then 
             if (ppm_reference_edge_limit == 1) then
                qzp(i,j,kc,QU    ) = Im(i,j,kc,3,2,QU)
                qzp(i,j,kc,QV    ) = Im(i,j,kc,3,2,QV)
             else
                qzp(i,j,kc,QU    ) = u
                qzp(i,j,kc,QV    ) = v
             endif
          else ! wave moving toward the interface
             qzp(i,j,kc,QU    ) = du
             qzp(i,j,kc,QV    ) = dv
          endif

          ! We may have already dealt with flattening in the parabola
          if (ppm_flatten_before_integrals == 0) then
             xi1 = ONE - flatn(i,j,k3d)
             xi = flatn(i,j,k3d)

             qzp(i,j,kc,QRHO  ) = xi1*rho  + xi*qzp(i,j,kc,QRHO  )
             qzp(i,j,kc,QW    ) = xi1*w    + xi*qzp(i,j,kc,QW    )
             qzp(i,j,kc,QU    ) = xi1*u    + xi*qzp(i,j,kc,QU    )
             qzp(i,j,kc,QV    ) = xi1*v    + xi*qzp(i,j,kc,QV    )
             qzp(i,j,kc,QREINT) = xi1*rhoe + xi*qzp(i,j,kc,QREINT)
             qzp(i,j,kc,QPRES ) = xi1*p    + xi*qzp(i,j,kc,QPRES )
          endif

          !--------------------------------------------------------------------
          ! This is all for qzm -- minus state on face kc
          !--------------------------------------------------------------------

          ! Note this is different from how we do 1D, 2D, and the
          ! x and y-faces in 3D, where the analogous thing would have
          ! been to find the minus state on face kc+1

          rho  = q(i,j,k3d-1,QRHO)

          cc   = c(i,j,k3d-1)
          csq  = cc**2
          Clag = rho*cc

          u    = q(i,j,k3d-1,QU)
          v    = q(i,j,k3d-1,QV)
          w    = q(i,j,k3d-1,QW)

          p    = q(i,j,k3d-1,QPRES)
          rhoe = q(i,j,k3d-1,QREINT)
          enth = ( (rhoe+p)/rho )/csq

          gam = gamc(i,j,k3d-1)

          game = q(i,j,k3d-1,QGAME)

          ! Set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. w + cc <= ZERO .and. &
                ppm_reference_edge_limit == 0) ) then
             rho_ref  = rho
             w_ref    = w

             p_ref    = p
             rhoe_ref = rhoe

             tau_ref  = ONE/rho

             gam_ref  = gam

             game_ref = game

          else
             ! This will be the fastest moving state to the right
             rho_ref  = Ip(i,j,km,3,3,QRHO)
             w_ref    = Ip(i,j,km,3,3,QW)

             p_ref    = Ip(i,j,km,3,3,QPRES)
             rhoe_ref = Ip(i,j,km,3,3,QREINT)

             tau_ref  = ONE/Ip(i,j,km,3,3,QRHO)

             gam_ref  = Ip_gc(i,j,km,3,3,1)

             game_ref = Ip(i,j,km,3,3,QGAME)
          endif

          ! For tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   = w_ref    - Ip(i,j,km,3,1,QW)
          dpm   = p_ref    - Ip(i,j,km,3,1,QPRES)

          drho  = rho_ref  - Ip(i,j,km,3,2,QRHO)
          dp    = p_ref    - Ip(i,j,km,3,2,QPRES)
          drhoe = rhoe_ref - Ip(i,j,km,3,2,QREINT)
          dtau  = tau_ref  - ONE/Ip(i,j,km,3,2,QRHO)

          dwp   = w_ref    - Ip(i,j,km,3,3,QW)
          dpp   = p_ref    - Ip(i,j,km,3,3,QPRES)

          ! If we are doing gravity tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_grav .eq. 1) then
             dwm = dwm - halfdt*Ip_g(i,j,km,3,1,igz)
             dwp = dwp - halfdt*Ip_g(i,j,km,3,3,igz)
          endif

          ! Optionally use the reference state in evaluating the
          ! eigenvectors
          if (ppm_reference_eigenvectors == 0) then
             rho_ev  = rho
             cc_ev   = cc
             csq_ev  = csq
             Clag_ev = Clag
             enth_ev = enth
             p_ev    = p
             tau_ev  = 1.0d0/rho
          else
             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             enth_ev = enth_ref
             p_ev    = p_ref
             tau_ev  = tau_ref
          endif
          
          if (ppm_tau_in_tracing == 0) then

             ! These are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)             
             alpham = HALF*(dpm/(rho_ev*cc_ev) - dwm)*rho_ev/cc_ev
             alphap = HALF*(dpp/(rho_ev*cc_ev) + dwp)*rho_ev/cc_ev
             alpha0r = drho - dp/csq_ev
             alpha0e = drhoe - dp*enth_ev

          else

             ! (tau, u, p, e) eigensystem
             ! or
             ! (tau, u, p, game) eigensystem

             ! This is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system
             de = (rhoe_ref/rho_ref - Ip(i,j,km,3,2,QREINT)/Ip(i,j,km,3,2,QRHO))
             dge = game_ref - Ip(i,j,km,3,2,QGAME)

             alpham = HALF*( dwm - dpm/Clag_ev)/Clag_ev
             alphap = HALF*(-dwp - dpp/Clag_ev)/Clag_ev
             alpha0r = dtau + dp/Clag_ev**2

             if (ppm_predict_gammae == 0) then
                alpha0e = de - dp*p_ev/Clag_ev**2
             else
                gfactor = (game + 1.0d0)*(game - gam)
                alpha0e = gfactor*dp/(tau_ev*Clag_ev**2) + dge
             endif

          endif 
             
          if (w-cc .gt. ZERO) then
             amleft = -alpham
          else if (w-cc .lt. ZERO) then
             amleft = ZERO
          else
             amleft = -HALF*alpham
          endif
          if (w+cc .gt. ZERO) then
             apleft = -alphap
          else if (w+cc .lt. ZERO) then
             apleft = ZERO
          else
             apleft = -HALF*alphap
          endif
          if (w .gt. ZERO) then
             azrleft = -alpha0r
             azeleft = -alpha0e
          else if (w .lt. ZERO) then
             azrleft = ZERO
             azeleft = ZERO
          else
             azrleft = -HALF*alpha0r
             azeleft = -HALF*alpha0e
          endif
             
          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          ! note that the a{mpz}left as defined above have the minus already
          if (ppm_tau_in_tracing == 0) then
             qzm(i,j,kc,QRHO  ) = rho_ref + apleft + amleft + azrleft
             qzm(i,j,kc,QW    ) = w_ref + (apleft - amleft)*cc_ev/rho_ev
             qzm(i,j,kc,QREINT) = rhoe_ref + (apleft + amleft)*enth_ev*csq_ev + azeleft
             qzm(i,j,kc,QPRES ) = p_ref + (apleft + amleft)*csq_ev
          else
             tau_s = tau_ref + apleft + amleft + azrleft
             qzm(i,j,kc,QRHO  ) = ONE/tau_s
             qzm(i,j,kc,QW    ) = w_ref + (amleft - apleft)*Clag_ev

             qzm(i,j,kc,QPRES ) = p_ref + (-apleft - amleft)*Clag_ev**2

             if (ppm_predict_gammae == 0) then
                e_s = rhoe_ref/rho_ref + (azeleft - p_ev*amleft - p_ev*apleft)
                qzm(i,j,kc,QREINT) = e_s/tau_s
             else
                qzm(i,j,kc,QGAME) = game_ref + gfactor*(amleft + apleft)/tau_ev + azeleft
                qzm(i,j,kc,QREINT) = qzm(i,j,kc,QPRES )/(qzm(i,j,kc,QGAME) - 1.0d0)
             endif

          endif

          ! Enforce small_*
          qzm(i,j,kc,QRHO ) = max(qzm(i,j,kc,QRHO ),small_dens)
          qzm(i,j,kc,QPRES) = max(qzm(i,j,kc,QPRES),small_pres)

          ! Transverse velocity
          du    = Ip(i,j,km,3,2,QU)
          dv    = Ip(i,j,km,3,2,QV)

          if (ppm_trace_grav .eq. 1) then
             du  = du  + halfdt*Ip_g(i,j,km,3,2,igx)
             dv  = dv  + halfdt*Ip_g(i,j,km,3,2,igy)
          endif

          if (w < ZERO) then
             if (ppm_reference_edge_limit == 1) then
                qzm(i,j,kc,QU    ) = Ip(i,j,km,3,2,QU)
                qzm(i,j,kc,QV    ) = Ip(i,j,km,3,2,QV)
             else
                qzm(i,j,kc,QU    ) = u
                qzm(i,j,kc,QV    ) = v
             endif
          else ! wave moving toward the interface
             qzm(i,j,kc,QU    ) = du
             qzm(i,j,kc,QV    ) = dv
          endif

          ! We may have already taken care of flattening in the parabola
          if (ppm_flatten_before_integrals == 0) then
             xi1 = ONE - flatn(i,j,k3d-1)
             xi = flatn(i,j,k3d-1)

             qzm(i,j,kc,QRHO  ) = xi1*rho  + xi*qzm(i,j,kc,QRHO  )
             qzm(i,j,kc,QW    ) = xi1*w    + xi*qzm(i,j,kc,QW    )
             qzm(i,j,kc,QU    ) = xi1*u    + xi*qzm(i,j,kc,QU    )
             qzm(i,j,kc,QV    ) = xi1*v    + xi*qzm(i,j,kc,QV    )
             qzm(i,j,kc,QREINT) = xi1*rhoe + xi*qzm(i,j,kc,QREINT)
             qzm(i,j,kc,QPRES ) = xi1*p    + xi*qzm(i,j,kc,QPRES )

          endif          

       end do
    end do
    !$OMP END PARALLEL DO

    !--------------------------------------------------------------------------
    ! passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    !$OMP parallel do private(n,w,i,j,ipassive,xi) IF(npassive .gt. 1)
    do ipassive = 1, npassive
         n = qpass_map(ipassive)
         do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1

                  ! Plus state on face kc
                  w = q(i,j,k3d,QW)

                  if (ppm_flatten_before_integrals == 0) then
                     xi = flatn(i,j,k3d)
                  else
                     xi = ONE
                  endif

                  if (w .gt. ZERO) then
                     qzp(i,j,kc,n) = q(i,j,k3d,n)
                  else if (w .lt. ZERO) then
                     qzp(i,j,kc,n) = q(i,j,k3d,n) + xi*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
                  else
                     qzp(i,j,kc,n) = q(i,j,k3d,n) + HALF*xi*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
                  endif

                  ! Minus state on face k
                  w = q(i,j,k3d-1,QW)
                  
                  if (ppm_flatten_before_integrals == 0) then
                     xi = flatn(i,j,k3d-1)
                  else
                     xi = ONE
                  endif

                  if (w .gt. ZERO) then
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n) + xi*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
                  else if (w .lt. ZERO) then
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n)
                  else
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n) + HALF*xi*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
                  endif

               enddo
         enddo
    enddo
    !$OMP end parallel do

  end subroutine tracez_ppm

end module trace_ppm_module
