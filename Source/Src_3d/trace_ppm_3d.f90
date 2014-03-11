module trace_ppm_module

  implicit none

  private

  public tracexy_ppm, tracez_ppm

contains

  subroutine tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         Ip,Im,Ip_g,Im_g,Ip_gc,Im_gc, &
                         qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                         gamc,gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3, &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QESGS, QPRES, QFA, QFS, nadv, &
         small_dens, small_pres, &
         ppm_type, ppm_reference, ppm_trace_grav, &
         ppm_tau_in_tracing, ppm_reference_eigenvectors, &
         ppm_reference_edge_limit, ppm_flatten_before_integrals

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
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

    double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2,gc_l3:gc_h3)

    double precision dx, dy, dt

    ! Local variables
    integer i, j
    integer n, iadv, ispec
    integer npassive,ipassive,qpass_map(QVAR)

    double precision cc, csq, Clag, rho, u, v, w, p, rhoe

    double precision drho, du, dv, dw, dp, drhoe, de, dtau
    double precision dup, dvp, dpp
    double precision dum, dvm, dpm

    double precision :: rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref, tau_ref
    double precision :: tau_s, e_s

    double precision :: cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref
    double precision :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev
    double precision :: gam

    double precision enth, alpham, alphap, alpha0r, alpha0e
    double precision alpha0u, alpha0v, alpha0w
    double precision apright, amright, azrright, azeright
    double precision azu1rght, azv1rght, azw1rght
    double precision apleft, amleft, azrleft, azeleft
    double precision azu1left, azv1left, azw1left

    double precision xi, xi1
    double precision halfdt

    integer, parameter :: igx = 1
    integer, parameter :: igy = 2
    integer, parameter :: igz = 3

    ! Group all the passively advected quantities together
    npassive = 0
    if (QESGS .gt. -1) then
       qpass_map(1) = QESGS
       npassive = 1
    endif
    do iadv = 1, nadv
       qpass_map(npassive + iadv) = QFA + iadv - 1
    enddo
    npassive = npassive + nadv
    if(QFS .gt. -1) then
       do ispec = 1, nspec+naux
          qpass_map(npassive + ispec) = QFS + ispec - 1
       enddo
       npassive = npassive + nspec + naux
    endif

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call bl_error("Error:: trace_ppm_3d.f90 :: tracexy_ppm")
    end if

    halfdt = 0.5d0 * dt


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
    !$OMP PRIVATE(rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref) &
    !$OMP PRIVATE(drho,dv,dw,dp,drhoe,de,dum,dpm,dup,dpp,alpham,alphap,alpha0r) &
    !$OMP PRIVATE(alpha0e,alpha0v,alpha0w,amright,apright,azrright,azeright,azv1rght,azw1rght) &
    !$OMP PRIVATE(amleft,apleft,azrleft,azeleft,azv1left,azw1left,xi,xi1) &
    !$OMP PRIVATE(cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref) &
    !$OMP PRIVATE(cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev) &
    !$OMP PRIVATE(gam, tau_ref, dtau, tau_s, e_s)

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


          !--------------------------------------------------------------------
          ! plus state on face i
          !--------------------------------------------------------------------

          if (i .ge. ilo1) then

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. u - cc >= 0.0d0 .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                u_ref    = u
                v_ref    = v
                w_ref    = w
                p_ref    = p
                rhoe_ref = rhoe
                tau_ref  = 1.d0/rho
                gam_ref  = gam
             else
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(i,j,kc,1,1,QRHO)
                u_ref    = Im(i,j,kc,1,1,QU)
                v_ref    = Im(i,j,kc,1,1,QV)
                w_ref    = Im(i,j,kc,1,1,QW)
                p_ref    = Im(i,j,kc,1,1,QPRES)
                rhoe_ref = Im(i,j,kc,1,1,QREINT)
                tau_ref  = 1.d0/Im(i,j,kc,1,1,QRHO)
                gam_ref  = Im_gc(i,j,kc,1,1,1)
             endif
   
             ! for tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)
   
             dum    = (u_ref    - Im(i,j,kc,1,1,QU))
             dpm    = (p_ref    - Im(i,j,kc,1,1,QPRES))
   
             drho  = (rho_ref  - Im(i,j,kc,1,2,QRHO))
             dv    = (v_ref    - Im(i,j,kc,1,2,QV))
             dw    = (w_ref    - Im(i,j,kc,1,2,QW))
             dp    = (p_ref    - Im(i,j,kc,1,2,QPRES))
             drhoe = (rhoe_ref - Im(i,j,kc,1,2,QREINT))
             dtau  = (tau_ref  - 1.d0/Im(i,j,kc,1,2,QRHO))

             dup    = (u_ref    - Im(i,j,kc,1,3,QU))
             dpp    = (p_ref    - Im(i,j,kc,1,3,QPRES))

             ! if we are doing gravity tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_grav .eq. 1) then
                dum = dum - halfdt*Im_g(i,j,kc,1,1,igx)
                dv  = dv  - halfdt*Im_g(i,j,kc,1,2,igy)
                dw  = dw  - halfdt*Im_g(i,j,kc,1,2,igz)
                dup = dup - halfdt*Im_g(i,j,kc,1,3,igx)
             endif
   
             ! optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
             endif


             if (ppm_tau_in_tracing == 0) then

                ! these are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is 
                ! simply (l . dq), where dq = qref - I(q)
   
                alpham = 0.5d0*(dpm/(rho_ev*cc_ev) - dum)*rho_ev/cc_ev
                alphap = 0.5d0*(dpp/(rho_ev*cc_ev) + dup)*rho_ev/cc_ev
                alpha0r = drho - dp/csq_ev
                alpha0e = drhoe - dp*enth_ev  ! note enth has a 1/c**2 in it
                alpha0v = dv
                alpha0w = dw
                
                if (u-cc .gt. 0.d0) then
                   amright = 0.d0
                else if (u-cc .lt. 0.d0) then
                   amright = -alpham
                else
                   amright = -0.5d0*alpham
                endif

                if (u+cc .gt. 0.d0) then
                   apright = 0.d0
                else if (u+cc .lt. 0.d0) then
                   apright = -alphap
                else
                   apright = -0.5d0*alphap
                endif

                if (u .gt. 0.d0) then
                   azrright = 0.d0
                   azeright = 0.d0
                   azv1rght = 0.d0
                   azw1rght = 0.d0
                else if (u .lt. 0.d0) then
                   azrright = -alpha0r
                   azeright = -alpha0e
                   azv1rght = -alpha0v
                   azw1rght = -alpha0w
                else
                   azrright = -0.5d0*alpha0r
                   azeright = -0.5d0*alpha0e
                   azv1rght = -0.5d0*alpha0v
                   azw1rght = -0.5d0*alpha0w
                endif
                
                ! the final interface states are just
                ! q_s = q_ref - sum(l . dq) r

                if (ppm_flatten_before_integrals == 0) then
                   xi1 = 1.0d0-flatn(i,j,k3d)
                   xi = flatn(i,j,k3d)
                else
                   xi1 = 0.0d0
                   xi = 1.0d0
                endif

                qxp(i,j,kc,QRHO  ) = xi1*rho  + xi*(rho_ref + apright + amright + azrright)
                qxp(i,j,kc,QU    ) = xi1*u    + xi*(u_ref + (apright - amright)*cc_ev/rho_ev)
                qxp(i,j,kc,QV    ) = xi1*v    + xi*(v_ref + azv1rght)
                qxp(i,j,kc,QW    ) = xi1*w    + xi*(w_ref + azw1rght)
                qxp(i,j,kc,QREINT) = xi1*rhoe + xi*(rhoe_ref + (apright + amright)*enth_ev*csq_ev + azeright)
                qxp(i,j,kc,QPRES ) = xi1*p    + xi*(p_ref + (apright + amright)*csq_ev)

                qxp(i,j,kc,QRHO ) = max(qxp(i,j,kc,QRHO ),small_dens)
                qxp(i,j,kc,QPRES) = max(qxp(i,j,kc,QPRES),small_pres)

             else
                ! (tau, u, p, e) eigensystem

                ! this is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_ref/rho_ref - Im(i,j,kc,1,2,QREINT)/Im(i,j,kc,1,2,QRHO))

                alpham = 0.5d0*( dum - dpm/Clag_ev)/Clag_ev
                alphap = 0.5d0*(-dup - dpp/Clag_ev)/Clag_ev
                alpha0r = dtau + dp/Clag_ev**2
                alpha0e = de - dp*p_ev/Clag_ev**2
                alpha0v = dv
                alpha0w = dw
                
                if (u-cc .gt. 0.d0) then
                   amright = 0.d0
                else if (u-cc .lt. 0.d0) then
                   amright = -alpham
                else
                   amright = -0.5d0*alpham
                endif

                if (u+cc .gt. 0.d0) then
                   apright = 0.d0
                else if (u+cc .lt. 0.d0) then
                   apright = -alphap
                else
                   apright = -0.5d0*alphap
                endif

                if (u .gt. 0.d0) then
                   azrright = 0.d0
                   azeright = 0.d0
                   azv1rght = 0.d0
                   azw1rght = 0.d0
                else if (u .lt. 0.d0) then
                   azrright = -alpha0r
                   azeright = -alpha0e
                   azv1rght = -alpha0v
                   azw1rght = -alpha0w
                else
                   azrright = -0.5d0*alpha0r
                   azeright = -0.5d0*alpha0e
                   azv1rght = -0.5d0*alpha0v
                   azw1rght = -0.5d0*alpha0w
                endif
                
                ! the final interface states are just
                ! q_s = q_ref - sum(l . dq) r

                if (ppm_flatten_before_integrals == 0) then
                   xi1 = 1.0d0-flatn(i,j,k3d)
                   xi = flatn(i,j,k3d)
                else
                   xi1 = 0.0d0
                   xi = 1.0d0
                endif

                tau_s = tau_ref + apright + amright + azrright
                qxp(i,j,kc,QRHO  ) = xi1*rho + xi/tau_s

                qxp(i,j,kc,QU    ) = xi1*u    + xi*(u_ref + (amright - apright)*Clag_ev)
                qxp(i,j,kc,QV    ) = xi1*v    + xi*(v_ref + azv1rght)
                qxp(i,j,kc,QW    ) = xi1*w    + xi*(w_ref + azw1rght)

                e_s = rhoe_ref/rho_ref + (azeright - p_ev*amright - p_ev*apright)
                qxp(i,j,kc,QREINT) = xi1*rhoe + xi*e_s/tau_s

                qxp(i,j,kc,QPRES ) = xi1*p    + xi*(p_ref + (-apright - amright)*Clag_ev**2)

                qxp(i,j,kc,QRHO ) = max(qxp(i,j,kc,QRHO ),small_dens)
                qxp(i,j,kc,QPRES) = max(qxp(i,j,kc,QPRES),small_pres)


             endif
                
          end if


          !--------------------------------------------------------------------
          ! minus state on face i + 1
          !--------------------------------------------------------------------
          if (i .le. ihi1) then

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. u + cc <= 0.0d0 .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc values
                rho_ref  = rho
                u_ref    = u
                v_ref    = v
                w_ref    = w
                p_ref    = p
                rhoe_ref = rhoe
                tau_ref  = 1.d0/rho
                gam_ref  = gam
             else
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,kc,1,3,QRHO)
                u_ref    = Ip(i,j,kc,1,3,QU)
                v_ref    = Ip(i,j,kc,1,3,QV)
                w_ref    = Ip(i,j,kc,1,3,QW)
                p_ref    = Ip(i,j,kc,1,3,QPRES)
                rhoe_ref = Ip(i,j,kc,1,3,QREINT)
                tau_ref  = 1.d0/Ip(i,j,kc,1,3,QRHO)
                gam_ref  = Ip_gc(i,j,kc,1,3,1)
             endif
   
             ! for tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)
   
             dum    = (u_ref    - Ip(i,j,kc,1,1,QU))
             dpm    = (p_ref    - Ip(i,j,kc,1,1,QPRES))
   
             drho  = (rho_ref  - Ip(i,j,kc,1,2,QRHO))
             dv    = (v_ref    - Ip(i,j,kc,1,2,QV))
             dw    = (w_ref    - Ip(i,j,kc,1,2,QW))
             dp    = (p_ref    - Ip(i,j,kc,1,2,QPRES))
             drhoe = (rhoe_ref - Ip(i,j,kc,1,2,QREINT))
             dtau  = (tau_ref  - 1.d0/Ip(i,j,kc,1,2,QRHO))

             dup    = (u_ref    - Ip(i,j,kc,1,3,QU))
             dpp    = (p_ref    - Ip(i,j,kc,1,3,QPRES))

             ! if we are doing gravity tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_grav .eq. 1) then
                dum = dum - halfdt*Ip_g(i,j,kc,1,1,igx)
                dv  = dv  - halfdt*Ip_g(i,j,kc,1,2,igy)
                dw  = dw  - halfdt*Ip_g(i,j,kc,1,2,igz)
                dup = dup - halfdt*Ip_g(i,j,kc,1,3,igx)
             endif

             ! optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
             endif

             if (ppm_tau_in_tracing == 0) then

                ! these are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is 
                ! simply (l . dq), where dq = qref - I(q)
                
                alpham = 0.5d0*(dpm/(rho_ev*cc_ev) - dum)*rho_ev/cc_ev
                alphap = 0.5d0*(dpp/(rho_ev*cc_ev) + dup)*rho_ev/cc_ev
                alpha0r = drho - dp/csq_ev
                alpha0e = drhoe - dp*enth_ev  ! enth has a 1/c**2 in it
                alpha0v = dv
                alpha0w = dw
                
                if (u-cc .gt. 0.d0) then
                   amleft = -alpham
                else if (u-cc .lt. 0.d0) then
                   amleft = 0.d0
                else
                   amleft = -0.5d0*alpham
                endif
                
                if (u+cc .gt. 0.d0) then
                   apleft = -alphap
                else if (u+cc .lt. 0.d0) then
                   apleft = 0.d0
                else
                   apleft = -0.5d0*alphap
                endif
                
                if (u .gt. 0.d0) then
                   azrleft = -alpha0r
                   azeleft = -alpha0e
                   azv1left = -alpha0v
                   azw1left = -alpha0w
                else if (u .lt. 0.d0) then
                   azrleft = 0.d0
                   azeleft = 0.d0
                   azv1left = 0.d0
                   azw1left = 0.d0
                else
                   azrleft = -0.5d0*alpha0r
                   azeleft = -0.5d0*alpha0e
                   azv1left = -0.5d0*alpha0v
                   azw1left = -0.5d0*alpha0w
                endif
                
                ! the final interface states are just
                ! q_s = q_ref - sum (l . dq) r

                if (ppm_flatten_before_integrals == 0) then
                   xi1 = 1.0d0 - flatn(i,j,k3d)
                   xi = flatn(i,j,k3d)
                else
                   xi1 = 0.0d0
                   xi = 1.0d0
                endif
                
                qxm(i+1,j,kc,QRHO  ) = xi1*rho  + xi*(rho_ref + apleft + amleft + azrleft)
                qxm(i+1,j,kc,QU    ) = xi1*u    + xi*(u_ref + (apleft - amleft)*cc_ev/rho_ev)
                qxm(i+1,j,kc,QV    ) = xi1*v    + xi*(v_ref + azv1left)
                qxm(i+1,j,kc,QW    ) = xi1*w    + xi*(w_ref + azw1left)
                qxm(i+1,j,kc,QREINT) = xi1*rhoe + xi*(rhoe_ref + (apleft + amleft)*enth_ev*csq_ev + azeleft)
                qxm(i+1,j,kc,QPRES ) = xi1*p    + xi*(p_ref + (apleft + amleft)*csq_ev)

                qxm(i+1,j,kc,QRHO  ) = max(qxm(i+1,j,kc,QRHO ),small_dens)
                qxm(i+1,j,kc,QPRES)  = max(qxm(i+1,j,kc,QPRES),small_pres)

             else
                ! (tau, u, p, e) eigensystem

                ! this is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_ref/rho_ref - Ip(i,j,kc,1,2,QREINT)/Ip(i,j,kc,1,2,QRHO))

                alpham = 0.5d0*( dum - dpm/Clag_ev)/Clag_ev
                alphap = 0.5d0*(-dup - dpp/Clag_ev)/Clag_ev
                alpha0r = dtau + dp/Clag_ev**2
                alpha0e = de - dp*p_ev/Clag_ev**2
                alpha0v = dv
                alpha0w = dw
                
                if (u-cc .gt. 0.d0) then
                   amleft = -alpham
                else if (u-cc .lt. 0.d0) then
                   amleft = 0.d0
                else
                   amleft = -0.5d0*alpham
                endif
                
                if (u+cc .gt. 0.d0) then
                   apleft = -alphap
                else if (u+cc .lt. 0.d0) then
                   apleft = 0.d0
                else
                   apleft = -0.5d0*alphap
                endif
                
                if (u .gt. 0.d0) then
                   azrleft = -alpha0r
                   azeleft = -alpha0e
                   azv1left = -alpha0v
                   azw1left = -alpha0w
                else if (u .lt. 0.d0) then
                   azrleft = 0.d0
                   azeleft = 0.d0
                   azv1left = 0.d0
                   azw1left = 0.d0
                else
                   azrleft = -0.5d0*alpha0r
                   azeleft = -0.5d0*alpha0e
                   azv1left = -0.5d0*alpha0v
                   azw1left = -0.5d0*alpha0w
                endif
                
                ! the final interface states are just
                ! q_s = q_ref - sum (l . dq) r

                if (ppm_flatten_before_integrals == 0) then
                   xi1 = 1.0d0 - flatn(i,j,k3d)
                   xi = flatn(i,j,k3d)
                else
                   xi1 = 0.0d0
                   xi = 1.0d0
                endif
                
                tau_s = tau_ref + (apleft + amleft + azrleft)
                qxm(i+1,j,kc,QRHO  ) = xi1*rho + xi/tau_s

                qxm(i+1,j,kc,QU    ) = xi1*u    + xi*(u_ref + (amleft - apleft)*Clag_ev)
                qxm(i+1,j,kc,QV    ) = xi1*v    + xi*(v_ref + azv1left)
                qxm(i+1,j,kc,QW    ) = xi1*w    + xi*(w_ref + azw1left)

                e_s = rhoe_ref/rho_ref + (azeleft - p_ev*amleft - p_ev*apleft)
                qxm(i+1,j,kc,QREINT) = xi1*rhoe + xi*e_s/tau_s

                qxm(i+1,j,kc,QPRES ) = xi1*p    + xi*(p_ref + (-apleft - amleft)*Clag_ev**2)

                qxm(i+1,j,kc,QRHO  ) = max(qxm(i+1,j,kc,QRHO ),small_dens)
                qxm(i+1,j,kc,QPRES)  = max(qxm(i+1,j,kc,QPRES),small_pres)


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
                xi = 1.0d0
             endif

             ! the flattening here is a little confusing.  What we              
             ! want to do is:                                                   
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

             if (u .gt. 0.d0) then
                qxp(i,j,kc,n) = q(i,j,k3d,n)
             else if (u .lt. 0.d0) then
                qxp(i,j,kc,n) = q(i,j,k3d,n) + xi*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else
                qxp(i,j,kc,n) = q(i,j,k3d,n) + 0.5d0*xi*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo

          ! Minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = 1.0d0
             endif

             if (u .gt. 0.d0) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) + xi*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else if (u .lt. 0.d0) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n)
             else
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) + 0.5d0*xi*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
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
    !$OMP PRIVATE(rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref) &
    !$OMP PRIVATE(drho,du,dw,dp,drhoe,de,dvm,dpm,dvp,dpp,alpham,alphap,alpha0r) &
    !$OMP PRIVATE(alpha0e,alpha0u,alpha0w,amright,apright,azrright,azeright,azu1rght,azw1rght,amleft) &
    !$OMP PRIVATE(apleft,azrleft,azeleft,azu1left,azw1left,xi,xi1) &
    !$OMP PRIVATE(cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref) &
    !$OMP PRIVATE(cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev) &
    !$OMP PRIVATE(gam, dtau, tau_ref, tau_s, e_s)

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


          !--------------------------------------------------------------------
          ! plus state on face j
          !--------------------------------------------------------------------

          if (j .ge. ilo2) then

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. v - cc >= 0.0d0 .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                u_ref    = u
                v_ref    = v
                w_ref    = w
                p_ref    = p
                rhoe_ref = rhoe
                tau_ref  = 1.d0/rho
                gam_ref  = gam
             else
                ! This will be the fastest moving state to the left
                rho_ref  = Im(i,j,kc,2,1,QRHO)
                u_ref    = Im(i,j,kc,2,1,QU)
                v_ref    = Im(i,j,kc,2,1,QV)
                w_ref    = Im(i,j,kc,2,1,QW)
                p_ref    = Im(i,j,kc,2,1,QPRES)
                rhoe_ref = Im(i,j,kc,2,1,QREINT)
                tau_ref  = 1.d0/Im(i,j,kc,2,1,QRHO)
                gam_ref  = Im_gc(i,j,kc,2,1,1)
             endif
   
             ! for tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm    = (v_ref    - Im(i,j,kc,2,1,QV))
             dpm    = (p_ref    - Im(i,j,kc,2,1,QPRES))
   
             drho  = (rho_ref  - Im(i,j,kc,2,2,QRHO))
             du    = (u_ref    - Im(i,j,kc,2,2,QU))
             dw    = (w_ref    - Im(i,j,kc,2,2,QW))
             dp    = (p_ref    - Im(i,j,kc,2,2,QPRES))
             drhoe = (rhoe_ref - Im(i,j,kc,2,2,QREINT))
             dtau  = (tau_ref  - 1.d0/Im(i,j,kc,2,2,QRHO))

             dvp    = (v_ref    - Im(i,j,kc,2,3,QV))
             dpp    = (p_ref    - Im(i,j,kc,2,3,QPRES))

             ! if we are doing gravity tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_grav .eq. 1) then
                dvm = dvm - halfdt*Im_g(i,j,kc,2,1,igy)
                du  = du  - halfdt*Im_g(i,j,kc,2,2,igx)
                dw  = dw  - halfdt*Im_g(i,j,kc,2,2,igz)
                dvp = dvp - halfdt*Im_g(i,j,kc,2,3,igy)
             endif

             ! optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
             endif


             if (ppm_tau_in_tracing == 0) then

                ! these are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This
                ! is simply (l . dq), where dq = qref - I(q)

                alpham = 0.5d0*(dpm/(rho_ev*cc_ev) - dvm)*rho_ev/cc_ev
                alphap = 0.5d0*(dpp/(rho_ev*cc_ev) + dvp)*rho_ev/cc_ev
                alpha0r = drho - dp/csq_ev
                alpha0e = drhoe - dp*enth_ev
                alpha0u = du
                alpha0w = dw
   
                if (v-cc .gt. 0.d0) then
                   amright = 0.d0
                else if (v-cc .lt. 0.d0) then
                   amright = -alpham
                else
                   amright = -0.5d0*alpham
                endif
                
                if (v+cc .gt. 0.d0) then
                   apright = 0.d0
                else if (v+cc .lt. 0.d0) then
                   apright = -alphap
                else
                   apright = -0.5d0*alphap
                endif
                
                if (v .gt. 0.d0) then
                   azrright = 0.d0
                   azeright = 0.d0
                   azu1rght = 0.d0
                   azw1rght = 0.d0
                else if (v .lt. 0.d0) then
                   azrright = -alpha0r
                   azeright = -alpha0e
                   azu1rght = -alpha0u
                   azw1rght = -alpha0w
                else
                   azrright = -0.5d0*alpha0r
                   azeright = -0.5d0*alpha0e
                   azu1rght = -0.5d0*alpha0u
                   azw1rght = -0.5d0*alpha0w
                endif
                
                ! the final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                if (ppm_flatten_before_integrals == 0) then
                   xi1 = 1.0d0 - flatn(i,j,k3d)
                   xi = flatn(i,j,k3d)
                else
                   xi1 = 0.0d0
                   xi = 1.0d0
                endif
                
                qyp(i,j,kc,QRHO  ) = xi1*rho  + xi*(rho_ref + apright + amright + azrright)
                qyp(i,j,kc,QV    ) = xi1*v    + xi*(v_ref + (apright - amright)*cc_ev/rho_ev)
                qyp(i,j,kc,QU    ) = xi1*u    + xi*(u_ref + azu1rght)
                qyp(i,j,kc,QW    ) = xi1*w    + xi*(w_ref + azw1rght)
                qyp(i,j,kc,QREINT) = xi1*rhoe + xi*(rhoe_ref + (apright + amright)*enth_ev*csq_ev + azeright)
                qyp(i,j,kc,QPRES ) = xi1*p    + xi*(p_ref + (apright + amright)*csq_ev)
                
                qyp(i,j,kc,QRHO ) = max(qyp(i,j,kc,QRHO ),small_dens)
                qyp(i,j,kc,QPRES) = max(qyp(i,j,kc,QPRES),small_pres)

             else
                ! (tau, u, p, e) eigensystem

                ! this is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_ref/rho_ref - Im(i,j,kc,2,2,QREINT)/Im(i,j,kc,2,2,QRHO))

                alpham = 0.5d0*( dvm - dpm/Clag_ev)/Clag_ev
                alphap = 0.5d0*(-dvp - dpp/Clag_ev)/Clag_ev
                alpha0r = dtau + dp/Clag_ev**2
                alpha0e = de - dp*p_ev/Clag_ev**2
                alpha0u = du
                alpha0w = dw
   
                if (v-cc .gt. 0.d0) then
                   amright = 0.d0
                else if (v-cc .lt. 0.d0) then
                   amright = -alpham
                else
                   amright = -0.5d0*alpham
                endif
                
                if (v+cc .gt. 0.d0) then
                   apright = 0.d0
                else if (v+cc .lt. 0.d0) then
                   apright = -alphap
                else
                   apright = -0.5d0*alphap
                endif
                
                if (v .gt. 0.d0) then
                   azrright = 0.d0
                   azeright = 0.d0
                   azu1rght = 0.d0
                   azw1rght = 0.d0
                else if (v .lt. 0.d0) then
                   azrright = -alpha0r
                   azeright = -alpha0e
                   azu1rght = -alpha0u
                   azw1rght = -alpha0w
                else
                   azrright = -0.5d0*alpha0r
                   azeright = -0.5d0*alpha0e
                   azu1rght = -0.5d0*alpha0u
                   azw1rght = -0.5d0*alpha0w
                endif
                
                ! the final interface states are just
                ! q_s = q_ref - sum (l . dq) r

                if (ppm_flatten_before_integrals == 0) then
                   xi1 = 1.0d0 - flatn(i,j,k3d)
                   xi = flatn(i,j,k3d)
                else
                   xi1 = 0.0d0
                   xi = 1.0d0
                endif
                
                tau_s = tau_ref + apright + amright + azrright                
                qyp(i,j,kc,QRHO  ) = xi1*rho  + xi/tau_s

                qyp(i,j,kc,QV    ) = xi1*v    + xi*(v_ref + (amright - apright)*Clag_ev)
                qyp(i,j,kc,QU    ) = xi1*u    + xi*(u_ref + azu1rght)
                qyp(i,j,kc,QW    ) = xi1*w    + xi*(w_ref + azw1rght)

                e_s = rhoe_ref/rho_ref + (azeright - p_ev*amright - p_ev*apright)
                qyp(i,j,kc,QREINT) = xi1*rhoe + xi*e_s/tau_s

                qyp(i,j,kc,QPRES ) = xi1*p    + xi*(p_ref + (-apright - amright)*Clag_ev**2)
                
                qyp(i,j,kc,QRHO ) = max(qyp(i,j,kc,QRHO ),small_dens)
                qyp(i,j,kc,QPRES) = max(qyp(i,j,kc,QPRES),small_pres)

             endif

          end if


          !--------------------------------------------------------------------
          ! minus state on face j+1
          !--------------------------------------------------------------------

          if (j .le. ihi2) then

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. v + cc <= 0.0d0 .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                u_ref    = u
                v_ref    = v
                w_ref    = w
                p_ref    = p
                rhoe_ref = rhoe
                tau_ref  = 1.d0/rho
                gam_ref  = gam
             else
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,kc,2,3,QRHO)
                u_ref    = Ip(i,j,kc,2,3,QU)
                v_ref    = Ip(i,j,kc,2,3,QV)
                w_ref    = Ip(i,j,kc,2,3,QW)
                p_ref    = Ip(i,j,kc,2,3,QPRES)
                rhoe_ref = Ip(i,j,kc,2,3,QREINT)
                tau_ref  = 1.d0/Ip(i,j,kc,2,3,QRHO)
                gam_ref  = Ip_gc(i,j,kc,2,3,1)
             endif

             ! for tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm    = (v_ref    - Ip(i,j,kc,2,1,QV))
             dpm    = (p_ref    - Ip(i,j,kc,2,1,QPRES))
   
             drho  = (rho_ref  - Ip(i,j,kc,2,2,QRHO))
             du    = (u_ref    - Ip(i,j,kc,2,2,QU))
             dw    = (w_ref    - Ip(i,j,kc,2,2,QW))
             dp    = (p_ref    - Ip(i,j,kc,2,2,QPRES))
             drhoe = (rhoe_ref - Ip(i,j,kc,2,2,QREINT))
             dtau  = (tau_ref  - 1.d0/Ip(i,j,kc,2,2,QRHO))

             dvp    = (v_ref    - Ip(i,j,kc,2,3,QV))
             dpp    = (p_ref    - Ip(i,j,kc,2,3,QPRES))

             ! if we are doing gravity tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_grav .eq. 1) then
                dvm = dvm - halfdt*Ip_g(i,j,kc,2,1,igy)
                du  = du  - halfdt*Ip_g(i,j,kc,2,2,igx)
                dw  = dw  - halfdt*Ip_g(i,j,kc,2,2,igz)
                dvp = dvp - halfdt*Ip_g(i,j,kc,2,3,igy)
             endif

             ! optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
             endif


             if (ppm_tau_in_tracing == 0) then

                ! these are analogous to the beta's from the original PPM
                ! paper.  This is simply (l . dq), where dq = qref - I(q)
 
                alpham = 0.5d0*(dpm/(rho_ev*cc_ev) - dvm)*rho_ev/cc_ev
                alphap = 0.5d0*(dpp/(rho_ev*cc_ev) + dvp)*rho_ev/cc_ev
                alpha0r = drho - dp/csq_ev
                alpha0e = drhoe - dp*enth_ev
                alpha0u = du
                alpha0w = dw

                if (v-cc .gt. 0.d0) then
                   amleft = -alpham
                else if (v-cc .lt. 0.d0) then
                   amleft = 0.d0
                else
                   amleft = -0.5d0*alpham
                endif
                
                if (v+cc .gt. 0.d0) then
                   apleft = -alphap
                else if (v+cc .lt. 0.d0) then
                   apleft = 0.d0
                else
                   apleft = -0.5d0*alphap
                endif
                
                if (v .gt. 0.d0) then
                   azrleft = -alpha0r
                   azeleft = -alpha0e
                   azu1left = -alpha0u
                   azw1left = -alpha0w
                else if (v .lt. 0.d0) then
                   azrleft = 0.d0
                   azeleft = 0.d0
                   azu1left = 0.d0
                   azw1left = 0.d0
                else
                   azrleft = -0.5d0*alpha0r
                   azeleft = -0.5d0*alpha0e
                   azu1left = -0.5d0*alpha0u
                   azw1left = -0.5d0*alpha0w
                endif
                
                ! the final interface states are just
                ! q_s = q_ref - sum (l . dq) r

                if (ppm_flatten_before_integrals == 0) then
                   xi1 = 1.0d0 - flatn(i,j,k3d)
                   xi = flatn(i,j,k3d)
                else
                   xi1 = 0.0d0
                   xi = 1.0d0
                endif
                
                qym(i,j+1,kc,QRHO  ) = xi1*rho  + xi*(rho_ref + apleft + amleft + azrleft)
                qym(i,j+1,kc,QV    ) = xi1*v    + xi*(v_ref + (apleft - amleft)*cc_ev/rho_ev)
                qym(i,j+1,kc,QU    ) = xi1*u    + xi*(u_ref + azu1left)
                qym(i,j+1,kc,QW    ) = xi1*w    + xi*(w_ref + azw1left)
                qym(i,j+1,kc,QREINT) = xi1*rhoe + xi*(rhoe_ref + (apleft + amleft)*enth_ev*csq_ev + azeleft)
                qym(i,j+1,kc,QPRES ) = xi1*p    + xi*(p_ref + (apleft + amleft)*csq_ev)
                
                qym(i,j+1,kc,QRHO ) = max(qym(i,j+1,kc,QRHO ),small_dens)
                qym(i,j+1,kc,QPRES) = max(qym(i,j+1,kc,QPRES),small_pres)

             else
                ! (tau, u, p, e) eigensystem

                ! this is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_ref/rho_ref - Ip(i,j,kc,2,2,QREINT)/Ip(i,j,kc,2,2,QRHO))

                alpham = 0.5d0*( dvm - dpm/Clag_ev)/Clag_ev
                alphap = 0.5d0*(-dvp - dpp/Clag_ev)/Clag_ev
                alpha0r = dtau + dp/Clag_ev**2
                alpha0e = de - dp*p_ev/Clag_ev**2
                alpha0u = du
                alpha0w = dw

                if (v-cc .gt. 0.d0) then
                   amleft = -alpham
                else if (v-cc .lt. 0.d0) then
                   amleft = 0.d0
                else
                   amleft = -0.5d0*alpham
                endif
                
                if (v+cc .gt. 0.d0) then
                   apleft = -alphap
                else if (v+cc .lt. 0.d0) then
                   apleft = 0.d0
                else
                   apleft = -0.5d0*alphap
                endif
                
                if (v .gt. 0.d0) then
                   azrleft = -alpha0r
                   azeleft = -alpha0e
                   azu1left = -alpha0u
                   azw1left = -alpha0w
                else if (v .lt. 0.d0) then
                   azrleft = 0.d0
                   azeleft = 0.d0
                   azu1left = 0.d0
                   azw1left = 0.d0
                else
                   azrleft = -0.5d0*alpha0r
                   azeleft = -0.5d0*alpha0e
                   azu1left = -0.5d0*alpha0u
                   azw1left = -0.5d0*alpha0w
                endif
                
                ! the final interface states are just
                ! q_s = q_ref - sum (l . dq) r

                if (ppm_flatten_before_integrals == 0) then
                   xi1 = 1.0d0 - flatn(i,j,k3d)
                   xi = flatn(i,j,k3d)
                else
                   xi1 = 0.0d0
                   xi = 1.0d0
                endif

                tau_s = tau_ref + apleft + amleft + azrleft
                qym(i,j+1,kc,QRHO  ) = xi1*rho  + xi/tau_s

                qym(i,j+1,kc,QV    ) = xi1*v    + xi*(v_ref + (amleft - apleft)*Clag_ev)
                qym(i,j+1,kc,QU    ) = xi1*u    + xi*(u_ref + azu1left)
                qym(i,j+1,kc,QW    ) = xi1*w    + xi*(w_ref + azw1left)

                e_s = rhoe_ref/rho_ref + (azeleft - p_ev*amleft - p_ev*apleft)
                qym(i,j+1,kc,QREINT) = xi1*rhoe + xi*e_s/tau_s

                qym(i,j+1,kc,QPRES ) = xi1*p    + xi*(p_ref + (-apleft - amleft)*Clag_ev**2)
                
                qym(i,j+1,kc,QRHO ) = max(qym(i,j+1,kc,QRHO ),small_dens)
                qym(i,j+1,kc,QPRES) = max(qym(i,j+1,kc,QPRES),small_pres)

             endif

          end if

       end do
    end do
    !$OMP END PARALLEL DO

    !--------------------------------------------------------------------------
    ! passively advected quantities
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
                xi = 1.0d0
             endif

             if (v .gt. 0.d0) then
                qyp(i,j,kc,n) = q(i,j,k3d,n)
             else if (v .lt. 0.d0) then
                qyp(i,j,kc,n) = q(i,j,k3d,n) + xi*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else
                qyp(i,j,kc,n) = q(i,j,k3d,n) + 0.5d0*xi*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo
          
          ! Minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,k3d,QV)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j,k3d)
             else
                xi = 1.0d0
             endif

             if (v .gt. 0.d0) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n) + xi*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else if (v .lt. 0.d0) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n)
             else
                qym(i,j+1,kc,n) = q(i,j,k3d,n) + 0.5d0*xi*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo
          
       enddo
    enddo

  end subroutine tracexy_ppm




  subroutine tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        Ip,Im,Ip_g,Im_g,Ip_gc,Im_gc, &
                        qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                        gamc,gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3, &
                        ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QESGS, QPRES, QFA, QFS, nadv, &
         small_dens, small_pres, &
         ppm_type, ppm_reference, ppm_trace_grav, &
         ppm_tau_in_tracing, ppm_reference_eigenvectors, &
         ppm_reference_edge_limit, ppm_flatten_before_integrals

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
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

    double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2,gc_l3:gc_h3)

    double precision dz, dt

    !     Local variables
    integer i, j
    integer n, iadv, ispec
    integer npassive,ipassive,qpass_map(QVAR)

    double precision cc, csq, Clag, rho, u, v, w, p, rhoe

    double precision drho, du, dv, dp, drhoe, de, dtau
    double precision dwp, dpp
    double precision dwm, dpm

    double precision :: rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref, tau_ref
    double precision :: tau_s, e_s

    double precision :: cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref
    double precision :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev
    double precision :: gam

    double precision enth, alpham, alphap, alpha0r, alpha0e
    double precision alpha0u, alpha0v
    double precision apright, amright, azrright, azeright
    double precision azu1rght, azv1rght
    double precision apleft, amleft, azrleft, azeleft
    double precision azu1left, azv1left

    double precision xi, xi1
    double precision halfdt

    integer, parameter :: igx = 1
    integer, parameter :: igy = 2
    integer, parameter :: igz = 3

    ! Group all the passively advected quantities together
    npassive = 0
    if (QESGS .gt. -1) then
       qpass_map(1) = QESGS
       npassive = 1
    endif
    do iadv = 1, nadv
       qpass_map(npassive + iadv) = QFA + iadv - 1
    enddo
    npassive = npassive + nadv
    if(QFS .gt. -1) then
       do ispec = 1, nspec+naux
          qpass_map(npassive + ispec) = QFS + ispec - 1
       enddo
       npassive = npassive + nspec + naux
    endif

    halfdt = 0.5d0 * dt

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
    !$OMP PRIVATE(rho_ref,u_ref,v_ref,w_ref,p_ref,rhoe_ref) &
    !$OMP PRIVATE(drho,du,dv,dp,drhoe,de,dwm,dpm,dwp,dpp,alpham,alphap,alpha0r,alpha0e) &
    !$OMP PRIVATE(alpha0u,alpha0v,amright,apright,azrright,azeright,azu1rght,azv1rght,amleft,apleft)&
    !$OMP PRIVATE(azrleft,azeleft,azu1left,azv1left,xi,xi1) &
    !$OMP PRIVATE(cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref) &
    !$OMP PRIVATE(cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev) &
    !$OMP PRIVATE(gam, dtau, tau_ref, tau_s, e_s)

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


          ! Set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. w - cc >= 0.0d0 .and. &
                ppm_reference_edge_limit == 0) ) then
             ! original Castro way -- cc value
             rho_ref  = rho
             u_ref    = u
             v_ref    = v
             w_ref    = w
             p_ref    = p
             rhoe_ref = rhoe
             tau_ref  = 1.d0/rho
             gam_ref  = gam
          else
             ! This will be the fastest moving state to the left
             rho_ref  = Im(i,j,kc,3,1,QRHO)
             u_ref    = Im(i,j,kc,3,1,QU)
             v_ref    = Im(i,j,kc,3,1,QV)
             w_ref    = Im(i,j,kc,3,1,QW)
             p_ref    = Im(i,j,kc,3,1,QPRES)
             rhoe_ref = Im(i,j,kc,3,1,QREINT)
             tau_ref  = 1.d0/Im(i,j,kc,3,1,QRHO)
             gam_ref  = Im_gc(i,j,kc,3,1,1)
          endif

          ! for tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm    = (w_ref    - Im(i,j,kc,3,1,QW))
          dpm    = (p_ref    - Im(i,j,kc,3,1,QPRES))

          drho  = (rho_ref  - Im(i,j,kc,3,2,QRHO))
          du    = (u_ref    - Im(i,j,kc,3,2,QU))
          dv    = (v_ref    - Im(i,j,kc,3,2,QV))
          dp    = (p_ref    - Im(i,j,kc,3,2,QPRES))
          drhoe = (rhoe_ref - Im(i,j,kc,3,2,QREINT))
          dtau  = (tau_ref  - 1.d0/Im(i,j,kc,3,2,QRHO))

          dwp    = (w_ref    - Im(i,j,kc,3,3,QW))
          dpp    = (p_ref    - Im(i,j,kc,3,3,QPRES))

          ! if we are doing gravity tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_grav .eq. 1) then
             dwm = dwm - halfdt*Im_g(i,j,kc,3,1,igz)
             du  = du  - halfdt*Im_g(i,j,kc,3,2,igx)
             dv  = dv  - halfdt*Im_g(i,j,kc,3,2,igy)
             dwp = dwp - halfdt*Im_g(i,j,kc,3,3,igz)
          endif

          ! optionally use the reference state in evaluating the
          ! eigenvectors
          if (ppm_reference_eigenvectors == 0) then
             rho_ev  = rho
             cc_ev   = cc
             csq_ev  = csq
             Clag_ev = Clag
             enth_ev = enth
             p_ev    = p
          else
             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             enth_ev = enth_ref
             p_ev    = p_ref
          endif

          if (ppm_tau_in_tracing == 0) then

             ! these are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)
             alpham = 0.5d0*(dpm/(rho_ev*cc_ev) - dwm)*rho_ev/cc_ev
             alphap = 0.5d0*(dpp/(rho_ev*cc_ev) + dwp)*rho_ev/cc_ev
             alpha0r = drho - dp/csq_ev
             alpha0e = drhoe - dp*enth_ev
             alpha0u = du
             alpha0v = dv

             if (w-cc .gt. 0.d0) then
                amright = 0.d0
             else if (w-cc .lt. 0.d0) then
                amright = -alpham
             else
                amright = -0.5d0*alpham
             endif
             if (w+cc .gt. 0.d0) then
                apright = 0.d0
             else if (w+cc .lt. 0.d0) then
                apright = -alphap
             else
                apright = -0.5d0*alphap
             endif
             if (w .gt. 0.d0) then
                azrright = 0.d0
                azeright = 0.d0
                azu1rght = 0.d0
                azv1rght = 0.d0
             else if (w .lt. 0.d0) then
                azrright = -alpha0r
                azeright = -alpha0e
                azu1rght = -alpha0u
                azv1rght = -alpha0v
             else
                azrright = -0.5d0*alpha0r
                azeright = -0.5d0*alpha0e
                azu1rght = -0.5d0*alpha0u
                azv1rght = -0.5d0*alpha0v
             endif

             ! the final interface states are just
             ! q_s = q_ref - sum (l . dq) r

             if (ppm_flatten_before_integrals == 0) then
                xi1 = 1.0d0 - flatn(i,j,k3d)
                xi = flatn(i,j,k3d)
             else
                xi1 = 0.0d0
                xi = 1.0d0
             endif
             
             qzp(i,j,kc,QRHO  ) = xi1*rho  + xi*(rho_ref + apright + amright + azrright)
             qzp(i,j,kc,QW    ) = xi1*w    + xi*(w_ref + (apright - amright)*cc_ev/rho_ev)
             qzp(i,j,kc,QU    ) = xi1*u    + xi*(u_ref + azu1rght)
             qzp(i,j,kc,QV    ) = xi1*v    + xi*(v_ref + azv1rght)
             qzp(i,j,kc,QREINT) = xi1*rhoe + xi*(rhoe_ref + (apright + amright)*enth_ev*csq_ev + azeright)
             qzp(i,j,kc,QPRES ) = xi1*p    + xi*(p_ref + (apright + amright)*csq_ev)
             
             qzp(i,j,kc,QRHO ) = max(qzp(i,j,kc,QRHO ),small_dens)
             qzp(i,j,kc,QPRES) = max(qzp(i,j,kc,QPRES),small_pres)

          else
             ! (tau, u, p, e) eigensystem
             
             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system
             de = (rhoe_ref/rho_ref - Im(i,j,kc,3,2,QREINT)/Im(i,j,kc,3,2,QRHO))

             ! these are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)
             alpham = 0.5d0*( dwm - dpm/Clag_ev)/Clag_ev
             alphap = 0.5d0*(-dwp - dpp/Clag_ev)/Clag_ev
             alpha0r = dtau + dp/Clag_ev**2
             alpha0e = de - dp*p_ev/Clag_ev**2
             alpha0u = du
             alpha0v = dv

             if (w-cc .gt. 0.d0) then
                amright = 0.d0
             else if (w-cc .lt. 0.d0) then
                amright = -alpham
             else
                amright = -0.5d0*alpham
             endif
             if (w+cc .gt. 0.d0) then
                apright = 0.d0
             else if (w+cc .lt. 0.d0) then
                apright = -alphap
             else
                apright = -0.5d0*alphap
             endif
             if (w .gt. 0.d0) then
                azrright = 0.d0
                azeright = 0.d0
                azu1rght = 0.d0
                azv1rght = 0.d0
             else if (w .lt. 0.d0) then
                azrright = -alpha0r
                azeright = -alpha0e
                azu1rght = -alpha0u
                azv1rght = -alpha0v
             else
                azrright = -0.5d0*alpha0r
                azeright = -0.5d0*alpha0e
                azu1rght = -0.5d0*alpha0u
                azv1rght = -0.5d0*alpha0v
             endif

             ! the final interface states are just
             ! q_s = q_ref - sum (l . dq) r

             if (ppm_flatten_before_integrals == 0) then
                xi1 = 1.0d0 - flatn(i,j,k3d)
                xi = flatn(i,j,k3d)
             else
                xi1 = 0.0d0
                xi = 1.0d0
             endif

             tau_s = tau_ref + apright + amright + azrright
             qzp(i,j,kc,QRHO  ) = xi1*rho  + xi/tau_s

             qzp(i,j,kc,QW    ) = xi1*w    + xi*(w_ref + (amright - apright)*Clag_ev)
             qzp(i,j,kc,QU    ) = xi1*u    + xi*(u_ref + azu1rght)
             qzp(i,j,kc,QV    ) = xi1*v    + xi*(v_ref + azv1rght)

             e_s = rhoe_ref/rho_ref + (azeright - p_ev*amright - p_ev*apright)           
             qzp(i,j,kc,QREINT) = xi1*rhoe + xi*e_s/tau_s

             qzp(i,j,kc,QPRES ) = xi1*p    + xi*(p_ref + (-apright - amright)*Clag_ev**2)
             
             qzp(i,j,kc,QRHO ) = max(qzp(i,j,kc,QRHO ),small_dens)
             qzp(i,j,kc,QPRES) = max(qzp(i,j,kc,QPRES),small_pres)

          endif


          !--------------------------------------------------------------------
          ! This is all for qzm -- minus state on face kc
          !--------------------------------------------------------------------

          ! note this is different from how we do 1D, 2D, and the
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


          ! Set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. w + cc <= 0.0d0 .and. &
                ppm_reference_edge_limit == 0) ) then
             rho_ref  = rho
             u_ref    = u
             v_ref    = v
             w_ref    = w
             p_ref    = p
             rhoe_ref = rhoe
             tau_ref  = 1.d0/rho
             gam_ref  = gam
          else
             ! This will be the fastest moving state to the right
             rho_ref  = Ip(i,j,km,3,3,QRHO)
             u_ref    = Ip(i,j,km,3,3,QU)
             v_ref    = Ip(i,j,km,3,3,QV)
             w_ref    = Ip(i,j,km,3,3,QW)
             p_ref    = Ip(i,j,km,3,3,QPRES)
             rhoe_ref = Ip(i,j,km,3,3,QREINT)
             tau_ref  = 1.d0/Ip(i,j,km,3,3,QRHO)
             gam_ref  = Ip_gc(i,j,km,3,3,1)
          endif

          ! for tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm    = (w_ref    - Ip(i,j,km,3,1,QW))
          dpm    = (p_ref    - Ip(i,j,km,3,1,QPRES))

          drho  = (rho_ref  - Ip(i,j,km,3,2,QRHO))
          du    = (u_ref    - Ip(i,j,km,3,2,QU))
          dv    = (v_ref    - Ip(i,j,km,3,2,QV))
          dp    = (p_ref    - Ip(i,j,km,3,2,QPRES))
          drhoe = (rhoe_ref - Ip(i,j,km,3,2,QREINT))
          dtau  = (tau_ref  - 1.d0/Ip(i,j,km,3,2,QRHO))

          dwp    = (w_ref    - Ip(i,j,km,3,3,QW))
          dpp    = (p_ref    - Ip(i,j,km,3,3,QPRES))

          ! if we are doing gravity tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_grav .eq. 1) then
             dwm = dwm - halfdt*Ip_g(i,j,km,3,1,igz)
             du  = du  - halfdt*Ip_g(i,j,km,3,2,igx)
             dv  = dv  - halfdt*Ip_g(i,j,km,3,2,igy)
             dwp = dwp - halfdt*Ip_g(i,j,km,3,3,igz)
          endif

          ! optionally use the reference state in evaluating the
          ! eigenvectors
          if (ppm_reference_eigenvectors == 0) then
             rho_ev  = rho
             cc_ev   = cc
             csq_ev  = csq
             Clag_ev = Clag
             enth_ev = enth
             p_ev    = p
          else
             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             enth_ev = enth_ref
             p_ev    = p_ref
          endif

          if (ppm_tau_in_tracing == 0) then

             ! these are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)
             alpham = 0.5d0*(dpm/(rho_ev*cc_ev) - dwm)*rho_ev/cc_ev
             alphap = 0.5d0*(dpp/(rho_ev*cc_ev) + dwp)*rho_ev/cc_ev
             alpha0r = drho - dp/csq_ev
             alpha0e = drhoe - dp*enth_ev
             alpha0u = du
             alpha0v = dv
             
             if (w-cc .gt. 0.d0) then
                amleft = -alpham
             else if (w-cc .lt. 0.d0) then
                amleft = 0.d0
             else
                amleft = -0.5d0*alpham
             endif
             if (w+cc .gt. 0.d0) then
                apleft = -alphap
             else if (w+cc .lt. 0.d0) then
                apleft = 0.d0
             else
                apleft = -0.5d0*alphap
             endif
             if (w .gt. 0.d0) then
                azrleft = -alpha0r
                azeleft = -alpha0e
                azu1left = -alpha0u
                azv1left = -alpha0v
             else if (w .lt. 0.d0) then
                azrleft = 0.d0
                azeleft = 0.d0
                azu1left = 0.d0
                azv1left = 0.d0
             else
                azrleft = -0.5d0*alpha0r
                azeleft = -0.5d0*alpha0e
                azu1left = -0.5d0*alpha0u
                azv1left = -0.5d0*alpha0v
             endif
             
             ! the final interface states are just
             ! q_s = q_ref - sum (l . dq) r

             if (ppm_flatten_before_integrals == 0) then
                xi1 = 1.0d0 - flatn(i,j,k3d-1)
                xi = flatn(i,j,k3d-1)
             else
                xi1 = 0.0d0
                xi = 1.0d0
             endif

             qzm(i,j,kc,QRHO  ) = xi1*rho  + xi*(rho_ref + apleft + amleft + azrleft)
             qzm(i,j,kc,QW    ) = xi1*w    + xi*(w_ref + (apleft - amleft)*cc_ev/rho_ev)
             qzm(i,j,kc,QU    ) = xi1*u    + xi*(u_ref + azu1left)
             qzm(i,j,kc,QV    ) = xi1*v    + xi*(v_ref + azv1left)
             qzm(i,j,kc,QREINT) = xi1*rhoe + xi*(rhoe_ref + (apleft + amleft)*enth_ev*csq_ev + azeleft)
             qzm(i,j,kc,QPRES ) = xi1*p    + xi*(p_ref + (apleft + amleft)*csq_ev)
             
             qzm(i,j,kc,QRHO ) = max(qzm(i,j,kc,QRHO ),small_dens)
             qzm(i,j,kc,QPRES) = max(qzm(i,j,kc,QPRES),small_pres)
             
          else
             ! (tau, u, p, e) eigensystem
             
             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system
             de = (rhoe_ref/rho_ref - Ip(i,j,km,3,2,QREINT)/Ip(i,j,km,3,2,QRHO))

             alpham = 0.5d0*( dwm - dpm/Clag_ev)/Clag_ev
             alphap = 0.5d0*(-dwp - dpp/Clag_ev)/Clag_ev
             alpha0r = dtau + dp/Clag_ev**2
             alpha0e = de - dp*p_ev/Clag_ev**2
             alpha0u = du
             alpha0v = dv
             
             if (w-cc .gt. 0.d0) then
                amleft = -alpham
             else if (w-cc .lt. 0.d0) then
                amleft = 0.d0
             else
                amleft = -0.5d0*alpham
             endif
             if (w+cc .gt. 0.d0) then
                apleft = -alphap
             else if (w+cc .lt. 0.d0) then
                apleft = 0.d0
             else
                apleft = -0.5d0*alphap
             endif
             if (w .gt. 0.d0) then
                azrleft = -alpha0r
                azeleft = -alpha0e
                azu1left = -alpha0u
                azv1left = -alpha0v
             else if (w .lt. 0.d0) then
                azrleft = 0.d0
                azeleft = 0.d0
                azu1left = 0.d0
                azv1left = 0.d0
             else
                azrleft = -0.5d0*alpha0r
                azeleft = -0.5d0*alpha0e
                azu1left = -0.5d0*alpha0u
                azv1left = -0.5d0*alpha0v
             endif
             
             ! the final interface states are just
             ! q_s = q_ref - sum (l . dq) r

             if (ppm_flatten_before_integrals == 0) then
                xi1 = 1.0d0 - flatn(i,j,k3d-1)
                xi = flatn(i,j,k3d-1)
             else
                xi1 = 0.0d0
                xi = 1.0d0
             endif

             tau_s = tau_ref + apleft + amleft + azrleft
             qzm(i,j,kc,QRHO  ) = xi1*rho  + xi/tau_s

             qzm(i,j,kc,QW    ) = xi1*w    + xi*(w_ref + (amleft - apleft)*Clag_ev)
             qzm(i,j,kc,QU    ) = xi1*u    + xi*(u_ref + azu1left)
             qzm(i,j,kc,QV    ) = xi1*v    + xi*(v_ref + azv1left)

             e_s = rhoe_ref/rho_ref + (azeleft - p_ev*amleft - p_ev*apleft)
             qzm(i,j,kc,QREINT) = xi1*rhoe + xi*e_s/tau_s

             qzm(i,j,kc,QPRES ) = xi1*p    + xi*(p_ref + (-apleft - amleft)*Clag_ev**2)
             
             qzm(i,j,kc,QRHO ) = max(qzm(i,j,kc,QRHO ),small_dens)
             qzm(i,j,kc,QPRES) = max(qzm(i,j,kc,QPRES),small_pres)

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
                     xi = 1.0d0
                  endif

                  if (w .gt. 0.d0) then
                     qzp(i,j,kc,n) = q(i,j,k3d,n)
                  else if (w .lt. 0.d0) then
                     qzp(i,j,kc,n) = q(i,j,k3d,n) + xi*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
                  else
                     qzp(i,j,kc,n) = q(i,j,k3d,n) + 0.5d0*xi*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
                  endif

                  ! Minus state on face k
                  w = q(i,j,k3d-1,QW)
                  
                  if (ppm_flatten_before_integrals == 0) then
                     xi = flatn(i,j,k3d-1)
                  else
                     xi = 1.0d0
                  endif

                  if (w .gt. 0.d0) then
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n) + xi*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
                  else if (w .lt. 0.d0) then
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n)
                  else
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n) + 0.5d0*xi*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
                  endif

               enddo
         enddo
    enddo
    !$OMP end parallel do

  end subroutine tracez_ppm

end module trace_ppm_module
