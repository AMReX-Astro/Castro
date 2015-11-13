module trace_ppm_module

  implicit none

  private

  public trace_ppm

contains

  subroutine trace_ppm(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                       dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                       qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       src,src_l1,src_l2,src_h1,src_h2, &
                       gamc,gc_l1,gc_l2,gc_h1,gc_h2, &
                       ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

    use network, only : nspec
    use eos_type_module
    use eos_module
    use bl_constants_module
    use meth_params_module, only : QVAR, QRHO, QU, QV, &
         QREINT, QPRES, QTEMP, QFS, QGAME, &
         small_dens, small_pres, &
         ppm_type, ppm_reference, ppm_trace_sources, ppm_temp_fix, &
         ppm_tau_in_tracing, ppm_reference_eigenvectors, ppm_reference_edge_limit, &
         ppm_flatten_before_integrals, ppm_predict_gammae, &
         npassive, qpass_map
    use ppm_module, only : ppm

    implicit none

    integer ilo1,ilo2,ihi1,ihi2
    integer qd_l1,qd_l2,qd_h1,qd_h2
    integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
    integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer src_l1,src_l2,src_h1,src_h2
    integer gc_l1,gc_l2,gc_h1,gc_h2

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
    
    double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)

    double precision  src(src_l1:src_h1,src_l2:src_h2,QVAR)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)

    double precision dx, dy, dt
    
    ! Local variables
    integer          :: i, j, iwave, idim
    integer          :: n, ipassive
    
    double precision :: dtdx, dtdy
    double precision :: cc, csq, Clag, rho, u, v, p, rhoe
    double precision :: drho, dp, drhoe, dtau
    double precision :: dup, dvp, dpp
    double precision :: dum, dvm, dpm

    double precision :: rho_ref, u_ref, v_ref, p_ref, rhoe_ref, tau_ref
    double precision :: tau_s, e_s, de, dge

    double precision :: cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref, game_ref, gfactor
    double precision :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev, tau_ev
    double precision :: gam, game
    
    double precision :: enth, alpham, alphap, alpha0r, alpha0e
    double precision :: apright, amright, azrright, azeright
    double precision :: azu1rght, azv1rght
    double precision :: apleft, amleft, azrleft, azeleft
    double precision :: azu1left, azv1left
    double precision :: sourcr,sourcp,source,courn,eta,dlogatmp

    double precision :: xi, xi1
    double precision :: halfdt

    double precision, allocatable :: Ip(:,:,:,:,:)
    double precision, allocatable :: Im(:,:,:,:,:)

    double precision, allocatable :: Ip_src(:,:,:,:,:)
    double precision, allocatable :: Im_src(:,:,:,:,:)

    ! gamma_c/1 on the interfaces
    double precision, allocatable :: Ip_gc(:,:,:,:,:)
    double precision, allocatable :: Im_gc(:,:,:,:,:)

    type (eos_t) :: eos_state

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
       call bl_error("Error:: ppm_2d.f90 :: trace_ppm")
    end if

    dtdx = dt/dx
    dtdy = dt/dy

    ! indices: (x, y, dimension, wave, variable)
    allocate(Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QVAR))
    allocate(Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QVAR))

    if (ppm_trace_sources .eq. 1) then
       allocate(Ip_src(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QVAR))
       allocate(Im_src(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QVAR))
    endif

    allocate(Ip_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,1))
    allocate(Im_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,1))

    halfdt = HALF * dt


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
    ! The indices are: Ip(i, j, dim, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.


    ! Compute Ip and Im -- this does the parabolic reconstruction,
    ! limiting, and returns the integral of each profile under
    ! each wave to each interface
    do n=1,QVAR
       call ppm(q(:,:,n),qd_l1,qd_l2,qd_h1,qd_h2, &
                q(:,:,QU:QV),c,qd_l1,qd_l2,qd_h1,qd_h2, &
                flatn, &
                Ip(:,:,:,:,n),Im(:,:,:,:,n), &
                ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    end do

    ! temperature-based PPM -- if desired, take the Ip(T)/Im(T)
    ! constructed above and use the EOS to overwrite Ip(p)/Im(p)
    if (ppm_temp_fix == 1) then
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1
             do idim = 1, 2
                do iwave = 1, 3
                   eos_state%rho   = Ip(i,j,idim,iwave,QRHO)
                   eos_state%T     = Ip(i,j,idim,iwave,QTEMP)
                   eos_state%xn(:) = Ip(i,j,idim,iwave,QFS:QFS-1+nspec)

                   call eos(eos_input_rt, eos_state)

                   Ip(i,j,idim,iwave,QPRES) = eos_state%p
                   Ip(i,j,idim,iwave,QREINT) = Ip(i,j,idim,iwave,QRHO)*eos_state%e
                   Ip_gc(i,j,idim,iwave,1) = eos_state%gam1

                   eos_state%rho   = Im(i,j,idim,iwave,QRHO)
                   eos_state%T     = Im(i,j,idim,iwave,QTEMP)
                   eos_state%xn(:) = Im(i,j,idim,iwave,QFS:QFS-1+nspec)

                   call eos(eos_input_rt, eos_state)

                   Im(i,j,idim,iwave,QPRES) = eos_state%p
                   Im(i,j,idim,iwave,QREINT) = Im(i,j,idim,iwave,QRHO)*eos_state%e
                   Im_gc(i,j,idim,iwave,1) = eos_state%gam1
                enddo
             enddo
          enddo
       enddo

    endif

    ! get an edge-based gam1 here if we didn't get it from the EOS
    ! call above (for ppm_temp_fix = 1)
    if (ppm_temp_fix /= 1) then
          call ppm(gamc(:,:),gc_l1,gc_l2,gc_h1,gc_h2, &
                   q(:,:,QU:QV),c,qd_l1,qd_l2,qd_h1,qd_h2, &
                   flatn, &
                   Ip_gc(:,:,:,:,1),Im_gc(:,:,:,:,1), &
                   ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
       endif

    ! if desired, do parabolic reconstruction of the momentum sources
    ! -- we'll use this for the force on the velocity
    if (ppm_trace_sources == 1) then
       do n = 1, QVAR
          call ppm(src(:,:,n),src_l1,src_l2,src_h1,src_h2, &
                   q(:,:,QU:QV),c,qd_l1,qd_l2,qd_h1,qd_h2, &
                   flatn, &
                   Ip_src(:,:,:,:,n),Im_src(:,:,:,:,n), &
                   ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
       enddo
    endif



    !-------------------------------------------------------------------------
    ! x-direction
    !-------------------------------------------------------------------------

    ! Trace to left and right edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          cc = c(i,j)
          csq = cc**2

          rho = q(i,j,QRHO)
          u = q(i,j,QU)
          v = q(i,j,QV)

          p = q(i,j,QPRES)
          rhoe = q(i,j,QREINT)
          enth = ( (rhoe+p)/rho )/csq

          Clag = rho*cc

          game = q(i,j,QGAME)

          gam = gamc(i,j)

          !-------------------------------------------------------------------
          ! plus state on face i
          !-------------------------------------------------------------------

          ! set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. u - cc >= ZERO .and. &
                ppm_reference_edge_limit == 0)) then
             ! original Castro way -- cc value
             rho_ref  = rho
             u_ref    = u

             p_ref    = p
             rhoe_ref = rhoe

             tau_ref  = ONE/rho

             gam_ref = gamc(i,j)

             game_ref = game

          else
             ! this will be the fastest moving state to the left --
             ! this is the method that Miller & Colella and Colella &
             ! Woodward use
             rho_ref  = Im(i,j,1,1,QRHO)
             u_ref    = Im(i,j,1,1,QU)

             p_ref    = Im(i,j,1,1,QPRES)
             rhoe_ref = Im(i,j,1,1,QREINT)

             tau_ref  = ONE/Im(i,j,1,1,QRHO)

             gam_ref  = Im_gc(i,j,1,1,1)

             game_ref  = Im(i,j,1,1,QGAME)
          endif

          ! for tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

          ! *m are the jumps carried by u-c
          ! *p are the jumps carried by u+c

          dum    = u_ref    - Im(i,j,1,1,QU)
          dpm    = p_ref    - Im(i,j,1,1,QPRES)

          drho  = rho_ref  - Im(i,j,1,2,QRHO)
          dp    = p_ref    - Im(i,j,1,2,QPRES)
          drhoe = rhoe_ref - Im(i,j,1,2,QREINT)
          dtau  = tau_ref  - ONE/Im(i,j,1,2,QRHO)

          dup    = u_ref    - Im(i,j,1,3,QU)
          dpp    = p_ref    - Im(i,j,1,3,QPRES)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dum = dum - halfdt*Im_src(i,j,1,1,QU)
             dup = dup - halfdt*Im_src(i,j,1,3,QU)
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

             ! these are analogous to the beta's from the original 
             ! PPM paper (except we work with rho instead of tau).  
             ! This is simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dpm/(rho_ev*cc_ev) - dum)*rho_ev/cc_ev
             alphap = HALF*(dpp/(rho_ev*cc_ev) + dup)*rho_ev/cc_ev
             alpha0r = drho - dp/csq_ev
             alpha0e = drhoe - dp*enth_ev  ! note enth has a 1/c**2 in it

          else

             ! (tau, u, p, e) eigensystem
             ! or
             ! (tau, u, p, game) eigensystem

             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system.

             ! we are dealing with e
             de = (rhoe_ref/rho_ref - Im(i,j,1,2,QREINT)/Im(i,j,1,2,QRHO))

             dge = game_ref - Im(i,j,1,2,QGAME)

             alpham = HALF*( dum - dpm/Clag_ev)/Clag_ev
             alphap = HALF*(-dup - dpp/Clag_ev)/Clag_ev
             alpha0r = dtau + dp/Clag_ev**2

             if (ppm_predict_gammae == 0) then
                alpha0e = de - dp*p_ev/Clag_ev**2
             else
                gfactor = (game - 1.0d0)*(game - gam)
                alpha0e = gfactor*dp/(tau_ev*Clag_ev**2) + dge
             endif

          endif ! which tracing method

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

          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          if (i .ge. ilo1) then

             if (ppm_tau_in_tracing == 0) then
                qxp(i,j,QRHO)   = rho_ref + apright + amright + azrright
                qxp(i,j,QU)     = u_ref + (apright - amright)*cc_ev/rho_ev
                qxp(i,j,QREINT) = rhoe_ref + (apright + amright)*enth_ev*csq_ev + azeright
                qxp(i,j,QPRES)  = p_ref + (apright + amright)*csq_ev
             else
                tau_s = tau_ref + apright + amright + azrright
                qxp(i,j,QRHO)   = ONE/tau_s

                qxp(i,j,QU)     = u_ref + (amright - apright)*Clag_ev

                qxp(i,j,QPRES)  = p_ref + (-apright - amright)*Clag_ev**2
                
                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_ref/rho_ref + (azeright - p_ev*amright -p_ev*apright)
                   qxp(i,j,QREINT) = e_s/tau_s
                else
                   qxp(i,j,QGAME) = game_ref + gfactor*(amright + apright)/tau_ev + azeright 
                   qxp(i,j,QREINT) = qxp(i,j,QPRES )/(qxp(i,j,QGAME) - 1.0d0)
                endif
             end if


             ! enforce small_*
             qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))
             qxp(i,j,QPRES) = max(qxp(i,j,QPRES), small_pres)

             ! transverse velocity -- there is no projection here, so
             ! we don't need a reference state.  We only care about 
             ! the state traced under the middle wave

             ! Recall that I already takes the limit of the parabola
             ! in the event that the wave is not moving toward the
             ! interface
             if (u > ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qxp(i,j,QV) = Im(i,j,1,2,QV)
                else
                   qxp(i,j,QV) = v
                endif
             else ! wave moving toward the interface
                qxp(i,j,QV) = Im(i,j,1,2,QV)
             endif

             if (ppm_trace_sources == 1) then
                qxp(i,j,QV)  = qxp(i,j,QV)  + halfdt*Im_src(i,j,1,2,QV)
             endif

             ! we may have done the flattening already in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE-flatn(i,j)
                xi = flatn(i,j)

                qxp(i,j,QRHO)   = xi1*rho  + xi*qxp(i,j,QRHO)
                qxp(i,j,QU)     = xi1*u    + xi*qxp(i,j,QU)
                qxp(i,j,QREINT) = xi1*rhoe + xi*qxp(i,j,QREINT)
                qxp(i,j,QPRES)  = xi1*p    + xi*qxp(i,j,QPRES) 
                qxp(i,j,QV)     = xi1*v    + xi*qxp(i,j,QV) 
             endif


          endif


          !-------------------------------------------------------------------
          ! minus state on face i+1
          !-------------------------------------------------------------------

          ! set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. u + cc <= ZERO .and. &
                ppm_reference_edge_limit == 0) ) then
             ! original Castro way -- cc values
             rho_ref  = rho
             u_ref    = u

             p_ref    = p
             rhoe_ref = rhoe

             tau_ref = ONE/rho

             gam_ref = gamc(i,j)

             game_ref = game

          else
             ! this will be the fastest moving state to the right
             rho_ref  = Ip(i,j,1,3,QRHO)
             u_ref    = Ip(i,j,1,3,QU)

             p_ref    = Ip(i,j,1,3,QPRES)
             rhoe_ref = Ip(i,j,1,3,QREINT)

             tau_ref  = ONE/Ip(i,j,1,3,QRHO)

             gam_ref    = Ip_gc(i,j,1,3,1)

             game_ref    = Ip(i,j,1,3,QGAME)
          endif

          ! for tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

          ! *m are the jumps carried by u-c
          ! *p are the jumps carried by u+c

          dum    = u_ref    - Ip(i,j,1,1,QU)
          dpm    = p_ref    - Ip(i,j,1,1,QPRES)
          
          drho  = rho_ref  - Ip(i,j,1,2,QRHO)
          dp    = p_ref    - Ip(i,j,1,2,QPRES)
          drhoe = rhoe_ref - Ip(i,j,1,2,QREINT)
          dtau  = tau_ref  - ONE/Ip(i,j,1,2,QRHO)
          
          dup    = u_ref    - Ip(i,j,1,3,QU)
          dpp    = p_ref    - Ip(i,j,1,3,QPRES)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dum = dum - halfdt*Ip_src(i,j,1,1,QU)
             dup = dup - halfdt*Ip_src(i,j,1,3,QU)
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

             ! these are analogous to the beta's from the original 
             ! PPM paper (except we work with rho instead of tau).
             ! This is simply (l . dq), where dq = qref - I(q)          
             alpham = HALF*(dpm/(rho_ev*cc_ev) - dum)*rho_ev/cc_ev
             alphap = HALF*(dpp/(rho_ev*cc_ev) + dup)*rho_ev/cc_ev
             alpha0r = drho - dp/csq_ev
             alpha0e = drhoe - dp*enth_ev  ! enth has a 1/c**2 in it

          else

             ! (tau, u, p, e) eigensystem
             ! or
             ! (tau, u, p, game) eigensystem

             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system.

             de = (rhoe_ref/rho_ref - Ip(i,j,1,2,QREINT)/Ip(i,j,1,2,QRHO))
             dge = game_ref - Ip(i,j,1,2,QGAME)

             alpham = HALF*( dum - dpm/Clag_ev)/Clag_ev
             alphap = HALF*(-dup - dpp/Clag_ev)/Clag_ev
             alpha0r = dtau + dp/Clag_ev**2

             if (ppm_predict_gammae == 0) then
                alpha0e = de - dp*p_ev/Clag_ev**2
             else
                gfactor = (game - 1.0d0)*(game - gam)
                alpha0e = gfactor*dp/(tau_ev*Clag_ev**2) + dge
             endif

          endif

             
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
          
          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          if (i .le. ihi1) then
             if (ppm_tau_in_tracing == 0) then

                qxm(i+1,j,QRHO)   = rho_ref + apleft + amleft + azrleft
                qxm(i+1,j,QU)     = u_ref + (apleft - amleft)*cc_ev/rho_ev  
                qxm(i+1,j,QREINT) = rhoe_ref + (apleft + amleft)*enth_ev*csq_ev + azeleft
                qxm(i+1,j,QPRES)  = p_ref + (apleft + amleft)*csq_ev
             else
                tau_s = tau_ref + (apleft + amleft + azrleft)
                qxm(i+1,j,QRHO)   = ONE/tau_s

                qxm(i+1,j,QU)     = u_ref + (amleft - apleft)*Clag_ev

                qxm(i+1,j,QPRES)  = p_ref + (-apleft - amleft)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_ref/rho_ref + (azeleft - p_ev*amleft -p_ev*apleft)
                   qxm(i+1,j,QREINT) = e_s/tau_s
                else
                   qxm(i+1,j,QGAME) = game_ref + gfactor*(amleft + apleft)/tau_ev + azeleft
                   qxm(i+1,j,QREINT) = qxm(i+1,j,QPRES )/(qxm(i+1,j,QGAME) - 1.0d0)
                endif

             end if


             ! enforce small_*
             qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
             qxm(i+1,j,QPRES) = max(qxm(i+1,j,QPRES), small_pres)

             ! transverse velocity -- there is no projection here, so
             ! we don't need a reference state.  We only care about 
             ! the state traced under the middle wave
             if (u < ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qxm(i+1,j,QV) = Ip(i,j,1,2,QV)
                else
                   qxm(i+1,j,QV) = v
                endif
             else ! wave moving toward interface
                qxm(i+1,j,QV) = Ip(i,j,1,2,QV)
             endif

             if (ppm_trace_sources == 1) then
                qxm(i+1,j,QV)  = qxm(i+1,j,QV)  + halfdt*Ip_src(i,j,1,2,QV)
             endif

             ! we may have already done the flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE-flatn(i,j)
                xi = flatn(i,j)

                qxm(i+1,j,QRHO)   = xi1*rho  + xi*qxm(i+1,j,QRHO)
                qxm(i+1,j,QU)     = xi1*u    + xi*qxm(i+1,j,QU)
                qxm(i+1,j,QV)     = xi1*v    + xi*qxm(i+1,j,QV)
                qxm(i+1,j,QREINT) = xi1*rhoe + xi*qxm(i+1,j,QREINT)
                qxm(i+1,j,QPRES)  = xi1*p    + xi*qxm(i+1,j,QPRES) 

             endif               

          endif

          !-------------------------------------------------------------------
          ! geometry source terms
          !-------------------------------------------------------------------

          if(dloga(i,j).ne.0)then
             courn = dtdx*(cc+abs(u))
             eta = (ONE-courn)/(cc*dt*abs(dloga(i,j)))
             dlogatmp = min(eta,ONE)*dloga(i,j)
             sourcr = -HALF*dt*rho*dlogatmp*u
             sourcp = sourcr*csq
             source = sourcp*enth

             if (i .le. ihi1) then
                qxm(i+1,j,QRHO) = qxm(i+1,j,QRHO) + sourcr
                qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
                qxm(i+1,j,QPRES) = qxm(i+1,j,QPRES) + sourcp
                qxm(i+1,j,QREINT) = qxm(i+1,j,QREINT) + source
             end if

             if (i .ge. ilo1) then
                qxp(i,j,QRHO) = qxp(i,j,QRHO) + sourcr
                qxp(i,j,QRHO) = max(qxp(i,j,QRHO),small_dens)
                qxp(i,j,QPRES) = qxp(i,j,QPRES) + sourcp
                qxp(i,j,QREINT) = qxp(i,j,QREINT) + source
             end if

          endif

       end do
    end do


    !-------------------------------------------------------------------------
    ! Now do the passively advected quantities
    !-------------------------------------------------------------------------
    
    ! We do all passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1
          
          ! plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,QU)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j)
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
                qxp(i,j,n) = q(i,j,n)    ! we might want to change this to 
                                         ! the limit of the parabola
             else if (u .lt. ZERO) then
                qxp(i,j,n) = q(i,j,n) + xi*(Im(i,j,1,2,n) - q(i,j,n))
             else
                qxp(i,j,n) = q(i,j,n) + HALF*xi*(Im(i,j,1,2,n) - q(i,j,n))
             endif
          enddo
          
          ! minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,QU)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j)
             else
                xi = ONE
             endif

             if (u .gt. ZERO) then
                qxm(i+1,j,n) = q(i,j,n) + xi*(Ip(i,j,1,2,n) - q(i,j,n))
             else if (u .lt. ZERO) then
                qxm(i+1,j,n) = q(i,j,n)
             else
                qxm(i+1,j,n) = q(i,j,n) + HALF*xi*(Ip(i,j,1,2,n) - q(i,j,n))
             endif
          enddo
          
       enddo
    enddo


    !-------------------------------------------------------------------------
    ! y-direction
    !-------------------------------------------------------------------------

    ! Trace to bottom and top edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1
          
          cc = c(i,j)
          csq = cc**2

          rho = q(i,j,QRHO)
          u = q(i,j,QU)
          v = q(i,j,QV)

          p = q(i,j,QPRES)
          rhoe = q(i,j,QREINT)
          enth = ( (rhoe+p)/rho )/csq

          Clag = rho*cc

          gam = gamc(i,j)

          game = q(i,j,QGAME)

          !------------------------------------------------------------------- 
          ! plus state on face j
          !-------------------------------------------------------------------

          ! set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. v - cc >= ZERO .and. &
                ppm_reference_edge_limit == 0)) then
             ! original Castro way -- cc value
             rho_ref  = rho
             v_ref    = v

             p_ref    = p
             rhoe_ref = rhoe

             tau_ref  = ONE/rho

             gam_ref = gamc(i,j)

             gam_ref = game
          else
             ! this will be the fastest moving state to the left
             rho_ref  = Im(i,j,2,1,QRHO)
             v_ref    = Im(i,j,2,1,QV)

             p_ref    = Im(i,j,2,1,QPRES)
             rhoe_ref = Im(i,j,2,1,QREINT)

             tau_ref  = ONE/Im(i,j,2,1,QRHO)

             gam_ref  = Im_gc(i,j,2,1,1)

             game_ref = Im(i,j,2,1,QGAME)
          endif
            
          ! for tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref
          
          ! *m are the jumps carried by v-c
          ! *p are the jumps carried by v+c

          dvm    = v_ref    - Im(i,j,2,1,QV)
          dpm    = p_ref    - Im(i,j,2,1,QPRES)
          
          drho  = rho_ref  - Im(i,j,2,2,QRHO)
          dp    = p_ref    - Im(i,j,2,2,QPRES)
          drhoe = rhoe_ref - Im(i,j,2,2,QREINT)
          dtau  = tau_ref  - ONE/Im(i,j,2,2,QRHO)

          dvp    = v_ref    - Im(i,j,2,3,QV)
          dpp    = p_ref    - Im(i,j,2,3,QPRES)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dvm = dvm - halfdt*Im_src(i,j,2,1,QV)
             dvp = dvp - halfdt*Im_src(i,j,2,3,QV)
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
             alpham = HALF*(dpm/(rho_ev*cc_ev) - dvm)*rho_ev/cc_ev
             alphap = HALF*(dpp/(rho_ev*cc_ev) + dvp)*rho_ev/cc_ev
             alpha0r = drho - dp/csq_ev
             alpha0e = drhoe - dp*enth_ev  ! enth has 1/c**2 in it
          
          else

             ! (tau, u, p, e) eigensystem
             ! or
             ! (tau, u, p, game) eigensystem

             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system.

             de = (rhoe_ref/rho_ref - Im(i,j,2,2,QREINT)/Im(i,j,2,2,QRHO))
             dge = game_ref - Im(i,j,2,2,QGAME)

             alpham = HALF*( dvm - dpm/Clag_ev)/Clag_ev
             alphap = HALF*(-dvp - dpp/Clag_ev)/Clag_ev
             alpha0r = dtau + dp/Clag_ev**2

             if (ppm_predict_gammae == 0) then
                alpha0e = de - dp*p_ev/Clag_ev**2
             else
                gfactor = (game - 1.0d0)*(game - gam)
                alpha0e = gfactor*dp/(tau_ev*Clag_ev**2) + dge
             endif
          
          endif

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
          
          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r          
          if (j .ge. ilo2) then
             if (ppm_tau_in_tracing == 0) then
                qyp(i,j,QRHO)   = rho_ref + apright + amright + azrright
                qyp(i,j,QV)     = v_ref + (apright - amright)*cc_ev/rho_ev
                qyp(i,j,QREINT) = rhoe_ref + (apright + amright)*enth_ev*csq_ev + azeright
                qyp(i,j,QPRES)  = p_ref + (apright + amright)*csq_ev

             else
                tau_s = tau_ref + apright + amright + azrright
                qyp(i,j,QRHO)   = ONE/tau_s

                qyp(i,j,QV)     = v_ref + (amright - apright)*Clag_ev

                qyp(i,j,QPRES)  = p_ref + (-apright - amright)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_ref/rho_ref + (azeright - p_ev*amright -p_ev*apright)
                   qyp(i,j,QREINT) = e_s/tau_s
                else
                   qyp(i,j,QGAME) = game_ref + gfactor*(amright + apright)/tau_ev + azeright
                   qyp(i,j,QREINT) = qyp(i,j,QPRES )/(qyp(i,j,QGAME) - 1.0d0)
                endif

             end if


             ! enforce small_*   
             qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))
             qyp(i,j,QPRES) = max(qyp(i,j,QPRES), small_pres)

             ! transverse velocity -- there is no projection here, so
             ! we don't need a reference state.  We only care about 
             ! the state traced under the middle wave
             if (v > ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qyp(i,j,QU)  = Im(i,j,2,2,QU)
                else
                   qyp(i,j,QU)  = u
                endif
             else ! wave moving toward the interface
                qyp(i,j,QU)     = Im(i,j,2,2,QU)
             endif

             if (ppm_trace_sources == 1) then
                qyp(i,j,QU)  = qyp(i,j,QU)  + halfdt*Im_src(i,j,2,2,QU)
             endif

             ! we may have already done the flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE-flatn(i,j)
                xi = flatn(i,j)

                qyp(i,j,QRHO)   = xi1*rho  + xi*qyp(i,j,QRHO)
                qyp(i,j,QV)     = xi1*v    + xi*qyp(i,j,QV)
                qyp(i,j,QU)     = xi1*u    + xi*qyp(i,j,QU)
                qyp(i,j,QREINT) = xi1*rhoe + xi*qyp(i,j,QREINT)
                qyp(i,j,QPRES)  = xi1*p    + xi*qyp(i,j,QPRES)
             endif
                
          endif

          !-------------------------------------------------------------------
          ! minus state on face j+1
          !-------------------------------------------------------------------

          ! set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. v + cc <= ZERO .and. &
                ppm_reference_edge_limit == 0) ) then
             ! original Castro way -- cc value
             rho_ref  = rho
             v_ref    = v

             p_ref    = p
             rhoe_ref = rhoe

             tau_ref  = ONE/rho

             gam_ref = gamc(i,j)

             game_ref = game

          else
             ! this will be the fastest moving state to the right
             rho_ref  = Ip(i,j,2,3,QRHO)
             v_ref    = Ip(i,j,2,3,QV)

             p_ref    = Ip(i,j,2,3,QPRES)
             rhoe_ref = Ip(i,j,2,3,QREINT)

             tau_ref  = ONE/Ip(i,j,2,3,QRHO)

             gam_ref    = Ip_gc(i,j,2,3,1)

             game_ref    = Ip(i,j,2,3,QGAME)
          endif

          ! for tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

          ! *m are the jumps carried by v-c
          ! *p are the jumps carried by v+c

          dvm    = v_ref    - Ip(i,j,2,1,QV)
          dpm    = p_ref    - Ip(i,j,2,1,QPRES)
          
          drho  = rho_ref  - Ip(i,j,2,2,QRHO)
          dp    = p_ref    - Ip(i,j,2,2,QPRES)
          drhoe = rhoe_ref - Ip(i,j,2,2,QREINT)
          dtau  = tau_ref  - ONE/Ip(i,j,2,2,QRHO)
          
          dvp    = v_ref    - Ip(i,j,2,3,QV)
          dpp    = p_ref    - Ip(i,j,2,3,QPRES)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dvm = dvm - halfdt*Ip_src(i,j,2,1,QV)
             dvp = dvp - halfdt*Ip_src(i,j,2,3,QV)
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
             alpham = HALF*(dpm/(rho_ev*cc_ev) - dvm)*rho_ev/cc_ev
             alphap = HALF*(dpp/(rho_ev*cc_ev) + dvp)*rho_ev/cc_ev
             alpha0r = drho - dp/csq_ev
             alpha0e = drhoe - dp*enth_ev

          else

             ! (tau, u, p, e) eigensystem
             ! or
             ! (tau, u, p, game) eigensystem

             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system.

             de = (rhoe_ref/rho_ref - Ip(i,j,2,2,QREINT)/Ip(i,j,2,2,QRHO))
             dge = game_ref - Ip(i,j,2,2,QGAME)

             alpham = HALF*( dvm - dpm/Clag_ev)/Clag_ev
             alphap = HALF*(-dvp - dpp/Clag_ev)/Clag_ev
             alpha0r = dtau + dp/Clag_ev**2

             if (ppm_predict_gammae == 0) then
                alpha0e = de - dp*p_ev/Clag_ev**2
             else
                gfactor = (game - 1.0d0)*(game - gam)
                alpha0e = gfactor*dp/(tau_ev*Clag_ev**2) + dge
             endif

          endif

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
          
          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r          
          if (j .le. ihi2) then
             if (ppm_tau_in_tracing == 0) then
                qym(i,j+1,QRHO)   = rho_ref + apleft + amleft + azrleft
                qym(i,j+1,QV)     = v_ref + (apleft - amleft)*cc_ev/rho_ev
                qym(i,j+1,QREINT) = rhoe_ref + (apleft + amleft)*enth_ev*csq_ev + azeleft
                qym(i,j+1,QPRES)  = p_ref + (apleft + amleft)*csq_ev
             else
                tau_s = tau_ref + apleft + amleft + azrleft
                qym(i,j+1,QRHO)   = ONE/tau_s
                qym(i,j+1,QV)     = v_ref + (amleft - apleft)*Clag_ev

                qym(i,j+1,QPRES)  = p_ref + (-apleft - amleft)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_ref/rho_ref + (azeleft - p_ev*amleft -p_ev*apleft)
                   qym(i,j+1,QREINT) = e_s/tau_s
                else
                   qym(i,j+1,QGAME) = game_ref + gfactor*(amleft + apleft)/tau_ev + azeleft 
                   qym(i,j+1,QREINT) = qym(i,j+1,QPRES )/(qym(i,j+1,QGAME) - 1.0d0)
                endif

             end if
             

             ! enforce small_*
             qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))
             qym(i,j+1,QPRES) = max(qym(i,j+1,QPRES), small_pres)

             ! transverse velocity -- there is no projection here, so
             ! we don't need a reference state.  We only care about 
             ! the state traced under the middle wave
             if (v < ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qym(i,j+1,QU) = Ip(i,j,2,2,QU)
                else
                   qym(i,j+1,QU) = u
                endif
             else
                qym(i,j+1,QU)    = Ip(i,j,2,2,QU)
             endif

             if (ppm_trace_sources == 1) then
                qym(i,j+1,QU)  = qym(i,j+1,QU)  + halfdt*Ip_src(i,j,2,2,QU)
             endif

             ! we may have already applied flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE-flatn(i,j)
                xi = flatn(i,j)

                qym(i,j+1,QRHO)   = xi1*rho  + xi*qym(i,j+1,QRHO)
                qym(i,j+1,QV)     = xi1*v    + xi*qym(i,j+1,QV)
                qym(i,j+1,QU)     = xi1*u    + xi*qym(i,j+1,QU)
                qym(i,j+1,QREINT) = xi1*rhoe + xi*qym(i,j+1,QREINT)
                qym(i,j+1,QPRES)  = xi1*p    + xi*qym(i,j+1,QPRES)
             endif
                
          endif

       end do
    end do
       
    
    !-------------------------------------------------------------------------
    ! Now do the passively advected quantities
    !-------------------------------------------------------------------------

    ! do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do i = ilo1-1, ihi1+1
          
          ! plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,QV)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j)
             else
                xi = ONE
             endif

             if (v .gt. ZERO) then
                qyp(i,j,n) = q(i,j,n)
             else if (v .lt. ZERO) then
                qyp(i,j,n) = q(i,j,n) + xi*(Im(i,j,2,2,n) - q(i,j,n))
             else
                qyp(i,j,n) = q(i,j,n) + HALF*xi*(Im(i,j,2,2,n) - q(i,j,n))
             endif
          enddo
          
          ! minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,QV)

             if (ppm_flatten_before_integrals == 0) then
                xi = flatn(i,j)
             else
                xi = ONE
             endif

             if (v .gt. ZERO) then
                qym(i,j+1,n) = q(i,j,n) + xi*(Ip(i,j,2,2,n) - q(i,j,n))
             else if (v .lt. ZERO) then
                qym(i,j+1,n) = q(i,j,n)
             else
                qym(i,j+1,n) = q(i,j,n) + HALF*xi*(Ip(i,j,2,2,n) - q(i,j,n))
             endif
          enddo
          
       enddo
    enddo
          
    
    deallocate(Ip,Im)
    if (ppm_trace_sources == 1) then
       deallocate(Ip_src,Im_src)
    endif

  end subroutine trace_ppm
  
end module trace_ppm_module
