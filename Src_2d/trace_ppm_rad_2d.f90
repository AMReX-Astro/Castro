module trace_ppm_rad_module

  implicit none

  private

  public trace_ppm_rad

contains

  subroutine trace_ppm_rad(lam, lam_l1, lam_l2, lam_h1, lam_h2, &
                           q,c,cg,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                           dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                           qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                           srcQ,src_l1,src_l2,src_h1,src_h2, &
                           ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

    use network, only : nspec
    use bl_constants_module
    use meth_params_module, only : QRHO, QU, QV, &
         QREINT, QPRES, QVAR, &
         small_dens, small_pres, &
         ppm_type, ppm_reference, ppm_trace_sources, ppm_temp_fix, &
         ppm_tau_in_tracing, ppm_reference_eigenvectors, ppm_reference_edge_limit, &
         ppm_flatten_before_integrals, ppm_predict_gammae, &
         npassive, qpass_map
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot
    use rad_params_module, only : ngroups
    use ppm_module, only : ppm

    implicit none

    integer ilo1,ilo2,ihi1,ihi2
    integer lam_l1,lam_l2,lam_h1,lam_h2
    integer qd_l1,qd_l2,qd_h1,qd_h2
    integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
    integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer src_l1,src_l2,src_h1,src_h2
    integer gc_l1,gc_l2,gc_h1,gc_h2

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision    cg(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)

    double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)

    double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)

    double precision dx, dy, dt

    ! Local variables
    integer i, j, g
    integer n, ipassive

    double precision, dimension(0:ngroups-1) :: er, der, alphar, sourcer, qrtmp,hr
    double precision, dimension(0:ngroups-1) :: lam0, lamp, lamm
    
    double precision dtdx, dtdy
    double precision cc, csq, Clag, rho, u, v, p, ptot, rhoe, enth, cgassq
    double precision drho, du, dv, drhoe, dptot
    double precision dup, dvp, dptotp
    double precision dum, dvm, dptotm

    double precision :: rho_ref, u_ref, v_ref, p_ref, rhoe_ref
    double precision :: ptot_ref, rhoe_g_ref
    double precision, dimension(0:ngroups-1) :: er_ref
    double precision :: cc_ref, csq_ref, Clag_ref, enth_ref

    double precision alpham, alphap, alpha0, alphae
    double precision alphau, alphav
    double precision sourcr,sourcp,source,courn,eta,dlogatmp

    double precision rhoe_g, h_g, alphae_g, drhoe_g

    double precision :: xi, xi1
    double precision :: halfdt

    double precision, allocatable :: Ip(:,:,:,:,:)
    double precision, allocatable :: Im(:,:,:,:,:)

    double precision, allocatable :: Ip_src(:,:,:,:,:)
    double precision, allocatable :: Im_src(:,:,:,:,:)

    double precision :: er_foo

    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
       call bl_error("Error:: ppm_2d.f90 :: trace_ppm")
    end if

    if (ppm_tau_in_tracing == 1) then
       call bl_error("ERROR: ppm_tau_in_tracing not implemented with radiation")
    endif

    if (ppm_predict_gammae == 1) then
       call bl_error("ERROR: ppm_predict_gammae not implemented with radiation")
    endif

    if (ppm_temp_fix > 0) then
       call bl_error("ERROR: ppm_temp_fix > 0 not implemented with radiation")
    endif

    if (ppm_reference_eigenvectors == 1) then
       call bl_error("ERROR: ppm_reference_eigenvectors not implemented with radiation")
    endif

    dtdx = dt/dx
    dtdy = dt/dy

    ! indices: (x, y, dimension, wave, variable)
    allocate(Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QRADVAR))
    allocate(Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QRADVAR))

    if (ppm_trace_sources == 1) then
       allocate(Ip_src(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QVAR))
       allocate(Im_src(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QVAR))
    endif

    !allocate(Ip_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,1))
    !allocate(Im_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,1))

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
    do n=1,QRADVAR
       call ppm(q(:,:,n),qd_l1,qd_l2,qd_h1,qd_h2, &
                q(:,:,QU:), c, qd_l1,qd_l2,qd_h1,qd_h2,&
                flatn, &
                Ip(:,:,:,:,n),Im(:,:,:,:,n), &
                ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    end do

    ! trace the gas gamma to the edge
    !if (ppm_temp_fix /= 1) then
    !      call ppm(gamc(:,:),gc_l1,gc_l2,gc_h1,gc_h2, &
    !               q(:,:,QU:),c,qd_l1,qd_l2,qd_h1,qd_h2, &
    !               flatn, &
    !               Ip_gc(:,:,:,:,1),Im_gc(:,:,:,:,1), &
    !               ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    !   endif

    ! if desired, do parabolic reconstruction of the source
    ! terms -- we'll use this for the force on the velocity
    if (ppm_trace_sources == 1) then
       do n = 1,QVAR
          call ppm(srcQ(:,:,n),src_l1,src_l2,src_h1,src_h2, &
                   q(:,:,QU:),c,qd_l1,qd_l2,qd_h1,qd_h2, &
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

          do g=0, ngroups-1
             lam0(g) = lam(i,j,g)
             lamp(g) = lam(i,j,g)
             lamm(g) = lam(i,j,g)
          end do

          ! cgassq is the gas soundspeed **2
          ! cc is the total soundspeed **2 (gas + radiation)
          cgassq = cg(i,j)**2
          cc = c(i,j)
          csq = cc**2

          rho = q(i,j,QRHO)
          u = q(i,j,QU)
          v = q(i,j,QV)
          p = q(i,j,QPRES)
          rhoe_g = q(i,j,QREINT)
          h_g = (p+rhoe_g) / rho

          ptot = q(i,j,qptot)
          rhoe = q(i,j,qreitot)
          enth = ( (rhoe+ptot)/rho )/csq

          er(:) = q(i,j,qrad:qradhi)
          hr(:) = (lam0+1.d0)*er/rho

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

             p_ref      = p
             rhoe_g_ref = rhoe_g

             ptot_ref = ptot
             rhoe_ref = rhoe

             er_ref(:) = er(:)
          else
             ! this will be the fastest moving state to the left --
             ! this is the method that Miller & Colella and Colella &
             ! Woodward use
             rho_ref  = Im(i,j,1,1,QRHO)
             u_ref    = Im(i,j,1,1,QU)

             p_ref    = Im(i,j,1,1,QPRES)
             rhoe_g_ref = Im(i,j,1,1,QREINT)

             ptot_ref = Im(i,j,1,1,qptot)
             rhoe_ref = Im(i,j,1,1,qreitot)

             er_ref(:) = Im(i,j,1,1,qrad:qradhi)
          endif

          ! *m are the jumps carried by u-c
          ! *p are the jumps carried by u+c

          dum    = u_ref    - Im(i,j,1,1,QU)
          dptotm = ptot_ref - Im(i,j,1,1,qptot)

          drho    = rho_ref    - Im(i,j,1,2,QRHO)
          dptot   = ptot_ref   - Im(i,j,1,2,qptot)
          drhoe   = rhoe_ref   - Im(i,j,1,2,qreitot)
          drhoe_g = rhoe_g_ref - Im(i,j,1,2,QREINT)
          der(:)  = er_ref(:)  - Im(i,j,1,2,qrad:qradhi)

          dup    = u_ref    - Im(i,j,1,3,QU)
          dptotp = ptot_ref - Im(i,j,1,3,qptot)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dum = dum - halfdt*Im_src(i,j,1,1,QU)
             dup = dup - halfdt*Im_src(i,j,1,3,QU)
          endif

          ! these are analogous to the beta's from the original
          ! PPM paper (except we work with rho instead of tau).
          ! This is simply (l . dq), where dq = qref - I(q)

          alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
          alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphar(:) = der(:) - dptot/csq*hr

          if (u-cc .gt. ZERO) then
             alpham = ZERO
          else if (u-cc .lt. ZERO) then
             alpham = -alpham
          else
             alpham = -HALF*alpham
          endif
          if (u+cc .gt. ZERO) then
             alphap = ZERO
          else if (u+cc .lt. ZERO) then
             alphap = -alphap
          else
             alphap = -HALF*alphap
          endif
          if (u .gt. ZERO) then
             alpha0 = ZERO
             alphae = ZERO
             alphae_g = ZERO
             alphar(:) = ZERO
          else if (u .lt. ZERO) then
             alpha0 = -alpha0
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else
             alpha0 = -HALF*alpha0
             alphae = -HALF*alphae
             alphae_g = -HALF*alphae_g
             alphar(:) = -HALF*alphar(:)
          endif

          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          if (i .ge. ilo1) then
             qxp(i,j,QRHO)   = rho_ref + alphap + alpham + alpha0
             qxp(i,j,QU)     = u_ref + (alphap - alpham)*cc/rho
             qxp(i,j,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g + alphae_g
             qxp(i,j,QPRES)  = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

             qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
             qxp(i,j,qrad:qradhi) = qrtmp

             qxp(i,j,qptot) = ptot_ref + (alphap + alpham)*csq
             qxp(i,j,qreitot) = qxp(i,j,QREINT) + sum(qrtmp)

             ! enforce small_*
             qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))

             do g=0, ngroups-1
                if (qxp(i,j,qrad+g) < ZERO) then
                   er_foo = - qxp(i,j,qrad+g)
                   qxp(i,j,qrad+g) = ZERO
                   qxp(i,j,qptot) = qxp(i,j,qptot) + lamp(g) * er_foo
                   qxp(i,j,qreitot) = qxp(i,j,qreitot) + er_foo
                end if
             end do

             if (qxp(i,j,QPRES) < ZERO) then
                qxp(i,j,QPRES) = p
             end if


             ! transverse velocity -- there is no projection here, so
             ! we don't need a reference state.  We only care about
             ! the state traced under the middle wave
             dv    = Im(i,j,1,2,QV)

             if (ppm_trace_sources == 1) then
                dv  = dv  + halfdt*Im_src(i,j,1,2,QV)
             endif

             ! Recall that I already takes the limit of the parabola
             ! in the event that the wave is not moving toward the
             ! interface
             if (u > ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qxp(i,j,QV)     = Im(i,j,1,2,QV)
                else
                   qxp(i,j,QV) = v
                endif
             else ! wave moving toward the interface
                qxp(i,j,QV) = dv
             endif


             ! we may have done the flattening already in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE-flatn(i,j)
                xi = flatn(i,j)

                qxp(i,j,QRHO)   = xi1*rho  + xi*qxp(i,j,QRHO)
                qxp(i,j,QU)     = xi1*u    + xi*qxp(i,j,QU)
                qxp(i,j,QREINT) = xi1*rhoe_g + xi*qxp(i,j,QREINT)
                qxp(i,j,QPRES)  = xi1*p    + xi*qxp(i,j,QPRES)
                qxp(i,j,QV)     = xi1*v    + xi*qxp(i,j,QV)

                qxp(i,j,qrad:qradhi) = xi1*er(:) + xi*qxp(i,j,qrad:qradhi)

                qxp(i,j,qptot)   = xi1*ptot + xi*qxp(i,j,qptot)
                qxp(i,j,qreitot) = xi1*rhoe + xi*qxp(i,j,qreitot)

             endif


          end if

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

             p_ref      = p
             rhoe_g_ref = rhoe_g

             ptot_ref = ptot
             rhoe_ref = rhoe

             er_ref(:) = er(:)
          else
             ! this will be the fastest moving state to the right
             rho_ref  = Ip(i,j,1,3,QRHO)
             u_ref    = Ip(i,j,1,3,QU)

             p_ref      = Ip(i,j,1,3,QPRES)
             rhoe_g_ref = Ip(i,j,1,3,QREINT)

             ptot_ref = Ip(i,j,1,3,qptot)
             rhoe_ref = Ip(i,j,1,3,qreitot)

             er_ref(:) = Ip(i,j,1,3,qrad:qradhi)
          endif

          !  *m are the jumps carried by u-c
          !  *p are the jumps carried by u+c

          dum    = u_ref    - Ip(i,j,1,1,QU)
          dptotm = ptot_ref - Ip(i,j,1,1,qptot)

          drho    = rho_ref    - Ip(i,j,1,2,QRHO)
          dptot   = ptot_ref   - Ip(i,j,1,2,qptot)
          drhoe   = rhoe_ref   - Ip(i,j,1,2,qreitot)
          drhoe_g = rhoe_g_ref - Ip(i,j,1,2,QREINT)
          der(:)  = er_ref(:)  - Ip(i,j,1,2,qrad:qradhi)

          dup    = u_ref    - Ip(i,j,1,3,QU)
          dptotp = ptot_ref - Ip(i,j,1,3,qptot)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dum = dum - halfdt*Ip_src(i,j,1,1,QU)
             dup = dup - halfdt*Ip_src(i,j,1,3,QU)
          endif

          ! these are analogous to the beta's from the original
          ! PPM paper (except we work with rho instead of tau).
          ! This is simply (l . dq), where dq = qref - I(q)
          alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
          alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphav = dv
          alphar(:) = der(:)- dptot/csq*hr

          if (u-cc .gt. ZERO) then
             alpham = -alpham
          else if (u-cc .lt. ZERO) then
             alpham = ZERO
          else
             alpham = -HALF*alpham
          endif
          if (u+cc .gt. ZERO) then
             alphap = -alphap
          else if (u+cc .lt. ZERO) then
             alphap = ZERO
          else
             alphap = -HALF*alphap
          endif
          if (u .gt. ZERO) then
             alpha0 = -alpha0
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else if (u .lt. ZERO) then
             alpha0 = ZERO
             alphae = ZERO
             alphae_g = ZERO
             alphar(:) = ZERO
          else
             alpha0 = -HALF*alpha0
             alphae = -HALF*alphae
             alphae_g = -HALF*alphae_g
             alphar(:) = -HALF*alphar(:)
          endif

          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          if (i .le. ihi1) then
             qxm(i+1,j,QRHO) = rho_ref + alphap + alpham + alpha0
             qxm(i+1,j,QU) = u_ref + (alphap - alpham)*cc/rho
             qxm(i+1,j,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g + alphae_g
             qxm(i+1,j,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))

             qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
             qxm(i+1,j,qrad:qradhi) = qrtmp

             qxm(i+1,j,qptot) = ptot_ref + (alphap + alpham)*csq
             qxm(i+1,j,qreitot) = qxm(i+1,j,QREINT) + sum(qrtmp)

             ! enforce small_*
             qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)

             do g=0, ngroups-1
                if (qxm(i+1,j,qrad+g) < ZERO) then
                   er_foo = - qxm(i+1,j,qrad+g)
                   qxm(i+1,j,qrad+g) = ZERO
                   qxm(i+1,j,qptot) = qxm(i+1,j,qptot) + lamm(g) * er_foo
                   qxm(i+1,j,qreitot) = qxm(i+1,j,qreitot) + er_foo
                end if
             end do

             if (qxm(i+1,j,QPRES) < ZERO) then
                qxm(i+1,j,QPRES) = p
             end if

             ! transverse velocity
             dv    = Ip(i,j,1,2,QV)

             if (ppm_trace_sources == 1) then
                dv  = dv  + halfdt*Ip_src(i,j,1,2,QV)
             endif

             if (u < ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qxm(i+1,j,QV) = Ip(i,j,1,2,QV)
                else
                   qxm(i+1,j,QV) = v
                endif
             else ! wave moving toward interface
                qxm(i+1,j,QV) = dv
             endif


             ! we may have already done the flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE-flatn(i,j)
                xi = flatn(i,j)

                qxm(i+1,j,QRHO)   = xi1*rho  + xi*qxm(i+1,j,QRHO)
                qxm(i+1,j,QU)     = xi1*u    + xi*qxm(i+1,j,QU)
                qxm(i+1,j,QV)     = xi1*v    + xi*qxm(i+1,j,QV)
                qxm(i+1,j,QREINT) = xi1*rhoe_g + xi*qxm(i+1,j,QREINT)
                qxm(i+1,j,QPRES)  = xi1*p    + xi*qxm(i+1,j,QPRES)

                qxm(i+1,j,qrad:qradhi) = xi1*er(:) + xi*qxm(i+1,j,qrad:qradhi)

                qxm(i+1,j,qptot)   = xi1*ptot + xi*qxm(i+1,j,qptot)
                qxm(i+1,j,qreitot) = xi1*rhoe + xi*qxm(i+1,j,qreitot)
             endif


          end if

          !-------------------------------------------------------------------
          ! geometry source terms
          !-------------------------------------------------------------------

          if(dloga(i,j).ne.0)then
             courn = dtdx*(cc+abs(u))
             eta = (ONE-courn)/(cc*dt*abs(dloga(i,j)))
             dlogatmp = min(eta,ONE)*dloga(i,j)
             sourcr = -HALF*dt*rho*dlogatmp*u
             sourcp = sourcr*cgassq
             source = sourcr*h_g
             sourcer(:) = -HALF*dt*dlogatmp*u*(lam0(:)+ONE)*er(:)
             if (i .le. ihi1) then
                qxm(i+1,j,QRHO) = qxm(i+1,j,QRHO) + sourcr
                qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
                qxm(i+1,j,QPRES) = qxm(i+1,j,QPRES) + sourcp
                qxm(i+1,j,QREINT) = qxm(i+1,j,QREINT) + source
                qxm(i+1,j,qrad:qradhi) = qxm(i+1,j,qrad:qradhi) + sourcer(:)
                !              qxm(i+1,j,qptot ) = sum(lamm(:)*qxm(i+1,j,qrad:qradhi)) + qxm(i+1,j,QPRES)
                qxm(i+1,j,qptot) = qxm(i+1,j,qptot) + sum(lamm(:)*sourcer(:)) + sourcp
                qxm(i+1,j,qreitot) = sum(qxm(i+1,j,qrad:qradhi))  + qxm(i+1,j,QREINT)
             end if
             if (i .ge. ilo1) then
                qxp(i,j,QRHO) = qxp(i,j,QRHO) + sourcr
                qxp(i,j,QRHO) = max(qxp(i,j,QRHO),small_dens)
                qxp(i,j,QPRES) = qxp(i,j,QPRES) + sourcp
                qxp(i,j,QREINT) = qxp(i,j,QREINT) + source
                qxp(i,j,qrad:qradhi) = qxp(i,j,qrad:qradhi) + sourcer(:)
                !              qxp(i,j,qptot ) = sum(lamp(:)*qxp(i,j,qrad:qradhi)) + qxp(i,j,QPRES)
                qxp(i,j,qptot) = qxp(i,j,qptot) + sum(lamp(:)*sourcer(:)) + sourcp
                qxp(i,j,qreitot) = sum(qxp(i,j,qrad:qradhi))  + qxp(i,j,QREINT)
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

          do g=0, ngroups-1
             lam0(g) = lam(i,j,g)
             lamp(g) = lam(i,j,g)
             lamm(g) = lam(i,j,g)
          end do

          ! cgassq is the gas soundspeed **2
          ! cc is the total soundspeed **2 (gas + radiation)
          cgassq = cg(i,j)**2
          cc = c(i,j)
          csq = cc**2

          rho = q(i,j,QRHO)
          u = q(i,j,QU)
          v = q(i,j,QV)
          p = q(i,j,QPRES)
          rhoe_g = q(i,j,QREINT)
          h_g = (p+rhoe_g) / rho
          ptot = q(i,j,qptot)
          rhoe = q(i,j,qreitot)
          enth = ( (rhoe+ptot)/rho )/csq

          er(:) = q(i,j,qrad:qradhi)
          hr(:) = (lam0+ONE)*er/rho

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

             p_ref      = p
             rhoe_g_ref = rhoe_g

             ptot_ref = ptot
             rhoe_ref = rhoe

             er_ref(:) = er(:)
          else
             ! this will be the fastest moving state to the left
             rho_ref  = Im(i,j,2,1,QRHO)
             v_ref    = Im(i,j,2,1,QV)

             p_ref = Im(i,j,2,1,QPRES)
             rhoe_g_ref = Im(i,j,2,1,QREINT)

             ptot_ref = Im(i,j,2,1,qptot)
             rhoe_ref = Im(i,j,2,1,qreitot)

             er_ref(:) = Im(i,j,2,1,qrad:qradhi)
          endif

          ! *m are the jumps carried by v-c
          ! *p are the jumps carried by v+c
          
          dvm    = v_ref    - Im(i,j,2,1,QV)
          dptotm = ptot_ref - Im(i,j,2,1,qptot)

          drho    = rho_ref    - Im(i,j,2,2,QRHO)
          dptot   = ptot_ref   - Im(i,j,2,2,qptot)
          drhoe   = rhoe_ref   - Im(i,j,2,2,qreitot)
          drhoe_g = rhoe_g_ref - Im(i,j,2,2,QREINT)
          der(:)  = er_ref(:)  - Im(i,j,2,2,qrad:qradhi)

          dvp    = v_ref    - Im(i,j,2,3,QV)
          dptotp = ptot_ref - Im(i,j,2,3,qptot)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dvm = dvm - halfdt*Im_src(i,j,2,1,QV)
             dvp = dvp - halfdt*Im_src(i,j,2,3,QV)
          endif

          ! these are analogous to the beta's from the original PPM
          ! paper (except we work with rho instead of tau).  This
          ! is simply (l . dq), where dq = qref - I(q)
          alpham = HALF*(dptotm/(rho*cc) - dvm)*rho/cc
          alphap = HALF*(dptotp/(rho*cc) + dvp)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphar(:) = der(:) - dptot/csq*hr

          if (v-cc .gt. ZERO) then
             alpham = ZERO
          else if (v-cc .lt. ZERO) then
             alpham = -alpham
          else
             alpham = -HALF*alpham
          endif
          if (v+cc .gt. ZERO) then
             alphap = ZERO
          else if (v+cc .lt. ZERO) then
             alphap = -alphap
          else
             alphap = -HALF*alphap
          endif
          if (v .gt. ZERO) then
             alpha0 = ZERO
             alphae = ZERO
             alphae_g = ZERO
             alphar(:) = ZERO
          else if (v .lt. ZERO) then
             alpha0 = -alpha0
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else
             alpha0 = -HALF*alpha0
             alphae = -HALF*alphae
             alphae_g = -HALF*alphae_g
             alphar(:) = -HALF*alphar(:)
          endif

          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          if (j .ge. ilo2) then
             qyp(i,j,QRHO) = rho_ref + alphap + alpham + alpha0
             qyp(i,j,QV) = v_ref + (alphap - alpham)*cc/rho
             qyp(i,j,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g + alphae_g
             qyp(i,j,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

             qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
             qyp(i,j,qrad:qradhi) = qrtmp

             qyp(i,j,qptot) = ptot_ref + (alphap + alpham)*csq
             qyp(i,j,qreitot) = qyp(i,j,QREINT) + sum(qrtmp)

             ! enforce small_*
             qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))

             do g=0, ngroups-1
                if (qyp(i,j,qrad+g) < ZERO) then
                   er_foo = - qyp(i,j,qrad+g)
                   qyp(i,j,qrad+g) = ZERO
                   qyp(i,j,qptot) = qyp(i,j,qptot) + lamp(g) * er_foo
                   qyp(i,j,qreitot) = qyp(i,j,qreitot) + er_foo
                end if
             end do

             if (qyp(i,j,QPRES) < ZERO) then
                qyp(i,j,QPRES) = p
             end if

             ! transverse velocity
             du    = Im(i,j,2,2,QU)

             if (ppm_trace_sources == 1) then
                du  = du  + halfdt*Im_src(i,j,2,2,QU)
             endif

             if (v > ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qyp(i,j,QU)     = Im(i,j,2,2,QU)
                else
                   qyp(i,j,QU)     = u
                endif
             else ! wave moving toward the interface
                qyp(i,j,QU)     = du
             endif


             ! we may have already done the flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE-flatn(i,j)
                xi = flatn(i,j)

                qyp(i,j,QRHO)   = xi1*rho  + xi*qyp(i,j,QRHO)
                qyp(i,j,QV)     = xi1*v    + xi*qyp(i,j,QV)
                qyp(i,j,QU)     = xi1*u    + xi*qyp(i,j,QU)
                qyp(i,j,QREINT) = xi1*rhoe_g + xi*qyp(i,j,QREINT)
                qyp(i,j,QPRES)  = xi1*p    + xi*qyp(i,j,QPRES)

                qyp(i,j,qrad:qradhi) = xi1*er(:) + xi*qyp(i,j,qrad:qradhi)

                qyp(i,j,qptot)   = xi1*ptot + xi*qyp(i,j,qptot)
                qyp(i,j,qreitot) = xi1*rhoe + xi*qyp(i,j,qreitot)

             endif

          end if

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

             p_ref      = p
             rhoe_g_ref = rhoe_g

             ptot_ref = ptot
             rhoe_ref = rhoe

             er_ref(:) = er(:)
          else
             ! this will be the fastest moving state to the right
             rho_ref = Ip(i,j,2,3,QRHO)
             v_ref   = Ip(i,j,2,3,QV)

             rhoe_g_ref = Ip(i,j,2,3,QREINT)
             p_ref      = Ip(i,j,2,3,QPRES)             

             ptot_ref = Ip(i,j,2,3,qptot)
             rhoe_ref = Ip(i,j,2,3,qreitot)

             er_ref(:) = Ip(i,j,2,3,qrad:qradhi)
          endif

          ! *m are the jumps carried by v-c
          ! *p are the jumps carried by v+c

          dvm    = v_ref    - Ip(i,j,2,1,QV)
          dptotm = ptot_ref - Ip(i,j,2,1,qptot)

          drho    = rho_ref    - Ip(i,j,2,2,QRHO)
          dptot   = ptot_ref   - Ip(i,j,2,2,qptot)
          drhoe   = rhoe_ref   - Ip(i,j,2,2,qreitot)
          drhoe_g = rhoe_g_ref - Ip(i,j,2,2,QREINT)
          der(:)  = er_ref(:)  - Ip(i,j,2,2,qrad:qradhi)

          dvp    = v_ref    - Ip(i,j,2,3,QV)
          dptotp = ptot_ref - Ip(i,j,2,3,qptot)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dvm = dvm - halfdt*Ip_src(i,j,2,1,QV)
             dvp = dvp - halfdt*Ip_src(i,j,2,3,QV)
          endif

          ! these are analogous to the beta's from the original PPM
          ! paper.  This is simply (l . dq), where dq = qref - I(q)
          alpham = HALF*(dptotm/(rho*cc) - dvm)*rho/cc
          alphap = HALF*(dptotp/(rho*cc) + dvp)*rho/cc
          alpha0 = drho - dptot/csq
          alphae = drhoe - dptot*enth
          alphae_g = drhoe_g - dptot/csq*h_g
          alphar(:) = der(:)- dptot/csq*hr

          if (v-cc .gt. ZERO) then
             alpham = -alpham
          else if (v-cc .lt. ZERO) then
             alpham = ZERO
          else
             alpham = -HALF*alpham
          endif
          if (v+cc .gt. ZERO) then
             alphap = -alphap
          else if (v+cc .lt. ZERO) then
             alphap = ZERO
          else
             alphap = -HALF*alphap
          endif
          if (v .gt. ZERO) then
             alpha0 = -alpha0
             alphae = -alphae
             alphae_g = -alphae_g
             alphar(:) = -alphar(:)
          else if (v .lt. ZERO) then
             alpha0 = ZERO
             alphae = ZERO
             alphae_g = ZERO
             alphar(:) = ZERO
          else
             alpha0 = -HALF*alpha0
             alphae = -HALF*alphae
             alphae_g = -HALF*alphae_g
             alphar(:) = -HALF*alphar(:)
          endif

          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          if (j .le. ihi2) then
             qym(i,j+1,QRHO) = rho_ref + alphap + alpham + alpha0
             qym(i,j+1,QV) = v_ref + (alphap - alpham)*cc/rho
             qym(i,j+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g + alphae_g
             qym(i,j+1,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))

             qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
             qym(i,j+1,qrad:qradhi) = qrtmp

             qym(i,j+1,qptot) = ptot_ref + (alphap + alpham)*csq
             qym(i,j+1,qreitot) = qym(i,j+1,QREINT) + sum(qrtmp)

             ! enforce small_*
             qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))

             do g=0, ngroups-1
                if (qym(i,j+1,qrad+g) < ZERO) then
                   er_foo = - qym(i,j+1,qrad+g)
                   qym(i,j+1,qrad+g) = ZERO
                   qym(i,j+1,qptot) = qym(i,j+1,qptot) + lamm(g) * er_foo
                   qym(i,j+1,qreitot) = qym(i,j+1,qreitot) + er_foo
                end if
             end do

             if (qym(i,j+1,QPRES) < ZERO) then
                qym(i,j+1,QPRES) = p
             end if

             ! transverse velocity
             du    =  Ip(i,j,2,2,QU)

             if (ppm_trace_sources == 1) then
                du  = du  + halfdt*Ip_src(i,j,2,2,QU)
             endif

             if (v < ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qym(i,j+1,QU) = Ip(i,j,2,2,QU)
                else
                   qym(i,j+1,QU) = u
                endif
             else
                qym(i,j+1,QU)     = du
             endif


             ! we may have already applied flattening in the parabola
             if (ppm_flatten_before_integrals == 0) then
                xi1 = ONE-flatn(i,j)
                xi = flatn(i,j)

                qym(i,j+1,QRHO)   = xi1*rho  + xi*qym(i,j+1,QRHO)
                qym(i,j+1,QV)     = xi1*v    + xi*qym(i,j+1,QV)
                qym(i,j+1,QU)     = xi1*u    + xi*qym(i,j+1,QU)
                qym(i,j+1,QREINT) = xi1*rhoe_g + xi*qym(i,j+1,QREINT)
                qym(i,j+1,QPRES)  = xi1*p    + xi*qym(i,j+1,QPRES)

                qym(i,j+1,qrad:qradhi) = xi1*er(:) + xi*qym(i,j+1,qrad:qradhi)

                qym(i,j+1,qptot)   = xi1*ptot + xi*qym(i,j+1,qptot)
                qym(i,j+1,qreitot) = xi1*rhoe + xi*qym(i,j+1,qreitot)

             endif

          end if

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

  end subroutine trace_ppm_rad

end module trace_ppm_rad_module
