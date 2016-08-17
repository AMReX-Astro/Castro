module trace_ppm_rad_module

  implicit none

  private

  public tracexy_ppm_rad, tracez_ppm_rad

contains

  subroutine tracexy_ppm_rad(lam, lam_lo, lam_hi, &
                             q, c, cg, flatn, qd_lo, qd_hi, &
                             Ip, Im, Ip_src, Im_src, &
                             qxm, qxp, qym, qyp, qs_lo, qs_hi, &
                             ilo1,ilo2,ihi1,ihi2,dt,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QPRES, &
         small_dens, small_pres, &
         ppm_type, ppm_reference, ppm_trace_sources, &
         ppm_reference_edge_limit, &
         npassive, qpass_map
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot
    use rad_params_module, only : ngroups
    use bl_constants_module

    implicit none

    integer, intent(in) :: lam_lo(3), lam_hi(3)
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qs_lo(3), qs_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: kc,k3d

    double precision, intent(in) :: lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)

    double precision, intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision, intent(in) ::     c(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision, intent(in) ::    cg(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision, intent(in) :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))

    double precision, intent(in) :: Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)
    double precision, intent(in) :: Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)

    double precision, intent(in) :: Ip_src(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)
    double precision, intent(in) :: Im_src(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)


    double precision, intent(inout) :: qxm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QRADVAR)
    double precision, intent(inout) :: qxp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QRADVAR)
    double precision, intent(inout) :: qym(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QRADVAR)
    double precision, intent(inout) :: qyp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QRADVAR)


    double precision, intent(in) :: dt

    ! Local variables
    integer :: i, j, g
    integer :: n, ipassive

    double precision :: hdt

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is 
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! rhoe : total specific internal energy (including radiation,
    !        if available)
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! htot : total specific enthalpy
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    double precision :: cc, csq, cgassq, Clag
    double precision :: rho, u, v, w, p, rhoe_g, h_g
    double precision :: ptot, rhoe, htot

    double precision :: drho, dptot, drhoe_g
    double precision :: drhoe
    double precision :: dup, dvp, dptotp
    double precision :: dum, dvm, dptotm

    double precision :: rho_ref, u_ref, v_ref, p_ref, rhoe_g_ref
    double precision :: tau_ref
    double precision :: ptot_ref, rhoe_ref, htot_ref


    double precision :: alpham, alphap, alpha0r, alpha0e_g, alpha0e


    double precision, dimension(0:ngroups-1) :: er, der, alphar, qrtmp,hr
    double precision, dimension(0:ngroups-1) :: lam0, lamp, lamm

    double precision, dimension(0:ngroups-1) :: er_ref


    double precision :: er_foo

    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call bl_error("Error:: RadHydro_3d.f90 :: tracexy_ppm_rad")
    end if

    hdt = HALF * dt


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

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          do g=0, ngroups-1
             lam0(g) = lam(i,j,k3d,g)
             lamp(g) = lam(i,j,k3d,g)
             lamm(g) = lam(i,j,k3d,g)
          end do

          rho = q(i,j,k3d,QRHO)

          ! cgassq is the gas soundspeed **2
          ! cc is the total soundspeed **2 (gas + radiation)
          cgassq = cg(i,j,k3d)**2
          cc = c(i,j,k3d)
          csq = cc**2

          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)

          p = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          h_g = ( (p+rhoe_g)/rho)/csq

          ptot = q(i,j,k3d,qptot)
          rhoe = q(i,j,k3d,qreitot)
          htot = ( (rhoe+ptot)/rho )/csq

          er(:) = q(i,j,k3d,qrad:qradhi)
          hr(:) = (lam0+ONE)*er/rho

          !-------------------------------------------------------------------
          ! plus state on face i
          !-------------------------------------------------------------------

          if (i .ge. ilo1) then

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. u - cc >= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                u_ref    = u

                p_ref    = p
                rhoe_g_ref = rhoe_g

                ptot_ref = ptot
                rhoe_ref = rhoe

                er_ref(:) = er(:)
             else
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(i,j,kc,1,1,QRHO)
                u_ref    = Im(i,j,kc,1,1,QU)

                p_ref    = Im(i,j,kc,1,1,QPRES)
                rhoe_g_ref = Im(i,j,kc,1,1,QREINT)

                ptot_ref = Im(i,j,kc,1,1,QPTOT)
                rhoe_ref = Im(i,j,kc,1,1,QREITOT)

                er_ref(:) = Im(i,j,kc,1,1,QRAD:QRADHI)
             endif

             rho_ref = max(rho_ref,small_dens)
             p_ref = max(p_ref,small_pres)

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c

             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)

             dum    = u_ref    - Im(i,j,kc,1,1,QU)
             dptotm = ptot_ref - Im(i,j,kc,1,1,qptot)

             drho    = rho_ref    - Im(i,j,kc,1,2,QRHO)
             dptot   = ptot_ref   - Im(i,j,kc,1,2,qptot)
             drhoe   = rhoe_ref   - Im(i,j,kc,1,2,qreitot)
             drhoe_g = rhoe_g_ref - Im(i,j,kc,1,2,QREINT)
             der(:)  = er_ref(:)  - Im(i,j,kc,1,2,qrad:qradhi)

             dup    = u_ref    - Im(i,j,kc,1,3,QU)
             dptotp = ptot_ref - Im(i,j,kc,1,3,qptot)

             ! If we are doing source term tracing, then add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_sources .eq. 1) then
                dum = dum - hdt*Im_src(i,j,kc,1,1,QU)
                dup = dup - hdt*Im_src(i,j,kc,1,3,QU)
             endif


             ! Optionally use the reference state in evaluating the
             ! eigenvectors -- NOT YET IMPLEMENTED


             alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
             alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dptot/csq
             alpha0e = drhoe - dptot*htot
             alpha0e_g = drhoe_g - dptot*h_g
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
                alpha0r = ZERO
                alpha0e = ZERO
                alpha0e_g = ZERO
                alphar(:) = ZERO
             else if (u .lt. ZERO) then
                alpha0r = -alpha0r
                alpha0e = -alpha0e
                alpha0e_g = -alpha0e_g
                alphar(:) = -alphar(:)
             else
                alpha0r = -HALF*alpha0r
                alpha0e = -HALF*alpha0e
                alpha0e_g = -HALF*alpha0e_g
                alphar(:) = -HALF*alphar(:)
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum(l . dq) r
             ! note that the a{mpz}right as defined above have the minus already

             qxp(i,j,kc,QRHO) = rho_ref + alphap + alpham + alpha0r
             qxp(i,j,kc,QU) = u_ref + (alphap - alpham)*cc/rho
             qxp(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
             qxp(i,j,kc,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

             qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
             qxp(i,j,kc,qrad:qradhi) = qrtmp

             qxp(i,j,kc,qptot) = ptot_ref + (alphap + alpham)*csq
             qxp(i,j,kc,qreitot) = qxp(i,j,kc,QREINT) + sum(qrtmp)

             ! Enforce small_*
             qxp(i,j,kc,QRHO) = max(qxp(i,j,kc,QRHO), small_dens)
             qxp(i,j,kc,QPRES) = max(qxp(i,j,kc,QPRES),small_pres)

             do g = 0, ngroups-1
                if (qxp(i,j,kc,qrad+g) < ZERO) then
                   er_foo = - qxp(i,j,kc,qrad+g)
                   qxp(i,j,kc,qrad+g) = ZERO
                   qxp(i,j,kc,qptot) = qxp(i,j,kc,qptot) + lamp(g) * er_foo
                   qxp(i,j,kc,qreitot) = qxp(i,j,kc,qreitot) + er_foo
                end if
             end do

             if (qxp(i,j,kc,QPRES) < ZERO) then
                qxp(i,j,kc,QPRES) = p
             end if

             ! Transverse velocities -- there's no projection here, so
             ! we don't need a reference state.  We only care about
             ! the state traced under the middle wave

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
                qxp(i,j,kc,QV) = Im(i,j,kc,1,2,QV)
                qxp(i,j,kc,QW) = Im(i,j,kc,1,2,QW)
             endif

             if (ppm_trace_sources .eq. 1) then
                qxp(i,j,kc,QV) = qxp(i,j,kc,QV) + hdt*Im_src(i,j,kc,1,2,QV)
                qxp(i,j,kc,QW) = qxp(i,j,kc,QW) + hdt*Im_src(i,j,kc,1,2,QW)
             endif

          endif


          !-------------------------------------------------------------------
          ! minus state on face i + 1
          !-------------------------------------------------------------------
          if (i .le. ihi1) then

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. u + cc <= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc values
                rho_ref  = rho
                u_ref    = u

                p_ref    = p
                rhoe_g_ref = rhoe_g

                ptot_ref = ptot
                rhoe_ref = rhoe

                er_ref(:) = er(:)
             else
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,kc,1,3,QRHO)
                u_ref    = Ip(i,j,kc,1,3,QU)

                p_ref    = Ip(i,j,kc,1,3,QPRES)
                rhoe_g_ref = Ip(i,j,kc,1,3,QREINT)

                ptot_ref = Ip(i,j,kc,1,3,QPTOT)
                rhoe_ref = Ip(i,j,kc,1,3,QREITOT)

                er_ref(:) = Ip(i,j,kc,1,3,QRAD:QRADHI)
             endif

             rho_ref = max(rho_ref,small_dens)
             p_ref = max(p_ref,small_pres)

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c

             dum    = u_ref    - Ip(i,j,kc,1,1,QU)
             dptotm = ptot_ref - Ip(i,j,kc,1,1,qptot)

             drho    = rho_ref    - Ip(i,j,kc,1,2,QRHO)
             dptot   = ptot_ref   - Ip(i,j,kc,1,2,qptot)
             drhoe   = rhoe_ref   - Ip(i,j,kc,1,2,qreitot)
             drhoe_g = rhoe_g_ref - Ip(i,j,kc,1,2,QREINT)
             der(:)  = er_ref(:)  - Ip(i,j,kc,1,2,qrad:qradhi)

             dup    = u_ref    - Ip(i,j,kc,1,3,QU)
             dptotp = ptot_ref - Ip(i,j,kc,1,3,qptot)

             ! If we are doing source term tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_sources .eq. 1) then
                dum = dum - hdt*Ip_src(i,j,kc,1,1,QU)
                dup = dup - hdt*Ip_src(i,j,kc,1,3,QU)
             endif


             ! Optionally use the reference state in evaluating the
             ! eigenvectors -- NOT YET IMPLEMENTED

             alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
             alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dptot/csq
             alpha0e = drhoe - dptot*htot
             alpha0e_g = drhoe_g - dptot*h_g
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
                alpha0r = -alpha0r
                alpha0e = -alpha0e
                alpha0e_g = -alpha0e_g
                alphar(:) = -alphar(:)
             else if (u .lt. ZERO) then
                alpha0r = ZERO
                alpha0e = ZERO
                alpha0e_g = ZERO
                alphar(:) = ZERO
             else
                alpha0r = -HALF*alpha0r
                alpha0e = -HALF*alpha0e
                alpha0e_g = -HALF*alpha0e_g
                alphar(:) = -HALF*alphar(:)
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}left as defined above have the minus already
             qxm(i+1,j,kc,QRHO) = rho_ref + alphap + alpham + alpha0r
             qxm(i+1,j,kc,QU) = u_ref + (alphap - alpham)*cc/rho
             qxm(i+1,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
             qxm(i+1,j,kc,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))

             qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
             qxm(i+1,j,kc,qrad:qradhi) = qrtmp

             qxm(i+1,j,kc,qptot) = ptot_ref + (alphap + alpham)*csq
             qxm(i+1,j,kc,qreitot) = qxm(i+1,j,kc,QREINT) + sum(qrtmp)

             ! Enforce small_*
             qxm(i+1,j,kc,QRHO)  = max(qxm(i+1,j,kc,QRHO),small_dens)
             qxm(i+1,j,kc,QPRES) = max(qxm(i+1,j,kc,QPRES),small_pres)

             do g=0,ngroups-1
                if (qxm(i+1,j,kc,qrad+g) < ZERO) then
                   er_foo = - qxm(i+1,j,kc,qrad+g)
                   qxm(i+1,j,kc,qrad+g) = ZERO
                   qxm(i+1,j,kc,qptot) = qxm(i+1,j,kc,qptot) + lamm(g) * er_foo
                   qxm(i+1,j,kc,qreitot) = qxm(i+1,j,kc,qreitot) + er_foo
                end if
             end do

             if (qxm(i+1,j,kc,QPRES) < ZERO) then
                qxm(i+1,j,kc,QPRES) = p
             end if

             ! transverse velocities
             if (u < ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qxm(i+1,j,kc,QV    ) = Ip(i,j,kc,1,2,QV)
                   qxm(i+1,j,kc,QW    ) = Ip(i,j,kc,1,2,QW)
                else
                   qxm(i+1,j,kc,QV    ) = v
                   qxm(i+1,j,kc,QW    ) = w
                endif
             else ! wave moving toward interface
                qxm(i+1,j,kc,QV    ) = Ip(i,j,kc,1,2,QV)
                qxm(i+1,j,kc,QW    ) = Ip(i,j,kc,1,2,QW)
             endif

             if (ppm_trace_sources .eq. 1) then
                qxm(i+1,j,kc,QV) = qxm(i+1,j,kc,QV) + hdt*Ip_src(i,j,kc,1,2,QV)
                qxm(i+1,j,kc,QW) = qxm(i+1,j,kc,QW) + hdt*Ip_src(i,j,kc,1,2,QW)
             endif

          end if

       end do
    end do

    !-------------------------------------------------------------------------
    ! passively advected quantities
    !-------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1

          ! Plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,QU)

             ! We have
             !
             ! q_l = q_ref - Proj{(q_ref - I)}
             !
             ! and Proj{} represents the characteristic projection.
             ! But for these, there is only 1-wave that matters, the u
             ! wave, so no projection is needed.  Since we are not
             ! projecting, the reference state doesn't matter

             if (u .gt. ZERO) then
                qxp(i,j,kc,n) = q(i,j,k3d,n)
             else if (u .lt. ZERO) then
                qxp(i,j,kc,n) = Im(i,j,kc,1,2,n)
             else
                qxp(i,j,kc,n) = q(i,j,k3d,n) + HALF*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo

          ! Minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)

             if (u .gt. ZERO) then
                qxm(i+1,j,kc,n) = Ip(i,j,kc,1,2,n)
             else if (u .lt. ZERO) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n)
             else
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) + HALF*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
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
             lam0(g) = lam(i,j,k3d,g)
             lamp(g) = lam(i,j,k3d,g)
             lamm(g) = lam(i,j,k3d,g)
          end do

          rho = q(i,j,k3d,QRHO)

          cgassq = cg(i,j,k3d)**2
          cc = c(i,j,k3d)
          csq = cc**2

          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)

          p = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          h_g = ( (p + rhoe_g)/rho)/csq

          ptot = q(i,j,k3d,qptot)
          rhoe = q(i,j,k3d,qreitot)
          htot = ( (rhoe+ptot)/rho )/csq

          er(:)= q(i,j,k3d,qrad:qradhi)
          hr(:) = (lam0+ONE)*er/rho

          !-------------------------------------------------------------------
          ! plus state on face j
          !-------------------------------------------------------------------

          if (j .ge. ilo2) then

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. v - cc >= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                v_ref    = v

                p_ref    = p
                rhoe_g_ref = rhoe_g

                ptot_ref    = ptot
                rhoe_ref = rhoe

                er_ref(:) = er(:)
             else
                ! This will be the fastest moving state to the left
                rho_ref  = Im(i,j,kc,2,1,QRHO)
                v_ref    = Im(i,j,kc,2,1,QV)

                p_ref    = Im(i,j,kc,2,1,QPRES)
                rhoe_g_ref = Im(i,j,kc,2,1,QREINT)

                ptot_ref    = Im(i,j,kc,2,1,QPTOT)
                rhoe_ref = Im(i,j,kc,2,1,QREITOT)

                er_ref(:) = Im(i,j,kc,2,1,QRAD:QRADHI)
             endif

             rho_ref = max(rho_ref,small_dens)
             p_ref = max(p_ref,small_pres)

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c

             dvm    = v_ref    - Im(i,j,kc,2,1,QV)
             dptotm = ptot_ref - Im(i,j,kc,2,1,qptot)

             drho    = rho_ref    - Im(i,j,kc,2,2,QRHO)
             dptot   = ptot_ref   - Im(i,j,kc,2,2,qptot)
             drhoe   = rhoe_ref   - Im(i,j,kc,2,2,qreitot)
             drhoe_g = rhoe_g_ref - Im(i,j,kc,2,2,QREINT)
             der(:)  = er_ref(:)  - Im(i,j,kc,2,2,qrad:qradhi)

             dvp    = v_ref    - Im(i,j,kc,2,3,QV)
             dptotp = ptot_ref - Im(i,j,kc,2,3,qptot)

             ! If we are doing source term tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_sources .eq. 1) then
                dvm = dvm - hdt*Im_src(i,j,kc,2,1,QV)
                dvp = dvp - hdt*Im_src(i,j,kc,2,3,QV)
             endif

             ! Optionally use the reference state in evaluating the
             ! eigenvectors -- NOT YET IMPLEMENTED

             alpham = HALF*(dptotm/(rho*cc) - dvm)*rho/cc
             alphap = HALF*(dptotp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dptot/csq
             alpha0e = drhoe - dptot*htot
             alpha0e_g = drhoe_g - dptot*h_g
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
                alpha0r = ZERO
                alpha0e = ZERO
                alpha0e_g = ZERO
                alphar(:) = ZERO
             else if (v .lt. ZERO) then
                alpha0r = -alpha0r
                alpha0e = -alpha0e
                alpha0e_g = -alpha0e_g
                alphar(:) = -alphar(:)
             else
                alpha0r = -HALF*alpha0r
                alpha0e = -HALF*alpha0e
                alpha0e_g = -HALF*alpha0e_g
                alphar(:) = -HALF*alphar(:)
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}right as defined above have the minus already
             qyp(i,j,kc,QRHO) = rho_ref + alphap + alpham + alpha0r
             qyp(i,j,kc,QV) = v_ref + (alphap - alpham)*cc/rho
             qyp(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
             qyp(i,j,kc,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

             qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
             qyp(i,j,kc,qrad:qradhi) = qrtmp

             qyp(i,j,kc,qptot) = ptot_ref + (alphap + alpham)*csq
             qyp(i,j,kc,qreitot) = qyp(i,j,kc,QREINT) + sum(qrtmp)

             ! Enforce small_*
             qyp(i,j,kc,QRHO ) = max(qyp(i,j,kc,QRHO ),small_dens)
             qyp(i,j,kc,QPRES) = max(qyp(i,j,kc,QPRES),small_pres)

             do g=0,ngroups-1
                if (qyp(i,j,kc,qrad+g) < ZERO) then
                   er_foo = - qyp(i,j,kc,qrad+g)
                   qyp(i,j,kc,qrad+g) = ZERO
                   qyp(i,j,kc,qptot) = qyp(i,j,kc,qptot) + lamp(g) * er_foo
                   qyp(i,j,kc,qreitot) = qyp(i,j,kc,qreitot) + er_foo
                end if
             end do

             if (qyp(i,j,kc,QPRES) < ZERO) then
                qyp(i,j,kc,QPRES) = p
             end if

             ! transverse velocities
             if (v > ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qyp(i,j,kc,QU    ) = Im(i,j,kc,2,2,QU)
                   qyp(i,j,kc,QW    ) = Im(i,j,kc,2,2,QW)
                else
                   qyp(i,j,kc,QU    ) = u
                   qyp(i,j,kc,QW    ) = w
                endif
             else ! wave moving toward the interface
                qyp(i,j,kc,QU    ) = Im(i,j,kc,2,2,QU)
                qyp(i,j,kc,QW    ) = Im(i,j,kc,2,2,QW)
             endif

             if (ppm_trace_sources .eq. 1) then
                qyp(i,j,kc,QU) = qyp(i,j,kc,QU) + hdt*Im_src(i,j,kc,2,2,QU)
                qyp(i,j,kc,QW) = qyp(i,j,kc,QW) + hdt*Im_src(i,j,kc,2,2,QW)
             endif

          end if

          !-------------------------------------------------------------------
          ! minus state on face j+1
          !-------------------------------------------------------------------

          if (j .le. ihi2) then

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. v + cc <= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value

                rho_ref  = rho
                v_ref    = v

                p_ref    = p
                rhoe_g_ref = rhoe_g

                ptot_ref = ptot
                rhoe_ref = rhoe

                er_ref(:) = er(:)
             else
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,kc,2,3,QRHO)
                v_ref    = Ip(i,j,kc,2,3,QV)

                p_ref    = Ip(i,j,kc,2,3,QPRES)
                rhoe_g_ref = Ip(i,j,kc,2,3,QREINT)

                ptot_ref    = Ip(i,j,kc,2,3,QPTOT)
                rhoe_ref = Ip(i,j,kc,2,3,QREITOT)

                er_ref(:) = Ip(i,j,kc,2,3,QRAD:QRADHI)
             endif

             rho_ref = max(rho_ref,small_dens)
             p_ref = max(p_ref,small_pres)

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c

             dvm    = v_ref    - Ip(i,j,kc,2,1,QV)
             dptotm = ptot_ref - Ip(i,j,kc,2,1,qptot)

             drho    = rho_ref    - Ip(i,j,kc,2,2,QRHO)
             dptot   = ptot_ref   - Ip(i,j,kc,2,2,qptot)
             drhoe   = rhoe_ref   - Ip(i,j,kc,2,2,qreitot)
             drhoe_g = rhoe_g_ref - Ip(i,j,kc,2,2,QREINT)
             der(:)  = er_ref(:)  - Ip(i,j,kc,2,2,qrad:qradhi)

             dvp    = v_ref    - Ip(i,j,kc,2,3,QV)
             dptotp = ptot_ref - Ip(i,j,kc,2,3,qptot)

             ! If we are doing source term tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_sources .eq. 1) then
                dvm = dvm - hdt*Ip_src(i,j,kc,2,1,QV)
                dvp = dvp - hdt*Ip_src(i,j,kc,2,3,QV)
             endif


             ! Optionally use the reference state in evaluating the
             ! eigenvectors -- NOT YET IMPLEMENTED

             alpham = HALF*(dptotm/(rho*cc) - dvm)*rho/cc
             alphap = HALF*(dptotp/(rho*cc) + dvp)*rho/cc
             alpha0e_g = drhoe_g - dptot*h_g
             alpha0r = drho - dptot/csq
             alpha0e = drhoe - dptot*htot
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
                alpha0r = -alpha0r
                alpha0e = -alpha0e
                alpha0e_g = -alpha0e_g
                alphar(:) = -alphar(:)
             else if (v .lt. ZERO) then
                alpha0r = ZERO
                alpha0e = ZERO
                alpha0e_g = ZERO
                alphar(:) = ZERO
             else
                alpha0r = -HALF*alpha0r
                alpha0e = -HALF*alpha0e
                alpha0e_g = -HALF*alpha0e_g
                alphar(:) = -HALF*alphar(:)
             endif

             qym(i,j+1,kc,QRHO) = rho_ref + alphap + alpham + alpha0r
             qym(i,j+1,kc,QV) = v_ref + (alphap - alpham)*cc/rho
             qym(i,j+1,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
             qym(i,j+1,kc,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))

             qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
             qym(i,j+1,kc,qrad:qradhi) = qrtmp

             qym(i,j+1,kc,qptot) = ptot_ref + (alphap + alpham)*csq
             qym(i,j+1,kc,qreitot) = qym(i,j+1,kc,QREINT) + sum(qrtmp)

             ! Enforce small_*
             qym(i,j+1,kc,QRHO ) = max(qym(i,j+1,kc,QRHO ),small_dens)
             qym(i,j+1,kc,QPRES) = max(qym(i,j+1,kc,QPRES),small_pres)

             do g=0,ngroups-1
                if (qym(i,j+1,kc,qrad+g) < ZERO) then
                   er_foo = - qym(i,j+1,kc,qrad+g)
                   qym(i,j+1,kc,qrad+g) = ZERO
                   qym(i,j+1,kc,qptot) = qym(i,j+1,kc,qptot) + lamm(g) * er_foo
                   qym(i,j+1,kc,qreitot) = qym(i,j+1,kc,qreitot) + er_foo
                end if
             end do

             if (qym(i,j+1,kc,QPRES) < ZERO) then
                qym(i,j+1,kc,QPRES) = p
             end if

             ! transverse velocities
             if (v < ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qym(i,j+1,kc,QU    ) = Ip(i,j,kc,2,2,QU)
                   qym(i,j+1,kc,QW    ) = Ip(i,j,kc,2,2,QW)
                else
                   qym(i,j+1,kc,QU    ) = u
                   qym(i,j+1,kc,QW    ) = w
                endif
             else ! wave is moving toward the interface
                qym(i,j+1,kc,QU    ) = Ip(i,j,kc,2,2,QU)
                qym(i,j+1,kc,QW    ) = Ip(i,j,kc,2,2,QW)
             endif

             if (ppm_trace_sources .eq. 1) then
                qym(i,j+1,kc,QU) = qym(i,j+1,kc,QU) + hdt*Ip_src(i,j,kc,2,2,QU)
                qym(i,j+1,kc,QW) = qym(i,j+1,kc,QW) + hdt*Ip_src(i,j,kc,2,2,QW)
             endif

          end if
       end do
    end do

    !-------------------------------------------------------------------------
    ! Passively advected quantities
    !-------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! Plus state on face j
       do j = ilo2, ihi2+1
          do i = ilo1-1, ihi1+1
             v = q(i,j,k3d,QV)

             if (v .gt. ZERO) then
                qyp(i,j,kc,n) = q(i,j,k3d,n)
             else if (v .lt. ZERO) then
                qyp(i,j,kc,n) = Im(i,j,kc,2,2,n)
             else
                qyp(i,j,kc,n) = q(i,j,k3d,n) + HALF*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo
       end do

       ! Minus state on face j+1
       do j = ilo2-1, ihi2
          do i = ilo1-1, ihi1+1
             v = q(i,j,k3d,QV)

             if (v .gt. ZERO) then
                qym(i,j+1,kc,n) = Ip(i,j,kc,2,2,n)
             else if (v .lt. ZERO) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n)
             else
                qym(i,j+1,kc,n) = q(i,j,k3d,n) + HALF*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo

       enddo
    enddo

  end subroutine tracexy_ppm_rad



  subroutine tracez_ppm_rad(lam, lam_lo, lam_hi, &
                            q, c, cg, flatn, qd_lo, qd_hi, &
                            Ip, Im, Ip_src, Im_src, &
                            qzm, qzp, qs_lo, qs_hi, &
                            ilo1,ilo2,ihi1,ihi2,dt,km,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, &
                                   small_dens, small_pres, &
                                   ppm_type, ppm_reference, ppm_trace_sources, &
                                   ppm_reference_edge_limit, &
                                   npassive, qpass_map
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot
    use rad_params_module, only : ngroups
    use bl_constants_module

    implicit none

    integer, intent(in) :: lam_lo(3), lam_hi(3)
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qs_lo(3), qs_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: km, kc, k3d

    double precision, intent(in) :: lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)

    double precision, intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision, intent(in) ::     c(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision, intent(in) ::    cg(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision, intent(in) :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))

    double precision, intent(in) :: Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)
    double precision, intent(in) :: Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)

    double precision, intent(in) :: Ip_src(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)
    double precision, intent(in) :: Im_src(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QRADVAR)


    double precision, intent(inout) :: qzm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QRADVAR)
    double precision, intent(inout) :: qzp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QRADVAR)


    double precision, intent(in) :: dt

    !     Local variables
    integer :: i, j, g
    integer :: n, ipassive

    double precision :: hdt

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is 
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! rhoe : total specific internal energy (including radiation,
    !        if available)
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! htot : total specific enthalpy
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    double precision :: cc, csq, cgassq, Clag
    double precision :: rho, u, v, w, p, rhoe_g, h_g
    double precision :: ptot, rhoe, htot

    double precision :: drho, dptot, drhoe_g
    double precision :: drhoe
    double precision :: dwp, dptotp
    double precision :: dwm, dptotm

    double precision :: rho_ref, w_ref, p_ref, rhoe_g_ref, h_g_ref
    double precision :: tau_ref
    double precision :: ptot_ref, rhoe_ref, htot_ref


    double precision :: alpham, alphap, alpha0r, alpha0e_g, alpha0e


    double precision, dimension(0:ngroups-1) :: er, der, alphar, qrtmp,hr
    double precision, dimension(0:ngroups-1) :: lam0, lamp, lamm

    double precision, dimension(0:ngroups-1) :: er_ref

    double precision :: er_foo

    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call bl_error("Error:: RadHydro_3d.f90 :: tracez_ppm_rad")
    end if

    hdt = HALF * dt


    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! Trace to left and right edges using upwind PPM
    !
    ! Note: in contrast to the above code for x and y, here the loop
    ! is over interfaces, not over cell-centers.


    !-------------------------------------------------------------------------
    ! construct qzp  -- plus state on face kc
    !-------------------------------------------------------------------------
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          do g=0, ngroups-1
             lam0(g) = lam(i,j,k3d,g)
             lamp(g) = lam(i,j,k3d,g)
             lamm(g) = lam(i,j,k3d,g)
          end do

          rho = q(i,j,k3d,QRHO)

          cgassq = cg(i,j,k3d)**2
          cc = c(i,j,k3d)
          csq = cc**2

          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)

          p = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          h_g = ( (p + rhoe_g)/rho)/csq

          ptot = q(i,j,k3d,qptot)
          rhoe = q(i,j,k3d,qreitot)
          htot = ( (rhoe+ptot)/rho )/csq

          er(:) = q(i,j,k3d,qrad:qradhi)
          hr(:) = (lam0+ONE)*er/rho

          ! Set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. w - cc >= ZERO .and. &
                ppm_reference_edge_limit == 0) ) then
             ! original Castro way -- cc value
             rho_ref  = rho
             w_ref    = w

             p_ref    = p
             rhoe_g_ref = rhoe_g

             ptot_ref = ptot
             rhoe_ref = rhoe

             er_ref(:) = er(:)
          else
             ! This will be the fastest moving state to the left
             rho_ref  = Im(i,j,kc,3,1,QRHO)
             w_ref    = Im(i,j,kc,3,1,QW)

             p_ref    = Im(i,j,kc,3,1,QPRES)
             rhoe_g_ref = Im(i,j,kc,3,1,QREINT)

             ptot_ref = Im(i,j,kc,3,1,QPTOT)
             rhoe_ref = Im(i,j,kc,3,1,QREITOT)

             er_ref(:) = Im(i,j,kc,3,1,QRAD:QRADHI)
          endif

          rho_ref = max(rho_ref,small_dens)
          p_ref = max(p_ref,small_pres)

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm    = w_ref    - Im(i,j,kc,3,1,QW)
          dptotm = ptot_ref - Im(i,j,kc,3,1,qptot)

          drho    = rho_ref    - Im(i,j,kc,3,2,QRHO)
          dptot   = ptot_ref   - Im(i,j,kc,3,2,qptot)
          drhoe   = rhoe_ref   - Im(i,j,kc,3,2,qreitot)
          drhoe_g = rhoe_g_ref - Im(i,j,kc,3,2,QREINT)
          der(:)  = er_ref(:)  - Im(i,j,kc,3,2,qrad:qradhi)

          dwp    = w_ref    - Im(i,j,kc,3,3,QW)
          dptotp = ptot_ref - Im(i,j,kc,3,3,qptot)

          ! If we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources .eq. 1) then
             dwm = dwm - hdt*Im_src(i,j,kc,3,1,QW)
             dwp = dwp - hdt*Im_src(i,j,kc,3,3,QW)
          endif

          ! Optionally use the reference state in evaluating the
          ! eigenvectors -- NOT YET IMPLEMENTED

          alpham = HALF*(dptotm/(rho*cc) - dwm)*rho/cc
          alphap = HALF*(dptotp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dptot/csq
          alpha0e = drhoe - dptot*htot
          alpha0e_g = drhoe_g - dptot*h_g
          alphar(:) = der(:) - dptot/csq*hr

          if (w-cc .gt. ZERO) then
             alpham = ZERO
          else if (w-cc .lt. ZERO) then
             alpham = -alpham
          else
             alpham = -HALF*alpham
          endif

          if (w+cc .gt. ZERO) then
             alphap = ZERO
          else if (w+cc .lt. ZERO) then
             alphap = -alphap
          else
             alphap = -HALF*alphap
          endif

          if (w .gt. ZERO) then
             alpha0r = ZERO
             alpha0e = ZERO
             alpha0e_g = ZERO
             alphar(:) = ZERO
          else if (w .lt. ZERO) then
             alpha0r = -alpha0r
             alpha0e = -alpha0e
             alpha0e_g = -alpha0e_g
             alphar(:) = -alphar(:)
          else
             alpha0r = -HALF*alpha0r
             alpha0e = -HALF*alpha0e
             alpha0e_g = -HALF*alpha0e_g
             alphar(:) = -HALF*alphar(:)
          endif

          qzp(i,j,kc,QRHO) = rho_ref + alphap + alpham + alpha0r
          qzp(i,j,kc,QW) = w_ref + (alphap - alpham)*cc/rho
          qzp(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
          qzp(i,j,kc,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

          qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
          qzp(i,j,kc,qrad:qradhi) = qrtmp

          qzp(i,j,kc,qptot) = ptot_ref + (alphap + alpham)*csq
          qzp(i,j,kc,qreitot) = qzp(i,j,kc,QREINT) + sum(qrtmp)

          ! Enforce small_*
          qzp(i,j,kc,QRHO ) = max(qzp(i,j,kc,QRHO ),small_dens)
          qzp(i,j,kc,QPRES) = max(qzp(i,j,kc,QPRES),small_pres)

          do g=0,ngroups-1
             if (qzp(i,j,kc,qrad+g) < ZERO) then
                er_foo = - qzp(i,j,kc,qrad+g)
                qzp(i,j,kc,qrad+g) = ZERO
                qzp(i,j,kc,qptot) = qzp(i,j,kc,qptot) + lamp(g) * er_foo
                qzp(i,j,kc,qreitot) = qzp(i,j,kc,qreitot) + er_foo
             end if
          end do

          if (qzp(i,j,kc,QPRES) < ZERO) then
             qzp(i,j,kc,QPRES) = p
          end if

          ! transverse velocities
          if (w > ZERO) then
             if (ppm_reference_edge_limit == 1) then
                qzp(i,j,kc,QU    ) = Im(i,j,kc,3,2,QU)
                qzp(i,j,kc,QV    ) = Im(i,j,kc,3,2,QV)
             else
                qzp(i,j,kc,QU    ) = u
                qzp(i,j,kc,QV    ) = v
             endif
          else ! wave moving toward the interface
             qzp(i,j,kc,QU    ) = Im(i,j,kc,3,2,QU)
             qzp(i,j,kc,QV    ) = Im(i,j,kc,3,2,QV)
          endif

          if (ppm_trace_sources .eq. 1) then
             qzp(i,j,kc,QU) = qzp(i,j,kc,QU) + hdt*Im_src(i,j,kc,3,2,QU)
             qzp(i,j,kc,QV) = qzp(i,j,kc,QV) + hdt*Im_src(i,j,kc,3,2,QV)
          endif


          !-------------------------------------------------------------------
          ! This is all for qzm -- minus state on face kc
          !-------------------------------------------------------------------

          ! Note this is different from how we do 1D, 2D, and the
          ! x and y-faces in 3D, where the analogous thing would have
          ! been to find the minus state on face kc+1

          do g=0, ngroups-1
             lam0(g) = lam(i,j,k3d-1,g)
             lamp(g) = lam(i,j,k3d-1,g)
             lamm(g) = lam(i,j,k3d-1,g)
          end do

          rho = q(i,j,k3d-1,QRHO)

          cgassq = cg(i,j,k3d-1)**2
          cc = c(i,j,k3d-1)
          csq = cc**2

          u = q(i,j,k3d-1,QU)
          v = q(i,j,k3d-1,QV)
          w = q(i,j,k3d-1,QW)

          p = q(i,j,k3d-1,QPRES)
          rhoe_g = q(i,j,k3d-1,QREINT)
          h_g = ( (p + rhoe_g)/rho)/csq

          ptot = q(i,j,k3d-1,qptot)
          rhoe = q(i,j,k3d-1,qreitot)
          htot = ( (rhoe+ptot)/rho )/csq

          er(:) = q(i,j,k3d-1,qrad:qradhi)
          hr(:) = (lam0+ONE)*er/rho


          ! Set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. w + cc <= ZERO .and. &
                ppm_reference_edge_limit == 0) ) then
             rho_ref  = rho
             w_ref    = w

             p_ref    = p
             rhoe_g_ref = rhoe_g

             ptot_ref    = ptot
             rhoe_ref = rhoe

             er_ref(:) = er(:)
          else
             ! This will be the fastest moving state to the right
             rho_ref  = Ip(i,j,km,3,3,QRHO)
             w_ref    = Ip(i,j,km,3,3,QW)

             p_ref    = Ip(i,j,km,3,3,QPRES)
             rhoe_g_ref = Ip(i,j,km,3,3,QREINT)

             ptot_ref    = Ip(i,j,km,3,3,QPTOT)
             rhoe_ref = Ip(i,j,km,3,3,QREITOT)

             er_ref(:) = Ip(i,j,km,3,3,QRAD:QRADHI)
          endif

          rho_ref = max(rho_ref,small_dens)
          p_ref = max(p_ref,small_pres)

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm    = w_ref    - Ip(i,j,km,3,1,QW)
          dptotm = ptot_ref - Ip(i,j,km,3,1,qptot)

          drho    = rho_ref    - Ip(i,j,km,3,2,QRHO)
          dptot   = ptot_ref   - Ip(i,j,km,3,2,qptot)
          drhoe   = rhoe_ref   - Ip(i,j,km,3,2,qreitot)
          drhoe_g = rhoe_g_ref - Ip(i,j,km,3,2,QREINT)
          der(:)  = er_ref(:)  - Ip(i,j,km,3,2,qrad:qradhi)

          dwp    = w_ref    - Ip(i,j,km,3,3,QW)
          dptotp = ptot_ref - Ip(i,j,km,3,3,qptot)

          ! If we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources .eq. 1) then
             dwm = dwm - hdt*Ip_src(i,j,km,3,1,QW)
             dwp = dwp - hdt*Ip_src(i,j,km,3,3,QW)
          endif

          ! Optionally use the reference state in evaluating the
          ! eigenvectors -- NOT YET IMPLEMENTED

          alpham = HALF*(dptotm/(rho*cc) - dwm)*rho/cc
          alphap = HALF*(dptotp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dptot/csq
          alpha0e = drhoe - dptot*htot
          alpha0e_g = drhoe_g - dptot*h_g
          alphar(:) = der(:) - dptot/csq*hr

          if (w-cc .gt. ZERO) then
             alpham = -alpham
          else if (w-cc .lt. ZERO) then
             alpham = ZERO
          else
             alpham = -HALF*alpham
          endif

          if (w+cc .gt. ZERO) then
             alphap = -alphap
          else if (w+cc .lt. ZERO) then
             alphap = ZERO
          else
             alphap = -HALF*alphap
          endif

          if (w .gt. ZERO) then
             alpha0r = -alpha0r
             alpha0e = -alpha0e
             alpha0e_g = -alpha0e_g
             alphar(:) = -alphar(:)
          else if (w .lt. ZERO) then
             alpha0r = ZERO
             alpha0e = ZERO
             alpha0e_g = ZERO
             alphar(:) = ZERO
          else
             alpha0r = -HALF*alpha0r
             alpha0e = -HALF*alpha0e
             alpha0e_g = -HALF*alpha0e_g
             alphar(:) = -HALF*alphar(:)
          endif

          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          ! note that the a{mpz}left as defined above have the minus already
          qzm(i,j,kc,QRHO) = rho_ref + alphap + alpham + alpha0r
          qzm(i,j,kc,QW) = w_ref + (alphap - alpham)*cc/rho
          qzm(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
          qzm(i,j,kc,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))

          qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
          qzm(i,j,kc,qrad:qradhi) = qrtmp

          qzm(i,j,kc,qptot) = ptot_ref + (alphap + alpham)*csq
          qzm(i,j,kc,qreitot) = qzm(i,j,kc,QREINT) + sum(qrtmp)

          ! Enforce small_*
          qzm(i,j,kc,QRHO ) = max(qzm(i,j,kc,QRHO ),small_dens)
          qzm(i,j,kc,QPRES) = max(qzm(i,j,kc,QPRES),small_pres)

          do g=0,ngroups-1
             if (qzm(i,j,kc,qrad+g) < ZERO) then
                er_foo = - qzm(i,j,kc,qrad+g)
                qzm(i,j,kc,qrad+g) = ZERO
                qzm(i,j,kc,qptot) = qzm(i,j,kc,qptot) + lamm(g) * er_foo
                qzm(i,j,kc,qreitot) = qzm(i,j,kc,qreitot) + er_foo
             end if
          end do

          if (qzm(i,j,kc,QPRES) < ZERO) then
             qzm(i,j,kc,QPRES) = p
          end if

          ! Transverse velocity
          if (w < ZERO) then
             if (ppm_reference_edge_limit == 1) then
                qzm(i,j,kc,QU    ) = Ip(i,j,km,3,2,QU)
                qzm(i,j,kc,QV    ) = Ip(i,j,km,3,2,QV)
             else
                qzm(i,j,kc,QU    ) = u
                qzm(i,j,kc,QV    ) = v
             endif
          else ! wave moving toward the interface
             qzm(i,j,kc,QU    ) = Ip(i,j,km,3,2,QU)
             qzm(i,j,kc,QV    ) = Ip(i,j,km,3,2,QV)
          endif

          if (ppm_trace_sources .eq. 1) then
             qzm(i,j,kc,QU) = qzm(i,j,kc,QU) + hdt*Ip_src(i,j,km,3,2,QU)
             qzm(i,j,kc,QV) = qzm(i,j,kc,QV) + hdt*Ip_src(i,j,km,3,2,QV)
          endif

       end do
    end do

    !-------------------------------------------------------------------------
    ! passively advected quantities
    !-------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! Plus state on face kc
             w = q(i,j,k3d,QW)

             if (w .gt. ZERO) then
                qzp(i,j,kc,n) = q(i,j,k3d,n)
             else if (w .lt. ZERO) then
                qzp(i,j,kc,n) = Im(i,j,kc,3,2,n)
             else
                qzp(i,j,kc,n) = q(i,j,k3d,n) + HALF*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
             endif

             ! Minus state on face k
             w = q(i,j,k3d-1,QW)

             if (w .gt. ZERO) then
                qzm(i,j,kc,n) = Ip(i,j,km,3,2,n)
             else if (w .lt. ZERO) then
                qzm(i,j,kc,n) = q(i,j,k3d-1,n)
             else
                qzm(i,j,kc,n) = q(i,j,k3d-1,n) + HALF*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
             endif

          enddo
       enddo
    enddo

  end subroutine tracez_ppm_rad

end module trace_ppm_rad_module
