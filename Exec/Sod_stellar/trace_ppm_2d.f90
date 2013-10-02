module trace_ppm_module

  implicit none

  private

  public trace_ppm

contains

  subroutine trace_ppm(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                       dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                       qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       grav,gv_l1,gv_l2,gv_h1,gv_h2, &
                       gamc,gc_l1,gc_l2,gc_h1,gc_h2, &
                       ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

    use network, only : nspec, naux
    use eos_type_module
    use eos_module

    use bl_constants_module, ONLY: ONE, ZERO, HALF

    use meth_params_module, only : iorder, QVAR, QRHO, QU, QV, &
         QREINT, QPRES, QFA, QFS, QFX, QTEMP, &
         nadv, small_dens, small_pres, small_temp,  &
         ppm_type, ppm_reference, ppm_trace_grav, ppm_temp_fix, &
         ppm_tau_in_tracing, ppm_reference_eigenvectors
    use ppm_module, only : ppm

    implicit none

    integer ilo1,ilo2,ihi1,ihi2
    integer qd_l1,qd_l2,qd_h1,qd_h2
    integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
    integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer gv_l1,gv_l2,gv_h1,gv_h2
    integer gc_l1,gc_l2,gc_h1,gc_h2

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
    
    double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)

    double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,2)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)

    double precision dx, dy, dt
    
    ! Local variables
    integer i, j, iwave, idim
    integer n, iadv
    integer ns, ispec, iaux
    
    double precision dtdx, dtdy
    double precision cc, csq, Clag, rho, u, v, p, rhoe, T
    double precision drho, du, dv, dp, drhoe, dT0, dtau
    double precision drhop, dup, dvp, dpp, drhoep, dTp, dtaup
    double precision drhom, dum, dvm, dpm, drhoem, dTm, dtaum

    double precision :: rho_ref, u_ref, v_ref, p_ref, rhoe_ref, T_ref, tau_ref
    double precision :: tau_s, e_s, de

    double precision :: cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref
    double precision :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev, T_ev, tau_ev
    double precision :: gam
    double precision :: p_r, p_T

    double precision enth, alpham, alphap, alpha0r, alpha0e
    double precision alpha0u, alpha0v
    double precision apright, amright, azrright, azeright
    double precision azu1rght, azv1rght
    double precision apleft, amleft, azrleft, azeleft
    double precision azu1left, azv1left
    double precision sourcr,sourcp,source,courn,eta,dlogatmp

    double precision :: xi, xi1
    double precision :: halfdt

    integer, parameter :: igx = 1
    integer, parameter :: igy = 2

    double precision, allocatable :: Ip(:,:,:,:,:)
    double precision, allocatable :: Im(:,:,:,:,:)

    double precision, allocatable :: Ip_g(:,:,:,:,:)
    double precision, allocatable :: Im_g(:,:,:,:,:)

    ! gamma_c/1 on the interfaces
    double precision, allocatable :: Ip_gc(:,:,:,:,:)
    double precision, allocatable :: Im_gc(:,:,:,:,:)

    double precision, allocatable :: tau(:,:)
    double precision, allocatable :: Ip_tau(:,:,:,:,:)
    double precision, allocatable :: Im_tau(:,:,:,:,:)

    double precision :: eval(3), beta(3), rvec(3,3), lvec(3,3), dq(3)

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

    if (ppm_trace_grav == 1) then
       allocate(Ip_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,2))
       allocate(Im_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,2))
    endif

    ! tau = 1/rho
    allocate(tau(qd_l1:qd_h1,qd_l2:qd_h2))

    do j = qd_l2, qd_h2
       do i = qd_l1, qd_h1
          if (q(i,j,QRHO) <= 0.0d0) then
             print *, 'error with rho'
          endif
          tau(i,j) = 1.0d0/q(i,j,QRHO)
       enddo
    enddo
    
    if (ppm_temp_fix == 3) then   
       allocate(Ip_tau(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,1))
       allocate(Im_tau(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,1))
    endif

    allocate(Ip_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,1))
    allocate(Im_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,1))

    halfdt = 0.5d0 * dt


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
                q(:,:,QU:),c,qd_l1,qd_l2,qd_h1,qd_h2, &
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
                         
                   call eos(eos_input_rt, eos_state, .false.)
                         
                   Ip(i,j,idim,iwave,QPRES) = eos_state%p
                   Ip(i,j,idim,iwave,QREINT) = Ip(i,j,idim,iwave,QRHO)*eos_state%e
                   Ip_gc(i,j,idim,iwave,1) = eos_state%gam1

                         
                   eos_state%rho   = Im(i,j,idim,iwave,QRHO)
                   eos_state%T     = Im(i,j,idim,iwave,QTEMP)
                   eos_state%xn(:) = Im(i,j,idim,iwave,QFS:QFS-1+nspec)
                         
                   call eos(eos_input_rt, eos_state, .false.)
                         
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
                q(:,:,QU:),c,qd_l1,qd_l2,qd_h1,qd_h2, &
                Ip_gc(:,:,:,:,1),Im_gc(:,:,:,:,1), &
                ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    endif


    if (ppm_temp_fix == 3) then
       call ppm(tau(:,:),qd_l1,qd_l2,qd_h1,qd_h2, &
                q(:,:,QU:),c,qd_l1,qd_l2,qd_h1,qd_h2, &
                Ip_tau(:,:,:,:,1),Im_tau(:,:,:,:,1), &
                ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    endif

    ! if desired, do parabolic reconstruction of the gravitational
    ! acceleration -- we'll use this for the force on the velocity
    if (ppm_trace_grav == 1) then
       do n = 1,2
          call ppm(grav(:,:,n),gv_l1,gv_l2,gv_h1,gv_h2, &
                   q(:,:,QU:),c,qd_l1,qd_l2,qd_h1,qd_h2, &
                   Ip_g(:,:,:,:,n),Im_g(:,:,:,:,n), &
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
          T = q(i,j,QTEMP)

          Clag = rho*cc

          gam = gamc(i,j)


          !-------------------------------------------------------------------
          ! plus state on face i
          !-------------------------------------------------------------------

          ! set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. u - cc >= 0.0d0)) then
             ! original Castro way -- cc value
             rho_ref  = rho
             u_ref    = u
             v_ref    = v

             p_ref    = p
             rhoe_ref = rhoe
             T_ref    = T

             tau_ref  = tau(i,j)

             gam_ref  = gamc(i,j)

          else
             ! this will be the fastest moving state to the left --
             ! this is the method that Miller & Colella and Colella &
             ! Woodward use
             rho_ref  = Im(i,j,1,1,QRHO)
             u_ref    = Im(i,j,1,1,QU)
             v_ref    = Im(i,j,1,1,QV)

             p_ref    = Im(i,j,1,1,QPRES)
             rhoe_ref = Im(i,j,1,1,QREINT)
             T_ref = Im(i,j,1,1,QTEMP)

             if (ppm_temp_fix == 3) then
                ! use the parabolic reconstruction of tau
                tau_ref = Im_tau(i,j,1,1,1)
             else
                ! use the parabolic reconstruction of rho
                tau_ref = 1.0d0/Im(i,j,1,1,QRHO)
             endif

             gam_ref = Im_gc(i,j,1,1,1)
          endif

          ! for tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

          ! *m are the jumps carried by u-c
          ! *p are the jumps carried by u+c

          dum    = (u_ref    - Im(i,j,1,1,QU))
          !dvm    = (v_ref    - Im(i,j,1,1,QV))
          dpm    = (p_ref    - Im(i,j,1,1,QPRES))
          drhoem = (rhoe_ref - Im(i,j,1,1,QREINT))
          dTm    = (T_ref    - Im(i,j,1,1,QTEMP)) 

          drho  = (rho_ref  - Im(i,j,1,2,QRHO))
          du    = (u_ref    - Im(i,j,1,2,QU))
          dv    = (v_ref    - Im(i,j,1,2,QV))
          dp    = (p_ref    - Im(i,j,1,2,QPRES))
          drhoe = (rhoe_ref - Im(i,j,1,2,QREINT))
          dT0   = (T_ref    - Im(i,j,1,2,QTEMP))

          dup    = (u_ref    - Im(i,j,1,3,QU))
          !dvp    = (v_ref    - Im(i,j,1,3,QV))
          dpp    = (p_ref    - Im(i,j,1,3,QPRES))
          drhoep = (rhoe_ref - Im(i,j,1,3,QREINT))
          dTp    = (T_ref    - Im(i,j,1,3,QTEMP))

          if (ppm_temp_fix < 3) then
             ! we are relying on tau as built from the reconstructed rho
             ! parabolas
             dtau  = (tau_ref  - 1.0d0/Im(i,j,1,2,QRHO))
          else
             ! we are directly reconstructing tau as parabola
             dtaum  = (tau_ref  - Im_tau(i,j,1,1,1))
             dtau  = (tau_ref  - Im_tau(i,j,1,2,1))
             dtaup  = (tau_ref  - Im_tau(i,j,1,3,1))
          endif

          ! if we are doing gravity tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the 
          ! trans_X routines
          if (ppm_trace_grav == 1) then
             dum = dum - halfdt*Im_g(i,j,1,1,igx)
             du  = du  - halfdt*Im_g(i,j,1,2,igx)
             dup = dup - halfdt*Im_g(i,j,1,3,igx)

             dv  = dv  - halfdt*Im_g(i,j,1,2,igy)
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
             T_ev    = T
             tau_ev  = tau(i,j)
          else
             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             enth_ev = enth_ref
             p_ev    = p_ref
             T_ev    = T_ref
             tau_ev  = tau_ref
          endif


          ! there are several options here on how to do the tracing
          ! for ppm_temp_fix < 3, we use:
          !
          !   rho, u, p  if ppm_tau_in_tracing = 0
          !   tau, u, p  if ppm_tau_in_tracing = 1
          !
          ! for ppm_temp_fix = 3, we use tau, u, T

          if (ppm_temp_fix < 3) then

             if (ppm_tau_in_tracing == 0) then

                ! these are analogous to the beta's from the original
                ! PPM paper (except we work with rho instead of tau).
                ! This is simply (l . dq), where dq = qref - I(q)

                alpham = 0.5d0*(dpm/(rho_ev*cc_ev) - dum)*rho_ev/cc_ev
                alphap = 0.5d0*(dpp/(rho_ev*cc_ev) + dup)*rho_ev/cc_ev
                alpha0r = drho - dp/csq_ev
                alpha0e = drhoe - dp*enth_ev  ! note enth has a 1/c**2 in it
                alpha0v = dv

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
                else if (u .lt. 0.d0) then
                   azrright = -alpha0r
                   azeright = -alpha0e
                   azv1rght = -alpha0v
                else
                   azrright = -0.5d0*alpha0r
                   azeright = -0.5d0*alpha0e
                   azv1rght = -0.5d0*alpha0v
                endif

                ! the final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                if (i .ge. ilo1) then
                
                   xi1 = 1.0d0-flatn(i,j)
                   xi = flatn(i,j)

                   qxp(i,j,QRHO)   = xi1*rho  + xi*(rho_ref + apright + amright + azrright)

                   qxp(i,j,QU)     = xi1*u    + xi*(u_ref + (apright - amright)*cc_ev/rho_ev)
                   qxp(i,j,QV)     = xi1*v    + xi*(v_ref + azv1rght)
                
                   qxp(i,j,QREINT) = xi1*rhoe + xi*(rhoe_ref + (apright + amright)*enth_ev*csq_ev + azeright)
                   qxp(i,j,QPRES)  = xi1*p    + xi*(p_ref + (apright + amright)*csq_ev)

                   qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))
                   qxp(i,j,QPRES) = max(qxp(i,j,QPRES), small_pres)
                end if

             else  
                ! (tau, u, p, e) eigensystem

                ! this is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system.

                ! we are dealing with e
                de = (rhoe_ref/rho_ref - Im(i,j,1,2,QREINT)/Im(i,j,1,2,QRHO))

                alpham = 0.5d0*( dum - dpm/Clag_ev)/Clag_ev
                alphap = 0.5d0*(-dup - dpp/Clag_ev)/Clag_ev
                alpha0r = dtau + dp/Clag_ev**2
                alpha0e = de - dp*p_ev/Clag_ev**2
                alpha0v = dv

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
                else if (u .lt. 0.d0) then
                   azrright = -alpha0r
                   azeright = -alpha0e
                   azv1rght = -alpha0v
                else
                   azrright = -0.5d0*alpha0r
                   azeright = -0.5d0*alpha0e
                   azv1rght = -0.5d0*alpha0v
                endif
                
                ! the final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                if (i .ge. ilo1) then

                   xi1 = 1.0d0-flatn(i,j)
                   xi = flatn(i,j)

                   tau_s = tau_ref + apright + amright + azrright
                   qxp(i,j,QRHO)   = xi1*rho + xi/tau_s
                   
                   qxp(i,j,QU)     = xi1*u    + xi*(u_ref + (amright - apright)*Clag_ev)
                   qxp(i,j,QV)     = xi1*v    + xi*(v_ref + azv1rght)
                   
                   e_s = rhoe_ref/rho_ref + (azeright - p_ev*amright -p_ev*apright)
                   qxp(i,j,QREINT) = xi1*rhoe + xi*e_s/tau_s

                   qxp(i,j,QPRES)  = xi1*p    + xi*(p_ref + (-apright - amright)*Clag_ev**2)
                   
                   qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))
                   qxp(i,j,QPRES) = max(qxp(i,j,QPRES), small_pres)                   
                endif

             endif
          else

             ! ppm_temp_fix = 3:
             ! (tau, u, T) eigensystem PPM 

             ! eos to get some thermodynamics
             eos_state%T = T_ev
             eos_state%rho = rho_ev
             eos_state%xn(:) = q(i,j,QFS:QFS-1+nspec)
             
             call eos(eos_input_rt, eos_state, .false.)

             p_r = eos_state%dpdr
             p_T = eos_state%dpdT

            
             ! construct eigenvectors and eigenvalues
             eval(1) = u - cc
             eval(2) = u
             eval(3) = u + cc


             ! compute the left eigenvectors                                             
             lvec(1,:) = [ HALF*p_r/(cc_ev*cc_ev),     HALF*tau_ev/cc_ev,   -HALF*tau_ev**2*p_T/(cc_ev*cc_ev) ]   ! u - c    
             lvec(2,:) = [ 1.0d0 - p_r/(cc_ev*cc_ev),  ZERO,               tau_ev**2*p_T/(cc_ev*cc_ev) ]   ! u        
             lvec(3,:) = [ HALF*p_r/(cc_ev*cc_ev),     -HALF*tau_ev/cc_ev,  -HALF*tau_ev**2*p_T/(cc_ev*cc_ev) ]   ! u + c    

             ! compute the right eigenvectors                                            
             rvec(1,:) = [ ONE,  cc_ev/rho_ev,  -(1.0d0/tau_ev**2)*(cc_ev*cc_ev - p_r)/p_T ]   ! u - c 
             rvec(2,:) = [ ONE,  ZERO,     (1.0d0/tau_ev**2)*p_r/p_T  ]   ! u     
             rvec(3,:) = [ ONE,  -cc_ev/rho_ev, -(1.0d0/tau_ev**2)*(cc_ev*cc_ev - p_r)/p_T ]   ! u + c   
             

             ! construct interface states of T, u, tau

             ! loop over waves and construct l . (qref - I(q))
             do iwave = 1, 3
                select case (iwave)
                case (1)
                   dq = [dtaum, dum, dTm]
                case (2)
                   dq = [dtau,  du,  dT0]
                case (3)
                   dq = [dtaup, dup, dTp]
                end select

                beta(iwave) = dot_product(lvec(iwave,:),dq(:))
             enddo

             if (i .ge. ilo1) then
                
                ! here we compute 
                !
                ! q_{l,r} = q_ref - sum { l . (q_ref - I(q)) r } 
                !
                ! where the sum is over the waves and
                ! limited to only include those moving toward the
                ! interface

                xi1 = 1.0d0-flatn(i,j)
                xi = flatn(i,j)
                
                ! tau (and density) state
                tau_s = tau_ref
                do iwave = 1, 3
                   if (eval(iwave) <= 0.0d0) tau_s = tau_s - beta(iwave)*rvec(iwave,1)
                enddo
                qxp(i,j,QRHO)  = xi1*rho + xi/tau_s

                ! u state
                qxp(i,j,QU)    = xi1*u + xi*u_ref 
                do iwave = 1, 3
                   if (eval(iwave) <= 0.0d0) qxp(i,j,QU) = qxp(i,j,QU) - xi*beta(iwave)*rvec(iwave,2)
                enddo

                ! T state
                qxp(i,j,QTEMP) = xi1*T + xi*T_ref 
                do iwave = 1, 3
                   if (eval(iwave) <= 0.0d0) qxp(i,j,QTEMP) = qxp(i,j,QTEMP) - xi*beta(iwave)*rvec(iwave,3)
                enddo
                qxp(i,j,QTEMP) = max(qxp(i,j,QTEMP), small_temp)

                ! transverse velocity
                if (u < 0.0d0) then
                   qxp(i,j,QV) = -dv
                else
                   qxp(i,j,QV) = 0.0d0
                endif
                
                qxp(i,j,QV) = qxp(i,j,QV) + xi1*v + xi*v_ref 

                ! limit
                qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))

             endif

          endif


          !-------------------------------------------------------------------
          ! minus state on face i+1
          !-------------------------------------------------------------------

          ! set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. u + cc <= 0.0d0) ) then
             ! original Castro way -- cc values
             rho_ref  = rho
             u_ref    = u
             v_ref    = v

             p_ref    = p
             rhoe_ref = rhoe

             T_ref    = T

             tau_ref  = tau(i,j)

             gam_ref  = gamc(i,j)

          else
             ! this will be the fastest moving state to the right
             rho_ref  = Ip(i,j,1,3,QRHO)
             u_ref    = Ip(i,j,1,3,QU)
             v_ref    = Ip(i,j,1,3,QV)

             p_ref    = Ip(i,j,1,3,QPRES)
             rhoe_ref = Ip(i,j,1,3,QREINT)

             T_ref    = Ip(i,j,1,3,QTEMP)

             if (ppm_temp_fix == 3) then
                ! use the parabolic reconstruction of tau
                tau_ref  = Ip_tau(i,j,1,3,1)
             else
                ! use the parabolic reconstruction of rho                
                tau_ref  = 1.0d0/Ip(i,j,1,3,QRHO)
             endif

             gam_ref = Ip_gc(i,j,1,3,1)
          endif

          ! for tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref          

          ! *m are the jumps carried by u-c
          ! *p are the jumps carried by u+c

          dum    = (u_ref    - Ip(i,j,1,1,QU))
          !dvm    = (v_ref    - Ip(i,j,1,1,QV))
          dpm    = (p_ref    - Ip(i,j,1,1,QPRES))
          drhoem = (rhoe_ref - Ip(i,j,1,1,QREINT))
          dTm    = (T_ref    - Ip(i,j,1,1,QTEMP))
          
          drho  = (rho_ref  - Ip(i,j,1,2,QRHO))
          du    = (u_ref    - Ip(i,j,1,2,QU))
          dv    = (v_ref    - Ip(i,j,1,2,QV))
          dp    = (p_ref    - Ip(i,j,1,2,QPRES))
          drhoe = (rhoe_ref - Ip(i,j,1,2,QREINT))
          dT0   = (T_ref    - Ip(i,j,1,2,QTEMP))
        
          dup    = (u_ref    - Ip(i,j,1,3,QU))
          !dvp    = (v_ref    - Ip(i,j,1,3,QV))
          dpp    = (p_ref    - Ip(i,j,1,3,QPRES))
          drhoep = (rhoe_ref - Ip(i,j,1,3,QREINT))
          dTp    = (T_ref    - Ip(i,j,1,3,QTEMP))

          if (ppm_temp_fix < 3) then
             ! we are relying on tau as built from the reconstructed rho
             ! parabolas
             dtau  = (tau_ref  - 1.0d0/Ip(i,j,1,2,QRHO))
          else
             dtaum  = (tau_ref  - Ip_tau(i,j,1,1,1))
             dtau  = (tau_ref  - Ip_tau(i,j,1,2,1))
             dtaup  = (tau_ref  - Ip_tau(i,j,1,3,1))
          endif

          ! if we are doing gravity tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the 
          ! trans_X routines
          if (ppm_trace_grav == 1) then
             dum = dum - halfdt*Ip_g(i,j,1,1,igx)
             du  = du  - halfdt*Ip_g(i,j,1,2,igx)
             dup = dup - halfdt*Ip_g(i,j,1,3,igx)

             dv  = dv  - halfdt*Ip_g(i,j,1,2,igy)
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
             T_ev    = T
             tau_ev  = tau(i,j)
          else
             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             enth_ev = enth_ref
             p_ev    = p_ref
             T_ev    = T_ref
             tau_ev  = tau_ref
          endif


          ! there are several options here on how to do the tracing
          ! for ppm_temp_fix < 3, we use:
          !
          !   rho, u, p  if ppm_tau_in_tracing = 0
          !   tau, u, p  if ppm_tau_in_tracing = 1
          !
          ! for ppm_temp_fix = 3, we use tau, u, T
          
          if (ppm_temp_fix < 3) then

             if (ppm_tau_in_tracing == 0) then

                ! these are analogous to the beta's from the original
                ! PPM paper (except we work with rho instead of tau).
                ! This is simply (l . dq), where dq = qref - I(q)
                alpham = 0.5d0*(dpm/(rho_ev*cc_ev) - dum)*rho_ev/cc_ev
                alphap = 0.5d0*(dpp/(rho_ev*cc_ev) + dup)*rho_ev/cc_ev
                alpha0r = drho - dp/csq_ev
                alpha0e = drhoe - dp*enth_ev   ! enth has a 1/c**2 in it
                alpha0v = dv
             
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
                else if (u .lt. 0.d0) then
                   azrleft = 0.d0
                   azeleft = 0.d0
                   azv1left = 0.d0
                else
                   azrleft = -0.5d0*alpha0r
                   azeleft = -0.5d0*alpha0e
                   azv1left = -0.5d0*alpha0v
                endif

                ! the final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                if (i .le. ihi1) then
                   xi1 = 1.0d0 - flatn(i,j)
                   xi = flatn(i,j)
                   
                   qxm(i+1,j,QRHO)   = xi1*rho  + xi*(rho_ref + apleft + amleft + azrleft)
                   qxm(i+1,j,QU)     = xi1*u    + xi*(u_ref + (apleft - amleft)*cc_ev/rho_ev)
                   qxm(i+1,j,QV)     = xi1*v    + xi*(v_ref + azv1left)
                   
                   qxm(i+1,j,QREINT) = xi1*rhoe + xi*(rhoe_ref + (apleft + amleft)*enth_ev*csq_ev + azeleft)
                   qxm(i+1,j,QPRES)  = xi1*p    + xi*(p_ref + (apleft + amleft)*csq_ev)
                   
                   qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
                   qxm(i+1,j,QPRES) = max(qxm(i+1,j,QPRES), small_pres)
                end if

             else
                ! (tau, u, p, e) eigensystem

                ! this is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system.

                de = (rhoe_ref/rho_ref - Ip(i,j,1,2,QREINT)/Ip(i,j,1,2,QRHO))

                alpham = 0.5d0*( dum - dpm/Clag_ev)/Clag_ev
                alphap = 0.5d0*(-dup - dpp/Clag_ev)/Clag_ev
                alpha0r = dtau + dp/Clag_ev**2
                alpha0e = de - dp*p_ev/Clag_ev**2
                alpha0v = dv

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
                else if (u .lt. 0.d0) then
                   azrleft = 0.d0
                   azeleft = 0.d0
                   azv1left = 0.d0
                else
                   azrleft = -0.5d0*alpha0r
                   azeleft = -0.5d0*alpha0e
                   azv1left = -0.5d0*alpha0v
                endif
             
                ! the final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                if (i .le. ihi1) then
                   
                   xi1 = 1.0d0 - flatn(i,j)
                   xi = flatn(i,j)
                
                   tau_s = tau_ref + (apleft + amleft + azrleft)
                   qxm(i+1,j,QRHO)   = xi1*rho  + xi/tau_s
                   
                   qxm(i+1,j,QU)     = xi1*u    + xi*(u_ref + (amleft - apleft)*Clag_ev)
                   qxm(i+1,j,QV)     = xi1*v    + xi*(v_ref + azv1left)

                   e_s = rhoe_ref/rho_ref + (azeleft - p_ev*amleft -p_ev*apleft)
                   qxm(i+1,j,QREINT) = xi1*rhoe + xi*e_s/tau_s

                   qxm(i+1,j,QPRES)  = xi1*p    + xi*(p_ref + (-apleft - amleft)*Clag_ev**2)
                   
                   qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
                   qxm(i+1,j,QPRES) = max(qxm(i+1,j,QPRES), small_pres)

                end if

             endif

          else

             ! T eigensystem

             ! eos to get some thermodynamics
             eos_state%T = T_ev
             eos_state%rho = rho_ev
             eos_state%xn(:) = q(i,j,QFS:QFS-1+nspec)

             call eos(eos_input_rt, eos_state, .false.)

             p_r = eos_state%dpdr
             p_T = eos_state%dpdT

            
             ! construct eigenvectors and eigenvalues
             eval(1) = u - cc
             eval(2) = u
             eval(3) = u + cc

             ! compute the left eigenvectors                                             
             lvec(1,:) = [ HALF*p_r/(cc_ev*cc_ev),     HALF*tau_ev/cc_ev,   -HALF*tau_ev**2*p_T/(cc_ev*cc_ev) ]   ! u - c    
             lvec(2,:) = [ 1.0d0 - p_r/(cc_ev*cc_ev),  ZERO,               tau_ev**2*p_T/(cc_ev*cc_ev) ]   ! u        
             lvec(3,:) = [ HALF*p_r/(cc_ev*cc_ev),     -HALF*tau_ev/cc_ev,  -HALF*tau_ev**2*p_T/(cc_ev*cc_ev) ]   ! u + c    

             ! compute the right eigenvectors                                            
             rvec(1,:) = [ ONE,  cc_ev/rho_ev,  -(1.0d0/tau_ev**2)*(cc_ev*cc_ev - p_r)/p_T ]   ! u - c 
             rvec(2,:) = [ ONE,  ZERO,     (1.0d0/tau_ev**2)*p_r/p_T  ]   ! u     
             rvec(3,:) = [ ONE,  -cc_ev/rho_ev, -(1.0d0/tau_ev**2)*(cc_ev*cc_ev - p_r)/p_T ]   ! u + c   
             

             ! construct interface states of T, u, tau

             ! loop over waves and construct l . (qref - I(q))
             do iwave = 1, 3
                select case (iwave)
                case (1)
                   dq = [dtaum, dum, dTm]
                case (2)
                   dq = [dtau,  du,  dT0]
                case (3)
                   dq = [dtaup, dup, dTp]
                end select

                beta(iwave) = dot_product(lvec(iwave,:),dq(:))
             enddo

             if (i .le. ihi1) then

                ! here we compute 
                !
                ! q_{l,r} = q_ref - sum { l . (q_ref - I(q)) r } 
                !
                ! where the sum is over the waves and
                ! limited to only include those moving toward the
                ! interface
                
                xi1 = 1.0d0-flatn(i,j)
                xi = flatn(i,j)

                ! tau (and density) state
                tau_s = tau_ref 
                do iwave = 1, 3
                   if (eval(iwave) >= 0.0d0) tau_s = tau_s - beta(iwave)*rvec(iwave,1)
                enddo
                qxm(i+1,j,QRHO) = xi1*rho + xi/tau_s


                ! u state
                qxm(i+1,j,QU)    = xi1*u + xi*u_ref 
                do iwave = 1, 3
                   if (eval(iwave) >= 0.0d0) qxm(i+1,j,QU) = qxm(i+1,j,QU) - xi*beta(iwave)*rvec(iwave,2)
                enddo

                ! T state
                qxm(i+1,j,QTEMP) = xi1*T + xi*T_ref 
                do iwave = 1, 3
                   if (eval(iwave) > 0.0d0) qxm(i+1,j,QTEMP) = qxm(i+1,j,QTEMP) - xi*beta(iwave)*rvec(iwave,3)
                enddo

                qxm(i+1,j,QTEMP) = max(qxm(i+1,j,QTEMP), small_temp)

                ! transverse velocity
                if (u > 0.0d0) then
                   qxm(i+1,j,QV) = -dv
                else
                   qxm(i+1,j,QV) = 0.0d0
                endif

                qxm(i+1,j,QV) = qxm(i+1,j,QV) + xi1*v + xi*v_ref 

                ! limit
                qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)

             endif

          endif


          !-------------------------------------------------------------------
          ! geometry source terms
          !-------------------------------------------------------------------

          if(dloga(i,j).ne.0)then
             courn = dtdx*(cc+abs(u))
             eta = (1.d0-courn)/(cc*dt*abs(dloga(i,j)))
             dlogatmp = min(eta,1.d0)*dloga(i,j)
             sourcr = -0.5d0*dt*rho*dlogatmp*u
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

    do iadv = 1, nadv
       n = QFA + iadv - 1
       do j = ilo2-1, ihi2+1
          
          ! plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,QU)
             if (u .gt. 0.d0) then
                qxp(i,j,n) = q(i,j,n)
             else if (u .lt. 0.d0) then
                qxp(i,j,n) = q(i,j,n) &
                     + flatn(i,j)*(Im(i,j,1,2,n) - q(i,j,n))
             else
                qxp(i,j,n) = q(i,j,n) &
                     + 0.5d0*flatn(i,j)*(Im(i,j,1,2,n) - q(i,j,n))
             endif
          enddo
          
          ! minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,QU)
             if (u .gt. 0.d0) then
                qxm(i+1,j,n) = q(i,j,n) &
                     + flatn(i,j)*(Ip(i,j,1,2,n) - q(i,j,n))
             else if (u .lt. 0.d0) then
                qxm(i+1,j,n) = q(i,j,n)
             else
                qxm(i+1,j,n) = q(i,j,n) &
                     + 0.5d0*flatn(i,j)*(Ip(i,j,1,2,n) - q(i,j,n))
             endif
          enddo
          
       enddo
    enddo


    ! species

    do ispec = 1, nspec
       ns = QFS + ispec - 1
       do j = ilo2-1, ihi2+1
          
          ! plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,QU)
             if (u .gt. 0.d0) then
                qxp(i,j,ns) = q(i,j,ns)
             else if (u .lt. 0.d0) then
                qxp(i,j,ns) = q(i,j,ns) &
                     + flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
             else
                qxp(i,j,ns) = q(i,j,ns) &
                     + 0.5d0*flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
             endif
          enddo
          
          ! minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,QU)
             if (u .gt. 0.d0) then
                qxm(i+1,j,ns) = q(i,j,ns) &
                     + flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
             else if (u .lt. 0.d0) then
                qxm(i+1,j,ns) = q(i,j,ns)
             else
                qxm(i+1,j,ns) = q(i,j,ns) &
                     + 0.5d0*flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
             endif
          enddo
          
       enddo
    enddo

    ! auxillary quantities

    do iaux = 1, naux
       ns = QFX + iaux - 1
       do j = ilo2-1, ihi2+1
          
          ! plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,QU)
             if (u .gt. 0.d0) then
                qxp(i,j,ns) = q(i,j,ns)
             else if (u .lt. 0.d0) then
                qxp(i,j,ns) = q(i,j,ns) &
                     + flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
             else
                qxp(i,j,ns) = q(i,j,ns) &
                     + 0.5d0*flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
             endif
          enddo

          ! minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,QU)
             if (u .gt. 0.d0) then
                qxm(i+1,j,ns) = q(i,j,ns) &
                     + flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
             else if (u .lt. 0.d0) then
                qxm(i+1,j,ns) = q(i,j,ns)
             else
                qxm(i+1,j,ns) = q(i,j,ns) &
                     + 0.5d0*flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
             endif
          enddo
          
       enddo
    enddo


    ! convert from T to p if we did the characteristic tracing in
    ! terms of T
    if (ppm_temp_fix == 3) then

       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             if (i .ge. ilo1) then
                ! now get the pressure and energy state via the EOS
                
                ! plus face
                eos_state%T     = qxp(i,j,QTEMP)
                eos_state%rho   = qxp(i,j,QRHO)
                eos_state%xn(:) = qxp(i,j,QFS:QFS-1+nspec)
                
                call eos(eos_input_rt, eos_state, .false.)
             
                qxp(i,j,QPRES) = eos_state%p
                qxp(i,j,QREINT) = qxp(i,j,QRHO)*eos_state%e
                
                qxp(i,j,QPRES) = max(qxp(i,j,QPRES), small_pres)

             endif

             if (i .le. ihi1) then
             
                ! minus face
                eos_state%T     = qxm(i+1,j,QTEMP)
                eos_state%rho   = qxm(i+1,j,QRHO)
                eos_state%xn(:) = qxm(i+1,j,QFS:QFS-1+nspec)

                call eos(eos_input_rt, eos_state, .false.)
                
                qxm(i+1,j,QPRES) = eos_state%p
                qxm(i+1,j,QREINT) = qxm(i+1,j,QRHO)*eos_state%e

                qxm(i+1,j,QPRES) = max(qxm(i+1,j,QPRES), small_pres)                
             endif

          enddo
       enddo

    endif

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

          !-------------------------------------------------------------------
          ! plus state on face j
          !-------------------------------------------------------------------

          ! set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. v - cc >= 0.0d0)) then
             ! original Castro way -- cc value
             rho_ref  = rho
             u_ref    = u
             v_ref    = v

             p_ref    = p
             rhoe_ref = rhoe

             tau_ref  = 1.0d0/rho

             gam_ref = gamc(i,j)
          else
             ! this will be the fastest moving state to the left
             rho_ref  = Im(i,j,2,1,QRHO)
             u_ref    = Im(i,j,2,1,QU)
             v_ref    = Im(i,j,2,1,QV)

             p_ref    = Im(i,j,2,1,QPRES)
             rhoe_ref = Im(i,j,2,1,QREINT)

             tau_ref  = 1.0d0/Im(i,j,2,1,QRHO)

             gam_ref  = Im_gc(i,j,2,1,1)
          endif

          ! for tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref
            
          ! *m are the jumps carried by v-c
          ! *p are the jumps carried by v+c

          !dum    = (u_ref    - Im(i,j,2,1,QU))
          dvm    = (v_ref    - Im(i,j,2,1,QV))
          dpm    = (p_ref    - Im(i,j,2,1,QPRES))
          drhoem = (rhoe_ref - Im(i,j,2,1,QREINT))
          
          drho  = (rho_ref  - Im(i,j,2,2,QRHO))
          du    = (u_ref    - Im(i,j,2,2,QU))
          !dv    = (v_ref    - Im(i,j,2,2,QV))
          dp    = (p_ref    - Im(i,j,2,2,QPRES))
          drhoe = (rhoe_ref - Im(i,j,2,2,QREINT))
          
          !dup    = (u_ref    - Im(i,j,2,3,QU))
          dvp    = (v_ref    - Im(i,j,2,3,QV))
          dpp    = (p_ref    - Im(i,j,2,3,QPRES))
          drhoep = (rhoe_ref - Im(i,j,2,3,QREINT))

          ! note: we do not implement temp_fix = 3 in the y-direction
          dtau  = (tau_ref  - 1.0d0/Im(i,j,2,2,QRHO))

          ! if we are doing gravity tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the 
          ! trans_X routines
          if (ppm_trace_grav == 1) then
             dvm = dvm - halfdt*Im_g(i,j,2,1,igy)
             du  = du  - halfdt*Im_g(i,j,2,2,igx)
             dvp = dvp - halfdt*Im_g(i,j,2,3,igy)
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
             alpha0e = drhoe - dp*enth_ev  ! enth has 1/c**2 in it
             alpha0u = du
          
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
             else if (v .lt. 0.d0) then
                azrright = -alpha0r
                azeright = -alpha0e
                azu1rght = -alpha0u
             else
                azrright = -0.5d0*alpha0r
                azeright = -0.5d0*alpha0e
                azu1rght = -0.5d0*alpha0u
             endif
             
             ! the final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             if (j .ge. ilo2) then
                xi1 = 1.0d0 - flatn(i,j)
                xi = flatn(i,j)
                
                qyp(i,j,QRHO)   = xi1*rho  + xi*(rho_ref + apright + amright + azrright)
                qyp(i,j,QV)     = xi1*v    + xi*(v_ref + (apright - amright)*cc_ev/rho_ev)
                qyp(i,j,QU)     = xi1*u    + xi*(u_ref + azu1rght)
                
                qyp(i,j,QREINT) = xi1*rhoe + xi*(rhoe_ref + (apright + amright)*enth_ev*csq_ev + azeright)
                qyp(i,j,QPRES)  = xi1*p    + xi*(p_ref + (apright + amright)*csq_ev)
                
                qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))
                qyp(i,j,QPRES) = max(qyp(i,j,QPRES), small_pres)
             end if
             
          else
             ! (tau, u, p, e) eigensystem

             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system.

             de = (rhoe_ref/rho_ref - Im(i,j,2,2,QREINT)/Im(i,j,2,2,QRHO))

             alpham = 0.5d0*( dvm - dpm/Clag_ev)/Clag_ev
             alphap = 0.5d0*(-dvp - dpp/Clag_ev)/Clag_ev
             alpha0r = dtau + dp/Clag_ev**2
             alpha0e = de - dp*p_ev/Clag_ev**2
             alpha0u = du
          
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
             else if (v .lt. 0.d0) then
                azrright = -alpha0r
                azeright = -alpha0e
                azu1rght = -alpha0u
             else
                azrright = -0.5d0*alpha0r
                azeright = -0.5d0*alpha0e
                azu1rght = -0.5d0*alpha0u
             endif

             ! the final interface states are just
             ! q_s = q_ref - sum (l . dq) r          
             if (j .ge. ilo2) then

                xi1 = 1.0d0 - flatn(i,j)
                xi = flatn(i,j)
                
                tau_s = tau_ref + apright + amright + azrright
                qyp(i,j,QRHO)   = xi1*rho  + xi/tau_s

                qyp(i,j,QV)     = xi1*v    + xi*(v_ref + (amright - apright)*Clag_ev)
                qyp(i,j,QU)     = xi1*u    + xi*(u_ref + azu1rght)

                e_s = rhoe_ref/rho_ref + (azeright - p_ev*amright -p_ev*apright)
                qyp(i,j,QREINT) = xi1*rhoe + xi*e_s/tau_s

                qyp(i,j,QPRES)  = xi1*p    + xi*(p_ref + (-apright - amright)*Clag_ev**2)
                
                qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))
                qyp(i,j,QPRES) = max(qyp(i,j,QPRES), small_pres)

             end if             

          endif

          !-------------------------------------------------------------------
          ! minus state on face j+1
          !-------------------------------------------------------------------

          ! set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. v + cc <= 0.0d0) ) then
             ! original Castro way -- cc value
             rho_ref  = rho
             u_ref    = u
             v_ref    = v

             p_ref    = p
             rhoe_ref = rhoe

             tau_ref  = 1.0d0/rho

             gam_ref = gamc(i,j)

          else
             ! this will be the fastest moving state to the right
             rho_ref  = Ip(i,j,2,3,QRHO)
             u_ref    = Ip(i,j,2,3,QU)
             v_ref    = Ip(i,j,2,3,QV)

             p_ref    = Ip(i,j,2,3,QPRES)
             rhoe_ref = Ip(i,j,2,3,QREINT)

             tau_ref  = 1.0d0/Ip(i,j,2,3,QRHO)

             gam_ref  = Ip_gc(i,j,2,3,1)
          endif

          ! for tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = ( (rhoe_ref+p_ref)/rho_ref )/csq_ref

          ! *m are the jumps carried by v-c
          ! *p are the jumps carried by v+c

          !dum    = (u_ref    - Ip(i,j,2,1,QU))
          dvm    = (v_ref    - Ip(i,j,2,1,QV))
          dpm    = (p_ref    - Ip(i,j,2,1,QPRES))
          drhoem = (rhoe_ref - Ip(i,j,2,1,QREINT))
          
          drho  = (rho_ref  - Ip(i,j,2,2,QRHO))
          du    = (u_ref    - Ip(i,j,2,2,QU))
          !dv    = (v_ref    - Ip(i,j,2,2,QV))
          dp    = (p_ref    - Ip(i,j,2,2,QPRES))
          drhoe = (rhoe_ref - Ip(i,j,2,2,QREINT))
          
          !dup    = (u_ref    - Ip(i,j,2,3,QU))
          dvp    = (v_ref    - Ip(i,j,2,3,QV))
          dpp    = (p_ref    - Ip(i,j,2,3,QPRES))
          drhoep = (rhoe_ref - Ip(i,j,2,3,QREINT))

          ! we are not implementing ppm_temp_fix = 3 in the y-direction
          dtau  = (tau_ref  - 1.0d0/Ip(i,j,2,2,QRHO))

          ! if we are doing gravity tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the 
          ! trans_X routines
          if (ppm_trace_grav == 1) then
             dvm = dvm - halfdt*Ip_g(i,j,2,1,igy)
             du  = du  - halfdt*Ip_g(i,j,2,2,igx)
             dvp = dvp - halfdt*Ip_g(i,j,2,3,igy)
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
             else if (v .lt. 0.d0) then
                azrleft = 0.d0
                azeleft = 0.d0
                azu1left = 0.d0
             else
                azrleft = -0.5d0*alpha0r
                azeleft = -0.5d0*alpha0e
                azu1left = -0.5d0*alpha0u
             endif

             ! the final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             if (j .le. ihi2) then
                xi1 = 1.0d0 - flatn(i,j)
                xi = flatn(i,j)
                
                qym(i,j+1,QRHO)   = xi1*rho  + xi*(rho_ref + apleft + amleft + azrleft)
                qym(i,j+1,QV)     = xi1*v    + xi*(v_ref + (apleft - amleft)*cc_ev/rho_ev)
                qym(i,j+1,QU)     = xi1*u    + xi*(u_ref + azu1left)
                
                qym(i,j+1,QREINT) = xi1*rhoe + xi*(rhoe_ref + (apleft + amleft)*enth_ev*csq_ev + azeleft)
                qym(i,j+1,QPRES)  = xi1*p    + xi*(p_ref + (apleft + amleft)*csq_ev)
                
                qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))
                qym(i,j+1,QPRES) = max(qym(i,j+1,QPRES), small_pres)
             end if

          else
             ! (tau, u, p, e) eigensystem

             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system.

             de = (rhoe_ref/rho_ref - Ip(i,j,2,2,QREINT)/Ip(i,j,2,2,QRHO))

             alpham = 0.5d0*( dvm - dpm/Clag_ev)/Clag_ev
             alphap = 0.5d0*(-dvp - dpp/Clag_ev)/Clag_ev
             alpha0r = dtau + dp/Clag_ev**2
             alpha0e = de - dp*p_ev/Clag_ev**2
             alpha0u = du
             
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
             else if (v .lt. 0.d0) then
                azrleft = 0.d0
                azeleft = 0.d0
                azu1left = 0.d0
             else
                azrleft = -0.5d0*alpha0r
                azeleft = -0.5d0*alpha0e
                azu1left = -0.5d0*alpha0u
             endif
             
             ! the final interface states are just
             ! q_s = q_ref - sum (l . dq) r          
             if (j .le. ihi2) then

                xi1 = 1.0d0 - flatn(i,j)
                xi = flatn(i,j)
                
                tau_s = tau_ref + apleft + amleft + azrleft
                qym(i,j+1,QRHO)   = xi1*rho  + xi/tau_s

                qym(i,j+1,QV)     = xi1*v    + xi*(v_ref + (amleft - apleft)*Clag_ev)
                qym(i,j+1,QU)     = xi1*u    + xi*(u_ref + azu1left)
             
                e_s = rhoe_ref/rho_ref + (azeleft - p_ev*amleft -p_ev*apleft)
                qym(i,j+1,QREINT) = xi1*rhoe + xi*e_s/tau_s                

                qym(i,j+1,QPRES)  = xi1*p    + xi*(p_ref + (-apleft - amleft)*Clag_ev**2)
                
                qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))
                qym(i,j+1,QPRES) = max(qym(i,j+1,QPRES), small_pres)

             end if

          endif

       end do
    end do


    !-------------------------------------------------------------------------
    ! Now do the passively advected quantities
    !-------------------------------------------------------------------------

    do iadv = 1, nadv
       n = QFA + iadv - 1
       do i = ilo1-1, ihi1+1
          
          ! plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,QV)
             if (v .gt. 0.d0) then
                qyp(i,j,n) = q(i,j,n)
             else if (v .lt. 0.d0) then
                qyp(i,j,n) = q(i,j,n) &
                     + flatn(i,j)*(Im(i,j,2,2,n) - q(i,j,n))
             else
                qyp(i,j,n) = q(i,j,n) &
                     + 0.5d0*flatn(i,j)*(Im(i,j,2,2,n) - q(i,j,n))
             endif
          enddo
          
          ! minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,QV)
             if (v .gt. 0.d0) then
                qym(i,j+1,n) = q(i,j,n) &
                     + flatn(i,j)*(Ip(i,j,2,2,n) - q(i,j,n))
             else if (v .lt. 0.d0) then
                qym(i,j+1,n) = q(i,j,n)
             else
                qym(i,j+1,n) = q(i,j,n) &
                     + 0.5d0*flatn(i,j)*(Ip(i,j,2,2,n) - q(i,j,n))
             endif
          enddo
          
       enddo
    enddo


    ! species

    do ispec = 1, nspec
       ns = QFS + ispec - 1
       do i = ilo1-1, ihi1+1
          
          ! plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,QV)
             if (v .gt. 0.d0) then
                qyp(i,j,ns) = q(i,j,ns)
             else if (v .lt. 0.d0) then
                qyp(i,j,ns) = q(i,j,ns) &
                     + flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
             else
                qyp(i,j,ns) = q(i,j,ns) &
                     + 0.5d0*flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
             endif
          enddo
          
          ! minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,QV)
             if (v .gt. 0.d0) then
                qym(i,j+1,ns) = q(i,j,ns) &
                     + flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
             else if (v .lt. 0.d0) then
                qym(i,j+1,ns) = q(i,j,ns)
             else
                qym(i,j+1,ns) = q(i,j,ns) &
                     + 0.5d0*flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
             endif
          enddo
          
       enddo
    enddo

    
    ! auxillary quantities
    
    do iaux = 1, naux
       ns = QFX + iaux - 1
       do i = ilo1-1, ihi1+1
          
          ! plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,QV)
             if (v .gt. 0.d0) then
                qyp(i,j,ns) = q(i,j,ns)
             else if (v .lt. 0.d0) then
                qyp(i,j,ns) = q(i,j,ns) &
                     + flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
             else
                qyp(i,j,ns) = q(i,j,ns) &
                     + 0.5d0*flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
             endif
          enddo
          
          ! minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,QV)
             if (v .gt. 0.d0) then
                qym(i,j+1,ns) = q(i,j,ns) &
                     + flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
             else if (v .lt. 0.d0) then
                qym(i,j+1,ns) = q(i,j,ns)
             else
                qym(i,j+1,ns) = q(i,j,ns) &
                     + 0.5d0*flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
             endif
          enddo
          
       enddo
    enddo
    
    deallocate(Ip,Im)
    if (ppm_trace_grav == 1) then
       deallocate(Ip_g,Im_g)
    endif

  end subroutine trace_ppm
  
end module trace_ppm_module
