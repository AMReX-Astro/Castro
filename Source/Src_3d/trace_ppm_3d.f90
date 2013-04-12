module trace_ppm_module

  implicit none

  private

  public tracexy_ppm, tracez_ppm

contains

  subroutine tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         Ip,Im, &
                         qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QESGS, QPRES, QFA, QFS, nadv, &
         ppm_type, ppm_reference, small_dens, small_pres

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer kc,k3d

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    double precision   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision dx, dy, dt

    ! Local variables
    integer i, j
    integer n, iadv, ispec
    integer npassive,ipassive,qpass_map(QVAR)

    double precision cc, csq, rho, u, v, w, p, rhoe

    double precision drho, du, dv, dw, dp, drhoe
    double precision dup, dvp, dpp
    double precision dum, dvm, dpm

    double precision enth, alpham, alphap, alpha0r, alpha0e
    double precision alpha0u, alpha0v, alpha0w
    double precision apright, amright, azrright, azeright
    double precision azu1rght, azv1rght, azw1rght
    double precision apleft, amleft, azrleft, azeleft
    double precision azu1left, azv1left, azw1left
    double precision rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref
    double precision xi, xi1

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

    !!!!!!!!!!!!!!!
    ! PPM CODE
    !!!!!!!!!!!!!!!

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).

    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.


    ! *********************************************************************************************
    ! x-direction
    ! *********************************************************************************************

    ! Trace to left and right edges using upwind PPM
    !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth) &
    !$OMP PRIVATE(rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref) &
    !$OMP PRIVATE(drho,dv,dw,dp,drhoe,dum,dpm,dup,dpp,alpham,alphap,alpha0r) &
    !$OMP PRIVATE(alpha0e,alpha0v,alpha0w,amright,apright,azrright,azeright,azv1rght,azw1rght) &
    !$OMP PRIVATE(amleft,apleft,azrleft,azeleft,azv1left,azw1left)
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)

          cc = c(i,j,k3d)
          csq = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! ******************************************************************************

          if (i .ge. ilo1) then

             ! Plus state on face i

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. u - cc >= 0.0d0) ) then
                rho_ref  = rho
                u_ref    = u
                v_ref    = v
                w_ref    = w
                p_ref    = p
                rhoe_ref = rhoe
   
             else
                ! This will be the fastest moving state to the left
                rho_ref  = Im(i,j,kc,1,1,QRHO)
                u_ref    = Im(i,j,kc,1,1,QU)
                v_ref    = Im(i,j,kc,1,1,QV)
                w_ref    = Im(i,j,kc,1,1,QW)
                p_ref    = Im(i,j,kc,1,1,QPRES)
                rhoe_ref = Im(i,j,kc,1,1,QREINT)
             endif
   
             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)
   
             dum    = flatn(i,j,k3d)*(u_ref    - Im(i,j,kc,1,1,QU))
             dpm    = flatn(i,j,k3d)*(p_ref    - Im(i,j,kc,1,1,QPRES))
   
             drho  = flatn(i,j,k3d)*(rho_ref  - Im(i,j,kc,1,2,QRHO))
             dv    = flatn(i,j,k3d)*(v_ref    - Im(i,j,kc,1,2,QV))
             dw    = flatn(i,j,k3d)*(w_ref    - Im(i,j,kc,1,2,QW))
             dp    = flatn(i,j,k3d)*(p_ref    - Im(i,j,kc,1,2,QPRES))
             drhoe = flatn(i,j,k3d)*(rhoe_ref - Im(i,j,kc,1,2,QREINT))

             dup    = flatn(i,j,k3d)*(u_ref    - Im(i,j,kc,1,3,QU))
             dpp    = flatn(i,j,k3d)*(p_ref    - Im(i,j,kc,1,3,QPRES))
   
             ! these are the beta's from the original PPM paper.  This is essentially
             ! (l . dq) r for each primitive quantity
   
             alpham = 0.5d0*(dpm/(rho*cc) - dum)*rho/cc
             alphap = 0.5d0*(dpp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
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

             xi1 = 1.0d0-flatn(i,j,k3d)
             xi = flatn(i,j,k3d)
             qxp(i,j,kc,QRHO  ) = xi1*rho  + xi* rho_ref + apright + amright + azrright
             qxp(i,j,kc,QU    ) = xi1*u    + xi*   u_ref + (apright - amright)*cc/rho
             qxp(i,j,kc,QV    ) = xi1*v    + xi*   v_ref + azv1rght
             qxp(i,j,kc,QW    ) = xi1*w    + xi*   w_ref + azw1rght
             qxp(i,j,kc,QREINT) = xi1*rhoe + xi*rhoe_ref + (apright + amright)*enth*csq + azeright
             qxp(i,j,kc,QPRES ) = xi1*p    + xi*   p_ref + (apright + amright)*csq

             qxp(i,j,kc,QRHO ) = max(qxp(i,j,kc,QRHO ),small_dens)
             qxp(i,j,kc,QPRES) = max(qxp(i,j,kc,QPRES),small_pres)
          end if

          ! ******************************************************************************

          if (i .le. ihi1) then

             ! Minus state on face i+1

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. u + cc <= 0.0d0) ) then
                rho_ref  = rho
                u_ref    = u
                v_ref    = v
                w_ref    = w
                p_ref    = p
                rhoe_ref = rhoe
             else
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,kc,1,3,QRHO)
                u_ref    = Ip(i,j,kc,1,3,QU)
                v_ref    = Ip(i,j,kc,1,3,QV)
                w_ref    = Ip(i,j,kc,1,3,QW)
                p_ref    = Ip(i,j,kc,1,3,QPRES)
             rhoe_ref = Ip(i,j,kc,1,3,QREINT)
             endif
   
             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)
   
             dum    = flatn(i,j,k3d)*(u_ref    - Ip(i,j,kc,1,1,QU))
             dpm    = flatn(i,j,k3d)*(p_ref    - Ip(i,j,kc,1,1,QPRES))
   
             drho  = flatn(i,j,k3d)*(rho_ref  - Ip(i,j,kc,1,2,QRHO))
             dv    = flatn(i,j,k3d)*(v_ref    - Ip(i,j,kc,1,2,QV))
             dw    = flatn(i,j,k3d)*(w_ref    - Ip(i,j,kc,1,2,QW))
             dp    = flatn(i,j,k3d)*(p_ref    - Ip(i,j,kc,1,2,QPRES))
             drhoe = flatn(i,j,k3d)*(rhoe_ref - Ip(i,j,kc,1,2,QREINT))

             dup    = flatn(i,j,k3d)*(u_ref    - Ip(i,j,kc,1,3,QU))
             dpp    = flatn(i,j,k3d)*(p_ref    - Ip(i,j,kc,1,3,QPRES))

             alpham = 0.5d0*(dpm/(rho*cc) - dum)*rho/cc
             alphap = 0.5d0*(dpp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
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

             xi1 = 1.0d0 - flatn(i,j,k3d)
             xi = flatn(i,j,k3d)
             qxm(i+1,j,kc,QRHO  ) = xi1*rho  + xi* rho_ref + apleft + amleft + azrleft
             qxm(i+1,j,kc,QU    ) = xi1*u    + xi*   u_ref + (apleft - amleft)*cc/rho
             qxm(i+1,j,kc,QV    ) = xi1*v    + xi*   v_ref + azv1left
             qxm(i+1,j,kc,QW    ) = xi1*w    + xi*   w_ref + azw1left
             qxm(i+1,j,kc,QREINT) = xi1*rhoe + xi*rhoe_ref + (apleft + amleft)*enth*csq + azeleft
             qxm(i+1,j,kc,QPRES ) = xi1*p    + xi*   p_ref + (apleft + amleft)*csq

             qxm(i+1,j,kc,QRHO  ) = max(qxm(i+1,j,kc,QRHO ),small_dens)
             qxm(i+1,j,kc,QPRES)  = max(qxm(i+1,j,kc,QPRES),small_pres)
          end if

       end do
    end do
    !$OMP END PARALLEL DO

    ! Do all of the passively advected quantities in one loop
    !$OMP parallel do private(ipassive,i,j,u,n) IF(npassive .gt. 1)
    do ipassive = 1, npassive
         n = qpass_map(ipassive)
         do j = ilo2-1, ihi2+1

               ! Plus state on face i
               do i = ilo1, ihi1+1
                  u = q(i,j,k3d,QU)
                  if (u .gt. 0.d0) then
                     qxp(i,j,kc,n) = q(i,j,k3d,n)
                  else if (u .lt. 0.d0) then
                     qxp(i,j,kc,n) = q(i,j,k3d,n) &
                          + flatn(i,j,k3d)*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
                  else
                     qxp(i,j,kc,n) = q(i,j,k3d,n) &
                          + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
                  endif
               enddo

               ! Minus state on face i+1
               do i = ilo1-1, ihi1
                  u = q(i,j,k3d,QU)
                  if (u .gt. 0.d0) then
                     qxm(i+1,j,kc,n) = q(i,j,k3d,n) &
                          + flatn(i,j,k3d)*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
                  else if (u .lt. 0.d0) then
                     qxm(i+1,j,kc,n) = q(i,j,k3d,n)
                  else
                     qxm(i+1,j,kc,n) = q(i,j,k3d,n) &
                          + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
                  endif
               enddo
        enddo
    enddo
    !$OMP end parallel do

    ! *********************************************************************************************
    ! y-direction
    ! *********************************************************************************************

    ! Trace to bottom and top edges using upwind PPM
    !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth) &
    !$OMP PRIVATE(rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref) &
    !$OMP PRIVATE(drho,du,dw,dp,drhoe,dvm,dpm,dvp,dpp,alpham,alphap,alpha0r) &
    !$OMP PRIVATE(alpha0e,alpha0u,alpha0w,amright,apright,azrright,azeright,azu1rght,azw1rght,amleft) &
    !$OMP PRIVATE(apleft,azrleft,azeleft,azu1left,azw1left)
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)

          cc = c(i,j,k3d)
          csq = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! ******************************************************************************

          if (j .ge. ilo2) then

             ! Plus state on face j

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. v - cc >= 0.0d0) ) then
                rho_ref  = rho
                u_ref    = u
                v_ref    = v
                w_ref    = w
                p_ref    = p
                rhoe_ref = rhoe
             else
                ! This will be the fastest moving state to the left
                rho_ref  = Im(i,j,kc,2,1,QRHO)
                u_ref    = Im(i,j,kc,2,1,QU)
                v_ref    = Im(i,j,kc,2,1,QV)
                w_ref    = Im(i,j,kc,2,1,QW)
                p_ref    = Im(i,j,kc,2,1,QPRES)
                rhoe_ref = Im(i,j,kc,2,1,QREINT)
             endif
   
             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm    = flatn(i,j,k3d)*(v_ref    - Im(i,j,kc,2,1,QV))
             dpm    = flatn(i,j,k3d)*(p_ref    - Im(i,j,kc,2,1,QPRES))
   
             drho  = flatn(i,j,k3d)*(rho_ref  - Im(i,j,kc,2,2,QRHO))
             du    = flatn(i,j,k3d)*(u_ref    - Im(i,j,kc,2,2,QU))
             dw    = flatn(i,j,k3d)*(w_ref    - Im(i,j,kc,2,2,QW))
             dp    = flatn(i,j,k3d)*(p_ref    - Im(i,j,kc,2,2,QPRES))
             drhoe = flatn(i,j,k3d)*(rhoe_ref - Im(i,j,kc,2,2,QREINT))

             dvp    = flatn(i,j,k3d)*(v_ref    - Im(i,j,kc,2,3,QV))
             dpp    = flatn(i,j,k3d)*(p_ref    - Im(i,j,kc,2,3,QPRES))
   
             alpham = 0.5d0*(dpm/(rho*cc) - dvm)*rho/cc
             alphap = 0.5d0*(dpp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
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

             xi1 = 1.0d0 - flatn(i,j,k3d)
             xi = flatn(i,j,k3d)
             qyp(i,j,kc,QRHO  ) = xi1*rho  + xi*rho_ref + apright + amright + azrright
             qyp(i,j,kc,QV    ) = xi1*v    + xi*  v_ref + (apright - amright)*cc/rho
             qyp(i,j,kc,QU    ) = xi1*u    + xi*  u_ref + azu1rght
             qyp(i,j,kc,QW    ) = xi1*w    + xi*  w_ref + azw1rght
             qyp(i,j,kc,QREINT) = xi1*rhoe + xi*rhoe_ref + (apright + amright)*enth*csq + azeright
             qyp(i,j,kc,QPRES ) = xi1*p    + xi*  p_ref + (apright + amright)*csq

             qyp(i,j,kc,QRHO ) = max(qyp(i,j,kc,QRHO ),small_dens)
             qyp(i,j,kc,QPRES) = max(qyp(i,j,kc,QPRES),small_pres)
          end if

          ! ******************************************************************************

          if (j .le. ihi2) then

             ! Minus state on face j+1 

             ! Set the reference state 
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. v + cc <= 0.0d0) ) then
                rho_ref  = rho
                u_ref    = u
                v_ref    = v
                w_ref    = w
                p_ref    = p
                rhoe_ref = rhoe
             else
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,kc,2,3,QRHO)
                u_ref    = Ip(i,j,kc,2,3,QU)
                v_ref    = Ip(i,j,kc,2,3,QV)
                w_ref    = Ip(i,j,kc,2,3,QW)
                p_ref    = Ip(i,j,kc,2,3,QPRES)
                rhoe_ref = Ip(i,j,kc,2,3,QREINT)
             endif

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c
   
             ! Note: for the transverse velocities, the jump is carried
             !       only by the v wave (the contact)

             dvm    = flatn(i,j,k3d)*(v_ref    - Ip(i,j,kc,2,1,QV))
             dpm    = flatn(i,j,k3d)*(p_ref    - Ip(i,j,kc,2,1,QPRES))
   
             drho  = flatn(i,j,k3d)*(rho_ref  - Ip(i,j,kc,2,2,QRHO))
             du    = flatn(i,j,k3d)*(u_ref    - Ip(i,j,kc,2,2,QU))
             dw    = flatn(i,j,k3d)*(w_ref    - Ip(i,j,kc,2,2,QW))
             dp    = flatn(i,j,k3d)*(p_ref    - Ip(i,j,kc,2,2,QPRES))
             drhoe = flatn(i,j,k3d)*(rhoe_ref - Ip(i,j,kc,2,2,QREINT))

             dvp    = flatn(i,j,k3d)*(v_ref    - Ip(i,j,kc,2,3,QV))
             dpp    = flatn(i,j,k3d)*(p_ref    - Ip(i,j,kc,2,3,QPRES))

             alpham = 0.5d0*(dpm/(rho*cc) - dvm)*rho/cc
             alphap = 0.5d0*(dpp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
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

             xi1 = 1.0d0 - flatn(i,j,k3d)
             xi = flatn(i,j,k3d)
             qym(i,j+1,kc,QRHO  ) = xi1*rho  + xi* rho_ref + apleft + amleft + azrleft
             qym(i,j+1,kc,QV    ) = xi1*v    + xi*   v_ref + (apleft - amleft)*cc/rho
             qym(i,j+1,kc,QU    ) = xi1*u    + xi*   u_ref + azu1left
             qym(i,j+1,kc,QW    ) = xi1*w    + xi*   w_ref + azw1left
             qym(i,j+1,kc,QREINT) = xi1*rhoe + xi*rhoe_ref + (apleft + amleft)*enth*csq + azeleft
             qym(i,j+1,kc,QPRES ) = xi1*p    + xi*   p_ref + (apleft + amleft)*csq

             qym(i,j+1,kc,QRHO ) = max(qym(i,j+1,kc,QRHO ),small_dens)
             qym(i,j+1,kc,QPRES) = max(qym(i,j+1,kc,QPRES),small_pres)
          end if

       end do
    end do
    !$OMP END PARALLEL DO

    ! Do all of the passively advected quantities in one loop
    !$OMP parallel do private(n,i,j,v,ipassive) IF(npassive .gt. 1)
    do ipassive = 1, npassive
         n = qpass_map(ipassive)
         do i = ilo1-1, ihi1+1

               ! Plus state on face j
               do j = ilo2, ihi2+1
                  v = q(i,j,k3d,QV)
                  if (v .gt. 0.d0) then
                     qyp(i,j,kc,n) = q(i,j,k3d,n)
                  else if (v .lt. 0.d0) then
                     qyp(i,j,kc,n) = q(i,j,k3d,n) &
                          + flatn(i,j,k3d)*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
                  else
                     qyp(i,j,kc,n) = q(i,j,k3d,n) &
                          + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
                  endif
               enddo

               ! Minus state on face j+1
               do j = ilo2-1, ihi2
                  v = q(i,j,k3d,QV)
                  if (v .gt. 0.d0) then
                     qym(i,j+1,kc,n) = q(i,j,k3d,n) &
                          + flatn(i,j,k3d)*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
                  else if (v .lt. 0.d0) then
                     qym(i,j+1,kc,n) = q(i,j,k3d,n)
                  else
                     qym(i,j+1,kc,n) = q(i,j,k3d,n) &
                          + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
                  endif
               enddo

         enddo
    enddo

  end subroutine tracexy_ppm

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        Ip,Im, &
                        qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QESGS, QPRES, QFA, QFS, nadv, &
         ppm_type, ppm_reference, small_dens, small_pres

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer km,kc,k3d

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    double precision   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    double precision qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision dz, dt

    !     Local variables
    integer i, j
    integer n, iadv, ispec

    double precision cc, csq, rho, u, v, w, p, rhoe
    double precision dwp, dpp
    double precision dwm, dpm

    double precision drho, du, dv, dp, drhoe
    double precision enth, alpham, alphap, alpha0r, alpha0e
    double precision alpha0u, alpha0v
    double precision apright, amright, azrright, azeright
    double precision azu1rght, azv1rght
    double precision apleft, amleft, azrleft, azeleft
    double precision azu1left, azv1left
    double precision rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref
    double precision xi, xi1

    integer npassive,ipassive,qpass_map(QVAR)

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
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call bl_error("Error:: trace_ppm_3d.f90 :: tracez_ppm")
    end if

    !!!!!!!!!!!!!!!
    ! PPM CODE
    !!!!!!!!!!!!!!!

    ! Trace to left and right edges using upwind PPM

    ! Note: in contrast to the above code for x and y, here the loop
    ! is over interfaces, not over cell-centers.

    !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth) &
    !$OMP PRIVATE(rho_ref,u_ref,v_ref,w_ref,p_ref,rhoe_ref) &
    !$OMP PRIVATE(drho,du,dv,dp,drhoe,dwm,dpm,dwp,dpp,alpham,alphap,alpha0r,alpha0e) &
    !$OMP PRIVATE(alpha0u,alpha0v,amright,apright,azrright,azeright,azu1rght,azv1rght,amleft,apleft)&
    !$OMP PRIVATE(azrleft,azeleft,azu1left,azv1left)
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          ! **************************************************************************
          ! This is all for qzp
          ! **************************************************************************

          rho  = q(i,j,k3d,QRHO)
          u    = q(i,j,k3d,QU)
          v    = q(i,j,k3d,QV)
          w    = q(i,j,k3d,QW)
          p    = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)

          cc   = c(i,j,k3d)
          csq  = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! Plus state on face kc

          ! Set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. w - cc >= 0.0d0) ) then
             rho_ref  = rho
             u_ref    = u
             v_ref    = v
             w_ref    = w
             p_ref    = p
             rhoe_ref = rhoe
          else
             ! This will be the fastest moving state to the left
             rho_ref  = Im(i,j,kc,3,1,QRHO)
             u_ref    = Im(i,j,kc,3,1,QU)
             v_ref    = Im(i,j,kc,3,1,QV)
             w_ref    = Im(i,j,kc,3,1,QW)
             p_ref    = Im(i,j,kc,3,1,QPRES)
             rhoe_ref = Im(i,j,kc,3,1,QREINT)
          endif

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm    = flatn(i,j,k3d)*(w_ref    - Im(i,j,kc,3,1,QW))
          dpm    = flatn(i,j,k3d)*(p_ref    - Im(i,j,kc,3,1,QPRES))

          drho  = flatn(i,j,k3d)*(rho_ref  - Im(i,j,kc,3,2,QRHO))
          du    = flatn(i,j,k3d)*(u_ref    - Im(i,j,kc,3,2,QU))
          dv    = flatn(i,j,k3d)*(v_ref    - Im(i,j,kc,3,2,QV))
          dp    = flatn(i,j,k3d)*(p_ref    - Im(i,j,kc,3,2,QPRES))
          drhoe = flatn(i,j,k3d)*(rhoe_ref - Im(i,j,kc,3,2,QREINT))

          dwp    = flatn(i,j,k3d)*(w_ref    - Im(i,j,kc,3,3,QW))
          dpp    = flatn(i,j,k3d)*(p_ref    - Im(i,j,kc,3,3,QPRES))

          alpham = 0.5d0*(dpm/(rho*cc) - dwm)*rho/cc
          alphap = 0.5d0*(dpp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
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

          xi1 = 1.0d0 - flatn(i,j,k3d)
          xi = flatn(i,j,k3d)
          qzp(i,j,kc,QRHO  ) = xi1*rho  + xi* rho_ref + apright + amright + azrright
          qzp(i,j,kc,QW    ) = xi1*w    + xi*   w_ref + (apright - amright)*cc/rho
          qzp(i,j,kc,QU    ) = xi1*u    + xi*   u_ref + azu1rght
          qzp(i,j,kc,QV    ) = xi1*v    + xi*   v_ref + azv1rght
          qzp(i,j,kc,QREINT) = xi1*rhoe + xi*rhoe_ref + (apright + amright)*enth*csq + azeright
          qzp(i,j,kc,QPRES ) = xi1*p    + xi*   p_ref + (apright + amright)*csq

          qzp(i,j,kc,QRHO ) = max(qzp(i,j,kc,QRHO ),small_dens)
          qzp(i,j,kc,QPRES) = max(qzp(i,j,kc,QPRES),small_pres)

          ! **************************************************************************
          ! This is all for qzm
          ! **************************************************************************

          ! Minus state on face kc

          ! note this is different from how we do 1D, 2D, and the
          ! x and y-faces in 3D, where the analogous thing would have
          ! been to find the minus state on face kc+1

          rho  = q(i,j,k3d-1,QRHO)
          u    = q(i,j,k3d-1,QU)
          v    = q(i,j,k3d-1,QV)
          w    = q(i,j,k3d-1,QW)
          p    = q(i,j,k3d-1,QPRES)
          rhoe = q(i,j,k3d-1,QREINT)

          cc   = c(i,j,k3d-1)
          csq  = cc**2
          enth = ( (rhoe+p)/rho )/csq

          ! Set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. w + cc <= 0.0d0) ) then
             rho_ref  = rho
             u_ref    = u
             v_ref    = v
             w_ref    = w
             p_ref    = p
             rhoe_ref = rhoe
          else
             ! This will be the fastest moving state to the right
             rho_ref  = Ip(i,j,km,3,3,QRHO)
             u_ref    = Ip(i,j,km,3,3,QU)
             v_ref    = Ip(i,j,km,3,3,QV)
             w_ref    = Ip(i,j,km,3,3,QW)
             p_ref    = Ip(i,j,km,3,3,QPRES)
             rhoe_ref = Ip(i,j,km,3,3,QREINT)
          endif

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm    = flatn(i,j,k3d-1)*(w_ref    - Ip(i,j,km,3,1,QW))
          dpm    = flatn(i,j,k3d-1)*(p_ref    - Ip(i,j,km,3,1,QPRES))

          drho  = flatn(i,j,k3d-1)*(rho_ref  - Ip(i,j,km,3,2,QRHO))
          du    = flatn(i,j,k3d-1)*(u_ref    - Ip(i,j,km,3,2,QU))
          dv    = flatn(i,j,k3d-1)*(v_ref    - Ip(i,j,km,3,2,QV))
          dp    = flatn(i,j,k3d-1)*(p_ref    - Ip(i,j,km,3,2,QPRES))
          drhoe = flatn(i,j,k3d-1)*(rhoe_ref - Ip(i,j,km,3,2,QREINT))

          dwp    = flatn(i,j,k3d-1)*(w_ref    - Ip(i,j,km,3,3,QW))
          dpp    = flatn(i,j,k3d-1)*(p_ref    - Ip(i,j,km,3,3,QPRES))

          alpham = 0.5d0*(dpm/(rho*cc) - dwm)*rho/cc
          alphap = 0.5d0*(dpp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
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

          xi1 = 1.0d0 - flatn(i,j,k3d-1)
          xi = flatn(i,j,k3d-1)
          qzm(i,j,kc,QRHO  ) = xi1*rho  + xi* rho_ref + apleft + amleft + azrleft
          qzm(i,j,kc,QW    ) = xi1*w    + xi*   w_ref + (apleft - amleft)*cc/rho
          qzm(i,j,kc,QU    ) = xi1*u    + xi*   u_ref + azu1left
          qzm(i,j,kc,QV    ) = xi1*v    + xi*   v_ref + azv1left
          qzm(i,j,kc,QREINT) = xi1*rhoe + xi*rhoe_ref + (apleft + amleft)*enth*csq + azeleft
          qzm(i,j,kc,QPRES ) = xi1*p    + xi*   p_ref + (apleft + amleft)*csq

          qzm(i,j,kc,QRHO ) = max(qzm(i,j,kc,QRHO ),small_dens)
          qzm(i,j,kc,QPRES) = max(qzm(i,j,kc,QPRES),small_pres)
       end do
    end do
    !$OMP END PARALLEL DO

    ! Do all of the passively advected quantities in one loop
    !$OMP parallel do private(n,w,i,j,ipassive) IF(npassive .gt. 1)
    do ipassive = 1, npassive
         n = qpass_map(ipassive)
         do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1

                  ! Plus state on face kc
                  w = q(i,j,k3d,QW)
                  if (w .gt. 0.d0) then
                     qzp(i,j,kc,n) = q(i,j,k3d,n)
                  else if (w .lt. 0.d0) then
                     qzp(i,j,kc,n) = q(i,j,k3d,n) &
                          + flatn(i,j,k3d)*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
                  else
                     qzp(i,j,kc,n) = q(i,j,k3d,n) &
                          + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
                  endif

                  ! Minus state on face k
                  w = q(i,j,k3d-1,QW)
                  if (w .gt. 0.d0) then
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n) &
                          + flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
                  else if (w .lt. 0.d0) then
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n)
                  else
                     qzm(i,j,kc,n) = q(i,j,k3d-1,n) &
                          + 0.5d0*flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
                  endif

               enddo
         enddo
    enddo
    !$OMP end parallel do

  end subroutine tracez_ppm

end module trace_ppm_module
