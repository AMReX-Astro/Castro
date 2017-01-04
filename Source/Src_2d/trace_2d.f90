module trace_module

  use bl_fort_module, only : rt => c_real
  implicit none

  private

  public trace

contains

  subroutine trace(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                   dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                   qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                   src,src_l1,src_l2,src_h1,src_h2, &
                   ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

    use meth_params_module, only : plm_iorder, QVAR, QRHO, QU, QV, &
                                   QREINT, QPRES, &
                                   npassive, qpass_map, small_dens, small_pres, ppm_type, use_pslope
    use slope_module, only : uslope, pslope, multid_slope
    use bl_constants_module

    use bl_fort_module, only : rt => c_real
    implicit none

    integer ilo1,ilo2,ihi1,ihi2
    integer qd_l1,qd_l2,qd_h1,qd_h2
    integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
    integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer src_l1,src_l2,src_h1,src_h2

    real(rt)         dx, dy, dt
    real(rt)             q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    real(rt)             c(qd_l1:qd_h1,qd_l2:qd_h2)
    real(rt)         flatn(qd_l1:qd_h1,qd_l2:qd_h2)
    real(rt)         dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
    real(rt)         qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    real(rt)         qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    real(rt)         qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    real(rt)         qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    real(rt)         src(src_l1:src_h1,src_l2:src_h2,QVAR)

    real(rt)        , allocatable :: dqx(:,:,:), dqy(:,:,:)

    ! Local variables
    integer i, j
    integer n, ipassive

    real(rt)         dtdx, dtdy
    real(rt)         cc, csq, rho, u, v, p, rhoe
    real(rt)         drho, du, dv, dp, drhoe
    
    real(rt)         enth, alpham, alphap, alpha0r, alpha0e
    real(rt)         alpha0u, alpha0v
    real(rt)         spzero
    real(rt)         apright, amright, azrright, azeright
    real(rt)         azu1rght, azv1rght
    real(rt)         apleft, amleft, azrleft, azeleft
    real(rt)         azu1left, azv1left
    real(rt)         acmprght, acmpleft, acmpbot, acmptop
    real(rt)         sourcr,sourcp,source,courn,eta,dlogatmp
    
    real(rt)         :: rho_ref, u_ref, v_ref, p_ref, rhoe_ref
    real(rt)         :: e(3)

    if (ppm_type .ne. 0) then
       print *,'Oops -- shouldnt be in trace with ppm_type != 0'
       call bl_error("Error:: trace_2d.f90")
    end if
    
    dtdx = dt/dx
    dtdy = dt/dy

    allocate(dqx(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR))
    allocate(dqy(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR))

    ! Compute slopes
    if (plm_iorder == 1) then
       dqx(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:QVAR) = ZERO
       dqy(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:QVAR) = ZERO

    elseif (plm_iorder == 2) then
       ! these are piecewise linear slopes.  The limiter is a 4th order
       ! limiter, but the overall method will be second order.
       call uslope(q, &
                   flatn, qd_l1, qd_l2, qd_h1, qd_h2, &
                   dqx  ,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                   ilo1,ilo2,ihi1,ihi2,QVAR,1)

       call uslope(q,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                   dqy,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                   ilo1,ilo2,ihi1,ihi2,QVAR,2)

       if (use_pslope .eq. 1) then
          call pslope(q(:,:,QPRES),q(:,:,QRHO),  &
                      flatn         , qd_l1, qd_l2, qd_h1, qd_h2, &
                      dqx(:,:,QPRES),qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                      src           ,src_l1,src_l2,src_h1,src_h2, &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,1)

          call pslope(q(:,:,QPRES),q(:,:,QRHO),  &
                      flatn         , qd_l1, qd_l2, qd_h1, qd_h2, &
                      dqy(:,:,QPRES),qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                      src           ,src_l1,src_l2,src_h1,src_h2, &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,2)
       endif

    elseif (plm_iorder == -2) then
       ! these are also piecewise linear, but it uses a multidimensional
       ! reconstruction based on the BDS advection method to construct
       ! the x- and y-slopes together
       do n = 1, QVAR
          call multid_slope(q(:,:,n), flatn, &
                            qd_l1, qd_l2, qd_h1, qd_h2, &
                            dqx(:,:,n), dqy(:,:,n), &
                            qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                            dx, dy, &
                            ilo1, ilo2, ihi1, ihi2)
       enddo

    else
       call bl_error("ERROR: invalid value of islope")
       
    endif

    !-------------------------------------------------------------------------
    ! x-direction
    !-------------------------------------------------------------------------
      
    ! Compute left and right traced states

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
          
          drho = dqx(i,j,QRHO)
          du = dqx(i,j,QU)
          dv = dqx(i,j,QV)
          dp = dqx(i,j,QPRES)
          drhoe = dqx(i,j,QREINT)

          alpham = HALF*(dp/(rho*cc) - du)*rho/cc
          alphap = HALF*(dp/(rho*cc) + du)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0v = dv

          e(1) = u-cc
          e(2) = u
          e(3) = u+cc

          ! construct the right state on the i-1/2 interface
                    
          rho_ref = rho - HALF*(ONE + dtdx*min(e(1),ZERO))*drho
          u_ref = u - HALF*(ONE + dtdx*min(e(1),ZERO))*du
          v_ref = v - HALF*(ONE + dtdx*min(e(1),ZERO))*dv
          p_ref = p - HALF*(ONE + dtdx*min(e(1),ZERO))*dp
          rhoe_ref = rhoe - HALF*(ONE + dtdx*min(e(1),ZERO))*drhoe

          ! this is -(1/2) ( 1 + dt/dx lambda) (l . dq) r
          apright = 0.25e0_rt*dtdx*(e(1) - e(3))*(ONE - sign(ONE,e(3)))*alphap
          amright = 0.25e0_rt*dtdx*(e(1) - e(1))*(ONE - sign(ONE,e(1)))*alpham
          
          azrright = 0.25e0*dtdx*(e(1)-e(2))*(ONE - sign(ONE,e(2)))*alpha0r
          azeright = 0.25e0*dtdx*(e(1)-e(2))*(ONE - sign(ONE,e(2)))*alpha0e
          azv1rght = 0.25e0*dtdx*(e(1)-e(2))*(ONE - sign(ONE,e(2)))*alpha0v
            
          if (i .ge. ilo1) then
             qxp(i,j,QRHO) = rho_ref + apright + amright + azrright
             qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))
             qxp(i,j,QU) = u_ref + (apright - amright)*cc/rho
             qxp(i,j,QV) = v_ref + azv1rght
             qxp(i,j,QPRES) = p_ref + (apright + amright)*csq
             qxp(i,j,QPRES) = max(qxp(i,j,QPRES),small_pres)
             qxp(i,j,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright
          end if


          ! construct the left state on the i+1/2 interface

          rho_ref = rho + HALF*(ONE - dtdx*max(e(3),ZERO))*drho
          u_ref = u + HALF*(ONE - dtdx*max(e(3),ZERO))*du
          v_ref = v + HALF*(ONE - dtdx*max(e(3),ZERO))*dv
          p_ref = p + HALF*(ONE - dtdx*max(e(3),ZERO))*dp
          rhoe_ref = rhoe + HALF*(ONE - dtdx*max(e(3),ZERO))*drhoe

          apleft = 0.25e0_rt*dtdx*(e(3) - e(3))*(ONE + sign(ONE,e(3)))*alphap
          amleft = 0.25e0_rt*dtdx*(e(3) - e(1))*(ONE + sign(ONE,e(1)))*alpham
          
          azrleft = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE,e(2)))*alpha0r
          azeleft = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE,e(2)))*alpha0e
          azv1left = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE,e(2)))*alpha0v
          
          if (i .le. ihi1) then
             qxm(i+1,j,QRHO) = rho_ref + apleft + amleft + azrleft
             qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
             qxm(i+1,j,QU) = u_ref + (apleft - amleft)*cc/rho
             qxm(i+1,j,QV) = v_ref + azv1left
             qxm(i+1,j,QPRES) = p_ref + (apleft + amleft)*csq
             qxm(i+1,j,QPRES) = max(qxm(i+1,j,QPRES), small_pres)
             qxm(i+1,j,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft
          end if
          

          ! geometry source terms
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
          
       enddo
    enddo

    ! We do all passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1

          ! Right state
          do i = ilo1, ihi1+1
             u = q(i,j,QU)
             if (u .gt. ZERO) then
                spzero = -ONE
             else
                spzero = u*dtdx
             endif
             acmprght = HALF*(-ONE - spzero )*dqx(i,j,n)
             qxp(i,j,n) = q(i,j,n) + acmprght
          enddo
          
          ! Left state
          do i = ilo1-1, ihi1
             u = q(i,j,QU)
             if (u .ge. ZERO) then
                spzero = u*dtdx
             else
                spzero = ONE
             endif
             acmpleft = HALF*(ONE - spzero )*dqx(i,j,n)
             qxm(i+1,j,n) = q(i,j,n) + acmpleft
          enddo

       enddo
    enddo



    !-------------------------------------------------------------------------
    ! y-direction
    !-------------------------------------------------------------------------
    
    ! Compute left and right traced states
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
          
          drho = dqy(i,j,QRHO)
          du = dqy(i,j,QU)
          dv = dqy(i,j,QV)
          dp = dqy(i,j,QPRES)
          drhoe = dqy(i,j,QREINT)
          
          alpham = HALF*(dp/(rho*cc) - dv)*rho/cc
          alphap = HALF*(dp/(rho*cc) + dv)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0u = du
          
          e(1) = v-cc
          e(2) = v
          e(3) = v+cc

          ! construct the right state on the j-1/2 interface
          rho_ref = rho - HALF*(ONE + dtdy*min(e(1),ZERO))*drho
          u_ref = u - HALF*(ONE + dtdy*min(e(1),ZERO))*du
          v_ref = v - HALF*(ONE + dtdy*min(e(1),ZERO))*dv
          p_ref = p - HALF*(ONE + dtdy*min(e(1),ZERO))*dp
          rhoe_ref = rhoe - HALF*(ONE + dtdy*min(e(1),ZERO))*drhoe

          apright = 0.25e0_rt*dtdy*(e(1) - e(3))*(ONE - sign(ONE,e(3)))*alphap
          amright = 0.25e0_rt*dtdy*(e(1) - e(1))*(ONE - sign(ONE,e(1)))*alpham
          
          azrright = 0.25e0*dtdy*(e(1)-e(2))*(ONE - sign(ONE,e(2)))*alpha0r
          azeright = 0.25e0*dtdy*(e(1)-e(2))*(ONE - sign(ONE,e(2)))*alpha0e
          azu1rght = 0.25e0*dtdy*(e(1)-e(2))*(ONE - sign(ONE,e(2)))*alpha0u
          
          if (j .ge. ilo2) then
             qyp(i,j,QRHO) = rho_ref + apright + amright + azrright
             qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))
             qyp(i,j,QV) = v_ref + (apright - amright)*cc/rho
             qyp(i,j,QU) = u_ref + azu1rght
             qyp(i,j,QPRES) = p_ref + (apright + amright)*csq
             qyp(i,j,QPRES) = max(qyp(i,j,QPRES), small_pres)
             qyp(i,j,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright
          end if


          ! construct the left state on the j+1/2 interface

          rho_ref = rho + HALF*(ONE - dtdy*max(e(3),ZERO))*drho
          u_ref = u + HALF*(ONE - dtdy*max(e(3),ZERO))*du
          v_ref = v + HALF*(ONE - dtdy*max(e(3),ZERO))*dv
          p_ref = p + HALF*(ONE - dtdy*max(e(3),ZERO))*dp
          rhoe_ref = rhoe + HALF*(ONE - dtdy*max(e(3),ZERO))*drhoe

          apleft = 0.25e0_rt*dtdy*(e(3) - e(3))*(ONE + sign(ONE,e(3)))*alphap
          amleft = 0.25e0_rt*dtdy*(e(3) - e(1))*(ONE + sign(ONE,e(1)))*alpham

          azrleft = 0.25e0_rt*dtdy*(e(3) - e(2))*(ONE + sign(ONE,e(2)))*alpha0r
          azeleft = 0.25e0_rt*dtdy*(e(3) - e(2))*(ONE + sign(ONE,e(2)))*alpha0e
          azu1left = 0.25e0_rt*dtdy*(e(3) - e(2))*(ONE + sign(ONE,e(2)))*alpha0u
          
          if (j .le. ihi2) then
             qym(i,j+1,QRHO) = rho_ref + apleft + amleft + azrleft
             qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))
             qym(i,j+1,QV) = v_ref + (apleft - amleft)*cc/rho
             qym(i,j+1,QU) = u_ref + azu1left
             qym(i,j+1,QPRES) = p_ref + (apleft + amleft)*csq
             qym(i,j+1,QPRES) = max(qym(i,j+1,QPRES), small_pres)
             qym(i,j+1,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft
          end if
          
       enddo
    enddo

    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do i = ilo1-1, ihi1+1
          
          ! Top state
          do j = ilo2, ihi2+1
             v = q(i,j,QV)
             if (v .gt. ZERO) then
                spzero = -ONE
             else
                spzero = v*dtdy
             endif
             acmptop = HALF*(-ONE - spzero )*dqy(i,j,n)
             qyp(i,j,n) = q(i,j,n) + acmptop
          enddo
          
          ! Bottom state
          do j = ilo2-1, ihi2
             v = q(i,j,QV)
             if (v .ge. ZERO) then
                spzero = v*dtdy
             else
                spzero = ONE
             endif
             acmpbot = HALF*(ONE - spzero )*dqy(i,j,n)
             qym(i,j+1,n) = q(i,j,n) + acmpbot
          enddo
          
       enddo
    enddo

    deallocate(dqx,dqy)
    
  end subroutine trace
  
end module trace_module
