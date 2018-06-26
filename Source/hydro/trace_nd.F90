module trace_module

  use amrex_error_module, only : amrex_error
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg

  implicit none

  private

  public tracexy, tracez

contains

  subroutine tracexy(q, q_lo, q_hi, &
                     qaux, qa_lo, qa_hi, &
                     dqx, dqy, dq_lo, dq_hi, &
                     qxm, qxp, qym, qyp, qpd_lo, qpd_hi, &
#if (BL_SPACEDIM < 3)
                     dloga, dloga_lo, dloga_hi, &
#endif
#if (BL_SPACEDIM == 1)
                     SrcQ, src_lo, Src_hi, &
#endif
                     ilo1, ilo2, ihi1, ihi2, domlo, domhi, &
                     dx, dt, kc, k3d)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, QVAR, QRHO, QU, QV, QW, QC, &
                                   QREINT, QPRES, &
                                   npassive, qpass_map, small_dens, small_pres, &
                                   ppm_type, fix_mass_flux
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow
    use amrex_fort_module, only : rt => amrex_real
    use ppm_module, only : ppm_reconstruct, ppm_int_profile

    implicit none

    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: dq_lo(3), dq_hi(3)
    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
#if (BL_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
#if (BL_SPACEDIM == 1)
    integer, intent(in) :: src_lo(3), src_hi(3)
#endif
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: kc, k3d
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) ::  dqx(dq_lo(1):dq_hi(1),dq_lo(2):dq_hi(2),dq_lo(3):dq_hi(3),NQ)
    real(rt), intent(in) ::  dqy(dq_lo(1):dq_hi(1),dq_lo(2):dq_hi(2),dq_lo(3):dq_hi(3),NQ)

    real(rt), intent(inout) :: qxm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: qxp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: qym(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: qyp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

#if (BL_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
#if (BL_SPACEDIM == 1)
    real(rt), intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)    
#endif
    real(rt), intent(in) :: dx(3), dt

    ! Local variables
    integer :: i, j, n, ipassive

    real(rt) :: dtdx, dtdy
    real(rt) :: cc, csq, rho, u, v, w, p, rhoe
    real(rt) :: drho, du, dv, dw, dp, drhoe

    real(rt) :: enth, alpham, alphap, alpha0r, alpha0e
    real(rt) :: alpha0u, alpha0v, alpha0w
    real(rt) :: apright, amright, azrright, azeright
    real(rt) :: azu1rght, azv1rght, azw1rght
    real(rt) :: apleft, amleft, azrleft, azeleft
    real(rt) :: azu1left, azv1left, azw1left
    real(rt) :: acmprght, acmpleft, acmpbot, acmptop
    real(rt) :: spzero
    real(rt) :: rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref
    real(rt) :: e(3)
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp
    logical :: fix_mass_flux_lo, fix_mass_flux_hi

    dtdx = dt/dx(1)
#if (BL_SPACEDIM >= 2)
    dtdy = dt/dx(2)
#endif

#ifndef AMREX_USE_CUDA
    if (ppm_type .ne. 0) then
       print *,'Oops -- shouldnt be in tracexy with ppm_type != 0'
       call amrex_error("Error:: trace_3d.f90 :: tracexy")
    end if
#endif

    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (ilo1 == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (ihi1 == domhi(1))


    !-----------------------------------------------------------------------
    ! x-direction
    !-----------------------------------------------------------------------

    ! Compute left and right traced states

    ! construct the right state on the i interface

    do j = ilo2-dg(2), ihi2+dg(2)
       do i = ilo1-1, ihi1+1

          cc = qaux(i,j,k3d,QC)
          csq = cc**2
          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)
          enth = (rhoe+p)/(rho*csq)

          drho = dqx(i,j,kc,QRHO)
          du = dqx(i,j,kc,QU)
          dv = dqx(i,j,kc,QV)
          dw = dqx(i,j,kc,QW)
          dp = dqx(i,j,kc,QPRES)
          drhoe = dqx(i,j,kc,QREINT)

          alpham = HALF*(dp/(rho*cc) - du)*(rho/cc)
          alphap = HALF*(dp/(rho*cc) + du)*(rho/cc)
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0v = dv
          alpha0w = dw

          e(1) = u - cc
          e(2) = u
          e(3) = u + cc

          rho_ref = rho - HALF*(ONE + dtdx*min(e(1),ZERO))*drho
          u_ref = u - HALF*(ONE + dtdx*min(e(1), ZERO))*du
          v_ref = v - HALF*(ONE + dtdx*min(e(1), ZERO))*dv
          w_ref = w - HALF*(ONE + dtdx*min(e(1), ZERO))*dw
          p_ref = p - HALF*(ONE + dtdx*min(e(1), ZERO))*dp
          rhoe_ref = rhoe - HALF*(ONE + dtdx*min(e(1),ZERO))*drhoe

          ! this is -(1/2) ( 1 + dt/dx lambda) (l . dq) r
          apright = 0.25e0_rt*dtdx*(e(1) - e(3))*(ONE - sign(ONE, e(3)))*alphap
          amright = 0.25e0_rt*dtdx*(e(1) - e(1))*(ONE - sign(ONE, e(1)))*alpham

          azrright = 0.25e0*dtdx*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0r
          azeright = 0.25e0*dtdx*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0e
          azv1rght = 0.25e0*dtdx*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0v
          azw1rght = 0.25e0*dtdx*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0w

          if (i .ge. ilo1) then
             qxp(i,j,kc,QRHO) = rho_ref + apright + amright + azrright
             qxp(i,j,kc,QRHO) = max(small_dens, qxp(i,j,kc,QRHO))
             qxp(i,j,kc,QU) = u_ref + (apright - amright)*cc/rho
             qxp(i,j,kc,QV) = v_ref + azv1rght
             qxp(i,j,kc,QW) = w_ref + azw1rght
             qxp(i,j,kc,QPRES) = p_ref + (apright + amright)*csq
             qxp(i,j,kc,QPRES) = max(qxp(i,j,kc,QPRES), small_pres)
             qxp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright

             ! add the source terms in 1-d, since we don't do this in
             ! the transverse routines
#if (BL_SPACEDIM == 1)
             qxp(i,j,kc,QRHO  ) = qxp(i,j,kc,QRHO  ) + HALF*dt*srcQ(i,j,k3d,QRHO)
             qxp(i,j,kc,QRHO  ) = max(small_dens, qxp(i,j,kc,QRHO))
             qxp(i,j,kc,QU    ) = qxp(i,j,kc,QU    ) + HALF*dt*srcQ(i,j,k3d,QU)
             qxp(i,j,kc,QREINT) = qxp(i,j,kc,QREINT) + HALF*dt*srcQ(i,j,k3d,QREINT)
             qxp(i,j,kc,QPRES ) = qxp(i,j,kc,QPRES ) + HALF*dt*srcQ(i,j,k3d,QPRES)
#endif
          end if

          ! now construct the left state on the i+1 interface

          rho_ref = rho + HALF*(ONE - dtdx*max(e(3), ZERO))*drho
          u_ref = u + HALF*(ONE - dtdx*max(e(3), ZERO))*du
          v_ref = v + HALF*(ONE - dtdx*max(e(3), ZERO))*dv
          w_ref = w + HALF*(ONE - dtdx*max(e(3), ZERO))*dw
          p_ref = p + HALF*(ONE - dtdx*max(e(3), ZERO))*dp
          rhoe_ref = rhoe + HALF*(ONE - dtdx*max(e(3), ZERO))*drhoe

          apleft = 0.25e0_rt*dtdx*(e(3) - e(3))*(ONE + sign(ONE, e(3)))*alphap
          amleft = 0.25e0_rt*dtdx*(e(3) - e(1))*(ONE + sign(ONE, e(1)))*alpham

          azrleft = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0r
          azeleft = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0e
          azv1left = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0v
          azw1left = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0w

          if (i .le. ihi1) then
             qxm(i+1,j,kc,QRHO) = rho_ref + apleft + amleft + azrleft
             qxm(i+1,j,kc,QRHO) = max(small_dens, qxm(i+1,j,kc,QRHO))
             qxm(i+1,j,kc,QU) = u_ref + (apleft - amleft)*cc/rho
             qxm(i+1,j,kc,QV) = v_ref + azv1left
             qxm(i+1,j,kc,QW) = w_ref + azw1left
             qxm(i+1,j,kc,QPRES) = p_ref + (apleft + amleft)*csq
             qxm(i+1,j,kc,QPRES) = max(qxm(i+1,j,kc,QPRES), small_pres)
             qxm(i+1,j,kc,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft

             ! add the source terms in 1-d, since we don't do this in
             ! the transverse routines
#if (BL_SPACEDIM == 1)
             qxm(i+1,j,kc,QRHO  ) = qxm(i+1,j,kc,QRHO  ) + HALF*dt*srcQ(i,j,k3d,QRHO)
             qxm(i+1,j,kc,QRHO  ) = max(small_dens, qxm(i+1,j,kc,QRHO))
             qxm(i+1,j,kc,QU    ) = qxm(i+1,j,kc,QU    ) + HALF*dt*srcQ(i,j,k3d,QU)
             qxm(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QREINT) + HALF*dt*srcQ(i,j,k3d,QREINT)
             qxm(i+1,j,kc,QPRES ) = qxm(i+1,j,kc,QPRES ) + HALF*dt*srcQ(i,j,k3d,QPRES)
#endif             
          endif

#if (BL_SPACEDIM < 3)
          ! geometry source terms
          if (dloga(i,j,k3d) /= ZERO) then
             courn = dtdx*(cc + abs(u))
             eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k3d)))
             dlogatmp = min(eta, ONE)*dloga(i,j,k3d)
             sourcr = -HALF*dt*rho*dlogatmp*u
             sourcp = sourcr*csq
             source = sourcp*enth
             if (i .le. ihi1) then
                qxm(i+1,j,kc,QRHO) = qxm(i+1,j,kc,QRHO) + sourcr
                qxm(i+1,j,kc,QRHO) = max(qxm(i+1,j,kc,QRHO),small_dens)
                qxm(i+1,j,kc,QPRES) = qxm(i+1,j,kc,QPRES) + sourcp
                qxm(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QREINT) + source
             end if
             if (i .ge. ilo1) then
                qxp(i,j,kc,QRHO) = qxp(i,j,kc,QRHO) + sourcr
                qxp(i,j,kc,QRHO) = max(qxp(i,j,kc,QRHO),small_dens)
                qxp(i,j,kc,QPRES) = qxp(i,j,kc,QPRES) + sourcp
                qxp(i,j,kc,QREINT) = qxp(i,j,kc,QREINT) + source
             end if
          endif
#endif

#if (BL_SPACEDIM == 1)
    ! Enforce constant mass flux rate if specified
    if (fix_mass_flux_lo) then
       qxm(ilo1,j,kc,QRHO  ) = q(domlo(1)-1,j,k3d,QRHO)
       qxm(ilo1,j,kc,QU    ) = q(domlo(1)-1,j,k3d,QU  )
       qxm(ilo1,j,kc,QPRES ) = q(domlo(1)-1,j,k3d,QPRES)
       qxm(ilo1,j,kc,QREINT) = q(domlo(1)-1,j,k3d,QREINT)
    end if

    ! Enforce constant mass flux rate if specified
    if (fix_mass_flux_hi) then
       qxp(ihi1+1,j,kc,QRHO  ) = q(domhi(1)+1,j,k3d,QRHO)
       qxp(ihi1+1,j,kc,QU    ) = q(domhi(1)+1,j,k3d,QU  )
       qxp(ihi1+1,j,kc,QPRES ) = q(domhi(1)+1,j,k3d,QPRES)
       qxp(ihi1+1,j,kc,QREINT) = q(domhi(1)+1,j,k3d,QREINT)
    end if
#endif

       enddo
    enddo

    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  But we already dealt with all 3 velocity
       ! components above, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       do j = ilo2-dg(2), ihi2+dg(2)

          ! Right state
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,QU)
             if (u .gt. ZERO) then
                spzero = -ONE
             else
                spzero = u*dtdx
             endif
             acmprght = HALF*(-ONE - spzero )*dqx(i,j,kc,n)
             qxp(i,j,kc,n) = q(i,j,k3d,n) + acmprght
          enddo

          ! Left state
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)
             if (u .ge. ZERO) then
                spzero = u*dtdx
             else
                spzero = ONE
             endif
             acmpleft = HALF*(ONE - spzero )*dqx(i,j,kc,n)
             qxm(i+1,j,kc,n) = q(i,j,k3d,n) + acmpleft
          enddo

#if (BL_SPACEDIM == 1)
       if (fix_mass_flux_hi) qxp(ihi1+1,j,kc,n) = q(ihi1+1,j,k3d,n)
       if (fix_mass_flux_lo) qxm(ilo1,j,kc,n) = q(ilo1-1,j,k3d,n)
#endif

       enddo
    enddo

#if (BL_SPACEDIM >= 2)
    !-----------------------------------------------------------------------
    ! y-direction
    !-----------------------------------------------------------------------

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          cc = qaux(i,j,k3d,QC)
          csq = cc**2
          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)
          enth = (rhoe+p)/(rho*csq)

          drho = dqy(i,j,kc,QRHO)
          du = dqy(i,j,kc,QU)
          dv = dqy(i,j,kc,QV)
          dw = dqy(i,j,kc,QW)
          dp = dqy(i,j,kc,QPRES)
          drhoe = dqy(i,j,kc,QREINT)

          alpham = HALF*(dp/(rho*cc) - dv)*(rho/cc)
          alphap = HALF*(dp/(rho*cc) + dv)*(rho/cc)
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0u = du
          alpha0w = dw

          e(1) = v - cc
          e(2) = v
          e(3) = v + cc

          ! construct the right state on the j-1/2 interface

          rho_ref = rho - HALF*(ONE + dtdy*min(e(1), ZERO))*drho
          u_ref = u - HALF*(ONE + dtdy*min(e(1), ZERO))*du
          v_ref = v - HALF*(ONE + dtdy*min(e(1), ZERO))*dv
          w_ref = w - HALF*(ONE + dtdy*min(e(1), ZERO))*dw
          p_ref = p - HALF*(ONE + dtdy*min(e(1), ZERO))*dp
          rhoe_ref = rhoe - HALF*(ONE + dtdy*min(e(1), ZERO))*drhoe

          apright = 0.25e0_rt*dtdy*(e(1) - e(3))*(ONE - sign(ONE, e(3)))*alphap
          amright = 0.25e0_rt*dtdy*(e(1) - e(1))*(ONE - sign(ONE, e(1)))*alpham

          azrright = 0.25e0*dtdy*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0r
          azeright = 0.25e0*dtdy*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0e
          azu1rght = 0.25e0*dtdy*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0u
          azw1rght = 0.25e0*dtdy*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0w

          if (j .ge. ilo2) then
             qyp(i,j,kc,QRHO) = rho_ref + apright + amright + azrright
             qyp(i,j,kc,QRHO) = max(small_dens, qyp(i,j,kc,QRHO))
             qyp(i,j,kc,QV) = v_ref + (apright - amright)*cc/rho
             qyp(i,j,kc,QU) = u_ref + azu1rght
             qyp(i,j,kc,QW) = w_ref + azw1rght
             qyp(i,j,kc,QPRES) = p_ref + (apright + amright)*csq
             qyp(i,j,kc,QPRES) = max(qyp(i,j,kc,QPRES), small_pres)
             qyp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright
          end if

          ! construct the left state on the j+1/2 interface

          rho_ref = rho + HALF*(ONE - dtdy*max(e(3), ZERO))*drho
          u_ref = u + HALF*(ONE - dtdy*max(e(3), ZERO))*du
          v_ref = v + HALF*(ONE - dtdy*max(e(3), ZERO))*dv
          w_ref = w + HALF*(ONE - dtdy*max(e(3), ZERO))*dw
          p_ref = p + HALF*(ONE - dtdy*max(e(3), ZERO))*dp
          rhoe_ref = rhoe + HALF*(ONE - dtdy*max(e(3), ZERO))*drhoe

          apleft = 0.25e0_rt*dtdy*(e(3) - e(3))*(ONE + sign(ONE, e(3)))*alphap
          amleft = 0.25e0_rt*dtdy*(e(3) - e(1))*(ONE + sign(ONE, e(1)))*alpham

          azrleft = 0.25e0_rt*dtdy*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0r
          azeleft = 0.25e0_rt*dtdy*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0e
          azu1left = 0.25e0_rt*dtdy*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0u
          azw1left = 0.25e0_rt*dtdy*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0w

          if (j .le. ihi2) then
             qym(i,j+1,kc,QRHO) = rho_ref + apleft + amleft + azrleft
             qym(i,j+1,kc,QRHO) = max(small_dens, qym(i,j+1,kc,QRHO))
             qym(i,j+1,kc,QV) = v_ref + (apleft - amleft)*cc/rho
             qym(i,j+1,kc,QU) = u_ref + azu1left
             qym(i,j+1,kc,QW) = w_ref + azw1left
             qym(i,j+1,kc,QPRES) = p_ref + (apleft + amleft)*csq
             qym(i,j+1,kc,QPRES) = max(qym(i,j+1,kc,QPRES), small_pres)
             qym(i,j+1,kc,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft
          endif

       enddo
    enddo

    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  But we already dealt with all 3 velocity
       ! components above, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       do i = ilo1-1, ihi1+1

          ! Top state
          do j = ilo2, ihi2+1
             v = q(i,j,k3d,QV)
             if (v .gt. ZERO) then
                spzero = -ONE
             else
                spzero = v*dtdy
             endif
             acmptop = HALF*(-ONE - spzero )*dqy(i,j,kc,n)
             qyp(i,j,kc,n) = q(i,j,k3d,n) + acmptop
          enddo

          ! Bottom state
          do j = ilo2-1, ihi2
             v = q(i,j,k3d,QV)
             if (v .ge. ZERO) then
                spzero = v*dtdy
             else
                spzero = ONE
             endif
             acmpbot = HALF*(ONE - spzero )*dqy(i,j,kc,n)
             qym(i,j+1,kc,n) = q(i,j,k3d,n) + acmpbot
          enddo

       enddo
    enddo
#endif

  end subroutine tracexy

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine tracez(q, q_lo, q_hi, &
                    qaux, qa_lo, qa_hi, &
                    dqz, dq_lo, dq_hi, &
                    qzm, qzp, qpd_lo, qpd_hi, &
                    ilo1, ilo2, ihi1, ihi2, domlo, domhi, &
                    dx, dt, km, kc, k3d)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, QVAR, QRHO, QU, QV, QW, QC, &
                                   QREINT, QPRES, &
                                   npassive, qpass_map, small_dens, small_pres, ppm_type
    use amrex_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: dq_lo(3), dq_hi(3)
    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: km, kc, k3d
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: dqz(dq_lo(1):dq_hi(1),dq_lo(2):dq_hi(2),dq_lo(3):dq_hi(3),NQ)
    real(rt), intent(inout) :: qzm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: qzp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: dx(3), dt

    ! Local variables
    integer :: i, j
    integer :: n, ipassive

    real(rt) :: dtdz
    real(rt) :: cc, csq, rho, u, v, w, p, rhoe

    real(rt) :: drho, du, dv, dw, dp, drhoe
    real(rt) :: enth, alpham, alphap, alpha0r, alpha0e
    real(rt) :: alpha0u, alpha0v
    real(rt) :: apright, amright, azrright, azeright
    real(rt) :: azu1rght, azv1rght
    real(rt) :: apleft, amleft, azrleft, azeleft
    real(rt) :: azu1left, azv1left
    real(rt) :: acmpbot, acmptop
    real(rt) :: spzero
    real(rt) :: rho_ref, u_ref, v_ref, w_ref, p_ref, rhoe_ref
    real(rt) :: e(3)

    dtdz = dt/dx(3)

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          cc = qaux(i,j,k3d,QC)
          csq = cc**2
          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)
          enth = (rhoe+p)/(rho*csq)

          drho = dqz(i,j,kc,QRHO)
          du = dqz(i,j,kc,QU)
          dv = dqz(i,j,kc,QV)
          dw = dqz(i,j,kc,QW)
          dp = dqz(i,j,kc,QPRES)
          drhoe = dqz(i,j,kc,QREINT)

          alpham = HALF*(dp/(rho*cc) - dw)*(rho/cc)
          alphap = HALF*(dp/(rho*cc) + dw)*(rho/cc)
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0u = du
          alpha0v = dv

          e(1) = w - cc
          e(2) = w
          e(3) = w + cc

          rho_ref = rho - HALF*(ONE + dtdz*min(e(1), ZERO))*drho
          u_ref = u - HALF*(ONE + dtdz*min(e(1), ZERO))*du
          v_ref = v - HALF*(ONE + dtdz*min(e(1), ZERO))*dv
          w_ref = w - HALF*(ONE + dtdz*min(e(1), ZERO))*dw
          p_ref = p - HALF*(ONE + dtdz*min(e(1), ZERO))*dp
          rhoe_ref = rhoe - HALF*(ONE + dtdz*min(e(1),ZERO))*drhoe

          apright = 0.25e0_rt*dtdz*(e(1) - e(3))*(ONE - sign(ONE, e(3)))*alphap
          amright = 0.25e0_rt*dtdz*(e(1) - e(1))*(ONE - sign(ONE, e(1)))*alpham

          azrright = 0.25e0*dtdz*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0r
          azeright = 0.25e0*dtdz*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0e
          azu1rght = 0.25e0*dtdz*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0u
          azv1rght = 0.25e0*dtdz*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0v

          qzp(i,j,kc,QRHO) = rho_ref + apright + amright + azrright
          qzp(i,j,kc,QRHO) = max(small_dens, qzp(i,j,kc,QRHO))
          qzp(i,j,kc,QW) = w_ref + (apright - amright)*(cc/rho)
          qzp(i,j,kc,QU) = u_ref + azu1rght
          qzp(i,j,kc,QV) = v_ref + azv1rght
          qzp(i,j,kc,QPRES) = p_ref + (apright + amright)*csq
          qzp(i,j,kc,QPRES) = max(qzp(i,j,kc,QPRES), small_pres)
          qzp(i,j,kc,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright

          ! repeat above with km (k3d-1) to get qzm at kc
          cc = qaux(i,j,k3d-1,QC)
          csq = cc**2
          rho = q(i,j,k3d-1,QRHO)
          u = q(i,j,k3d-1,QU)
          v = q(i,j,k3d-1,QV)
          w = q(i,j,k3d-1,QW)
          p = q(i,j,k3d-1,QPRES)
          rhoe = q(i,j,k3d-1,QREINT)
          enth = (rhoe+p)/(rho*csq)

          drho = dqz(i,j,km,QRHO)
          du = dqz(i,j,km,QU)
          dv = dqz(i,j,km,QV)
          dw = dqz(i,j,km,QW)
          dp = dqz(i,j,km,QPRES)
          drhoe = dqz(i,j,km,QREINT)

          alpham = HALF*(dp/(rho*cc) - dw)*(rho/cc)
          alphap = HALF*(dp/(rho*cc) + dw)*(rho/cc)
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0u = du
          alpha0v = dv

          e(1) = w - cc
          e(2) = w
          e(3) = w + cc

          rho_ref = rho + HALF*(ONE - dtdz*max(e(3), ZERO))*drho
          u_ref = u + HALF*(ONE - dtdz*max(e(3), ZERO))*du
          v_ref = v + HALF*(ONE - dtdz*max(e(3), ZERO))*dv
          w_ref = w + HALF*(ONE - dtdz*max(e(3), ZERO))*dw
          p_ref = p + HALF*(ONE - dtdz*max(e(3), ZERO))*dp
          rhoe_ref = rhoe + HALF*(ONE - dtdz*max(e(3), ZERO))*drhoe

          apleft = 0.25e0_rt*dtdz*(e(3) - e(3))*(ONE + sign(ONE, e(3)))*alphap
          amleft = 0.25e0_rt*dtdz*(e(3) - e(1))*(ONE + sign(ONE, e(1)))*alpham

          azrleft = 0.25e0_rt*dtdz*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0r
          azeleft = 0.25e0_rt*dtdz*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0e
          azu1left = 0.25e0_rt*dtdz*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0u
          azv1left = 0.25e0_rt*dtdz*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0v

          qzm(i,j,kc,QRHO) = rho_ref + apleft + amleft + azrleft
          qzm(i,j,kc,QRHO) = max(small_dens, qzm(i,j,kc,QRHO))
          qzm(i,j,kc,QW) = w_ref + (apleft - amleft)*(cc/rho)
          qzm(i,j,kc,QU) = u_ref + azu1left
          qzm(i,j,kc,QV) = v_ref + azv1left
          qzm(i,j,kc,QPRES) = p_ref + (apleft + amleft)*csq
          qzm(i,j,kc,QPRES) = max(qzm(i,j,kc,QPRES), small_pres)
          qzm(i,j,kc,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft

       enddo
    enddo

    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  But we already dealt with all 3 velocity
       ! components above, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! Top state
             w = q(i,j,k3d,QW)
             if (w .gt. ZERO) then
                spzero = -ONE
             else
                spzero = w*dtdz
             endif
             acmptop = HALF*(-ONE - spzero )*dqz(i,j,kc,n)
             qzp(i,j,kc,n) = q(i,j,k3d,n) + acmptop

             ! Bottom state
             w = q(i,j,k3d-1,QW)
             if (w .ge. ZERO) then
                spzero = w*dtdz
             else
                spzero = ONE
             endif
             acmpbot = HALF*(ONE - spzero )*dqz(i,j,km,n)
             qzm(i,j,kc,n) = q(i,j,k3d-1,n) + acmpbot
          enddo
       enddo
    enddo

  end subroutine tracez

end module trace_module
