module trace_plm_module

  use amrex_error_module, only : amrex_error
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg

  implicit none

  private

  public trace_plm

contains

  subroutine trace_plm(idir, q, q_lo, q_hi, &
                       qaux, qa_lo, qa_hi, &
                       dq, dq_lo, dq_hi, &
                       qm, qp, qpd_lo, qpd_hi, &
#if (AMREX_SPACEDIM < 3)
                       dloga, dloga_lo, dloga_hi, &
#endif
#if (AMREX_SPACEDIM == 1)
                       SrcQ, src_lo, Src_hi, &
#endif
                       lo, hi, domlo, domhi, &
                       dx, dt)

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

    integer, intent(in) :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: dq_lo(3), dq_hi(3)
    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
#if (AMREX_SPACEDIM == 1)
    integer, intent(in) :: src_lo(3), src_hi(3)
#endif
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) ::  dq(dq_lo(1):dq_hi(1),dq_lo(2):dq_hi(2),dq_lo(3):dq_hi(3),NQ)

    real(rt), intent(inout) :: qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
#if (AMREX_SPACEDIM == 1)
    real(rt), intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
#endif
    real(rt), intent(in) :: dx(3), dt

    ! Local variables
    integer :: i, j, k, n, ipassive

    real(rt) :: dtdx
    real(rt) :: cc, csq, rho, un, ut, utt, p, rhoe
    real(rt) :: drho, dun, dut, dutt, dp, drhoe

    real(rt) :: enth, alpham, alphap, alpha0r, alpha0e
    real(rt) :: alpha0un, alpha0ut, alpha0utt
    real(rt) :: apright, amright, azrright, azeright
    real(rt) :: azun1rght, azut1rght, azutt1rght
    real(rt) :: apleft, amleft, azrleft, azeleft
    real(rt) :: azun1left, azut1left, azutt1left
    real(rt) :: acmprght, acmpleft, acmpbot, acmptop
    real(rt) :: spzero
    real(rt) :: rho_ref, un_ref, ut_ref, utt_ref, p_ref, rhoe_ref
    real(rt) :: e(3)
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp
    logical :: fix_mass_flux_lo, fix_mass_flux_hi

    integer :: ix, iy, iz, QUN, QUT, QUTT

    dtdx = dt/dx(idir)

#ifndef AMREX_USE_CUDA
    if (ppm_type .ne. 0) then
       print *,'Oops -- shouldnt be in tracexy with ppm_type != 0'
       call amrex_error("Error:: trace_3d.f90 :: tracexy")
    end if
#endif

    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (lo(1) == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (hi(1) == domhi(1))

    if (idir == 1) then
       ix = 1
       iy = 0
       iz = 0
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       ix = 0
       iy = 1
       iz = 0
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       ix = 0
       iy = 0
       iz = 1
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! Compute left and right traced states

    ! construct the right state on the i interface

    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-1, hi(1)+1

             cc = qaux(i,j,k,QC)
             csq = cc**2
             rho = q(i,j,k,QRHO)
             un = q(i,j,k,QUN)
             ut = q(i,j,k,QUT)
             utt = q(i,j,k,QUTT)
             p = q(i,j,k,QPRES)
             rhoe = q(i,j,k,QREINT)
             enth = (rhoe+p)/(rho*csq)

             drho = dq(i,j,k,QRHO)
             dun = dq(i,j,k,QUN)
             dut = dq(i,j,k,QUT)
             dutt = dq(i,j,k,QUTT)
             dp = dq(i,j,k,QPRES)
             drhoe = dq(i,j,k,QREINT)

             alpham = HALF*(dp/(rho*cc) - dun)*(rho/cc)
             alphap = HALF*(dp/(rho*cc) + dun)*(rho/cc)
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
             alpha0ut = dut
             alpha0utt = dutt

             e(1) = un - cc
             e(2) = un
             e(3) = un + cc

             rho_ref = rho - HALF*(ONE + dtdx*min(e(1),ZERO))*drho
             un_ref = un - HALF*(ONE + dtdx*min(e(1), ZERO))*dun
             ut_ref = ut - HALF*(ONE + dtdx*min(e(1), ZERO))*dut
             utt_ref = utt - HALF*(ONE + dtdx*min(e(1), ZERO))*dutt
             p_ref = p - HALF*(ONE + dtdx*min(e(1), ZERO))*dp
             rhoe_ref = rhoe - HALF*(ONE + dtdx*min(e(1),ZERO))*drhoe

             ! this is -(1/2) ( 1 + dt/dx lambda) (l . dq) r
             apright = 0.25e0_rt*dtdx*(e(1) - e(3))*(ONE - sign(ONE, e(3)))*alphap
             amright = 0.25e0_rt*dtdx*(e(1) - e(1))*(ONE - sign(ONE, e(1)))*alpham

             azrright = 0.25e0*dtdx*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0r
             azeright = 0.25e0*dtdx*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0e
             azut1rght = 0.25e0*dtdx*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0ut
             azutt1rght = 0.25e0*dtdx*(e(1) - e(2))*(ONE - sign(ONE, e(2)))*alpha0utt

             if ((idir == 1 .and. i >= lo(1)) .or. &
                 (idir == 2 .and. j >= lo(2)) .or. &
                 (idir == 3 .and. k >= lo(3))) then

                qp(i,j,k,QRHO) = rho_ref + apright + amright + azrright
                qp(i,j,k,QRHO) = max(small_dens, qp(i,j,k,QRHO))
                qp(i,j,k,QUN) = un_ref + (apright - amright)*cc/rho
                qp(i,j,k,QUT) = ut_ref + azut1rght
                qp(i,j,k,QUTT) = utt_ref + azutt1rght
                qp(i,j,k,QPRES) = p_ref + (apright + amright)*csq
                qp(i,j,k,QPRES) = max(qp(i,j,k,QPRES), small_pres)
                qp(i,j,k,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright

                ! add the source terms in 1-d, since we don't do this in
                ! the transverse routines
#if (AMREX_SPACEDIM == 1)
                qp(i,j,k,QRHO  ) = qp(i,j,k,QRHO  ) + HALF*dt*srcQ(i,j,k,QRHO)
                qp(i,j,k,QRHO  ) = max(small_dens, qp(i,j,k,QRHO))
                qp(i,j,k,QUN   ) = qp(i,j,k,QUN    ) + HALF*dt*srcQ(i,j,k,QUN)
                qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + HALF*dt*srcQ(i,j,k,QREINT)
                qp(i,j,k,QPRES ) = qp(i,j,k,QPRES ) + HALF*dt*srcQ(i,j,k,QPRES)
#endif
             end if

             ! now construct the left state on the i+1 interface

             rho_ref = rho + HALF*(ONE - dtdx*max(e(3), ZERO))*drho
             un_ref = un + HALF*(ONE - dtdx*max(e(3), ZERO))*dun
             ut_ref = ut + HALF*(ONE - dtdx*max(e(3), ZERO))*dut
             utt_ref = utt + HALF*(ONE - dtdx*max(e(3), ZERO))*dutt
             p_ref = p + HALF*(ONE - dtdx*max(e(3), ZERO))*dp
             rhoe_ref = rhoe + HALF*(ONE - dtdx*max(e(3), ZERO))*drhoe

             apleft = 0.25e0_rt*dtdx*(e(3) - e(3))*(ONE + sign(ONE, e(3)))*alphap
             amleft = 0.25e0_rt*dtdx*(e(3) - e(1))*(ONE + sign(ONE, e(1)))*alpham

             azrleft = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0r
             azeleft = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0e
             azut1left = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0ut
             azutt1left = 0.25e0_rt*dtdx*(e(3) - e(2))*(ONE + sign(ONE, e(2)))*alpha0utt

             if ((idir == 1 .and. i <= hi(1)) .or. &
                 (idir == 2 .and. j <= hi(2)) .or. &
                 (idir == 3 .and. k <= hi(3))) then

                qm(i+ix,j+iy,k+iz,QRHO) = rho_ref + apleft + amleft + azrleft
                qm(i+ix,j+iy,k+iz,QRHO) = max(small_dens, qm(i+ix,j+iy,k+iz,QRHO))
                qm(i+ix,j+iy,k+iz,QUN) = un_ref + (apleft - amleft)*cc/rho
                qm(i+ix,j+iy,k+iz,QUT) = ut_ref + azut1left
                qm(i+ix,j+iy,k+iz,QUTT) = utt_ref + azutt1left
                qm(i+ix,j+iy,k+iz,QPRES) = p_ref + (apleft + amleft)*csq
                qm(i+ix,j+iy,k+iz,QPRES) = max(qm(i+ix,j+iy,k+iz,QPRES), small_pres)
                qm(i+ix,j+iy,k+iz,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft

                ! add the source terms in 1-d, since we don't do this in
                ! the transverse routines
#if (AMREX_SPACEDIM == 1)
                qm(i+ix,j+iy,k+iz,QRHO  ) = qm(i+ix,j+iy,k+iz,QRHO  ) + HALF*dt*srcQ(i,j,k,QRHO)
                qm(i+ix,j+iy,k+iz,QRHO  ) = max(small_dens, qm(i+ix,j+iy,k+iz,QRHO))
                qm(i+ix,j+iy,k+iz,QUN   ) = qm(i+ix,j+iy,k+iz,QUN    ) + HALF*dt*srcQ(i,j,k,QUN)
                qm(i+ix,j+iy,k+iz,QREINT) = qm(i+ix,j+iy,k+iz,QREINT) + HALF*dt*srcQ(i,j,k,QREINT)
                qm(i+ix,j+iy,k+iz,QPRES ) = qm(i+ix,j+iy,k+iz,QPRES ) + HALF*dt*srcQ(i,j,k,QPRES)
#endif
             endif

#if (AMREX_SPACEDIM < 3)
             ! geometry source terms
             if (dloga(i,j,k) /= ZERO) then
                courn = dtdx*(cc + abs(un))
                eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k)))
                dlogatmp = min(eta, ONE)*dloga(i,j,k)
                sourcr = -HALF*dt*rho*dlogatmp*un
                sourcp = sourcr*csq
                source = sourcp*enth

                if ((idir == 1 .and. i <= hi(1)) .or. &
                    (idir == 2 .and. j <= hi(2))) then
                   qm(i+ix,j+iy,k+iz,QRHO) = qm(i+ix,j+iy,k+iz,QRHO) + sourcr
                   qm(i+ix,j+iy,k+iz,QRHO) = max(qm(i+ix,j+iy,k+iz,QRHO),small_dens)
                   qm(i+ix,j+iy,k+iz,QPRES) = qm(i+ix,j+iy,k+iz,QPRES) + sourcp
                   qm(i+ix,j+iy,k+iz,QREINT) = qm(i+ix,j+iy,k+iz,QREINT) + source
                end if
                if ((idir == 1 .and. i >= lo(1)) .or. &
                    (idir == 2 .and. j >= lo(2))) then
                   qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr
                   qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO),small_dens)
                   qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                end if
             end if
#endif

#if (AMREX_SPACEDIM == 1)
             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_lo) then
                qm(lo(1),j,k,QRHO  ) = q(domlo(1)-1,j,k,QRHO)
                qm(lo(1),j,k,QUN   ) = q(domlo(1)-1,j,k,QUN )
                qm(lo(1),j,k,QPRES ) = q(domlo(1)-1,j,k,QPRES)
                qm(lo(1),j,k,QREINT) = q(domlo(1)-1,j,k,QREINT)
             end if

             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_hi) then
                qp(hi(1)+1,j,k,QRHO  ) = q(domhi(1)+1,j,k,QRHO)
                qp(hi(1)+1,j,k,QUN    ) = q(domhi(1)+1,j,k,QUN )
                qp(hi(1)+1,j,k,QPRES ) = q(domhi(1)+1,j,k,QPRES)
                qp(hi(1)+1,j,k,QREINT) = q(domhi(1)+1,j,k,QREINT)
             end if
#endif

          end do
       end do
    end do

    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  But we already dealt with all 3 velocity
       ! components above, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2)-dg(2), hi(2)+dg(2)
             do i = lo(1)-1, hi(1)+1

                ! Right state
                if ((idir == 1 .and. i >= lo(1)) .or. &
                    (idir == 2 .and. j >= lo(2)) .or. &
                    (idir == 3 .and. k >= lo(3))) then

                   un = q(i,j,k,QUN)
                   if (un .gt. ZERO) then
                      spzero = -ONE
                   else
                      spzero = un*dtdx
                   endif
                   acmprght = HALF*(-ONE - spzero )*dq(i,j,k,n)
                   qp(i,j,k,n) = q(i,j,k,n) + acmprght
                endif
                
                ! Left state
                if ((idir == 1 .and. i <= hi(1)) .or. &
                    (idir == 2 .and. j <= hi(2)) .or. &
                    (idir == 3 .and. k <= hi(3))) then

                   un = q(i,j,k,QUN)
                   if (un .ge. ZERO) then
                      spzero = un*dtdx
                   else
                      spzero = ONE
                   endif
                   acmpleft = HALF*(ONE - spzero )*dq(i,j,k,n)
                   qm(i+ix,j+iy,k+iz,n) = q(i,j,k,n) + acmpleft
                endif
             enddo

#if (AMREX_SPACEDIM == 1)
             if (fix_mass_flux_hi) qp(hi(1)+1,j,k,n) = q(hi(1)+1,j,k,n)
             if (fix_mass_flux_lo) qm(lo(1),j,k,n) = q(lo(1)-1,j,k,n)
#endif

          end do
       end do

    end do

  end subroutine trace_plm

end module trace_plm_module
