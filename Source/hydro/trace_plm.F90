module trace_plm_module

  use amrex_error_module, only : amrex_error
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg

  implicit none

  private

  public trace_plm

contains

  subroutine trace_plm(lo, hi, &
                       idir, q, q_lo, q_hi, &
                       qaux, qa_lo, qa_hi, &
                       dq, dq_lo, dq_hi, &
                       qm, qp, qpd_lo, qpd_hi, &
#if (AMREX_SPACEDIM < 3)
                       dloga, dloga_lo, dloga_hi, &
#endif
                       SrcQ, src_lo, Src_hi, &
                       vlo, vhi, domlo, domhi, &
                       dx, dt)

    ! here, lo and hi are the range we loop over -- this can include ghost cells
    ! vlo and vhi are the bounds of the valid box (no ghost cells)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, QVAR, QRHO, QU, QV, QW, QC, &
                                   QREINT, QPRES, &
                                   npassive, qpass_map, small_dens, small_pres, &
                                   ppm_type, fix_mass_flux
    use amrex_constants_module
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: dq_lo(3), dq_hi(3)
    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) ::  dq(dq_lo(1):dq_hi(1),dq_lo(2):dq_hi(2),dq_lo(3):dq_hi(3),NQ,AMREX_SPACEDIM)

    real(rt), intent(inout) :: qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    real(rt), intent(in) :: dx(3), dt

    ! Local variables
    integer :: i, j, k, n, ipassive

    real(rt) :: dtdx
    real(rt) :: cc, csq, rho, un, ut, utt, p, rhoe
    real(rt) :: drho, dun, dut, dutt, dp, drhoe

    real(rt) :: enth, alpham, alphap, alpha0r, alpha0e
    real(rt) :: alpha0ut, alpha0utt
    real(rt) :: apright, amright, azrright, azeright
    real(rt) :: azut1rght, azutt1rght
    real(rt) :: apleft, amleft, azrleft, azeleft
    real(rt) :: azut1left, azutt1left
    real(rt) :: acmprght, acmpleft
    real(rt) :: spzero
    real(rt) :: rho_ref, un_ref, ut_ref, utt_ref, p_ref, rhoe_ref
    real(rt) :: e(3)
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp
    logical :: fix_mass_flux_lo, fix_mass_flux_hi

    integer :: QUN, QUT, QUTT
    real(rt) :: ref_fac, trace_fac1, trace_fac2, trace_fac3

    dtdx = dt/dx(idir)

#ifndef AMREX_USE_CUDA
    if (ppm_type .ne. 0) then
       print *,'Oops -- shouldnt be in tracexy with ppm_type != 0'
       call amrex_error("Error:: trace_3d.f90 :: tracexy")
    end if
#endif

    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (vlo(1) == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (vhi(1) == domhi(1))

    if (idir == 1) then
       QUN = QU
       QUT = QV
       QUTT = QW
    else if (idir == 2) then
       QUN = QV
       QUT = QW
       QUTT = QU
    else if (idir == 3) then
       QUN = QW
       QUT = QU
       QUTT = QV
    endif

    ! Compute left and right traced states

    ! construct the right state on the i interface

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             cc = qaux(i,j,k,QC)
             csq = cc**2
             rho = q(i,j,k,QRHO)
             un = q(i,j,k,QUN)
             ut = q(i,j,k,QUT)
             utt = q(i,j,k,QUTT)
             p = q(i,j,k,QPRES)
             rhoe = q(i,j,k,QREINT)
             enth = (rhoe+p)/(rho*csq)

             drho = dq(i,j,k,QRHO,idir)
             dun = dq(i,j,k,QUN,idir)
             dut = dq(i,j,k,QUT,idir)
             dutt = dq(i,j,k,QUTT,idir)
             dp = dq(i,j,k,QPRES,idir)
             drhoe = dq(i,j,k,QREINT,idir)

             alpham = HALF*(dp/(rho*cc) - dun)*(rho/cc)
             alphap = HALF*(dp/(rho*cc) + dun)*(rho/cc)
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth
             alpha0ut = dut
             alpha0utt = dutt

             e(1) = un - cc
             e(2) = un
             e(3) = un + cc

             ref_fac = HALF*(ONE + dtdx*min(e(1), ZERO))
             rho_ref = rho - ref_fac*drho
             un_ref = un - ref_fac*dun
             ut_ref = ut - ref_fac*dut
             utt_ref = utt - ref_fac*dutt
             p_ref = p - ref_fac*dp
             rhoe_ref = rhoe - ref_fac*drhoe

             ! this is -(1/2) ( 1 + dt/dx lambda) (l . dq) r
             trace_fac1 = ZERO !  FOURTH*dtdx*(e(1) - e(1))*(ONE - sign(ONE, e(1)))
             trace_fac2 = FOURTH*dtdx*(e(1) - e(2))*(ONE - sign(ONE, e(2)))
             trace_fac3 = FOURTH*dtdx*(e(1) - e(3))*(ONE - sign(ONE, e(3)))

             apright = trace_fac3*alphap
             amright = trace_fac1*alpham

             azrright = trace_fac2*alpha0r
             azeright = trace_fac2*alpha0e
             azut1rght = trace_fac2*alpha0ut
             azutt1rght = trace_fac2*alpha0utt

             if ((idir == 1 .and. i >= vlo(1)) .or. &
                 (idir == 2 .and. j >= vlo(2)) .or. &
                 (idir == 3 .and. k >= vlo(3))) then

                qp(i,j,k,QRHO) = max(small_dens, rho_ref + apright + amright + azrright)
                qp(i,j,k,QUN) = un_ref + (apright - amright)*cc/rho
                qp(i,j,k,QUT) = ut_ref + azut1rght
                qp(i,j,k,QUTT) = utt_ref + azutt1rght
                qp(i,j,k,QPRES) = max(small_pres, p_ref + (apright + amright)*csq)
                qp(i,j,k,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright

                ! add the source terms 
                qp(i,j,k,QRHO  ) = qp(i,j,k,QRHO  ) + HALF*dt*srcQ(i,j,k,QRHO)
                qp(i,j,k,QRHO  ) = max(small_dens, qp(i,j,k,QRHO))
                qp(i,j,k,QUN   ) = qp(i,j,k,QUN   ) + HALF*dt*srcQ(i,j,k,QUN)
                qp(i,j,k,QUT   ) = qp(i,j,k,QUT   ) + HALF*dt*srcQ(i,j,k,QUT)
                qp(i,j,k,QUTT  ) = qp(i,j,k,QUTT  ) + HALF*dt*srcQ(i,j,k,QUTT)
                qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + HALF*dt*srcQ(i,j,k,QREINT)
                qp(i,j,k,QPRES ) = qp(i,j,k,QPRES ) + HALF*dt*srcQ(i,j,k,QPRES)

             end if

             ! now construct the left state on the i+1 interface

             ref_fac = HALF*(ONE - dtdx*max(e(3), ZERO))
             rho_ref = rho + ref_fac*drho
             un_ref = un + ref_fac*dun
             ut_ref = ut + ref_fac*dut
             utt_ref = utt + ref_fac*dutt
             p_ref = p + ref_fac*dp
             rhoe_ref = rhoe + ref_fac*drhoe

             trace_fac1 = FOURTH*dtdx*(e(3) - e(1))*(ONE + sign(ONE, e(1)))
             trace_fac2 = FOURTH*dtdx*(e(3) - e(2))*(ONE + sign(ONE, e(2)))
             trace_fac3 = ZERO !  FOURTH*dtdx*(e(3) - e(3))*(ONE + sign(ONE, e(3)))

             apleft = trace_fac3*alphap
             amleft = trace_fac1*alpham

             azrleft = trace_fac2*alpha0r
             azeleft = trace_fac2*alpha0e
             azut1left = trace_fac2*alpha0ut
             azutt1left = trace_fac2*alpha0utt

             if (idir == 1 .and. i <= vhi(1)) then
                qm(i+1,j,k,QRHO) = max(small_dens, rho_ref + apleft + amleft + azrleft)
                qm(i+1,j,k,QUN) = un_ref + (apleft - amleft)*cc/rho
                qm(i+1,j,k,QUT) = ut_ref + azut1left
                qm(i+1,j,k,QUTT) = utt_ref + azutt1left
                qm(i+1,j,k,QPRES) = max(small_pres, p_ref + (apleft + amleft)*csq)
                qm(i+1,j,k,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft

                ! add the source terms
                qm(i+1,j,k,QRHO) = max(small_dens, qm(i+1,j,k,QRHO) + HALF*dt*srcQ(i,j,k,QRHO))
                qm(i+1,j,k,QUN) = qm(i+1,j,k,QUN) + HALF*dt*srcQ(i,j,k,QUN)
                qm(i+1,j,k,QUT) = qm(i+1,j,k,QUT) + HALF*dt*srcQ(i,j,k,QUT)
                qm(i+1,j,k,QUTT) = qm(i+1,j,k,QUTT) + HALF*dt*srcQ(i,j,k,QUTT)
                qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + HALF*dt*srcQ(i,j,k,QREINT)
                qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES ) + HALF*dt*srcQ(i,j,k,QPRES)

             else if (idir == 2 .and. j <= vhi(2)) then
                qm(i,j+1,k,QRHO) = max(small_dens, rho_ref + apleft + amleft + azrleft)
                qm(i,j+1,k,QUN) = un_ref + (apleft - amleft)*cc/rho
                qm(i,j+1,k,QUT) = ut_ref + azut1left
                qm(i,j+1,k,QUTT) = utt_ref + azutt1left
                qm(i,j+1,k,QPRES) = max(small_pres, p_ref + (apleft + amleft)*csq)
                qm(i,j+1,k,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft

                ! add the source terms
                qm(i,j+1,k,QRHO) = max(small_dens, qm(i,j+1,k,QRHO) + HALF*dt*srcQ(i,j,k,QRHO))
                qm(i,j+1,k,QUN) = qm(i,j+1,k,QUN) + HALF*dt*srcQ(i,j,k,QUN)
                qm(i,j+1,k,QUT) = qm(i,j+1,k,QUT) + HALF*dt*srcQ(i,j,k,QUT)
                qm(i,j+1,k,QUTT) = qm(i,j+1,k,QUTT) + HALF*dt*srcQ(i,j,k,QUTT)
                qm(i,j+1,k,QREINT) = qm(i,j+1,k,QREINT) + HALF*dt*srcQ(i,j,k,QREINT)
                qm(i,j+1,k,QPRES) = qm(i,j+1,k,QPRES ) + HALF*dt*srcQ(i,j,k,QPRES)

             else if (idir == 3 .and. k <= vhi(3)) then
                qm(i,j,k+1,QRHO) = max(small_dens, rho_ref + apleft + amleft + azrleft)
                qm(i,j,k+1,QUN) = un_ref + (apleft - amleft)*cc/rho
                qm(i,j,k+1,QUT) = ut_ref + azut1left
                qm(i,j,k+1,QUTT) = utt_ref + azutt1left
                qm(i,j,k+1,QPRES) = max(small_pres, p_ref + (apleft + amleft)*csq)
                qm(i,j,k+1,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft

                ! add the source terms
                qm(i,j,k+1,QRHO) = max(small_dens, qm(i,j,k+1,QRHO) + HALF*dt*srcQ(i,j,k,QRHO))
                qm(i,j,k+1,QUN) = qm(i,j,k+1,QUN) + HALF*dt*srcQ(i,j,k,QUN)
                qm(i,j,k+1,QUT) = qm(i,j,k+1,QUT) + HALF*dt*srcQ(i,j,k,QUT)
                qm(i,j,k+1,QUTT) = qm(i,j,k+1,QUTT) + HALF*dt*srcQ(i,j,k,QUTT)
                qm(i,j,k+1,QREINT) = qm(i,j,k+1,QREINT) + HALF*dt*srcQ(i,j,k,QREINT)
                qm(i,j,k+1,QPRES) = qm(i,j,k+1,QPRES ) + HALF*dt*srcQ(i,j,k,QPRES)

             endif

#if (AMREX_SPACEDIM < 3)
             ! geometry source terms -- these only apply to the x-states
             if (idir == 1 .and. dloga(i,j,k) /= ZERO) then
                courn = dtdx*(cc + abs(un))
                eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k)))
                dlogatmp = min(eta, ONE)*dloga(i,j,k)
                sourcr = -HALF*dt*rho*dlogatmp*un
                sourcp = sourcr*csq
                source = sourcp*enth

                if (i <= vhi(1)) then
                   qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr
                   qm(i+1,j,k,QRHO) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                end if
                if (i >= vlo(1)) then
                   qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr
                   qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                end if
             end if
#endif

#if (AMREX_SPACEDIM == 1)
             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_lo) then
                qm(vlo(1),j,k,QRHO  ) = q(domlo(1)-1,j,k,QRHO)
                qm(vlo(1),j,k,QUN   ) = q(domlo(1)-1,j,k,QUN )
                qm(vlo(1),j,k,QPRES ) = q(domlo(1)-1,j,k,QPRES)
                qm(vlo(1),j,k,QREINT) = q(domlo(1)-1,j,k,QREINT)
             end if

             ! Enforce constant mass flux rate if specified
             if (fix_mass_flux_hi) then
                qp(vhi(1)+1,j,k,QRHO  ) = q(domhi(1)+1,j,k,QRHO)
                qp(vhi(1)+1,j,k,QUN    ) = q(domhi(1)+1,j,k,QUN )
                qp(vhi(1)+1,j,k,QPRES ) = q(domhi(1)+1,j,k,QPRES)
                qp(vhi(1)+1,j,k,QREINT) = q(domhi(1)+1,j,k,QREINT)
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

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Right state
                if ((idir == 1 .and. i >= vlo(1)) .or. &
                    (idir == 2 .and. j >= vlo(2)) .or. &
                    (idir == 3 .and. k >= vlo(3))) then

                   un = q(i,j,k,QUN)
                   spzero = merge(-ONE, un*dtdx, un >= ZERO)
                   acmprght = HALF*(-ONE - spzero)*dq(i,j,k,n,idir)
                   qp(i,j,k,n) = q(i,j,k,n) + acmprght + HALF*dt*srcQ(i,j,k,n)
                endif

                ! Left state
                un = q(i,j,k,QUN)
                spzero = merge(un*dtdx, ONE, un >= ZERO)
                acmpleft = HALF*(ONE - spzero )*dq(i,j,k,n,idir)

                if (idir == 1 .and. i <= vhi(1)) then
                   qm(i+1,j,k,n) = q(i,j,k,n) + acmpleft + HALF*dt*srcQ(i,j,k,n)
                else if (idir == 2 .and. j <= vhi(2)) then
                   qm(i,j+1,k,n) = q(i,j,k,n) + acmpleft + HALF*dt*srcQ(i,j,k,n)
                else if (idir == 3 .and. k <= vhi(3)) then
                   qm(i,j,k+1,n) = q(i,j,k,n) + acmpleft + HALF*dt*srcQ(i,j,k,n)
                endif
             enddo

#if (AMREX_SPACEDIM == 1)
             if (fix_mass_flux_hi) qp(vhi(1)+1,j,k,n) = q(vhi(1)+1,j,k,n)
             if (fix_mass_flux_lo) qm(vlo(1),j,k,n) = q(vlo(1)-1,j,k,n)
#endif

          end do
       end do

    end do

  end subroutine trace_plm

end module trace_plm_module
