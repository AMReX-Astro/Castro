! advection routines in support of method of lines integration
!
module fourth_order_hydro

  implicit none

contains

  subroutine avisc(lo, hi, &
                   q, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   dx, avis, a_lo, a_hi, idir) bind(C, name="avisc")

    ! this computes the *face-centered* artifical viscosity using the
    ! 4th order expression from McCorquodale & Colella (Eq. 35)

    use meth_params_module, only: QU, QV, QW, QC, NQ, NQAUX
    use amrex_constants_module, only: HALF, FOURTH, ONE, ZERO
    use prob_params_module, only: dg
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: a_lo(3), a_hi(3)
    integer, intent(in), value :: idir
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(inout) :: avis(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    real(rt) :: coeff, cmin

    real(rt), parameter :: beta = 0.3_rt

    integer :: i, j, k

    real(rt) :: dxinv, dyinv, dzinv

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    do k = lo(3), hi(3)+dg(3)
       do j = lo(2), hi(2)+dg(2)
          do i = lo(1), hi(1)+1

             if (idir == 1) then

                ! normal direction
                avis(i,j,k) = (q(i,j,k,QU) - q(i-1,j,k,QU)) * dxinv
#if BL_SPACEDIM >= 2
                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i,j+1,k,QV) - q(i,j-1,k,QV) + &
                     q(i-1,j+1,k,QV) - q(i-1,j-1,k,QV)) * dyinv
#endif
#if BL_SPACEDIM >= 3
                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i,j,k+1,QW) - q(i,j,k-1,QW) + &
                     q(i-1,j,k+1,QW) - q(i-1,j,k-1,QW)) * dzinv
#endif

                cmin = min(qaux(i,j,k,QC), qaux(i-1,j,k,QC))

             else if (idir == 2) then

                ! normal direction
                avis(i,j,k) = (q(i,j,k,QV) - q(i,j-1,k,QV)) * dyinv

                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i+1,j,k,QU) - q(i-1,j,k,QU) + &
                     q(i+1,j-1,k,QU) - q(i-1,j-1,k,QU)) * dxinv

#if BL_SPACEDIM >= 3
                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i,j,k+1,QW) - q(i,j,k-1,QW) + &
                     q(i,j-1,k+1,QW) - q(i,j-1,k-1,QW)) * dzinv
#endif

                cmin = min(qaux(i,j,k,QC), qaux(i,j-1,k,QC))

             else

                ! normal direction
                avis(i,j,k) = (q(i,j,k,QW) - q(i,j,k-1,QW)) * dzinv

                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i,j+1,k,QV) - q(i,j-1,k,QV) + &
                     q(i,j+1,k-1,QV) - q(i,j-1,k-1,QV)) * dyinv

                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i+1,j,k,QU) - q(i-1,j,k,QU) + &
                     q(i+1,j,k-1,QU) - q(i-1,j,k-1,QU)) * dxinv

                cmin = min(qaux(i,j,k,QC), qaux(i,j,k-1,QC))

             endif

             ! MC Eq. 36
             coeff = min(ONE, (dx(idir)*avis(i,j,k))**2/(beta * cmin**2))

             if (avis(i,j,k) < ZERO) then
                avis(i,j,k) = dx(idir)*avis(i,j,k)*coeff
             else
                avis(i,j,k) = ZERO
             endif

          enddo
       enddo
    enddo

  end subroutine avisc

#ifdef DIFFUSION
  subroutine add_diffusive_flux(lo, hi, &
                                q, q_lo, q_hi, ncomp, temp_comp, &
                                qint, qi_lo, qi_hi, &
                                F, F_lo, F_hi, &
                                dx, idir, is_avg) bind(C, name="add_diffusive_flux")

    ! add the diffusive flux to the energy fluxes
    !

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NQ, NVAR, &
                                   URHO, &
                                   UEDEN, UEINT, UTEMP, &
                                   QRHO, QREINT, QFS
    use eos_type_module, only : eos_t, eos_input_re
    use conductivity_module, only : conducteos
    use network, only : nspec
    use castro_error_module, only : castro_error

    integer, intent(in), value :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qi_lo(3), qi_hi(3)
    integer, intent(in) :: F_lo(3), F_hi(3)
    integer, intent(in), value :: ncomp, temp_comp

    real(rt), intent(in) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), ncomp)
    real(rt), intent(in) :: qint(qi_lo(1):qi_hi(1), qi_lo(2):qi_hi(2), qi_lo(3):qi_hi(3), NQ)
    real(rt), intent(out) :: F(F_lo(1):F_hi(1), F_lo(2):F_hi(2), F_lo(3):F_hi(3), NVAR)
    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dx(3)
    integer, intent(in), value :: is_avg

    integer :: i, j, k

    type(eos_t) :: eos_state
    real(rt) :: dTdx

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state % rho = qint(i,j,k,QRHO)
             eos_state % T = q(i,j,k,temp_comp)   ! initial guess
             eos_state % e = qint(i,j,k,QREINT) / qint(i,j,k,QRHO)
             eos_state % xn(:) = qint(i,j,k,QFS:QFS-1+nspec)

             call conducteos(eos_input_re, eos_state)

             if (idir == 1) then

                if (is_avg == 0) then
                   ! we are working with the cell-center state
                   dTdx = (-q(i+1,j,k,temp_comp) + 27*q(i,j,k,temp_comp) - &
                        27*q(i-1,j,k,temp_comp) + q(i-2,j,k,temp_comp))/(24.0_rt * dx(1))

                else
                   ! we are working with the cell-average state
                   dTdx = (-q(i+1,j,k,temp_comp) + 15*q(i,j,k,temp_comp) - &
                        15*q(i-1,j,k,temp_comp) + q(i-2,j,k,temp_comp))/(12.0_rt * dx(1))
                end if

             else if (idir == 2) then

                if (is_avg == 0) then
                   ! we are working with the cell-center state
                   dTdx = (-q(i,j+1,k,temp_comp) + 27*q(i,j,k,temp_comp) - &
                        27*q(i,j-1,k,temp_comp) + q(i,j-2,k,temp_comp))/(24.0_rt * dx(2))

                else
                   ! we are working with the cell-average state
                   dTdx = (-q(i,j+1,k,temp_comp) + 15*q(i,j,k,temp_comp) - &
                        15*q(i,j-1,k,temp_comp) + q(i,j-2,k,temp_comp))/(12.0_rt * dx(2))
                end if

             else

                if (is_avg == 0) then
                   ! we are working with the cell-center state
                   dTdx = (-q(i,j,k+1,temp_comp) + 27*q(i,j,k,temp_comp) - &
                        27*q(i,j,k-1,temp_comp) + q(i,j,k-2,temp_comp))/(24.0_rt * dx(3))

                else
                   ! we are working with the cell-average state
                   dTdx = (-q(i,j,k+1,temp_comp) + 15*q(i,j,k,temp_comp) - &
                        15*q(i,j,k-1,temp_comp) + q(i,j,k-2,temp_comp))/(12.0_rt * dx(3))
                end if

             endif

             F(i,j,k,UEINT) = F(i,j,k,UEINT) - eos_state % conductivity * dTdx
             F(i,j,k,UEDEN) = F(i,j,k,UEDEN) - eos_state % conductivity * dTdx

          end do
       end do
    end do

  end subroutine add_diffusive_flux
#endif

end module fourth_order_hydro
