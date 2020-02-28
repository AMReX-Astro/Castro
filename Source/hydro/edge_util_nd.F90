module edge_util_module

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg

  implicit none

contains


  subroutine reset_edge_state_thermo(lo, hi, &
                                     qedge, qd_lo, qd_hi) &
                                     bind(C, name="reset_edge_state_thermo")

    use amrex_constants_module, only : ZERO, ONE, HALF

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, &
                                   QPRES, QREINT, QFS, QFX, &
#ifdef RADIATION
                                   QPTOT, &
#endif
                                   small_pres, small_temp, &
                                   transverse_use_eos, transverse_reset_rhoe

    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    real(rt), intent(inout) :: qedge(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    logical :: reset
    type (eos_t) :: eos_state
    integer :: ii, jj, kk
    real(rt) :: old_p_state

    !$gpu

    do kk = lo(3), hi(3)
       do jj = lo(2), hi(2)
          do ii = lo(1), hi(1)

             reset = .false.

#ifdef RADIATION
             old_p_state = qedge(ii,jj,kk,QPRES)
#endif

             if (transverse_reset_rhoe == 1) then
                ! if we are still negative, then we need to reset
                if (qedge(ii,jj,kk,QREINT) < ZERO) then
                   reset = .true.

                   eos_state % rho = qedge(ii,jj,kk,QRHO)
                   eos_state % T = small_temp
                   eos_state % xn(:) = qedge(ii,jj,kk,QFS:QFS-1+nspec)
                   eos_state % aux(:) = qedge(ii,jj,kk,QFX:QFX-1+naux)

                   call eos(eos_input_rt, eos_state)

                   qedge(ii,jj,kk,QREINT) = qedge(ii,jj,kk,QRHO)*eos_state % e
                   qedge(ii,jj,kk,QPRES) = eos_state % p
                endif

             end if

             if (transverse_use_eos == 1) then
                eos_state % rho = qedge(ii,jj,kk,QRHO)
                eos_state % e   = qedge(ii,jj,kk,QREINT) / qedge(ii,jj,kk,QRHO)
                eos_state % T   = small_temp
                eos_state % xn  = qedge(ii,jj,kk,QFS:QFS+nspec-1)
                eos_state % aux = qedge(ii,jj,kk,QFX:QFX+naux-1)

                call eos(eos_input_re, eos_state)

                qedge(ii,jj,kk,QREINT) = eos_state % e * eos_state % rho
                qedge(ii,jj,kk,QPRES) = max(eos_state % p, small_pres)
             end if

#ifdef RADIATION
             ! correct the total pressure (gas + radiation) with any
             ! change to the gas pressure state
             qedge(ii,jj,kk,QPTOT) = qedge(ii,jj,kk,QPTOT) + (qedge(ii,jj,kk,QPRES) - old_p_state)
#endif

          end do
       end do
    end do

  end subroutine reset_edge_state_thermo


  subroutine edge_state_temp_to_pres(lo, hi, &
                                     qm, qm_lo, qm_hi, &
                                     qp, qp_lo, qp_hi) bind(C, name="edge_state_temp_to_pres")

    use castro_error_module
    use meth_params_module, only : NQ, NVAR, NGDNV, GDPRES, &
                                   UTEMP, UMX, &
                                   plm_well_balanced, QPRES, &
                                   QTEMP, QFS, QFX, QREINT, QRHO, QU, QV, QW, &
                                   first_order_hydro, hybrid_riemann, &
                                   ppm_temp_fix, const_grav
    use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
    use slope_module, only : uslope
    use amrex_fort_module, only : rt => amrex_real
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use network, only : nspec, naux
    use prob_params_module, only : dg, coord_type, Symmetry, physbc_lo, physbc_hi

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1), qm_lo(2):qm_hi(2), qm_lo(3):qm_hi(3), NQ)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), NQ)

    integer :: i, j, k
    type (eos_t) :: eos_state

    !$gpu

    ! use T to define p

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! we just got the extremes corresponding to a particular cell-center, but now
             ! we need to assign them to interfaces

             eos_state%rho    = qp(i,j,k,QRHO)
             eos_state%T      = qp(i,j,k,QTEMP)
             eos_state%xn(:)  = qp(i,j,k,QFS:QFS-1+nspec)
             eos_state%aux(:) = qp(i,j,k,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

             qp(i,j,k,QPRES) = eos_state%p
             qp(i,j,k,QREINT) = qp(i,j,k,QRHO)*eos_state%e
             ! should we try to do something about Gamma_! on interface?

             eos_state%rho    = qm(i,j,k,QRHO)
             eos_state%T      = qm(i,j,k,QTEMP)
             eos_state%xn(:)  = qm(i,j,k,QFS:QFS-1+nspec)
             eos_state%aux(:) = qm(i,j,k,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

             qm(i,j,k,QPRES) = eos_state%p
             qm(i,j,k,QREINT) = qm(i,j,k,QRHO)*eos_state%e
             ! should we try to do something about Gamma_! on interface?

          end do
       end do
    end do

  end subroutine edge_state_temp_to_pres

end module edge_util_module
