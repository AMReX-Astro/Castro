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

end module edge_util_module
