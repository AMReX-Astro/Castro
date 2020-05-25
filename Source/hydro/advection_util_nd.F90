module advection_util_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_enforce_minimum_density(lo, hi, &
                                        state, s_lo, s_hi, &
                                        verbose) bind(c,name='ca_enforce_minimum_density')

    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEINT, UEDEN, UFS, small_temp, small_dens, npassive, upass_map
    use amrex_constants_module, only: ZERO
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif
    use amrex_fort_module, only: rt => amrex_real
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_module, only: eos
    use castro_util_module, only: position ! function
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only: linear_to_hybrid ! function
    use meth_params_module, only: UMR, UMP
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer,  intent(in   ), value :: verbose

    ! Local variables
    integer  :: i, j, k

    integer          :: n, ipassive
    type (eos_t)     :: eos_state

#ifdef HYBRID_MOMENTUM
    real(rt)         :: loc(3)
#endif

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,URHO) < small_dens) then

#ifndef AMREX_USE_CUDA
                if (verbose .gt. 0) then
                   print *,'   '
                   if (state(i,j,k,URHO) < ZERO) then
                      print *,'>>> RESETTING NEG.  DENSITY AT ', i, j, k
                   else
                      print *,'>>> RESETTING SMALL DENSITY AT ', i, j, k
                   endif
                   print *,'>>> FROM ', state(i,j,k,URHO), ' TO ', small_dens
                   print *,'>>> IN GRID ', lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
                   print *,'   '
                end if
#endif

                do ipassive = 1, npassive
                   n = upass_map(ipassive)
                   state(i,j,k,n) = state(i,j,k,n) * (small_dens / state(i,j,k,URHO))
                end do

                eos_state % rho = small_dens
                eos_state % T   = small_temp
                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) / small_dens
                eos_state % aux = state(i,j,k,UFS:UFS+naux-1) / small_dens

                call eos(eos_input_rt, eos_state)

                state(i,j,k,URHO ) = eos_state % rho
                state(i,j,k,UTEMP) = eos_state % T

                state(i,j,k,UMX  ) = ZERO
                state(i,j,k,UMY  ) = ZERO
                state(i,j,k,UMZ  ) = ZERO

                state(i,j,k,UEINT) = eos_state % rho * eos_state % e
                state(i,j,k,UEDEN) = state(i,j,k,UEINT)

#ifdef HYBRID_MOMENTUM
                loc = position(i, j, k)
                state(i,j,k,UMR:UMP) = linear_to_hybrid(loc, state(i,j,k,UMX:UMZ))
#endif

             end if

          end do
       end do
    end do

  end subroutine ca_enforce_minimum_density

end module advection_util_module
