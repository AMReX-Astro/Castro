module advection_util_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_enforce_minimum_density(lo, hi, &
                                        state, s_lo, s_hi, &
                                        verbose) bind(c,name='ca_enforce_minimum_density')

    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, small_dens, density_reset_method
    use amrex_constants_module, only: ZERO
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif
    use amrex_fort_module, only: rt => amrex_real
    use reduction_module, only: reduce_min

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer,  intent(in   ), value :: verbose

    ! Local variables
    integer  :: i, j, k
    integer  :: ii, jj, kk
    integer  :: i_set, j_set, k_set
    real(rt) :: max_dens, old_rho
    real(rt) :: uold(NVAR), unew(NVAR)
    integer  :: num_positive_zones

    !$gpu

    max_dens = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,URHO) .eq. ZERO) then

#ifndef AMREX_USE_GPU
                print *,'DENSITY EXACTLY ZERO AT CELL ', i, j, k
                print *,'  in grid ',lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
                call castro_error("Error :: ca_enforce_minimum_density")
#endif

             else if (state(i,j,k,URHO) < small_dens) then

                old_rho = state(i,j,k,URHO)

                if (density_reset_method == 1) then

                   ! Reset to the characteristics of the adjacent state with the highest density.

                   max_dens = state(i,j,k,URHO)
                   i_set = i
                   j_set = j
                   k_set = k
                   do kk = -1, 1
                      do jj = -1, 1
                         do ii = -1, 1

                            if (i+ii >= s_lo(1) .and. j+jj >= s_lo(2) .and. k+kk >= s_lo(3) .and. &
                                i+ii <= s_hi(1) .and. j+jj <= s_hi(2) .and. k+kk <= s_hi(3)) then

                               if (state(i+ii,j+jj,k+kk,URHO) .gt. max_dens) then

                                  i_set = i+ii
                                  j_set = j+jj
                                  k_set = k+kk
                                  max_dens = state(i_set,j_set,k_set,URHO)

                               end if

                            end if

                         end do
                      end do
                   end do

                   if (max_dens < small_dens) then

                      ! We could not find any nearby zones with sufficient density.

                      uold = state(i,j,k,:)
                      call reset_to_small_state(uold, [i, j, k], s_lo, s_hi, verbose)
                      state(i,j,k,:) = uold

                   else

                      uold = state(i,j,k,:)
                      unew = state(i_set,j_set,k_set,:)

                      call reset_to_zone_state(uold, unew, [i, j, k], s_lo, s_hi, verbose)

                      state(i,j,k,:) = uold

                   endif

                else if (density_reset_method == 2) then

                   ! Reset to the average of adjacent zones. The median is independently calculated for each variable.

                   num_positive_zones = 0
                   unew(:) = ZERO

                   do kk = -1, 1
                      do jj = -1, 1
                         do ii = -1, 1

                            if (i+ii >= s_lo(1) .and. j+jj >= s_lo(2) .and. k+kk >= s_lo(3) .and. &
                                i+ii <= s_hi(1) .and. j+jj <= s_hi(2) .and. k+kk <= s_hi(3)) then

                               if (state(i+ii,j+jj,k+kk,URHO) .ge. small_dens) then

                                  unew(:) = unew(:) + state(i+ii,j+jj,k+kk,:)
                                  num_positive_zones = num_positive_zones + 1

                               end if

                            end if

                         end do
                      end do
                   end do

                   if (num_positive_zones == 0) then

                      ! We could not find any nearby zones with sufficient density.

                      uold = state(i,j,k,:)
                      call reset_to_small_state(uold, [i, j, k], s_lo, s_hi, verbose)
                      state(i,j,k,:) = uold

                   else

                      uold = state(i,j,k,:)
                      unew(:) = unew(:) / num_positive_zones

                      call reset_to_zone_state(uold, unew, [i, j, k], s_lo, s_hi, verbose)

                      state(i,j,k,:) = uold

                   endif

#ifndef AMREX_USE_CUDA
                else

                   call castro_error("Unknown density_reset_method in subroutine ca_enforce_minimum_density.")
#endif
                endif

             end if

          end do
       end do
    end do

  end subroutine ca_enforce_minimum_density


  subroutine reset_to_small_state(state, idx, lo, hi, verbose)
    ! If no neighboring zones are above small_dens, our only recourse
    ! is to set the density equal to small_dens, and the temperature
    ! equal to small_temp. We set the velocities to zero,
    ! though any choice here would be arbitrary.
    !

    use amrex_constants_module, only: ZERO
    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEINT, UEDEN, UFS, small_temp, small_dens, npassive, upass_map
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_module, only: eos
    use castro_util_module, only: position ! function
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only: linear_to_hybrid ! function
    use meth_params_module, only: UMR, UMP
#endif

    use amrex_fort_module, only: rt => amrex_real
    implicit none

    real(rt)         :: state(NVAR)
    integer          :: idx(3), lo(3), hi(3), verbose

    integer          :: n, ipassive
    type (eos_t)     :: eos_state

#ifdef HYBRID_MOMENTUM
    real(rt)         :: loc(3)
#endif

    !$gpu

#ifndef AMREX_USE_CUDA
    if (verbose .gt. 0) then
       print *,'   '
       if (state(URHO) < ZERO) then
          print *,'>>> RESETTING NEG.  DENSITY AT ', idx(1), idx(2), idx(3)
       else
          print *,'>>> RESETTING SMALL DENSITY AT ', idx(1), idx(2), idx(3)
       endif
       print *,'>>> FROM ', state(URHO), ' TO ', small_dens
       print *,'>>> IN GRID ', lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
       print *,'   '
    end if
#endif

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       state(n) = state(n) * (small_dens / state(URHO))
    end do

    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = state(UFS:UFS+nspec-1) / small_dens
    eos_state % aux = state(UFS:UFS+naux-1) / small_dens

    call eos(eos_input_rt, eos_state)

    state(URHO ) = eos_state % rho
    state(UTEMP) = eos_state % T

    state(UMX  ) = ZERO
    state(UMY  ) = ZERO
    state(UMZ  ) = ZERO

    state(UEINT) = eos_state % rho * eos_state % e
    state(UEDEN) = state(UEINT)

#ifdef HYBRID_MOMENTUM
    loc = position(idx(1),idx(2),idx(3))
    state(UMR:UMP) = linear_to_hybrid(loc, state(UMX:UMZ))
#endif

  end subroutine reset_to_small_state



  subroutine reset_to_zone_state(state, input_state, idx, lo, hi, verbose)

    use amrex_constants_module, only: ZERO
    use meth_params_module, only: NVAR, URHO
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt) :: state(NVAR), input_state(NVAR)
    integer  :: idx(3), lo(3), hi(3), verbose

    !$gpu

#ifndef AMREX_USE_CUDA
    if (verbose .gt. 0) then
       if (state(URHO) < ZERO) then
          print *,'   '
          print *,'>>> RESETTING NEG.  DENSITY AT ',idx(1), idx(2), idx(3)
          print *,'>>> FROM ', state(URHO) ,' TO ', input_state(URHO)
          print *,'>>> IN GRID ', lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
          print *,'   '
       else
          print *,'   '
          print *,'>>> RESETTING SMALL DENSITY AT ', idx(1), idx(2), idx(3)
          print *,'>>> FROM ', state(URHO), ' TO ', input_state(URHO)
          print *,'>>> IN GRID ', lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
          print *,'   '
       end if
    end if
#endif

    state(:) = input_state(:)

  end subroutine reset_to_zone_state

end module advection_util_module
