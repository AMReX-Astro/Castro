module timestep_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  public

contains

#ifdef REACTIONS

  subroutine ca_estdt_burning(lo, hi, &
                              snew, sn_lo, sn_hi, &
                              rnew, rn_lo, rn_hi, &
                              dx, dt) &
                              bind(C, name="ca_estdt_burning")
    ! Reactions-limited timestep
    !
    ! .. note::
    !    Binds to C function ``ca_estdt_burning``

    use amrex_constants_module, only: HALF, ONE
    use network, only: nspec, naux, aion
    use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, UFS, dtnuc_e, dtnuc_X, dtnuc_X_threshold
    use prob_params_module, only : dim
#if naux > 0
    use meth_params_module, only : UFX
#endif
    use actual_rhs_module, only: actual_rhs
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use react_util_module, only: okay_to_burn_type ! function
    use burn_type_module, only : burn_t, net_ienuc, burn_to_eos, eos_to_burn, neqs
    use temperature_integration_module, only: self_heat
    use amrex_fort_module, only : rt => amrex_real
    use extern_probin_module, only: small_x
    use reduction_module, only: reduce_min

    implicit none

    integer,  intent(in) :: sn_lo(3), sn_hi(3)
    integer,  intent(in) :: rn_lo(3), rn_hi(3)
    integer,  intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: snew(sn_lo(1):sn_hi(1),sn_lo(2):sn_hi(2),sn_lo(3):sn_hi(3),NVAR)
    real(rt), intent(in) :: rnew(rn_lo(1):rn_hi(1),rn_lo(2):rn_hi(2),rn_lo(3):rn_hi(3),nspec+3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: dt

    real(rt)      :: e, X(nspec), dedt, dXdt(nspec)
    integer       :: i, j, k
    integer       :: n

    type (burn_t) :: state_new
    real(rt) :: ydot(neqs)
    type (eos_t)  :: eos_state
    real(rt)      :: rhoninv
    real(rt) :: dt_tmp

    ! Set a floor on the minimum size of a derivative. This floor
    ! is small enough such that it will result in no timestep limiting.

    real(rt), parameter :: derivative_floor = 1.e-50_rt

    !$gpu

    ! We want to limit the timestep so that it is not larger than
    ! dtnuc_e * (e / (de/dt)).  If the timestep factor dtnuc is
    ! equal to 1, this says that we don't want the
    ! internal energy to change by any more than its current
    ! magnitude in the next timestep.
    !
    ! If dtnuc is less than one, it controls the fraction we will
    ! allow the internal energy to change in this timestep due to
    ! nuclear burning, provided that the last timestep's burning is a
    ! good estimate for the current timestep's burning.
    !
    ! We do the same thing for the species, using a timestep
    ! limiter dtnuc_X * (X_k / (dX_k/dt)). To prevent changes
    ! due to trace isotopes that we probably are not interested in,
    ! only apply the limiter to species with an abundance greater
    ! than a user-specified threshold.
    !
    ! To estimate de/dt and dX/dt, we are going to call the RHS of the
    ! burner given the current state data. We need to do an EOS
    ! call before we do the RHS call so that we have accurate
    ! values for the thermodynamic data like abar, zbar, etc.
    ! But we will call in (rho, T) mode, which is inexpensive.

    if (dtnuc_e > 1.e199_rt .and. dtnuc_X > 1.e199_rt) return

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoninv = ONE / snew(i,j,k,URHO)

             state_new % rho = snew(i,j,k,URHO)
             state_new % T   = snew(i,j,k,UTEMP)
             state_new % e   = snew(i,j,k,UEINT) * rhoninv
             state_new % xn  = snew(i,j,k,UFS:UFS+nspec-1) * rhoninv
#if naux > 0
             state_new % aux = snew(i,j,k,UFX:UFX+naux-1) * rhoninv
#endif

             if (.not. okay_to_burn_type(state_new)) cycle

             e    = state_new % e
             X    = max(state_new % xn, small_x)

             call burn_to_eos(state_new, eos_state)
             call eos(eos_input_rt, eos_state)
             call eos_to_burn(eos_state, state_new)

             state_new % dx = minval(dx(1:dim))

#ifndef SIMPLIFIED_SDC
             state_new % self_heat = self_heat
#else
             state_new % self_heat = .true.
#endif
             call actual_rhs(state_new, ydot)

             dedt = ydot(net_ienuc)
             dXdt = ydot(1:nspec) * aion

             ! Apply a floor to the derivatives. This ensures that we don't
             ! divide by zero; it also gives us a quick method to disable
             ! the timestep limiting, because the floor is small enough
             ! that the implied timestep will be very large, and thus
             ! ignored compared to other limiters.

             dedt = max(abs(dedt), derivative_floor)

             do n = 1, nspec
                if (X(n) .ge. dtnuc_X_threshold) then
                   dXdt(n) = max(abs(dXdt(n)), derivative_floor)
                else
                   dXdt(n) = derivative_floor
                end if
             end do

             dt_tmp = min(dtnuc_e * e / dedt, dtnuc_X * minval(X / dXdt))

             call reduce_min(dt, dt_tmp)

          enddo
       enddo
    enddo

  end subroutine ca_estdt_burning
#endif

end module timestep_module
