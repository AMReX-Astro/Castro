module react_util_module

  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec_evolve
  implicit none

  ! indices for working with the intermediate "w" state for the SDC Jacobian
  integer, parameter :: iwrho = 0
  integer, parameter :: iwfs = 1
  integer, parameter :: iwT = nspec_evolve + 1

contains



  pure function okay_to_burn(state) result(burn_flag)

    use meth_params_module, only : NVAR, URHO, UTEMP, &
                                   react_T_min, react_T_max, react_rho_min, react_rho_max
    implicit none

    real(rt), intent(in) :: state(NVAR)

    logical :: burn_flag

    !$gpu

    burn_flag = .true.

    if (state(UTEMP) < react_T_min .or. state(UTEMP) > react_T_max .or. &
        state(URHO) < react_rho_min .or. state(URHO) > react_rho_max) then
       burn_flag = .false.
    end if

    return

  end function okay_to_burn


  subroutine single_zone_react_source(state, R, i, j, k, burn_state)

    use burn_type_module, only : burn_t, net_ienuc, neqs, burn_to_eos, eos_to_burn
    use network, only : nspec, nspec_evolve, aion
    use eos_module, only : eos
    use eos_type_module, only: eos_t, eos_input_re, eos_get_small_temp

    use meth_params_module, only : NVAR, URHO, UTEMP, UEDEN, UEINT, UFS, UFX
    use amrex_constants_module, only : ZERO, HALF, ONE
    use actual_rhs_module
    use extern_probin_module, only : SMALL_X_SAFE, MAX_TEMP

    implicit none

    integer, intent(in) :: i, j, k
    real(rt), intent(in) :: state(NVAR)
    real(rt), intent(out) :: R(NVAR)
    type(burn_t), intent(inout) :: burn_state
    real(rt) :: ydot(neqs)

    type(eos_t) :: eos_state
    real(rt) :: rhoInv
    integer :: n
    real(rt) :: small_temp

    !$gpu

    rhoInv = ONE / state(URHO)

    burn_state % rho = state(URHO)
    burn_state % T   = state(UTEMP)
    burn_state % e = state(UEINT) * rhoInv

    do n = 1, nspec
       burn_state % xn(n) = state(UFS+n-1) * rhoInv
    enddo

    burn_state % xn(1:nspec_evolve) = max(min(burn_state % xn(1:nspec_evolve), ONE), SMALL_X_SAFE) 

#if NAUX_NET > 0
    do n = 1, naux
       burn_state % aux(n) = state(UFX+n-1) * rhoInv
    enddo
#endif

    ! Ensure that the temperature going in is consistent with the internal energy.
    call burn_to_eos(burn_state, eos_state)
    call eos(eos_input_re, eos_state)
    call eos_to_burn(eos_state, burn_state)

    call eos_get_small_temp(small_temp)
    burn_state % T = min(MAX_TEMP, max(burn_state % T, small_temp))

    burn_state % i = i
    burn_state % j = j
    burn_state % k = k

    burn_state % self_heat = .false.

    call actual_rhs(burn_state, ydot)

    ! store the instantaneous R
    R(:) = ZERO

    ! species rates come back in terms of molar fractions
    R(UFS:UFS-1+nspec_evolve) = &
         state(URHO) * aion(1:nspec_evolve) * ydot(1:nspec_evolve)

    R(UEDEN) = state(URHO) * ydot(net_ienuc)
    R(UEINT) = state(URHO) * ydot(net_ienuc)

  end subroutine single_zone_react_source


  subroutine single_zone_jac(state, burn_state, dRdw)

    ! we assume that we are coming in with a valid burn_state, e.g., as called
    ! from single_zone_react_source

    use burn_type_module, only : burn_t, net_ienuc, net_itemp, neqs
    use network, only : nspec, nspec_evolve, aion, aion_inv
    use meth_params_module, only : NVAR, URHO, &
                                   sdc_use_analytic_jac
    use amrex_constants_module, only : ZERO, HALF, ONE
    use actual_rhs_module
    use numerical_jac_module
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    implicit none

    real(rt), intent(in) :: state(NVAR)
    type(burn_t), intent(inout) :: burn_state
    real(rt), intent(out) :: dRdw(0:nspec_evolve+1, 0:nspec_evolve+1)

    integer :: m, n
    type(burn_t) :: burn_state_pert

    ! for computing a numerical derivative
    real(rt) :: eps = 1.e-8_rt
    real(rt) :: jac(neqs, neqs)
    real(rt) :: ydot(neqs), ydot_pert(neqs)

#ifdef SIMPLIFIED_SDC
#ifndef AMREX_USE_GPU
    call castro_error("we shouldn't be here with the simplified SDC method (USE_SIMPLIFIED_SDC=TRUE)")
#endif
#else

    call actual_rhs(burn_state, ydot)

    if (sdc_use_analytic_jac == 0) then
       ! note the numerical Jacobian will be returned in terms of X
       call numerical_jac(burn_state, jac)
    else
       call actual_jac(burn_state, jac)

       ! The Jacobian from the nets is in terms of dYdot/dY, but we want
       ! it was dXdot/dX, so convert here.
       do n = 1, nspec_evolve
          jac(n,:) = jac(n,:) * aion(n)
          jac(:,n) = jac(:,n) * aion_inv(n)
       enddo

    endif
#endif

    ! at this point, our Jacobian should be entirely in terms of X,
    ! not Y.  Let's now fix the rhs terms themselves to be in terms of
    ! dX/dt and not dY/dt.
    ydot(1:nspec_evolve) = ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    ! Our jacobian, dR/dw has the form:
    !
    !  /      0                  0                  0                       0          \
    !  | d(rho X1dot)/drho  d(rho X1dot)/dX1   d(rho X1dit)/dX2   ...  d(rho X1dot)/dT |
    !  | d(rho X2dot)/drho  d(rho X2dot)/dX1   d(rho X2dot)/dX2   ...  d(rho X2dot)/dT |
    !  |   ...                                                                         |
    !  \ d(rho Edot)/drho   d(rho Edot)/dX1    d(rho Edot)/dX2    ...  d(rho Edot)/dT  /

    dRdw(:,:) = ZERO

    ! now perturb density and call the RHS to compute the derivative wrt rho
    ! species rates come back in terms of molar fractions
    call copy_burn_t(burn_state_pert, burn_state)
    burn_state_pert % rho = burn_state % rho * (ONE + eps)

    burn_state_pert % i = burn_state % i
    burn_state_pert % j = burn_state % j
    burn_state_pert % k = burn_state % k

    call actual_rhs(burn_state_pert, ydot_pert)

    ! make the rates dX/dt and not dY/dt
    ydot_pert(1:nspec_evolve) = ydot_pert(1:nspec_evolve) * aion(1:nspec_evolve)

    ! fill the column of dRdw corresponding to the derivative
    ! with respect to rho
    do m = 1, nspec_evolve
       ! d( d(rho X_m)/dt)/drho
       dRdw(m, iwrho) = ydot(m) + &
            state(URHO) * (ydot_pert(m) - ydot(m))/(eps * burn_state % rho)
    enddo

    ! d( d(rho E)/dt)/drho
    dRdw(nspec_evolve+1, iwrho) = ydot(net_ienuc) + &
         state(URHO) * (ydot_pert(net_ienuc) - ydot(net_ienuc))/(eps * burn_state % rho)

    ! fill the columns of dRdw corresponding to each derivative
    ! with respect to species mass fraction
    do n = 1, nspec_evolve
       dRdw(0, iwfs-1+n) = ZERO  ! density source

       do m = 1, nspec_evolve
          ! d( d(rho X_m)/dt)/dX_n
          dRdw(m, iwfs-1+n) = state(URHO) * jac(m, n)
       enddo

       ! d( d(rho E)/dt)/dX_n
       dRdw(nspec_evolve+1, iwfs-1+n) = state(URHO) * jac(net_ienuc, n)

    enddo

    ! now fill the column corresponding to derivatives with respect to
    ! temperature -- this column is iwT
    dRdw(0, iwT) = ZERO

    ! d( d(rho X_m)/dt)/dT
    do m = 1, nspec_evolve
       dRdw(m, iwT) = state(URHO) * jac(m, net_itemp)
    enddo

    ! d( d(rho E)/dt)/dT
    dRdw(nspec_evolve+1, iwT) = state(URHO) * jac(net_ienuc, net_itemp)


  end subroutine single_zone_jac

end module react_util_module
