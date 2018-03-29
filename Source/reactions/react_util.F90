module react_util_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine single_zone_react_source(state, R, i, j, k)

    use burn_type_module, only : burn_t
    use network, only : NSPEC
    use eos_module, only : eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only : NVAR, URHO, UTEMP, UEDEN, UEINT, UMX, UMZ, UFS, UFX, &
                                   dual_energy_eta3
    use bl_constants_module, only : ZERO, HALF, ONE
    use actual_rhs_module

    implicit none

    integer, intent(in) :: i, j, k
    real(rt), intent(in) :: state(NVAR)
    real(rt), intent(out) :: R(NVAR)

    type(burn_t) :: burn_state
    type(eos_t) :: eos_state
    real(rt) :: rhoInv, rho_e_K
    integer :: n

    rhoInv = ONE / state(URHO)

    burn_state % rho = state(URHO)
    burn_state % T   = state(UTEMP)

    rho_e_K = state(UEDEN) - HALF * rhoInv * sum(state(UMX:UMZ)**2)

    ! Dual energy formalism: switch between e and (E - K) depending on (E - K) / E.

    if ( rho_e_K / state(UEDEN) .gt. dual_energy_eta3 .and. rho_e_K .gt. ZERO ) then
       burn_state % e = rho_E_K * rhoInv
    else
       burn_state % e = state(UEINT) * rhoInv
    endif

    do n = 1, nspec
       burn_state % xn(n) = state(UFS+n-1) * rhoInv
    enddo

#if naux > 0
    do n = 1, naux
       burn_state % aux(n) = state(UFX+n-1) * rhoInv
    enddo
#endif

    ! Ensure that the temperature going in is consistent with the internal energy.
    call burn_to_eos(burn_state, eos_state)
    call eos(eos_input_re, eos_state)
    call eos_to_burn(eos_state, burn_state)

    burn_state % i = i
    burn_state % j = j
    burn_state % k = k

    call actual_rhs(burn_state)

    ! store the instantaneous R
    R(:) = ZERO

    ! species rates come back in terms of molar fractions
    R(UFS:UFS-1+nspec_evolve) = &
         state(URHO) * aion(1:nspec_evolve) * burn_state % ydot(1:nspec_evolve)

    R(UEDEN) = state(URHO) * burn_state % ydot(net_ienuc)
    R(UEINT) = state(URHO) * burn_state % ydot(net_ienuc)

  end subroutine single_zone_react_source


  subroutine single_zone_jac(burn_state, dRdw)

    ! we assume that we are coming in with a valid burn_state, e.g., as called
    ! from single_zone_react_source

    use burn_type_module, only : burn_t
    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UTEMP, UEDEN, UEINT, UMX, UMZ, UFS, UFX, &
                                   dual_energy_eta3
    use bl_constants_module, only : ZERO, HALF, ONE
    use actual_jac_module

    implicit none

    type(burn_t), intent(in) :: burn_state
    real(rt), intent(out) :: dRdw(nspec+2, nspec+2)


    call actual_jac(burn_state)

    ! store the instantaneous R
    dRdw(:,:) = ZERO

    ! now perturb density and call the RHS to compute the derivative wrt rho

    ! species rates come back in terms of molar fractions
    R(UFS:UFS-1+nspec_evolve) = &
         state(URHO) * aion(1:nspec_evolve) * burn_state % ydot(1:nspec_evolve)

    R(UEDEN) = state(URHO) * burn_state % ydot(net_ienuc)
    R(UEINT) = state(URHO) * burn_state % ydot(net_ienuc)

  end subroutine single_zone_react_source

end module react_util_module
