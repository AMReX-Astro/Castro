module sdc_vode_rhs_module

  use amrex_fort_module, only : rt => amrex_real
  use castro_error_module, only : castro_error

  implicit none

contains

  subroutine f_rhs(t, U, dUdt, rpar)
    ! this is the righthand side for the ODE system that we will use
    ! with VODE

    use meth_params_module, only : nvar, URHO, UFS, UEDEN, UTEMP, UMX, UMZ, UEINT
    use burn_type_module
    use vode_rpar_indices
    use react_util_module
    use cuvode_parameters_module

    implicit none

    real(rt), intent(in) :: t
    real(rt), intent(in) :: U(0:VODE_NEQS-1)
    real(rt), intent(out) :: dUdt(0:VODE_NEQS-1)
    real(rt), intent(inout) :: rpar(n_rpar_comps)

    real(rt) :: U_full(nvar),  R_full(nvar)
    real(rt) :: R_react(0:VODE_NEQS-1), C_react(0:VODE_NEQS-1)
    type(burn_t) :: burn_state

    ! evaluate R

    ! we are not solving the momentum equations
    ! create a full state -- we need this for some interfaces
    U_full(URHO) = U(0)
    U_full(UFS:UFS-1+nspec_evolve) = U(1:nspec_evolve)
    U_full(UEINT) = U(nspec_evolve+1)
    U_full(UEDEN) = rpar(irp_evar)

    U_full(UMX:UMZ) = rpar(irp_mom:irp_mom+2)
    U_full(UFS+nspec_evolve:UFS-1+nspec) = rpar(irp_spec:irp_spec-1+(nspec-nspec_evolve))

    ! initialize the temperature -- a better value will be found when we do the EOS
    ! call in single_zone_react_source
    U_full(UTEMP) = rpar(irp_temp)

    call single_zone_react_source(U_full, R_full, 0,0,0, burn_state)

    ! update our temperature for next time
    rpar(irp_temp) = burn_state % T

    R_react(0) = R_full(URHO)
    R_react(1:nspec_evolve) = R_full(UFS:UFS-1+nspec_evolve)
    R_react(nspec_evolve+1) = R_full(UEINT)

    ! C comes in through rpar
    C_react(:) = rpar(irp_f_source:irp_f_source-1+nspec_evolve+2)

    ! create the RHS
    dUdt(:) = R_react(:) + C_react(:)

  end subroutine f_rhs

  subroutine jac(time, U, ml, mu, Jacobian, nrowpd, rpar)
    ! this is the Jacobian function for use with VODE

    use amrex_constants_module, only : ZERO, ONE
    use meth_params_module, only : nvar, URHO, UFS, UEDEN, UTEMP, UMX, UMZ, UEINT
    use burn_type_module
    use eos_type_module, only : eos_t, eos_input_re
    use eos_composition_module, only : eos_comp_t, composition_derivatives
    use eos_module, only : eos
    use vode_rpar_indices
    use react_util_module
    use cuvode_parameters_module

    implicit none

    integer   , intent(IN   ) :: ml, mu, nrowpd
    real(rt), intent(INOUT) :: U(0:VODE_NEQS-1), rpar(n_rpar_comps), time
    real(rt), intent(  OUT) :: Jacobian(0:VODE_NEQS-1, 0:VODE_NEQS-1)

    type(burn_t) :: burn_state
    type(eos_t) :: eos_state
    type(eos_comp_t) :: eos_comp

    real(rt) :: U_full(nvar),  R_full(nvar)

    real(rt) :: denom
    real(rt) :: dRdw(0:nspec_evolve+1, 0:nspec_evolve+1), dwdU(0:nspec_evolve+1, 0:nspec_evolve+1)

    integer :: m

    ! we are not solving the momentum equations
    ! create a full state -- we need this for some interfaces
    U_full(URHO) = U(0)
    U_full(UFS:UFS-1+nspec_evolve) = U(1:nspec_evolve)
    U_full(UEINT) = U(nspec_evolve+1)
    U_full(UEDEN) = rpar(irp_evar)

    U_full(UMX:UMZ) = rpar(irp_mom:irp_mom+2)
    U_full(UFS+nspec_evolve:UFS-1+nspec) = rpar(irp_spec:irp_spec-1+(nspec-nspec_evolve))

    ! compute the temperature and species derivatives --
    ! maybe this should be done using the burn_state
    ! returned by single_zone_react_source, since it is
    ! more consistent T from e
    eos_state % rho = U_full(URHO)
    eos_state % T = rpar(irp_temp)   ! initial guess
    eos_state % xn(:) = U_full(UFS:UFS-1+nspec)/U_full(URHO)
    eos_state % e = U_full(UEINT)/U_full(URHO)  !(U_full(UEDEN) - HALF*sum(U_full(UMX:UMZ))/U_full(URHO))/U_full(URHO)

    call eos(eos_input_re, eos_state)

    U_full(UTEMP) = eos_state % T

    call single_zone_react_source(U_full, R_full, 0,0,0, burn_state)

    call single_zone_jac(U_full, burn_state, dRdw)

    ! construct dwdU
    dwdU(:, :) = ZERO

    ! the density row
    dwdU(0, 0) = ONE

    ! the X_k rows
    do m = 1, nspec_evolve
       dwdU(m,0) = -U(m)/U(0)**2
       dwdU(m,m) = ONE/U(0)
    enddo

    call composition_derivatives(eos_state, eos_comp)

    ! now the T row -- this depends on whether we are evolving (rho E) or (rho e)
    denom = ONE/(eos_state % rho * eos_state % dedT)
    dwdU(nspec_evolve+1,0) = denom*(sum(eos_state % xn(1:nspec_evolve) * eos_comp % dedX(1:nspec_evolve)) - &
                                    eos_state % rho * eos_state % dedr - eos_state % e)

    do m = 1, nspec_evolve
       dwdU(nspec_evolve+1,m) = -denom * eos_comp % dedX(m)
    enddo

    dwdU(nspec_evolve+1, nspec_evolve+1) = denom

    ! construct the Jacobian -- we can get most of the
    ! terms from the network itself, but we do not rely on
    ! it having derivative wrt density
    Jacobian(:,:) = matmul(dRdw, dwdU)

  end subroutine jac

end module sdc_vode_rhs_module
