module sdc_vode_rhs_module

  use amrex_fort_module, only : rt => amrex_real
  use castro_error_module, only : castro_error
  use cuvode_types_module, only : dvode_t

  implicit none

contains

  subroutine f_rhs(t, vode_state, dUdt)
    ! this is the righthand side for the ODE system that we will use
    ! with VODE

    use meth_params_module, only : nvar, URHO, UFS, UEDEN, UTEMP, UMX, UMZ, UEINT
    use burn_type_module
    use vode_rpar_indices
    use react_util_module
    use cuvode_parameters_module

    implicit none

    real(rt), intent(in) :: t
    type(dvode_t), intent(inout) :: vode_state
    real(rt), intent(out) :: dUdt(0:VODE_NEQS-1)

    real(rt) :: U_full(nvar),  R_full(nvar)
    real(rt) :: R_react(0:VODE_NEQS-1), C_react(0:VODE_NEQS-1)
    type(burn_t) :: burn_state

    ! evaluate R

    ! we are not solving the momentum equations
    ! create a full state -- we need this for some interfaces
    ! note: vode_state % y(:) is 1-based
    U_full(URHO) = vode_state % y(1)
    U_full(UFS:UFS-1+nspec) = vode_state % y(2:nspec+1)
    U_full(UEINT) = vode_state % y(nspec+2)
    U_full(UEDEN) = vode_state % rpar(irp_evar)

    U_full(UMX:UMZ) = vode_state % rpar(irp_mom:irp_mom+2)

    ! initialize the temperature -- a better value will be found when we do the EOS
    ! call in single_zone_react_source
    U_full(UTEMP) = vode_state % rpar(irp_temp)

    call single_zone_react_source(U_full, R_full, 0,0,0, burn_state)

    ! update our temperature for next time
    vode_state % rpar(irp_temp) = burn_state % T

    R_react(0) = R_full(URHO)
    R_react(1:nspec) = R_full(UFS:UFS-1+nspec)
    R_react(nspec+1) = R_full(UEINT)

    ! C comes in through rpar
    C_react(:) = vode_state % rpar(irp_f_source:irp_f_source-1+nspec+2)

    ! create the RHS
    dUdt(:) = R_react(:) + C_react(:)

  end subroutine f_rhs

  subroutine jac(time, vode_state, ml, mu, Jacobian, nrowpd)
    ! this is the Jacobian function for use with VODE

    use amrex_constants_module, only : ZERO, ONE
    use meth_params_module, only : nvar, URHO, UFS, UEDEN, UTEMP, UMX, UMZ, UEINT
    use burn_type_module
    use eos_type_module, only : eos_t, eos_input_re
    use eos_composition_module, only : eos_xderivs_t, composition_derivatives
    use eos_module, only : eos
    use vode_rpar_indices
    use react_util_module
    use cuvode_parameters_module

    implicit none

    integer   , intent(IN   ) :: ml, mu, nrowpd
    real(rt), intent(IN) :: time
    type(dvode_t), intent(inout) :: vode_state
    real(rt), intent(  OUT) :: Jacobian(0:VODE_NEQS-1, 0:VODE_NEQS-1)

    type(burn_t) :: burn_state
    type(eos_t) :: eos_state
    type(eos_xderivs_t) :: eos_xderivs

    real(rt) :: U_full(nvar),  R_full(nvar)

    real(rt) :: denom
    real(rt) :: dRdw(0:nspec+1, 0:nspec+1), dwdU(0:nspec+1, 0:nspec+1)

    integer :: m

    ! we are not solving the momentum equations
    ! create a full state -- we need this for some interfaces
    U_full(URHO) = vode_state % y(1)
    U_full(UFS:UFS-1+nspec_evolve) = vode_state % y(2:nspec+1)
    U_full(UEINT) = vode_state % y(nspec+2)
    U_full(UEDEN) = vode_state % rpar(irp_evar)

    U_full(UMX:UMZ) = vode_state % rpar(irp_mom:irp_mom+2)

    ! compute the temperature and species derivatives --
    ! maybe this should be done using the burn_state
    ! returned by single_zone_react_source, since it is
    ! more consistent T from e
    eos_state % rho = U_full(URHO)
    eos_state % T = vode_state % rpar(irp_temp)   ! initial guess
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
    do m = 1, nspec
       dwdU(m,0) = -vode_state % y(m+1)/ vode_state % y(1)**2
       dwdU(m,m) = ONE/vode_state % y(1)
    enddo

    call composition_derivatives(eos_state, eos_xderivs)

    ! now the T row -- this depends on whether we are evolving (rho E) or (rho e)
    denom = ONE/(eos_state % rho * eos_state % dedT)
    dwdU(nspec+1,0) = denom*(sum(eos_state % xn(1:nspec) * eos_xderivs % dedX(1:nspec)) - &
                                    eos_state % rho * eos_state % dedr - eos_state % e)

    do m = 1, nspec
       dwdU(nspec+1,m) = -denom * eos_xderivs % dedX(m)
    enddo

    dwdU(nspec+1, nspec+1) = denom

    ! construct the Jacobian -- we can get most of the
    ! terms from the network itself, but we do not rely on
    ! it having derivative wrt density
    Jacobian(:,:) = matmul(dRdw, dwdU)

  end subroutine jac

end module sdc_vode_rhs_module
