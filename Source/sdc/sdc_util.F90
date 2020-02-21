module sdc_util

  use amrex_fort_module, only : rt => amrex_real
  use castro_error_module, only : castro_error

  implicit none

  ! error codes
  integer, parameter :: NEWTON_SUCCESS = 0
  integer, parameter :: SINGULAR_MATRIX = -1
  integer, parameter :: CONVERGENCE_FAILURE = -2

  ! solvers
  integer, parameter :: NEWTON_SOLVE = 1
  integer, parameter :: VODE_SOLVE = 2
  integer, parameter :: HYBRID_SOLVE = 3

contains

#ifdef REACTIONS
  subroutine sdc_solve(dt_m, U_old, U_new, C, sdc_iteration)
    ! this is the main interface to solving the discretized nonlinear
    ! reaction update.  It either directly calls the Newton method or first
    ! tries VODE and then does the Newton update.

    use meth_params_module, only : NVAR, sdc_solver, UTEMP, URHO, UFS
    use amrex_constants_module, only : ZERO, HALF, ONE
    use burn_type_module, only : burn_t
    use react_util_module
    use network, only : nspec, nspec_evolve

    implicit none

    real(rt), intent(in) :: dt_m
    real(rt), intent(in) :: U_old(NVAR)
    real(rt), intent(inout) :: U_new(NVAR)
    real(rt), intent(in) :: C(NVAR)
    integer, intent(in) :: sdc_iteration

    integer :: ierr

    real(rt) :: err_out

    ! for debugging
    real(rt) :: U_orig(NVAR)

    U_orig(:) = U_old(:)

    if (sdc_solver == NEWTON_SOLVE) then
       ! we are going to assume we already have a good guess for the
       ! solving in U_new and just pass the solve onto the main Newton
       ! solve
       call sdc_newton_subdivide(dt_m, U_old, U_new, C, sdc_iteration, err_out, ierr)

       ! failing?
       if (ierr /= NEWTON_SUCCESS) then
          print *, "Newton convergence failure"
          print *, "convergence failure, error = ", err_out
          print *, "density: ", U_old(URHO)
          print *, "(old) temperature: ", U_old(UTEMP)
          print *, "mass fractions: ", U_old(UFS:UFS-1+nspec)/U_old(URHO)
          call castro_error("Newton subcycling failed in sdc_solve")
       end if

    else if (sdc_solver == VODE_SOLVE) then
       ! use VODE to do the solution
       call sdc_vode_solve(dt_m, U_old, U_new, C, sdc_iteration)

    else if (sdc_solver == HYBRID_SOLVE) then
       ! if it is the first iteration, we will use VODE to predict
       ! the solution.  Otherwise, we will use Newton.
       if (sdc_iteration == 0) then
          call sdc_vode_solve(dt_m, U_old, U_new, C, sdc_iteration)
       endif

       ! now U_new is the update that VODE predicts, so we will use
       ! that as the initial guess to the Newton solve
       call sdc_newton_subdivide(dt_m, U_old, U_new, C, sdc_iteration, err_out, ierr)

       ! failing?
       if (ierr /= NEWTON_SUCCESS) then
          call castro_error("Newton failure in sdc_solve")
       end if
    end if

  end subroutine sdc_solve

  subroutine sdc_newton_subdivide(dt_m, U_old, U_new, C, sdc_iteration, err_out, ierr)
    ! This is the driver for solving the nonlinear update for the
    ! reacting/advecting system using Newton's method.  It attempts to
    ! do the solution for the full dt_m requested, but if it fails,
    ! will subdivide the domain until it converges or reaches our
    ! limit on the number of subintervals.
    use meth_params_module, only : NVAR, URHO, UFS
    use amrex_constants_module, only : ZERO, HALF, ONE
    use network, only : nspec, nspec_evolve
    use extern_probin_module, only : small_x

    implicit none

    real(rt), intent(in) :: dt_m
    real(rt), intent(in) :: U_old(NVAR)
    real(rt), intent(inout) :: U_new(NVAR)
    real(rt), intent(in) :: C(NVAR)
    integer, intent(in) :: sdc_iteration
    real(rt), intent(out) :: err_out
    integer, intent(inout) :: ierr

    integer :: isub, nsub
    real(rt) :: dt_sub
    real(rt) :: U_begin(NVAR)
    integer :: n
    real(rt) :: sum_rhoX

    integer, parameter :: MAX_NSUB = 64

    ! subdivide the timestep and do multiple Newtons.  We come in here
    ! with an initial guess for the new solution stored in U_new.
    ! That only really makes sense for the case where we have 1
    ! substep.  Otherwise, we should just use the old time solution.
    nsub = 1
    ierr = CONVERGENCE_FAILURE
    U_begin(:) = U_old(:)
    do while (nsub < MAX_NSUB .and. ierr /= NEWTON_SUCCESS)
       if (nsub > 1) then
          U_new(:) = U_old(:)
       end if
       dt_sub = dt_m / nsub
       do isub = 1, nsub

          ! normalize species
          do n = 1, nspec
             U_begin(UFS-1+n) = max(small_x, U_begin(UFS-1+n))
          end do

          sum_rhoX = sum(U_begin(UFS:UFS-1+nspec))
          U_begin(UFS:UFS-1+nspec) = U_begin(UFS:UFS-1+nspec) * U_begin(URHO)/sum_rhoX

          call sdc_newton_solve(dt_sub, U_begin, U_new, C, sdc_iteration, err_out, ierr)
          U_begin(:) = U_new(:)

       end do
       nsub = nsub * 2
    end do

  end subroutine sdc_newton_subdivide

  subroutine sdc_newton_solve(dt_m, U_old, U_new, C, sdc_iteration, err_out, ierr)
    ! the purpose of this function is to solve the system
    ! U - dt R(U) = U_old + dt C using a Newton solve.
    !
    ! here, U_new should come in as a guess for the new U and will be
    ! returned with the value that satisfies the nonlinear function

    use meth_params_module, only : NVAR, UEDEN, UEINT, URHO, UFS, UMX, UMZ, UTEMP, &
                                   sdc_order, &
                                   sdc_solver_tol_dens, sdc_solver_tol_spec, sdc_solver_tol_ener, &
                                   sdc_solver_atol, &
                                   sdc_solver_relax_factor, &
                                   sdc_solve_for_rhoe
    use amrex_constants_module, only : ZERO, HALF, ONE
    use burn_type_module, only : burn_t
    use react_util_module
    use network, only : nspec, nspec_evolve
    use vode_rpar_indices
    use extern_probin_module, only : small_x
#if INTEGRATOR == 0
    use linpack_module
#endif
    implicit none

    real(rt), intent(in) :: dt_m
    real(rt), intent(in) :: U_old(NVAR)
    real(rt), intent(inout) :: U_new(NVAR)
    real(rt), intent(in) :: C(NVAR)
    integer, intent(in) :: sdc_iteration
    real(rt), intent(out) :: err_out
    integer, intent(out) :: ierr

    real(rt) :: Jac(0:nspec_evolve+1, 0:nspec_evolve+1)
    real(rt) :: w(0:nspec_evolve+1)

    real(rt) :: rpar(n_rpar_comps)

    integer :: ipvt(nspec_evolve+2)
    integer :: info

    logical :: converged

    real(rt) :: tol_dens, tol_spec, tol_ener, relax_fac
    real(rt) :: eps_tot(0:nspec_evolve+1)

    ! we will do the implicit update of only the terms that have reactive sources
    !
    !   0               : rho
    !   1:nspec_evolve  : species
    !   nspec_evolve+1  : (rho E) or (rho e)

    real(rt) :: U_react(0:nspec_evolve+1), f_source(0:nspec_evolve+1)
    real(rt) :: dU_react(0:nspec_evolve+1), f(0:nspec_evolve+1), f_rhs(0:nspec_evolve+1)

    real(rt) :: err, eta

    integer, parameter :: MAX_ITER = 100
    integer :: iter

    integer :: max_newton_iter

    real(rt) :: xn(nspec)
    integer :: k

    ierr = NEWTON_SUCCESS

    ! the tolerance we are solving to may depend on the iteration
    relax_fac = sdc_solver_relax_factor**(sdc_order - sdc_iteration - 1)
    tol_dens = sdc_solver_tol_dens * relax_fac
    tol_spec = sdc_solver_tol_spec * relax_fac
    tol_ener = sdc_solver_tol_ener * relax_fac

    ! update the momenta for this zone -- they don't react
    U_new(UMX:UMZ) = U_old(UMX:UMZ) + dt_m * C(UMX:UMZ)

    ! update the non-reacting species
    U_new(UFS+nspec_evolve:UFS-1+nspec) = U_old(UFS+nspec_evolve:UFS-1+nspec) + &
         dt_m * C(UFS+nspec_evolve:UFS-1+nspec)

    ! now only save the subset that participates in the nonlinear
    ! solve -- note: we include the old state in f_source

    ! load rpar

    ! for the Jacobian solve, we are solving
    !   f(U) = U - dt R(U) - U_old - dt C = 0
    ! we define f_source = U_old + dt C so we are solving
    !   f(U) = U - dt R(U) - f_source = 0

    f_source(0) = U_old(URHO) + dt_m * C(URHO)
    f_source(1:nspec_evolve) = U_old(UFS:UFS-1+nspec_evolve) + dt_m * C(UFS:UFS-1+nspec_evolve)
    if (sdc_solve_for_rhoe == 1) then
       f_source(nspec_evolve+1) = U_old(UEINT) + dt_m * C(UEINT)
    else
       f_source(nspec_evolve+1) = U_old(UEDEN) + dt_m * C(UEDEN)
    endif

    rpar(irp_f_source:irp_f_source-1+nspec_evolve+2) = f_source(:)
    rpar(irp_dt) = dt_m
    rpar(irp_mom:irp_mom-1+3) = U_new(UMX:UMZ)

    ! temperature will be used as an initial guess in the EOS
    rpar(irp_temp) = U_old(UTEMP)

    ! we should be able to do an update for this somehow?
    if (sdc_solve_for_rhoe == 1) then
       rpar(irp_evar) = U_new(UEDEN)
    else
       rpar(irp_evar) = U_new(UEINT)
    endif

    rpar(irp_spec:irp_spec-1+(nspec-nspec_evolve)) = &
         U_new(UFS+nspec_evolve:UFS-1+nspec)

    ! store the subset for the nonlinear solve
    ! We use an initial guess if possible
    U_react(0) = U_new(URHO)
    U_react(1:nspec_evolve) = U_new(UFS:UFS-1+nspec_evolve)
    if (sdc_solve_for_rhoe == 1) then
       U_react(nspec_evolve+1) = U_new(UEINT)
    else
       U_react(nspec_evolve+1) = U_new(UEDEN)
    endif

#if (INTEGRATOR == 0)

    ! do a simple Newton solve

    ! iterative loop
    iter = 0
    max_newton_iter = MAX_ITER

    err = 1.e30_rt
    converged = .false.
    do while (.not. converged .and. iter < max_newton_iter)

       call f_sdc_jac(nspec_evolve+2, U_react, f, Jac, nspec_evolve+2, info, rpar)

       ! solve the linear system: Jac dU_react = -f
       call dgefa(Jac, ipvt, info)
       if (info /= 0) then
          ierr = SINGULAR_MATRIX
          return
       endif

       f_rhs(:) = -f(:)

       call dgesl(Jac, ipvt, f_rhs)

       dU_react(:) = f_rhs(:)

       ! how much of dU_react should we apply?
       eta = ONE
       dU_react(:) = eta * dU_react(:)

       U_react(:) = U_react(:) + dU_react(:)

       ! we still need to normalize here
       ! xn(1:nspec_evolve) = U_react(1:nspec_evolve)/U_react(0)
       ! xn(nspec_evolve+1:nspec) = rpar(irp_spec:irp_spec-1+(nspec-nspec_evolve))/U_react(0)
 
       ! do k = 1, nspec
       !    xn(k) = max(small_x, xn(k))
       ! end do
       ! xn(:) = xn(:)/sum(xn)

       ! U_react(1:nspec_evolve) = U_react(0) * xn(1:nspec_evolve)
       ! rpar(irp_spec:irp_spec-1+(nspec-nspec_evolve)) = U_react(0) * xn(nspec_evolve+1:nspec)

       eps_tot(0) = tol_dens * abs(U_react(0)) + sdc_solver_atol
       ! for species, atol is the mass fraction limit, so we multiply by density to get a partial density limit
       eps_tot(1:nspec_evolve) = tol_spec * abs(U_react(1:nspec_evolve)) + sdc_solver_atol * abs(U_react(0))
       eps_tot(nspec_evolve+1) = tol_ener * abs(U_react(nspec_evolve+1)) + sdc_solver_atol

       ! compute the norm of the weighted error, where the weights are 1/eps_tot
       err = sqrt(sum((dU_react/eps_tot)**2)/(nspec_evolve+2))

       if (err < ONE) then
          converged = .true.
       endif

       iter = iter + 1
    enddo

    err_out = err

    if (.not. converged) then
       !print *, "dens: ", U_react(0), dU_react(0), eps_tot(0), abs(dU_react(0))/eps_tot(0)
       !do n = 1, nspec_evolve
       !   print *, "spec: ", n, U_react(n), dU_react(n), eps_tot(n), abs(dU_react(n))/eps_tot(n)
       !end do
       !print *, "enuc: ", U_react(nspec_evolve+1), dU_react(nspec_evolve+1), eps_tot(nspec_evolve+1), dU_react(nspec_evolve+1)/eps_tot(nspec_evolve+1)
       ierr = CONVERGENCE_FAILURE
       return
    endif

#endif

    ! update the full U_new
    ! if we updated total energy, then correct internal, or vice versa
    U_new(URHO) = U_react(0)
    U_new(UFS:UFS-1+nspec_evolve) = U_react(1:nspec_evolve)
    if (sdc_solve_for_rhoe == 1) then
       U_new(UEINT) = U_react(nspec_evolve+1)
       U_new(UEDEN) = U_new(UEINT) + HALF*sum(U_new(UMX:UMZ)**2)/U_new(URHO)
    else
       U_new(UEDEN) = U_react(nspec_evolve+1)
       U_new(UEINT) = U_new(UEDEN) - HALF*sum(U_new(UMX:UMZ)**2)/U_new(URHO)
    endif

  end subroutine sdc_newton_solve


  subroutine sdc_vode_solve(dt_m, U_old, U_new, C, sdc_iteration)
    ! the purpose of this function is to solve the system the
    ! approximate system dU/dt = R + C using the VODE ODE integrator.
    ! the solution we get here will then be used as the initial guess
    ! to the Newton solve on the real system.

    use meth_params_module, only : NVAR, UEDEN, UEINT, URHO, UFS, UMX, UMZ, UTEMP, &
                                   sdc_order, &
                                   sdc_solver_tol_dens, sdc_solver_tol_spec, sdc_solver_tol_ener, &
                                   sdc_solver_atol, &
                                   sdc_solver_relax_factor, &
                                   sdc_solve_for_rhoe, sdc_use_analytic_jac
    use amrex_constants_module, only : ZERO, HALF, ONE
    use burn_type_module, only : burn_t
    use react_util_module
    use network, only : nspec, nspec_evolve
    use vode_rpar_indices
    use cuvode_parameters_module
    use cuvode_types_module, only : dvode_t, rwork_t
    use extern_probin_module, only : use_jacobian_caching
    use cuvode_module, only: dvode

    implicit none

    real(rt), intent(in) :: dt_m
    real(rt), intent(in) :: U_old(NVAR)
    real(rt), intent(inout) :: U_new(NVAR)
    real(rt), intent(in) :: C(NVAR)
    integer, intent(in) :: sdc_iteration

    type(dvode_t) :: dvode_state
    type(rwork_t) :: rwork
    integer :: iwork(VODE_LIW)

    integer :: imode

    real(rt) :: time
    real(rt) :: tol_dens, tol_spec, tol_ener, relax_fac
    real(rt) :: rtol(0:nspec_evolve+1), atol(0:nspec_evolve+1)

    ! we will do the implicit update of only the terms that have reactive sources
    !
    !   0               : rho
    !   1:nspec_evolve  : species
    !   nspec_evolve+1  : (rho E) or (rho e)

    real(rt) :: C_react(0:nspec_evolve+1)

#if (INTEGRATOR == 0)

    ! the tolerance we are solving to may depend on the iteration
    relax_fac = sdc_solver_relax_factor**(sdc_order - sdc_iteration - 1)
    tol_dens = sdc_solver_tol_dens * relax_fac
    tol_spec = sdc_solver_tol_spec * relax_fac
    tol_ener = sdc_solver_tol_ener * relax_fac

    ! update the momenta for this zone -- they don't react
    U_new(UMX:UMZ) = U_old(UMX:UMZ) + dt_m * C(UMX:UMZ)

    ! update the non-reacting species
    U_new(UFS+nspec_evolve:UFS-1+nspec) = U_old(UFS+nspec_evolve:UFS-1+nspec) + &
         dt_m * C(UFS+nspec_evolve:UFS-1+nspec)

    ! now only save the subset that participates in the nonlinear
    ! solve -- note: we include the old state in f_source

    ! load rpar

    ! if we are solving the system as an ODE, then we
    ! are solving
    !    dU/dt = R(U) + C
    ! so we simply pass in C
    C_react(0) = C(URHO)
    C_react(1:nspec_evolve) = C(UFS:UFS-1+nspec_evolve)
    C_react(nspec_evolve+1) = C(UEINT)

    dvode_state % rpar(irp_f_source:irp_f_source-1+nspec_evolve+2) = C_react(:)
    dvode_state % rpar(irp_dt) = dt_m
    dvode_state % rpar(irp_mom:irp_mom-1+3) = U_new(UMX:UMZ)

    ! temperature will be used as an initial guess in the EOS
    dvode_state % rpar(irp_temp) = U_old(UTEMP)

    ! we are always solving for rhoe with the VODE predict
    dvode_state % rpar(irp_evar) = U_new(UEDEN)

    dvode_state % rpar(irp_spec:irp_spec-1+(nspec-nspec_evolve)) = &
         U_new(UFS+nspec_evolve:UFS-1+nspec)


    ! store the subset for the nonlinear solve.  We only consider (rho
    ! e), not (rho E).  This is because at present we do not have a
    ! method of updating the velocities during the multistep
    ! integration

    ! Also note that the dvode_state is 1-based, but we'll access it
    ! as 0-based in our implementation of the RHS routine

    dvode_state % y(1) = U_old(URHO)
    dvode_state % y(2:nspec_evolve+1) = U_old(UFS:UFS-1+nspec_evolve)
    dvode_state % y(nspec_evolve+2) = U_old(UEINT)

    dvode_state % istate = 1

    iwork(:) = 0

    ! set the maximum number of steps allowed -- the VODE default is 500
    iwork(6) = 25000

    rwork % CONDOPT = ZERO
    rwork % YH = ZERO
    rwork % WM = ZERO
    rwork % EWT = ZERO
    rwork % SAVF = ZERO
    rwork % ACOR = ZERO

    dvode_state % T = ZERO
    dvode_state % TOUT = dt_m

    if (sdc_use_analytic_jac == 1) then
       imode = MF_ANALYTIC_JAC_CACHED
    else
       imode = MF_NUMERICAL_JAC_CACHED
    endif

    if (.not. use_jacobian_caching) then
       imode = -imode
    endif

    ! relative tolerances
    rtol(0) = tol_dens
    rtol(1:nspec_evolve) = tol_spec
    rtol(nspec_evolve+1) = tol_ener

    ! absolute tolerances
    atol(0) = sdc_solver_atol * U_old(URHO)
    atol(1:nspec_evolve) = sdc_solver_atol * U_old(URHO)   ! this way, atol is the minimum x
    if (sdc_solve_for_rhoe == 1) then
       atol(nspec_evolve+1) = sdc_solver_atol * U_old(UEINT)
    else
       atol(nspec_evolve+1) = sdc_solver_atol * U_old(UEDEN)
    endif

    dvode_state % atol(:) = atol(:)
    dvode_state % rtol(:) = rtol(:)

    call dvode(dvode_state, rwork, iwork, ITASK, IOPT, imode)

    if (dvode_state % istate < 0) then
       print *, "VODE error, istate = ", dvode_state % istate
       call castro_error("vode termination poorly")
    endif

    ! update the full U_new
    U_new(URHO) = dvode_state % y(1)
    U_new(UFS:UFS-1+nspec_evolve) = dvode_state % y(2:nspec_evolve+1)
    U_new(UEINT) = dvode_state % y(nspec_evolve+2)
    U_new(UEDEN) = U_new(UEINT) + HALF*sum(U_new(UMX:UMZ)**2)/U_new(URHO)

    ! keep our temperature guess
    U_new(UTEMP) = dvode_state % rpar(irp_temp)

#endif

  end subroutine sdc_vode_solve

  subroutine f_sdc_jac(neq, U, f, Jac, ldjac, iflag, rpar)
    ! this is used with the Newton solve and returns f and the Jacobian

    use vode_rpar_indices
    use meth_params_module, only : nvar, URHO, UFS, UEINT, UEDEN, UMX, UMZ, UTEMP, &
         sdc_solve_for_rhoe
    use network, only : nspec, nspec_evolve
    use burn_type_module
    use react_util_module
    use eos_type_module, only : eos_t, eos_input_re
    use eos_module, only : eos
    use eos_composition_module, only : eos_xderivs_t, composition_derivatives
    use amrex_constants_module, only : ZERO, HALF, ONE
    use vode_rpar_indices
    use extern_probin_module, only : small_x

    ! this computes the function we need to zero for the SDC update
    implicit none

    integer,intent(in) :: neq, ldjac
    real(rt), intent(in)  :: U(0:neq-1)
    real(rt), intent(out) :: f(0:neq-1)
    real(rt), intent(out) :: Jac(0:ldjac-1,0:neq-1)
    integer, intent(inout) :: iflag  !! leave this untouched
    real(rt), intent(inout) :: rpar(n_rpar_comps)

    real(rt) :: U_full(nvar),  R_full(nvar)
    real(rt) :: R_react(0:neq-1), f_source(0:neq-1)
    type(burn_t) :: burn_state
    type(eos_t) :: eos_state
    type(eos_xderivs_t) :: eos_xderivs
    real(rt) :: dt_m

    real(rt) :: denom
    real(rt) :: dRdw(0:nspec_evolve+1, 0:nspec_evolve+1), dwdU(0:nspec_evolve+1, 0:nspec_evolve+1)
    integer :: m, k
    real(rt) :: sum_rhoX

    ! we are not solving the momentum equations
    ! create a full state -- we need this for some interfaces
    U_full(URHO) = U(0)
    U_full(UFS:UFS-1+nspec_evolve) = U(1:nspec_evolve)
    if (sdc_solve_for_rhoe == 1) then
       U_full(UEINT) = U(nspec_evolve+1)
       U_full(UEDEN) = rpar(irp_evar)
    else
       U_full(UEDEN) = U(nspec_evolve+1)
       U_full(UEINT) = rpar(irp_evar)
    endif

    U_full(UMX:UMZ) = rpar(irp_mom:irp_mom+2)
    U_full(UFS+nspec_evolve:UFS-1+nspec) = rpar(irp_spec:irp_spec-1+(nspec-nspec_evolve))

    ! normalize the species
    do k = 1, nspec
       U_full(UFS-1+k) = max(small_x, U_full(UFS-1+k))
    end do

    sum_rhoX = sum(U_full(UFS:UFS-1+nspec))
    U_full(UFS:UFS-1+nspec) = U_full(UFS:UFS-1+nspec) * U_full(URHO)/sum_rhoX

    ! unpack rpar
    dt_m = rpar(irp_dt)
    f_source(:) = rpar(irp_f_source:irp_f_source-1+nspec_evolve+2)

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

    ! store the subset of R used in the Jacobian
    R_react(0) = R_full(URHO)
    R_react(1:nspec_evolve) = R_full(UFS:UFS-1+nspec_evolve)
    if (sdc_solve_for_rhoe == 1) then
       R_react(nspec_evolve+1) = R_full(UEINT)
    else
       R_react(nspec_evolve+1) = R_full(UEDEN)
    endif

    f(:) = U(:) - dt_m * R_react(:) - f_source(:)

    ! get dRdw -- this may do a numerical approxiation or use the
    ! network's analytic Jac
    call single_zone_jac(U_full, burn_state, dRdw)

    ! construct dwdU
    dwdU(:, :) = ZERO

    ! the density row
    dwdU(iwrho, 0) = ONE

    ! the X_k rows
    do m = 1, nspec_evolve
       dwdU(iwfs-1+m,0) = -U(m)/U(0)**2
       dwdU(iwfs-1+m,m) = ONE/U(0)
    enddo

    call composition_derivatives(eos_state, eos_xderivs)

    ! now the T row -- this depends on whether we are evolving (rho E) or (rho e)
    denom = ONE/(eos_state % rho * eos_state % dedT)
    if (sdc_solve_for_rhoe == 1) then
       dwdU(iwT,0) = denom*(sum(eos_state % xn(1:nspec_evolve) * eos_xderivs % dedX(1:nspec_evolve)) - &
                                       eos_state % rho * eos_state % dedr - eos_state % e)
    else
       dwdU(iwT,0) = denom*(sum(eos_state % xn(1:nspec_evolve) * eos_xderivs % dedX(1:nspec_evolve)) - &
                                       eos_state % rho * eos_state % dedr - eos_state % e - &
                                       HALF*sum(U_full(UMX:UMZ)**2)/eos_state % rho**2)
    endif

    do m = 1, nspec_evolve
       dwdU(iwT,m) = -denom * eos_xderivs % dedX(m)
    enddo

    dwdU(iwT, nspec_evolve+1) = denom

    ! construct the Jacobian -- we can get most of the
    ! terms from the network itself, but we do not rely on
    ! it having derivative wrt density
    Jac(:, :) = ZERO
    do m = 0, nspec_evolve+1
       Jac(m, m) = ONE
    enddo

    Jac(:,:) = Jac(:,:) - dt_m * matmul(dRdw, dwdU)

  end subroutine f_sdc_jac
#endif

#ifdef REACTIONS
  subroutine ca_sdc_compute_C4_lobatto(lo, hi, &
                                       dt_m, dt, &
                                       A_m, Amlo, Amhi, &
                                       A_0_old, A0lo, A0hi, &
                                       A_1_old, A1lo, A1hi, &
                                       A_2_old, A2lo, A2hi, &
                                       R_0_old, R0lo, R0hi, &
                                       R_1_old, R1lo, R1hi, &
                                       R_2_old, R2lo, R2hi, &
                                       C, Clo, Chi, &
                                       m_start) bind(C, name="ca_sdc_compute_C4_lobatto")

    ! compute the 'C' term for the 4th-order solve with reactions
    ! note: this 'C' is cell-averages

    use meth_params_module, only : NVAR
    use amrex_constants_module, only : ONE, HALF, TWO, FIVE, EIGHT

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: Amlo(3), Amhi(3)
    integer, intent(in) :: A0lo(3), A0hi(3)
    integer, intent(in) :: A1lo(3), A1hi(3)
    integer, intent(in) :: A2lo(3), A2hi(3)
    integer, intent(in) :: R0lo(3), R0hi(3)
    integer, intent(in) :: R1lo(3), R1hi(3)
    integer, intent(in) :: R2lo(3), R2hi(3)
    integer, intent(in) :: Clo(3), Chi(3)
    integer, intent(in) :: m_start

    real(rt), intent(in) :: dt_m, dt
    real(rt), intent(in) :: A_m(Amlo(1):Amhi(1), Amlo(2):Amhi(2), Amlo(3):Amhi(3), NVAR)
    real(rt), intent(in) :: A_0_old(A0lo(1):A0hi(1), A0lo(2):A0hi(2), A0lo(3):A0hi(3), NVAR)
    real(rt), intent(in) :: A_1_old(A1lo(1):A1hi(1), A1lo(2):A1hi(2), A1lo(3):A1hi(3), NVAR)
    real(rt), intent(in) :: A_2_old(A2lo(1):A2hi(1), A2lo(2):A2hi(2), A2lo(3):A2hi(3), NVAR)
    real(rt), intent(in) :: R_0_old(R0lo(1):R0hi(1), R0lo(2):R0hi(2), R0lo(3):R0hi(3), NVAR)
    real(rt), intent(in) :: R_1_old(R1lo(1):R1hi(1), R1lo(2):R1hi(2), R1lo(3):R1hi(3), NVAR)
    real(rt), intent(in) :: R_2_old(R2lo(1):R2hi(1), R2lo(2):R2hi(2), R2lo(3):R2hi(3), NVAR)
    real(rt), intent(out) :: C(Clo(1):Chi(1), Clo(2):Chi(2), Clo(3):Chi(3), NVAR)

    integer :: i, j, k
    real(rt) :: integral(NVAR)

    ! Gauss-Lobatto (Simpsons)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (m_start == 0) then
                ! compute the integral from [t_m, t_{m+1}], normalized by dt_m
                integral(:) = ONE/12.0_rt * (FIVE*(A_0_old(i,j,k,:) + R_0_old(i,j,k,:)) + &
                                             EIGHT*(A_1_old(i,j,k,:) + R_1_old(i,j,k,:)) - &
                                             (A_2_old(i,j,k,:) + R_2_old(i,j,k,:)))

                C(i,j,k,:) = (A_m(i,j,k,:) - A_0_old(i,j,k,:)) - R_1_old(i,j,k,:) + integral

             else if (m_start == 1) then
                ! compute the integral from [t_m, t_{m+1}], normalized by dt_m
                integral(:) = ONE/12.0_rt * (-(A_0_old(i,j,k,:) + R_0_old(i,j,k,:)) + &
                                             EIGHT*(A_1_old(i,j,k,:) + R_1_old(i,j,k,:)) + &
                                             FIVE*(A_2_old(i,j,k,:) + R_2_old(i,j,k,:)))

                C(i,j,k,:) = (A_m(i,j,k,:) - A_1_old(i,j,k,:)) - R_2_old(i,j,k,:) + integral

             else
                call castro_error("error in ca_sdc_compute_C4 -- should not be here")
             endif

          enddo
       enddo
    enddo

  end subroutine ca_sdc_compute_C4_lobatto


  subroutine ca_sdc_compute_C4_radau(lo, hi, &
                                     dt_m, dt, &
                                     A_m, Amlo, Amhi, &
                                     A_0_old, A0lo, A0hi, &
                                     A_1_old, A1lo, A1hi, &
                                     A_2_old, A2lo, A2hi, &
                                     A_3_old, A3lo, A3hi, &
                                     R_0_old, R0lo, R0hi, &
                                     R_1_old, R1lo, R1hi, &
                                     R_2_old, R2lo, R2hi, &
                                     R_3_old, R3lo, R3hi, &
                                     C, Clo, Chi, &
                                     m_start) bind(C, name="ca_sdc_compute_C4_radau")
    ! compute the 'C' term for the 4th-order solve with reactions
    ! note: this 'C' is cell-averages

    use meth_params_module, only : NVAR
    use amrex_constants_module, only : ONE, HALF, TWO, FIVE, EIGHT

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: Amlo(3), Amhi(3)
    integer, intent(in) :: A0lo(3), A0hi(3)
    integer, intent(in) :: A1lo(3), A1hi(3)
    integer, intent(in) :: A2lo(3), A2hi(3)
    integer, intent(in) :: A3lo(3), A3hi(3)
    integer, intent(in) :: R0lo(3), R0hi(3)
    integer, intent(in) :: R1lo(3), R1hi(3)
    integer, intent(in) :: R2lo(3), R2hi(3)
    integer, intent(in) :: R3lo(3), R3hi(3)
    integer, intent(in) :: Clo(3), Chi(3)
    integer, intent(in) :: m_start

    real(rt), intent(in) :: dt_m, dt
    real(rt), intent(in) :: A_m(Amlo(1):Amhi(1), Amlo(2):Amhi(2), Amlo(3):Amhi(3), NVAR)
    real(rt), intent(in) :: A_0_old(A0lo(1):A0hi(1), A0lo(2):A0hi(2), A0lo(3):A0hi(3), NVAR)
    real(rt), intent(in) :: A_1_old(A1lo(1):A1hi(1), A1lo(2):A1hi(2), A1lo(3):A1hi(3), NVAR)
    real(rt), intent(in) :: A_2_old(A2lo(1):A2hi(1), A2lo(2):A2hi(2), A2lo(3):A2hi(3), NVAR)
    real(rt), intent(in) :: A_3_old(A3lo(1):A3hi(1), A3lo(2):A3hi(2), A3lo(3):A3hi(3), NVAR)
    real(rt), intent(in) :: R_0_old(R0lo(1):R0hi(1), R0lo(2):R0hi(2), R0lo(3):R0hi(3), NVAR)
    real(rt), intent(in) :: R_1_old(R1lo(1):R1hi(1), R1lo(2):R1hi(2), R1lo(3):R1hi(3), NVAR)
    real(rt), intent(in) :: R_2_old(R2lo(1):R2hi(1), R2lo(2):R2hi(2), R2lo(3):R2hi(3), NVAR)
    real(rt), intent(in) :: R_3_old(R3lo(1):R3hi(1), R3lo(2):R3hi(2), R3lo(3):R3hi(3), NVAR)
    real(rt), intent(out) :: C(Clo(1):Chi(1), Clo(2):Chi(2), Clo(3):Chi(3), NVAR)

    integer :: i, j, k
    real(rt) :: integral(NVAR)

    ! Gauss-Lobatto (Simpsons)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (m_start == 0) then
                ! compute the integral from [t_m, t_{m+1}], normalized by dt_m
                integral(:) = (dt/dt_m) * (ONE/1800.0_rt) * &
                     ((-35.0_rt*sqrt(6.0_rt) + 440.0_rt)*(A_1_old(i,j,k,:) + R_1_old(i,j,k,:)) + &
                      (-169.0_rt*sqrt(6.0_rt) + 296.0_rt)*(A_2_old(i,j,k,:) + R_2_old(i,j,k,:)) + &
                      (-16.0_rt + 24.0_rt*sqrt(6.0_rt))*(A_3_old(i,j,k,:) + R_3_old(i,j,k,:)))

                C(i,j,k,:) = (A_m(i,j,k,:) - A_0_old(i,j,k,:)) - R_1_old(i,j,k,:) + integral

             else if (m_start == 1) then
                ! compute the integral from [t_m, t_{m+1}], normalized by dt_m
                integral(:) = (dt/dt_m) * (ONE/150.0_rt) * &
                     ((-12.0_rt + 17.0_rt*sqrt(6.0_rt))*(A_1_old(i,j,k,:) + R_1_old(i,j,k,:)) + &
                      (12.0_rt + 17.0_rt*sqrt(6.0_rt))*(A_2_old(i,j,k,:) + R_2_old(i,j,k,:)) + &
                      (-4.0_rt*sqrt(6.0_rt))*(A_3_old(i,j,k,:) + R_3_old(i,j,k,:)))

                C(i,j,k,:) = (A_m(i,j,k,:) - A_1_old(i,j,k,:)) - R_2_old(i,j,k,:) + integral


             else if (m_start == 2) then
                ! compute the integral from [t_m, t_{m+1}], normalized by dt_m
                integral(:) = (dt/dt_m) * (ONE/600.0_rt) * &
                     ((168.0_rt - 73.0_rt*sqrt(6.0_rt))*(A_1_old(i,j,k,:) + R_1_old(i,j,k,:)) + &
                     (120.0_rt + 5.0_rt*sqrt(6.0_rt))*(A_2_old(i,j,k,:) + R_2_old(i,j,k,:)) + &
                     (72.0_rt + 8.0_rt*sqrt(6.0_rt))*(A_3_old(i,j,k,:) + R_3_old(i,j,k,:)))

                C(i,j,k,:) = (A_m(i,j,k,:) - A_2_old(i,j,k,:)) - R_3_old(i,j,k,:) + integral

             else
                call castro_error("error in ca_sdc_compute_C4 -- should not be here")
             endif

          end do
       end do
    end do

  end subroutine ca_sdc_compute_C4_radau


  subroutine ca_sdc_compute_C2_lobatto(lo, hi, dt_m, dt, &
                                       A_m, Amlo, Amhi, &
                                       A_0_old, A0lo, A0hi, &
                                       A_1_old, A1lo, A1hi, &
                                       R_0_old, R0lo, R0hi, &
                                       R_1_old, R1lo, R1hi, &
                                       C, Clo, Chi, &
                                       m_start) bind(C, name="ca_sdc_compute_C2_lobatto")
    ! compute the source term C for the 2nd order Lobatto update

    ! Here, dt_m is the timestep between time-nodes m and m+1

    use meth_params_module, only : NVAR
    use amrex_constants_module, only : ZERO, HALF
    use burn_type_module, only : burn_t
    use network, only : nspec, nspec_evolve
    use react_util_module

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dt_m, dt
    integer, intent(in) :: Amlo(3), Amhi(3)
    integer, intent(in) :: A0lo(3), A0hi(3)
    integer, intent(in) :: A1lo(3), A1hi(3)
    integer, intent(in) :: R0lo(3), R0hi(3)
    integer, intent(in) :: R1lo(3), R1hi(3)
    integer, intent(in) :: Clo(3), Chi(3)
    integer, intent(in) :: m_start

    real(rt), intent(in) :: A_m(Amlo(1):Amhi(1), Amlo(2):Amhi(2), Amlo(3):Amhi(3), NVAR)
    real(rt), intent(in) :: A_0_old(A0lo(1):A0hi(1), A0lo(2):A0hi(2), A0lo(3):A0hi(3), NVAR)
    real(rt), intent(in) :: A_1_old(A1lo(1):A1hi(1), A1lo(2):A1hi(2), A1lo(3):A1hi(3), NVAR)

    real(rt), intent(in) :: R_0_old(R0lo(1):R0hi(1), R0lo(2):R0hi(2), R0lo(3):R0hi(3), NVAR)
    real(rt), intent(in) :: R_1_old(R1lo(1):R1hi(1), R1lo(2):R1hi(2), R1lo(3):R1hi(3), NVAR)

    real(rt), intent(out) :: C(Clo(1):Chi(1), Clo(2):Chi(2), Clo(3):Chi(3), NVAR)

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! construct the source term to the update for 2nd order
             ! Lobatto, there is no advective correction, and we have
             ! C = - R(U^{m+1,k}) + I_m^{m+1}/dt
             C(i,j,k,:) = -R_1_old(i,j,k,:) + &
                  HALF * (A_0_old(i,j,k,:) + A_1_old(i,j,k,:)) + &
                  HALF * (R_0_old(i,j,k,:) + R_1_old(i,j,k,:))
          end do
       end do
    end do

  end subroutine ca_sdc_compute_C2_lobatto


  subroutine ca_sdc_compute_C2_radau(lo, hi, dt_m, dt, &
                                     A_m, Amlo, Amhi, &
                                     A_0_old, A0lo, A0hi, &
                                     A_1_old, A1lo, A1hi, &
                                     A_2_old, A2lo, A2hi, &
                                     R_0_old, R0lo, R0hi, &
                                     R_1_old, R1lo, R1hi, &
                                     R_2_old, R2lo, R2hi, &
                                     C, Clo, Chi, &
                                     m_start) bind(C, name="ca_sdc_compute_C2_radau")
    ! compute the source term C for the 2nd order Radau update

    ! Here, dt_m is the timestep between time-nodes m and m+1

    use meth_params_module, only : NVAR
    use amrex_constants_module, only : ZERO, HALF, ONE, FIVE
    use burn_type_module, only : burn_t
    use network, only : nspec, nspec_evolve
    use react_util_module

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dt_m, dt
    integer, intent(in) :: Amlo(3), Amhi(3)
    integer, intent(in) :: A0lo(3), A0hi(3)
    integer, intent(in) :: A1lo(3), A1hi(3)
    integer, intent(in) :: A2lo(3), A2hi(3)
    integer, intent(in) :: R0lo(3), R0hi(3)
    integer, intent(in) :: R1lo(3), R1hi(3)
    integer, intent(in) :: R2lo(3), R2hi(3)
    integer, intent(in) :: Clo(3), Chi(3)
    integer, intent(in) :: m_start

    real(rt), intent(in) :: A_m(Amlo(1):Amhi(1), Amlo(2):Amhi(2), Amlo(3):Amhi(3), NVAR)
    real(rt), intent(in) :: A_0_old(A0lo(1):A0hi(1), A0lo(2):A0hi(2), A0lo(3):A0hi(3), NVAR)
    real(rt), intent(in) :: A_1_old(A1lo(1):A1hi(1), A1lo(2):A1hi(2), A1lo(3):A1hi(3), NVAR)
    real(rt), intent(in) :: A_2_old(A2lo(1):A2hi(1), A2lo(2):A2hi(2), A2lo(3):A2hi(3), NVAR)

    real(rt), intent(in) :: R_0_old(R0lo(1):R0hi(1), R0lo(2):R0hi(2), R0lo(3):R0hi(3), NVAR)
    real(rt), intent(in) :: R_1_old(R1lo(1):R1hi(1), R1lo(2):R1hi(2), R1lo(3):R1hi(3), NVAR)
    real(rt), intent(in) :: R_2_old(R2lo(1):R2hi(1), R2lo(2):R2hi(2), R2lo(3):R2hi(3), NVAR)

    real(rt), intent(out) :: C(Clo(1):Chi(1), Clo(2):Chi(2), Clo(3):Chi(3), NVAR)

    integer :: i, j, k

    ! construct the source term to the update for 2nd order
    ! Radau

    if (m_start == 0) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                C(i,j,k,:) = -R_1_old(i,j,k,:) + &
                     (A_m(i,j,k,:) - A_0_old(i,j,k,:)) + &
                     (dt/dt_m) * (ONE/12.0_rt) * &
                     FIVE*(A_1_old(i,j,k,:) + R_1_old(i,j,k,:)) - &
                          (A_2_old(i,j,k,:) + R_2_old(i,j,k,:))
             end do
          end do
       end do

    else if (m_start == 1) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                C(i,j,k,:) = -R_2_old(i,j,k,:) + &
                     (A_m(i,j,k,:) - A_1_old(i,j,k,:)) + &
                     (dt/dt_m) * (ONE/3.0_rt) * &
                     (A_1_old(i,j,k,:) + R_1_old(i,j,k,:)) + &
                     (A_2_old(i,j,k,:) + R_2_old(i,j,k,:))
             end do
          end do
       end do
    end if

  end subroutine ca_sdc_compute_C2_radau


  subroutine ca_sdc_update_o2(lo, hi, dt_m, &
                              k_m, kmlo, kmhi, &
                              k_n, knlo, knhi, &
                              A_m, Amlo, Amhi, &
                              R_m_old, Rmlo, Rmhi, &
                              C, Clo, Chi, &
                              sdc_iteration, &
                              m_start) bind(C, name="ca_sdc_update_o2")
    ! update k_m to k_n via advection -- this is a second-order accurate update

    ! Here, dt_m is the timestep between time-nodes m and m+1

    use meth_params_module, only : NVAR
    use amrex_constants_module, only : ZERO, HALF
    use burn_type_module, only : burn_t
    use network, only : nspec, nspec_evolve
    use react_util_module

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dt_m
    integer, intent(in) :: kmlo(3), kmhi(3)
    integer, intent(in) :: knlo(3), knhi(3)
    integer, intent(in) :: Amlo(3), Amhi(3)
    integer, intent(in) :: Rmlo(3), Rmhi(3)
    integer, intent(in) :: Clo(3), Chi(3)
    integer, intent(in) :: sdc_iteration, m_start


    real(rt), intent(in) :: k_m(kmlo(1):kmhi(1), kmlo(2):kmhi(2), kmlo(3):kmhi(3), NVAR)
    real(rt), intent(inout) :: k_n(knlo(1):knhi(1), knlo(2):knhi(2), knlo(3):knhi(3), NVAR)

    real(rt), intent(in) :: A_m(Amlo(1):Amhi(1), Amlo(2):Amhi(2), Amlo(3):Amhi(3), NVAR)
    real(rt), intent(in) :: R_m_old(Rmlo(1):Rmhi(1), Rmlo(2):Rmhi(2), Rmlo(3):Rmhi(3), NVAR)
    real(rt), intent(in) :: C(Clo(1):Chi(1), Clo(2):Chi(2), Clo(3):Chi(3), NVAR)

    integer :: i, j, k

    type(burn_t) :: burn_state

    real(rt) :: U_old(NVAR), U_new(NVAR), R_full(NVAR), C_zone(NVAR)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             U_old(:) = k_m(i,j,k,:)
             C_zone(:) = C(i,j,k,:)

             ! only burn if we are within the temperature and density
             ! limits for burning
             if (.not. okay_to_burn(U_old)) then
                R_full(:) = ZERO

             else

                ! this is the full state -- this will be updated as we
                ! solve the nonlinear system.  We want to start with a
                ! good initial guess.  For later iterations, we should
                ! begin with the result from the previous iteration.  For
                ! the first iteration, let's try to extrapolate forward
                ! in time.
                if (sdc_iteration == 0) then
                   U_new(:) = U_old(:) + dt_m * A_m(i,j,k,:) + dt_m * R_m_old(i,j,k,:)
                else
                   U_new(:) = k_n(i,j,k,:)
                endif

                call sdc_solve(dt_m, U_old, U_new, C, sdc_iteration)

                ! we solved our system to some tolerance, but let's be sure we are conservative by
                ! reevaluating the reactions and then doing the full step update
                call single_zone_react_source(U_new, R_full, i, j, k, burn_state)

             end if

             U_new(:) = U_old(:) + dt_m * R_full(:) + dt_m * C_zone(:)

             ! copy back to k_n
             k_n(i,j,k,:) = U_new(:)

          enddo
       enddo
    enddo

  end subroutine ca_sdc_update_o2


  subroutine ca_sdc_update_centers_o4(lo, hi, dt_m, &
                                      U_old, Uo_lo, Uo_hi, &
                                      U_new, Un_lo, Un_hi, &
                                      C, C_lo, C_hi, &
                                      sdc_iteration) &
                                      bind(C, name="ca_sdc_update_centers_o4")
    ! update U_old to U_new on cell-centers.  This is an implicit
    ! solve because of reactions.  Here U_old corresponds to time node
    ! m and U_new is node m+1.  dt_m is the timestep between m and
    ! m+1

    use meth_params_module, only : NVAR
    use react_util_module, only : okay_to_burn

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dt_m
    integer, intent(in) :: Uo_lo(3), Uo_hi(3)
    integer, intent(in) :: Un_lo(3), Un_hi(3)
    integer, intent(in) :: C_lo(3), C_hi(3)
    integer, intent(in) :: sdc_iteration

    real(rt), intent(in) :: U_old(Uo_lo(1):Uo_hi(1), Uo_lo(2):Uo_hi(2), Uo_lo(3):Uo_hi(3), NVAR)
    real(rt), intent(out) :: U_new(Un_lo(1):Un_hi(1), Un_lo(2):Un_hi(2), Un_lo(3):Un_hi(3), NVAR)
    real(rt), intent(in) :: C(C_lo(1):C_hi(1), C_lo(2):C_hi(2), C_lo(3):C_hi(3), NVAR)

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! we come in with U_new being a guess for the updated solution
             if (okay_to_burn(U_old(i,j,k,:))) then
                call sdc_solve(dt_m, U_old(i,j,k,:), U_new(i,j,k,:), &
                               C(i,j,k,:), sdc_iteration)
             else
                ! no reactions, so it is a straightforward update
                U_new(i,j,k,:) = U_old(i,j,k,:) + dt_m * C(i,j,k,:)

             end if
          enddo
       enddo
    enddo

  end subroutine ca_sdc_update_centers_o4


  subroutine ca_sdc_conservative_update(lo, hi, dt_m, &
                                        U_old, Uo_lo, Uo_hi, &
                                        U_new, Un_lo, Un_hi, &
                                        C, C_lo, C_hi, &
                                        R_new, R_lo, R_hi) &
                                        bind(C, name="ca_sdc_conservative_update")
    ! given <U>_old, <R>_new, and <C>, compute <U>_new

    use meth_params_module, only : NVAR

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dt_m
    integer, intent(in) :: Uo_lo(3), Uo_hi(3)
    integer, intent(in) :: Un_lo(3), Un_hi(3)
    integer, intent(in) :: C_lo(3), C_hi(3)
    integer, intent(in) :: R_lo(3), R_hi(3)

    real(rt), intent(in) :: U_old(Uo_lo(1):Uo_hi(1), Uo_lo(2):Uo_hi(2), Uo_lo(3):Uo_hi(3), NVAR)
    real(rt), intent(out) :: U_new(Un_lo(1):Un_hi(1), Un_lo(2):Un_hi(2), Un_lo(3):Un_hi(3), NVAR)
    real(rt), intent(in) :: C(C_lo(1):C_hi(1), C_lo(2):C_hi(2), C_lo(3):C_hi(3), NVAR)
    real(rt), intent(in) :: R_new(R_lo(1):R_hi(1), R_lo(2):R_hi(2), R_lo(3):R_hi(3), NVAR)

    integer :: i, j, k

    ! now consider the reacting system
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             U_new(i,j,k,:) = U_old(i,j,k,:) + dt_m * R_new(i,j,k,:) + dt_m * C(i,j,k,:)
          enddo
       enddo
    enddo

  end subroutine ca_sdc_conservative_update


  subroutine ca_instantaneous_react(lo, hi, &
                                    state, s_lo, s_hi, &
                                    R_source, r_lo, r_hi) &
                                    bind(C, name="ca_instantaneous_react")

    use amrex_constants_module, only : ZERO
    use burn_type_module
    use meth_params_module, only : NVAR, NQ, NQAUX
    use network, only : nspec, nspec_evolve
    use react_util_module


    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: r_lo(3), r_hi(3)

    real(rt), intent(in) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)
    real(rt), intent(inout) :: R_source(r_lo(1):r_hi(1), r_lo(2):r_hi(2), r_lo(3):r_hi(3), NVAR)

    integer :: i, j, k
    type(burn_t) :: burn_state

    ! convert from cons to prim -- show this be here or in C++-land?
    ! or should I do things like we do in burn_state and convert it manually?
    ! (in that case, I am not sure if I can assume UTEMP is defined)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (okay_to_burn(state(i,j,k,:))) then
                call single_zone_react_source(state(i,j,k,:), R_source(i,j,k,:), i,j,k, burn_state)
             else
                R_source(i,j,k,:) = ZERO
             end if
          enddo
       enddo
    enddo

  end subroutine ca_instantaneous_react

#endif

end module sdc_util
