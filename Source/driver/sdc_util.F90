module rpar_sdc_module

  use network, only : nspec, nspec_evolve
  implicit none

  integer, parameter :: irp_f_source = 0  ! nspec_evolve + 2 components
  integer, parameter :: irp_dt = irp_f_source + nspec_evolve + 2
  integer, parameter :: irp_mom = irp_dt + 1    ! 3 components
  integer, parameter :: irp_evar = irp_mom + 3
  integer, parameter :: irp_spec = irp_evar + 1  ! nspec - nspec_evolve components

  integer, parameter :: n_rpar = nspec_evolve + 7 + (nspec - nspec_evolve)

end module rpar_sdc_module


module sdc_util

  use amrex_fort_module, only : rt => amrex_real
  use amrex_error_module, only : amrex_error

  implicit none

contains

  subroutine sdc_solve(dt_m, U_old, U_new, C, sdc_iteration)

    ! the purpose of this function is to solve the system
    ! U - dt R(U) = U_old + dt C

    ! here, U_new should come in as a guess for the new U (for
    ! sdc_solver = 1) and will be returned with the value that
    ! satisfies the nonlinear function

    use meth_params_module, only : NVAR, UEDEN, UEINT, URHO, UFS, UMX, UMZ, UTEMP, &
         sdc_solver, sdc_solver_tol, sdc_solve_for_rhoe
    use amrex_constants_module, only : ZERO, HALF, ONE
    use burn_type_module, only : burn_t
    use network, only : nspec, nspec_evolve
    use react_util_module
    use rpar_sdc_module
    use extern_probin_module, only : SMALL_X_SAFE

    implicit none

    real(rt), intent(in) :: dt_m
    real(rt), intent(in) :: U_old(NVAR)
    real(rt), intent(inout) :: U_new(NVAR)
    real(rt), intent(in) :: C(NVAR)
    integer, intent(in) :: sdc_iteration

    real(rt) :: Jac(0:nspec_evolve+1, 0:nspec_evolve+1)
    real(rt) :: w(0:nspec_evolve+1)

    real(rt) :: rpar(0:n_rpar-1)
    integer :: ipar

    integer :: ipvt(nspec_evolve+2)
    integer :: info, istate, iopt

    logical :: converged

    integer, parameter :: lrw = 22 + 9*(nspec_evolve+2) + 2*(nspec_evolve+2)**2
    integer, parameter :: liw = 30 + nspec_evolve + 2

    real(rt) :: rwork(lrw)
    integer :: iwork(liw)
    real(rt) :: time

    ! we will do the implicit update of only the terms that have reactive sources
    !
    !   0               : rho
    !   1:nspec_evolve  : species
    !   nspec_evolve+1  : (rho E) or (rho e)

    real(rt) :: U_react(0:nspec_evolve+1), f_source(0:nspec_evolve+1), R_react(0:nspec_evolve+1), C_react(0:nspec_evolve+1)
    real(rt) :: dU_react(0:nspec_evolve+1), f(0:nspec_evolve+1), f_rhs(0:nspec_evolve+1)

    integer :: m, n

    real(rt) :: err
    real(rt), parameter :: tol = 1.e-5_rt
    integer, parameter :: MAX_ITER = 100
    integer :: iter

    integer, parameter :: NEWTON_SOLVE = 1
    integer, parameter :: VODE_SOLVE = 2
    integer :: solver

    if (sdc_solver == 1) then
       solver = NEWTON_SOLVE
    else if (sdc_solver == 2) then
       solver = VODE_SOLVE
    else if (sdc_solver == 3) then
       if (sdc_iteration == 0) then
          solver = VODE_SOLVE
       else
          solver = NEWTON_SOLVE
       endif
    else
       call amrex_error("invalid sdc_solver")
    endif

    ! update the momenta for this zone -- they don't react
    U_new(UMX:UMZ) = U_old(UMX:UMZ) + dt_m * C(UMX:UMZ)

    ! update the non-reacting species
    U_new(UFS+nspec_evolve:UFS-1+nspec) = U_old(UFS+nspec_evolve:UFS-1+nspec) + &
         dt_m * C(UFS+nspec_evolve:UFS-1+nspec)

    ! now only save the subset that participates in the nonlinear
    ! solve -- note: we include the old state in f_source

    ! load rpar
    if (solver == NEWTON_SOLVE) then

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
    else

       ! if we are solving the system as an ODE, then we
       ! are solving
       !    dU/dt = R(U) + C
       ! so we simply pass in C
       C_react(0) = C(URHO)
       C_react(1:nspec_evolve) = C(UFS:UFS-1+nspec_evolve)
       C_react(nspec_evolve+1) = C(UEINT)

       rpar(irp_f_source:irp_f_source-1+nspec_evolve+2) = C_react(:)
       rpar(irp_dt) = dt_m
       rpar(irp_mom:irp_mom-1+3) = U_new(UMX:UMZ)
    endif

    ! we should be able to do an update for this somehow?
    if (sdc_solve_for_rhoe == 1) then
       rpar(irp_evar) = U_new(UEDEN)
    else
       rpar(irp_evar) = U_new(UEINT)
    endif

    rpar(irp_spec:irp_spec-1+(nspec-nspec_evolve)) = &
         U_new(UFS+nspec_evolve:UFS-1+nspec)

    ! store the subset for the nonlinear solve
    if (sdc_solver == NEWTON_SOLVE) then

       ! Newton solve -- we use an initial guess if possible
       U_react(0) = U_new(URHO)
       U_react(1:nspec_evolve) = U_new(UFS:UFS-1+nspec_evolve)
       if (sdc_solve_for_rhoe == 1) then
          U_react(nspec_evolve+1) = U_new(UEINT)
       else
          U_react(nspec_evolve+1) = U_new(UEDEN)
       endif
    else

       ! VODE ODE solve -- we only consider (rho e), not (rho E)
       U_react(0) = U_old(URHO)
       U_react(1:nspec_evolve) = U_old(UFS:UFS-1+nspec_evolve)
       U_react(nspec_evolve+1) = U_old(UEINT)
    endif

    if (sdc_solver == NEWTON_SOLVE) then
       ! do a simple Newton solve

       err = 1.e30_rt

       ! iterative loop
       iter = 0
       converged = .false.
       do while (.not. converged .and. iter < MAX_ITER)

          call f_sdc_jac(nspec_evolve+2, U_react, f, Jac, nspec_evolve+2, info, rpar)

          ! solve the linear system: Jac dU_react = -f
          call dgefa(Jac, nspec_evolve+2, nspec_evolve+2, ipvt, info)
          if (info /= 0) then
             call amrex_error("singular matrix")
          endif

          f_rhs(:) = -f(:)
          call dgesl(Jac, nspec_evolve+2, nspec_evolve+2, ipvt, f_rhs, 0)

          dU_react(:) = f_rhs(:)

          U_react(:) = U_react(:) + dU_react(:)

          ! construct the norm of the correction -- only worry about
          ! species here, and use some protection against divide by 0
          w(:) = abs(dU_react(:)/(U_react(:) + SMALL_X_SAFE))

          err = sqrt(sum(w(1:nspec_evolve)**2))

          if (err < tol) then
             converged = .true.
          endif

          iter = iter + 1
       enddo

       if (.not. converged) then
          call amrex_error("did not converge in SDC")
       endif

    else if (sdc_solver == VODE_SOLVE) then

       ! use VODE to do the solve

       istate = 1
       iopt = 0
       iwork(:) = 0
       rwork(:) = ZERO
       time = ZERO
       call dvode(f_ode, nspec_evolve+2, U_react, time, dt_m, &
                  1, sdc_solver_tol, 1.e-100_rt, &
                  1, istate, iopt, rwork, lrw, iwork, liw, jac_ode, 22, rpar, ipar)

       if (istate < 0) then
          call amrex_error("vode termination poorly, istate = ", istate)
       endif

    endif

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

  end subroutine sdc_solve


  subroutine f_ode(n, t, U, dUdt, rpar, ipar)

    ! this is the righthand side for the ODE system that we will use
    ! with VODE

    use meth_params_module, only : nvar, URHO, UFS, UEDEN, UTEMP, UMX, UMZ, UEINT
    use burn_type_module
    use rpar_sdc_module
    use react_util_module

    implicit none

    integer, intent(in) :: n
    real(rt), intent(in) :: t
    real(rt), intent(in) :: U(0:n-1)
    real(rt), intent(out) :: dUdt(0:n-1)
    real(rt), intent(in) :: rpar(0:n_rpar-1)
    integer, intent(in) :: ipar

    real(rt) :: U_full(nvar),  R_full(nvar)
    real(rt) :: R_react(0:n-1), C_react(0:n-1)
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

    ! we'll get the temperature here, so just initialize it to something
    U_full(UTEMP) = 1.e4_rt

    call single_zone_react_source(U_full, R_full, 0,0,0, burn_state)

    R_react(0) = R_full(URHO)
    R_react(1:nspec_evolve) = R_full(UFS:UFS-1+nspec_evolve)
    R_react(nspec_evolve+1) = R_full(UEINT)

    ! C comes in through rpar
    C_react(:) = rpar(irp_f_source:irp_f_source-1+nspec_evolve+2)

    ! create the RHS
    dUdt(:) = R_react(:) + C_react(:)

  end subroutine f_ode

  subroutine jac_ode(neq, time, U, ml, mu, pd, nrowpd, rpar, ipar)

    use rpar_sdc_module
    implicit none

    integer   , intent(IN   ) :: neq, ml, mu, nrowpd, ipar
    real(rt), intent(INOUT) :: U(neq), rpar(n_rpar), time
    real(rt), intent(  OUT) :: pd(neq,neq)

    ! this is a stub -- we are using a numerical Jacobian at the moment

  end subroutine jac_ode

  subroutine f_sdc_jac(n, U, f, Jac, ldjac, iflag, rpar)

    use rpar_sdc_module
    use meth_params_module, only : nvar, URHO, UFS, UEINT, UEDEN, UMX, UMZ, UTEMP, sdc_solve_for_rhoe
    use network, only : nspec, nspec_evolve
    use burn_type_module
    use react_util_module
    use eos_type_module, only : eos_t, eos_input_re
    use eos_module, only : eos
    use amrex_constants_module, only : ZERO, HALF, ONE

    ! this computes the function we need to zero for the SDC update
    implicit none

    integer,intent(in) :: n, ldjac
    real(rt), intent(in)  :: U(0:n-1)
    real(rt), intent(out) :: f(0:n-1)
    real(rt), intent(out) :: Jac(0:ldjac-1,0:n-1)
    integer, intent(inout) :: iflag  !! leave this untouched
    real(rt), intent(in) :: rpar(0:n_rpar-1)

    real(rt) :: U_full(nvar),  R_full(nvar)
    real(rt) :: R_react(0:n-1), f_source(0:n-1)
    type(burn_t) :: burn_state
    type(eos_t) :: eos_state
    real(rt) :: dt_m

    real(rt) :: denom
    real(rt) :: dRdw(0:nspec_evolve+1, 0:nspec_evolve+1), dwdU(0:nspec_evolve+1, 0:nspec_evolve+1)
    integer :: m

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

    ! unpack rpar
    dt_m = rpar(irp_dt)
    f_source(:) = rpar(irp_f_source:irp_f_source-1+nspec_evolve+2)

    ! compute the temperature and species derivatives --
    ! maybe this should be done using the burn_state
    ! returned by single_zone_react_source, since it is
    ! more consistent T from e
    eos_state % rho = U_full(URHO)
    eos_state % T = 1.e6_rt   ! initial guess
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

    ! get dRdw
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

    ! now the T row -- this depends on whether we are evolving (rho E) or (rho e)
    denom = ONE/(eos_state % rho * eos_state % dedT)
    if (sdc_solve_for_rhoe == 1) then
       dwdU(nspec_evolve+1,0) = denom*(sum(eos_state % xn(1:nspec_evolve) * eos_state % dedX(1:nspec_evolve)) - &
                                       eos_state % rho * eos_state % dedr - eos_state % e)
    else
       dwdU(nspec_evolve+1,0) = denom*(sum(eos_state % xn(1:nspec_evolve) * eos_state % dedX(1:nspec_evolve)) - &
                                       eos_state % rho * eos_state % dedr - eos_state % e + &
                                       HALF*sum(U_full(UMX:UMZ)**2)/eos_state % rho**2)
    endif

    do m = 1, nspec_evolve
       dwdU(nspec_evolve+1,m) = -denom * eos_state % dedX(m)
    enddo

    dwdU(nspec_evolve+1, nspec_evolve+1) = denom

    ! construct the Jacobian -- we can get most of the
    ! terms from the network itself, but we do not rely on
    ! it having derivative wrt density
    Jac(:, :) = ZERO
    do m = 0, nspec_evolve+1
       Jac(m, m) = ONE
    enddo

    Jac(:,:) = Jac(:,:) - dt_m * matmul(dRdw, dwdU)

    f(:) = U(:) - dt_m * R_react(:) - f_source(:)

  end subroutine f_sdc_jac


  subroutine ca_sdc_update_advection_o2(lo, hi, dt_m, &
                                        k_m, kmlo, kmhi, &
                                        k_n, knlo, knhi, &
                                        A_m, Amlo, Amhi, &
                                        A_0_old, A0lo, A0hi, &
                                        A_1_old, A1lo, A1hi, &
                                        m_start) bind(C, name="ca_sdc_update_advection_o2")

    ! update k_m to k_n via advection -- this is a second-order accurate update

    use meth_params_module, only : NVAR
    use amrex_constants_module, only : HALF

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dt_m
    integer, intent(in) :: kmlo(3), kmhi(3)
    integer, intent(in) :: knlo(3), knhi(3)
    integer, intent(in) :: Amlo(3), Amhi(3)
    integer, intent(in) :: A0lo(3), A0hi(3)
    integer, intent(in) :: A1lo(3), A1hi(3)
    integer, intent(in) :: m_start


    real(rt), intent(in) :: k_m(kmlo(1):kmhi(1), kmlo(2):kmhi(2), kmlo(3):kmhi(3), NVAR)
    real(rt), intent(inout) :: k_n(knlo(1):knhi(1), knlo(2):knhi(2), knlo(3):knhi(3), NVAR)

    real(rt), intent(in) :: A_m(Amlo(1):Amhi(1), Amlo(2):Amhi(2), Amlo(3):Amhi(3), NVAR)
    real(rt), intent(in) :: A_0_old(A0lo(1):A0hi(1), A0lo(2):A0hi(2), A0lo(3):A0hi(3), NVAR)
    real(rt), intent(in) :: A_1_old(A1lo(1):A1hi(1), A1lo(2):A1hi(2), A1lo(3):A1hi(3), NVAR)

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             k_n(i,j,k,:) = k_m(i,j,k,:) + HALF * dt_m * (A_0_old(i,j,k,:) + A_1_old(i,j,k,:))
          enddo
       enddo
    enddo

  end subroutine ca_sdc_update_advection_o2


  subroutine ca_sdc_update_advection_o4(lo, hi, dt, &
                                        k_m, kmlo, kmhi, &
                                        k_n, knlo, knhi, &
                                        A_m, Amlo, Amhi, &
                                        A_0_old, A0lo, A0hi, &
                                        A_1_old, A1lo, A1hi, &
                                        A_2_old, A2lo, A2hi, &
                                        m_start) bind(C, name="ca_sdc_update_advection_o4")

    ! update k_m to k_n via advection -- this is a second-order accurate update
    ! dt is the total timestep from n to n+1

    use meth_params_module, only : NVAR
    use amrex_constants_module, only : HALF, TWO, FIVE, EIGHT

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dt
    integer, intent(in) :: kmlo(3), kmhi(3)
    integer, intent(in) :: knlo(3), knhi(3)
    integer, intent(in) :: Amlo(3), Amhi(3)
    integer, intent(in) :: A0lo(3), A0hi(3)
    integer, intent(in) :: A1lo(3), A1hi(3)
    integer, intent(in) :: A2lo(3), A2hi(3)
    integer, intent(in) :: m_start


    real(rt), intent(in) :: k_m(kmlo(1):kmhi(1), kmlo(2):kmhi(2), kmlo(3):kmhi(3), NVAR)
    real(rt), intent(inout) :: k_n(knlo(1):knhi(1), knlo(2):knhi(2), knlo(3):knhi(3), NVAR)

    real(rt), intent(in) :: A_m(Amlo(1):Amhi(1), Amlo(2):Amhi(2), Amlo(3):Amhi(3), NVAR)
    real(rt), intent(in) :: A_0_old(A0lo(1):A0hi(1), A0lo(2):A0hi(2), A0lo(3):A0hi(3), NVAR)
    real(rt), intent(in) :: A_1_old(A1lo(1):A1hi(1), A1lo(2):A1hi(2), A1lo(3):A1hi(3), NVAR)
    real(rt), intent(in) :: A_2_old(A2lo(1):A2hi(1), A2lo(2):A2hi(2), A2lo(3):A2hi(3), NVAR)

    integer :: i, j, k
    real(rt) :: dt_m


    dt_m = HALF * dt

    if (m_start == 0) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                k_n(i,j,k,:) = k_m(i,j,k,:) + &
                     dt_m * (A_m(i,j,k,:) - A_0_old(i,j,k,:)) + &
                     dt/24.0_rt * (FIVE*A_0_old(i,j,k,:) + EIGHT*A_1_old(i,j,k,:) - A_2_old(i,j,k,:))
             enddo
          enddo
       enddo

    else if (m_start == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                k_n(i,j,k,:) = k_m(i,j,k,:) + &
                     dt_m * (A_m(i,j,k,:) - A_1_old(i,j,k,:)) + &
                     dt/24.0_rt * (-A_0_old(i,j,k,:) + EIGHT*A_1_old(i,j,k,:) + FIVE*A_2_old(i,j,k,:))
             enddo
          enddo
       enddo

    else
       call amrex_error("error in ca_sdc_update_advection_o4 -- shouldn't be here")
    endif

  end subroutine ca_sdc_update_advection_o4


#ifdef REACTIONS
  subroutine ca_sdc_compute_C4(lo, hi, &
                               A_m, Amlo, Amhi, &
                               A_0_old, A0lo, A0hi, &
                               A_1_old, A1lo, A1hi, &
                               A_2_old, A2lo, A2hi, &
                               R_0_old, R0lo, R0hi, &
                               R_1_old, R1lo, R1hi, &
                               R_2_old, R2lo, R2hi, &
                               C, Clo, Chi, &
                               m_start) bind(C, name="ca_sdc_compute_C4")

    ! compute the 'C' term for the 4th-order solve with reactions

    use meth_params_module, only : NVAR
    use amrex_constants_module, only : HALF, TWO, FIVE, EIGHT, TWELFTH

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

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! compute the integral (without the dt)
             if (m_start == 0) then
                integral(:) = TWELFTH * (FIVE*(A_0_old(i,j,k,:) + R_0_old(i,j,k,:)) + &
                                         EIGHT*(A_1_old(i,j,k,:) + R_1_old(i,j,k,:)) - &
                                         (A_2_old(i,j,k,:) + R_2_old(i,j,k,:)))

                C(i,j,k,:) = (A_m(i,j,k,:) - A_0_old(i,j,k,:)) - R_1_old(i,j,k,:) + integral

             else if (m_start == 1) then
                integral(:) = TWELFTH * (-(A_0_old(i,j,k,:) + R_0_old(i,j,k,:)) + &
                                         EIGHT*(A_1_old(i,j,k,:) + R_1_old(i,j,k,:)) + &
                                         FIVE*(A_2_old(i,j,k,:) + R_2_old(i,j,k,:)))

                C(i,j,k,:) = (A_m(i,j,k,:) - A_1_old(i,j,k,:)) - R_2_old(i,j,k,:) + integral

             else
                call amrex_error("error in ca_sdc_compute_C4 -- shouldn't be here")
             endif

          enddo
       enddo
    enddo

  end subroutine ca_sdc_compute_C4


  subroutine ca_sdc_update_o2(lo, hi, dt_m, &
                              k_m, kmlo, kmhi, &
                              k_n, knlo, knhi, &
                              A_m, Amlo, Amhi, &
                              A_0_old, A0lo, A0hi, &
                              A_1_old, A1lo, A1hi, &
                              R_0_old, R0lo, R0hi, &
                              R_1_old, R1lo, R1hi, &
                              sdc_iteration, &
                              m_start) bind(C, name="ca_sdc_update_o2")

    ! update k_m to k_n via advection -- this is a second-order accurate update

    use meth_params_module, only : NVAR
    use amrex_constants_module, only : HALF
    use burn_type_module, only : burn_t
    use network, only : nspec, nspec_evolve
    use react_util_module

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dt_m
    integer, intent(in) :: kmlo(3), kmhi(3)
    integer, intent(in) :: knlo(3), knhi(3)
    integer, intent(in) :: Amlo(3), Amhi(3)
    integer, intent(in) :: A0lo(3), A0hi(3)
    integer, intent(in) :: A1lo(3), A1hi(3)
    integer, intent(in) :: R0lo(3), R0hi(3)
    integer, intent(in) :: R1lo(3), R1hi(3)
    integer, intent(in) :: sdc_iteration, m_start


    real(rt), intent(in) :: k_m(kmlo(1):kmhi(1), kmlo(2):kmhi(2), kmlo(3):kmhi(3), NVAR)
    real(rt), intent(inout) :: k_n(knlo(1):knhi(1), knlo(2):knhi(2), knlo(3):knhi(3), NVAR)

    real(rt), intent(in) :: A_m(Amlo(1):Amhi(1), Amlo(2):Amhi(2), Amlo(3):Amhi(3), NVAR)
    real(rt), intent(in) :: A_0_old(A0lo(1):A0hi(1), A0lo(2):A0hi(2), A0lo(3):A0hi(3), NVAR)
    real(rt), intent(in) :: A_1_old(A1lo(1):A1hi(1), A1lo(2):A1hi(2), A1lo(3):A1hi(3), NVAR)

    real(rt), intent(in) :: R_0_old(R0lo(1):R0hi(1), R0lo(2):R0hi(2), R0lo(3):R0hi(3), NVAR)
    real(rt), intent(in) :: R_1_old(R1lo(1):R1hi(1), R1lo(2):R1hi(2), R1lo(3):R1hi(3), NVAR)

    integer :: i, j, k

    type(burn_t) :: burn_state

    real(rt) :: U_old(NVAR), U_new(NVAR), C(NVAR), R_full(NVAR)


    ! now consider the reacting system
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             U_old(:) = k_m(i,j,k,:)

             ! construct the source term to the update
             ! for 2nd order, there is no advective correction, and we have
             ! C = - R(U^{m+1,k}) + I_m^{m+1}/dt
             C(:) = -R_1_old(i,j,k,:) + &
                  HALF * (A_0_old(i,j,k,:) + A_1_old(i,j,k,:)) + &
                  HALF * (R_0_old(i,j,k,:) + R_1_old(i,j,k,:))

             ! this is the full state -- this will be updated as we
             ! solve the nonlinear system.  We want to start with a
             ! good initial guess.  For later iterations, we should
             ! begin with the result from the previous iteration.  For
             ! the first iteration, let's try to extrapolate forward
             ! in time.
             if (sdc_iteration == 0) then
                U_new(:) = U_old(:) + dt_m * A_m(i,j,k,:) + dt_m * R_0_old(i,j,k,:)
             else
                U_new(:) = k_n(i,j,k,:)
             endif

             call sdc_solve(dt_m, U_old, U_new, C, sdc_iteration)

             ! we solved our system to some tolerance, but let's be sure we are conservative by
             ! reevaluating the reactions and then doing the full step update
             call single_zone_react_source(U_new, R_full, i, j, k, burn_state)

             ! redo the update of the momenta to reduce accumulation of roundoff
             U_new(:) = U_old(:) + dt_m * R_full(:) + dt_m * C(:)

             ! copy back to k_n
             k_n(i,j,k,:) = U_new(:)

          enddo
       enddo
    enddo

  end subroutine ca_sdc_update_o2

  subroutine ca_instantaneous_react(lo, hi, &
                                    state, s_lo, s_hi, &
                                    R_source, r_lo, r_hi) &
                                    bind(C, name="ca_instantaneous_react")

    use amrex_constants_module, only : ZERO
    use burn_type_module
    use meth_params_module, only : NVAR, NQ, NQAUX, QFS, QRHO, QTEMP, UFS, UEDEN, UEINT
    use network, only : nspec, nspec_evolve, aion
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
             call single_zone_react_source(state(i,j,k,:), R_source(i,j,k,:), i,j,k, burn_state)
          enddo
       enddo
    enddo

  end subroutine ca_instantaneous_react

  subroutine ca_store_reaction_state(lo, hi, &
                                     R_old, r_lo, r_hi, &
                                     state, s_lo, s_hi, &
                                     R_store, rs_lo, rs_hi) &
                                     bind(C, name="ca_store_reaction_state")

    use meth_params_module, only : NVAR, URHO, UEDEN, UFS
    use network, only : nspec

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: r_lo(3), r_hi(3)
    integer, intent(in) :: rs_lo(3), rs_hi(3)

    real(rt), intent(in) :: R_old(r_lo(1):r_hi(1), r_lo(2):r_hi(2), r_lo(3):r_hi(3), NVAR)
    real(rt), intent(in) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)
    real(rt), intent(inout) :: R_store(rs_lo(1):rs_hi(1), rs_lo(2):rs_hi(2), rs_lo(3):rs_hi(3), nspec+2)

    integer :: i, j, k

    ! copy the data from the last node's reactive source to the state data

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! for R_store we use the indices defined in Castro_setup.cpp for
             ! Reactions_Type
             R_store(i,j,k,1:nspec) = R_old(i,j,k,UFS:UFS-1+nspec)/state(i,j,k,URHO)
             R_store(i,j,k,nspec+1) = R_old(i,j,k,UEDEN)/state(i,j,k,URHO)
             R_store(i,j,k,nspec+2) = R_old(i,j,k,UEDEN)
          enddo
       enddo
    enddo

  end subroutine ca_store_reaction_state

#endif

end module sdc_util
