module sdc_util

  use amrex_fort_module, only : rt => amrex_real
  use amrex_error_module, only : amrex_error

  implicit none

contains

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
  subroutine ca_sdc_update_o2(lo, hi, dt_m, &
                              k_m, kmlo, kmhi, &
                              k_n, knlo, knhi, &
                              A_m, Amlo, Amhi, &
                              A_0_old, A0lo, A0hi, &
                              A_1_old, A1lo, A1hi, &
                              R_0_old, R0lo, R0hi, &
                              R_1_old, R1lo, R1hi, &
                              m_start) bind(C, name="ca_sdc_update_o2")

    ! update k_m to k_n via advection -- this is a second-order accurate update

    use meth_params_module, only : NVAR, UEDEN, UEINT, URHO, UFS, UMX, UMZ, UTEMP
    use amrex_constants_module, only : ZERO, HALF, ONE
    use burn_type_module, only : burn_t
    use eos_type_module, only : eos_t, eos_input_re
    use eos_module
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
    integer, intent(in) :: m_start


    real(rt), intent(in) :: k_m(kmlo(1):kmhi(1), kmlo(2):kmhi(2), kmlo(3):kmhi(3), NVAR)
    real(rt), intent(inout) :: k_n(knlo(1):knhi(1), knlo(2):knhi(2), knlo(3):knhi(3), NVAR)

    real(rt), intent(in) :: A_m(Amlo(1):Amhi(1), Amlo(2):Amhi(2), Amlo(3):Amhi(3), NVAR)
    real(rt), intent(in) :: A_0_old(A0lo(1):A0hi(1), A0lo(2):A0hi(2), A0lo(3):A0hi(3), NVAR)
    real(rt), intent(in) :: A_1_old(A1lo(1):A1hi(1), A1lo(2):A1hi(2), A1lo(3):A1hi(3), NVAR)

    real(rt), intent(in) :: R_0_old(R0lo(1):R0hi(1), R0lo(2):R0hi(2), R0lo(3):R0hi(3), NVAR)
    real(rt), intent(in) :: R_1_old(R1lo(1):R1hi(1), R1lo(2):R1hi(2), R1lo(3):R1hi(3), NVAR)

    integer :: i, j, k

    type(burn_t) :: burn_state
    type(eos_t) :: eos_state

    real(rt) :: err
    real(rt), parameter :: tol = 1.e-8_rt

    real(rt) :: U_new(NVAR), C(NVAR), R_full(NVAR)

    ! we will do the implicit update of only the terms that have reactive sources
    !
    !   0               : rho
    !   1:nspec_evolve  : species
    !   nspec_evolve+1  : rho E

    real(rt) :: U_react(0:nspec_evolve+1), C_react(0:nspec_evolve+1), R_react(0:nspec_evolve+1)
    real(rt) :: dU_react(0:nspec_evolve+1), f(0:nspec_evolve+1)

    integer :: m, n
    real(rt) :: Jac(0:nspec_evolve+1, 0:nspec_evolve+1)
    real(rt) :: dRdw(0:nspec_evolve+1, 0:nspec_evolve+1), dwdU(0:nspec_evolve+1, 0:nspec_evolve+1)

    real(rt) :: denom

    integer :: ipvt(nspec_evolve+2)
    integer :: info

    ! now consider the reacting system
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! this is the full state -- this will be updates as we
             ! solve the nonlinear system
             U_new(:) = k_m(i,j,k,:)

             ! construct the source term to the update
             ! for 2nd order, there is no advective correction, and we have
             ! C = U^{m,(k+1)} - dt * R(U^{m+1,k}) + I_m^{m+1}
             C(:) = U_new(:) - dt_m * R_1_old(i,j,k,:) + &
                  HALF * dt_m * (A_0_old(i,j,k,:) + A_1_old(i,j,k,:)) + &
                  HALF * dt_m * (R_0_old(i,j,k,:) + R_1_old(i,j,k,:))

             ! update the momenta for this zone -- this never gets updated again
             U_new(UMX:UMZ) = C(UMX:UMZ)

             ! update the non-reacting species
             U_new(UFS-1+nspec_evolve:UFS-1+nspec) = C(UFS-1+nspec_evolve:UFS-1+nspec)

             ! now only save the subset that participates in the nonlinear solve
             C_react(0) = C(URHO)
             C_react(1:nspec_evolve) = C(UFS:UFS-1+nspec_evolve)
             C_react(nspec_evolve+1) = C(UEDEN)  ! need to consider which energy

             err = 1.e30_rt

             ! iterative loop
             do while (err > tol)

                ! compute the temperature
                eos_state % rho = U_new(URHO)
                eos_state % T = 1.e6_rt   ! initial guess
                eos_state % xn(:) = U_new(UFS:UFS-1+nspec)/U_new(URHO)
                eos_state % e = (U_new(UEDEN) - HALF*sum(U_new(UMX:UMZ))/U_new(URHO))/U_new(URHO)

                call eos(eos_input_re, eos_state)

                U_new(UTEMP) = eos_state % T

                ! get R for the new guess
                call single_zone_react_source(U_new, R_full, i,j,k, burn_state)

                ! store the subset for the nonlinear solve
                U_react(0) = U_new(URHO)
                U_react(1:nspec_evolve) = U_new(UFS:UFS-1+nspec_evolve)
                U_react(nspec_evolve+1) = U_new(UEDEN)  ! we have a choice of which energy variable to update

                R_react(0) = R_full(URHO)
                R_react(1:nspec_evolve) = R_full(UFS:UFS-1+nspec_evolve)
                R_react(nspec_evolve+1) = R_full(UEDEN)

                ! get dRdw
                call single_zone_jac(U_new, burn_state, dRdw)

                ! construct dwdU
                dwdU(:, :) = ZERO

                ! the density row
                dwdU(0, 0) = ONE

                ! the X_k rows
                do n = 1, nspec_evolve
                   dwdU(n,0) = -U_react(n)/U_react(0)**2
                   dwdU(n,n) = ONE/U_react(0)
                enddo

                ! now the T row
                denom = ONE/(eos_state % rho * eos_state % dedT)
                dwdU(nspec_evolve+1,0) = denom*(sum(eos_state % xn(:) * eos_state % dedX(:)) - &
                                                eos_state % rho * eos_state % dedr - eos_state % e + &
                                                HALF*sum(U_new(UMX:UMZ)**2)/eos_state % rho)
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

                ! compute the RHS of the linear system, f
                f(:) = -U_react(:) + dt_m * R_react(:) + C_react(:)

                ! solve the linear system: Jac dU_react = f
                call dgefa(Jac, nspec_evolve+2, nspec_evolve+2, ipvt, info)
                if (info /= 0) then
                   call amrex_error("singular matrix")
                endif

                call dgesl(Jac, nspec_evolve+2, nspec_evolve+2, ipvt, f, 0)

                dU_react(:) = f(:)

                ! correct the full state
                U_new(URHO) = U_new(URHO) + dU_react(0)
                U_new(UFS:UFS-1+nspec_evolve) = U_new(UFS:UFS-1+nspec_evolve) + dU_react(1:nspec_evolve)
                U_new(UEDEN) = U_new(UEDEN) + dU_react(nspec_evolve+1)

                ! if we updated total energy, then correct internal, or vice versa
                U_new(UEINT) = U_new(UEDEN) - HALF*(sum(U_new(UMX:UMZ)**2)/U_new(URHO))

                ! construct the norm of the correction
                err = sum(dU_react**2)/sum(U_react**2)
             enddo

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
