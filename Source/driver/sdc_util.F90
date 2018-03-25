module sdc_util

  use amrex_fort_module, only : rt => amrex_real

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
    use bl_constants_module, only : HALF

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
    use bl_constants_module, only : HALF, TWO, FIVE, EIGHT

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
       call bl_error("error in ca_sdc_update_advection_o4 -- shouldn't be here")
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

    use meth_params_module, only : NVAR, UEDEN, URHO, UFS
    use bl_constants_module, only : HALF
    use burn_type_module, only : burn_t
    use network, only : nspec

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

    real(rt) :: err, tol
    real(rt) :: C(nvar), R_full(nvar)
    real(rt) :: U_react(nspec+2), C_react(nspec+2), U_new(nspec+2), R_react(nspec+2)


    ! for those variables without reactive sources, we can do the
    ! explicit update -- we'll do that here for everything, and then
    ! overwrite the reacting ones
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             k_n(i,j,k,:) = k_m(i,j,k,:) + HALF * dt_m * (A_0_old(i,j,k,:) + A_1_old(i,j,k,:))
          enddo
       enddo
    enddo

    ! now consider the reacting system
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             U_react(1:nspec) = k_m(i,j,k,UFS:UFS-1+nspec)
             U_react(nspec+1) = k_m(i,j,k,URHO)

             ! we have a choice of which energy variable to update
             U_react(nspec+2) = k_m(i,j,k,UEDEN)

             ! construct the source term to the update
             ! for 2nd order, there is no advective correction, and we have
             ! C = U^{m,(k+1)} - dt * R(U^{m+1,k}) + I_m^{m+1}
             C(:) = U_react(:) - dt_m * R_1_old(i,j,k,:) + &
                  HALF * dt_m * (A_0_old(i,j,k,:) + A_1_old(i,j,k,:)) + &
                  HALF * dt_m * (R_0_old(i,j,k,:) + R_1_old(i,j,k,:))

             ! now only save the subset
             C_react(1:nspec) = C(1:nspec)
             C_react(nspec+1) = C(URHO)
             C_react(nspec+2) = C(UEDEN)  ! need to consider which energy

             ! set the initial guess
             U_new(:) = U_react(:)

             ! iterative loop
             do while (err > tol)

                ! get R for the new guess
                call single_zone_react_source(U_new, R_full, i,j,k)
                R_react(URHO) = R_full(URHO)
                R_react(UFS:UFS-1+nspec) = R_full(UFS:UFS-1+nspec)
                R_react(UEDEN) = R_full(UEDEN)

                ! construct the Jacobian -- we can get most of the
                ! terms from the network itself, but we do not rely on
                ! it having derivative wrt density
                do m = 1, nspec+2
                   J(m, m) = ONE
                enddo

                ! solve the linear system

                ! correct

                ! construct the norm

             enddo

             ! if we updated total energy, then correct internal, or vice versa

             ! copy back to k_n

          enddo
       enddo
    enddo

  end subroutine ca_sdc_update_o2

  subroutine ca_instantaneous_react(lo, hi, &
                                    state, s_lo, s_hi, &
                                    R_source, r_lo, r_hi) &
                                    bind(C, name="ca_instantaneous_react")

    use bl_constants_module, only : ZERO
    use burn_type_module
    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : NVAR, NQ, NQAUX, QFS, QRHO, QTEMP, UFS, UEDEN, UEINT
    use network, only : nspec, aion
    use react_util_module


    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: r_lo(3), r_hi(3)

    real(rt), intent(in) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)
    real(rt), intent(inout) :: R_source(r_lo(1):r_hi(1), r_lo(2):r_hi(2), r_lo(3):r_hi(3), NVAR)

    integer :: i, j, k

    ! convert from cons to prim -- show this be here or in C++-land?
    ! or should I do things like we do in burn_state and convert it manually?
    ! (in that case, I am not sure if I can assume UTEMP is defined)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             call single_zone_react_source(state(i,j,k,:), R_source(i,j,k,:), i,j,k)
          enddo
       enddo
    enddo

  end subroutine ca_instantaneous_react

#endif

end module sdc_util
