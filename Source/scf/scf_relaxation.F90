module scf_relaxation_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  ! Internal data for SCF relaxation

  ! Position of points A and B relative to system center
  real(rt), save :: scf_r_A(3), scf_r_B(3)

contains

  subroutine scf_setup_relaxation(d_A, axis_A, d_B, axis_B) bind(C, name='scf_setup_relaxation')

    use prob_params_Module, only: center

    implicit none

    real(rt), intent(in), value :: d_A, d_B
    integer,  intent(in), value :: axis_A, axis_B

    ! We need to fix two points to uniquely determine an equilibrium
    ! configuration for a rotating star. We are provided the distance
    ! from the center, and the axis we're measuring on.

    scf_r_A(:) = center
    scf_r_B(:) = center

    ! Note that the sign is somewhat arbitrary here.

    scf_r_A(axis_A) = scf_r_A(axis_A) + d_A
    scf_r_B(axis_B) = scf_r_B(axis_B) + d_B

  end subroutine scf_setup_relaxation




  ! Calculate the maximum allowable enthalpy on the domain.

  subroutine scf_calculate_target_h_max(lo, hi, &
                                        state, s_lo, s_hi, &
                                        maximum_density, &
                                        target_h_max) bind(C, name='scf_calculate_target_h_max')

    use meth_params_module, only: NVAR, URHO, UTEMP, UFS
    use network, only: nspec
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: target_h_max
    real(rt), intent(in   ), value :: maximum_density

    integer  :: i, j, k

    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho = maximum_density
             eos_state%T   = state(i,j,k,UTEMP)
             eos_state%xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)
             eos_state%aux  = state(i,j,k,UFX:UFX+naux-1) / state(i,j,k,URHO)

             call eos(eos_input_rt, eos_state)

             target_h_max = max(target_h_max, eos_state%h)

          enddo
       enddo
    enddo

  end subroutine scf_calculate_target_h_max



  ! Calculate the phi and psi factors that go into
  ! updating the rotation frequency.

  subroutine scf_update_for_omegasq(lo, hi, &
                                    phi, phi_lo, phi_hi, &
                                    psi, psi_lo, psi_hi, &
                                    dx, &
                                    phi_A, psi_A, phi_B, psi_B) bind(C, name='scf_update_for_omegasq')

    use amrex_constants_module, only: ZERO, HALF, ONE
    use meth_params_module, only: NVAR
    use prob_params_module, only: problo

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: psi_lo(3), psi_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ) :: psi(psi_lo(1):psi_hi(1),psi_lo(2):psi_hi(2),psi_lo(3):psi_hi(3))
    real(rt), intent(inout) :: phi_A, phi_B, psi_A, psi_B

    integer  :: i, j, k
    real(rt) :: r(3), rr(3), scale

    ! The below assumes we are rotating on the z-axis.

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1)

             ! Do a trilinear interpolation to find the contribution from
             ! this grid point. Limit so that only the nearest zone centers
             ! can participate. This implies that the maximum allowable
             ! distance from the target location is 0.5 * dx.

             rr = abs(r - scf_r_A) / dx

             if (any(rr > ONE)) then
                scale = ZERO
             else
                scale = (ONE - rr(1)) * (ONE - rr(2)) * (ONE - rr(3))
             end if

             phi_A = phi_A + scale * phi(i,j,k)
             psi_A = psi_A + scale * psi(i,j,k)

          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1)

             rr = abs(r - scf_r_B) / dx

             if (any(rr > ONE)) then
                scale = ZERO
             else
                scale = (ONE - rr(1)) * (ONE - rr(2)) * (ONE - rr(3))
             end if

             phi_B = phi_B + scale * phi(i,j,k)
             psi_B = psi_B + scale * psi(i,j,k)

          enddo
       enddo
    enddo

  end subroutine scf_update_for_omegasq



  subroutine scf_get_bernoulli_const(lo, hi, &
                                     phi, phi_lo, phi_hi, &
                                     phi_rot, phr_lo, phr_hi, &
                                     dx, bernoulli) bind(C, name='scf_get_bernoulli_const')

    use amrex_constants_module, only: ZERO, HALF, ONE
    use meth_params_module, only: NVAR
    use prob_params_module, only: problo

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: phr_lo(3), phr_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ) :: phi_rot(phr_lo(1):phr_hi(1),phr_lo(2):phr_hi(2),phr_lo(3):phr_hi(3))
    real(rt), intent(inout) :: bernoulli

    integer  :: i, j, k
    real(rt) :: r(3), rr(3), scale

    ! The below assumes we are rotating on the z-axis.

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1)

             rr = abs(r - scf_r_A) / dx

             if (any(rr > ONE)) then
                scale = ZERO
             else
                scale = (ONE - rr(1)) * (ONE - rr(2)) * (ONE - rr(3))
             end if

             bernoulli = bernoulli + scale * (phi(i,j,k) + phi_rot(i,j,k))

          enddo
       enddo
    enddo

  end subroutine scf_get_bernoulli_const



  subroutine scf_construct_enthalpy(lo, hi, &
                                    phi, phi_lo, phi_hi, &
                                    phi_rot, phr_lo, phr_hi, &
                                    enthalpy, h_lo, h_hi, &
                                    dx, bernoulli) bind(C, name='scf_construct_enthalpy')

    use meth_params_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: phr_lo(3), phr_hi(3)
    integer,  intent(in   ) :: h_lo(3), h_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ) :: phi_rot(phr_lo(1):phr_hi(1),phr_lo(2):phr_hi(2),phr_lo(3):phr_hi(3))
    real(rt), intent(inout) :: enthalpy(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3))
    real(rt), intent(in   ), value :: bernoulli

    integer  :: i, j, k

    ! The Bernoulli equation says that energy is conserved:
    ! enthalpy + gravitational potential + rotational potential = const
    ! We already have the constant, so our goal is to construct the enthalpy field.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             enthalpy(i,j,k) = bernoulli - phi(i,j,k) - phi_rot(i,j,k)

          enddo
       enddo
    enddo

  end subroutine scf_construct_enthalpy



  subroutine scf_update_density(lo, hi, &
                                state, s_lo, s_hi, &
                                enthalpy, h_lo, h_hi, &
                                dx, actual_rho_max, &
                                actual_h_max, target_h_max, &
                                Linf_norm) bind(C, name='scf_update_density')

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: NVAR, URHO, UTEMP, UMX, UMZ, UEDEN, UEINT, UFS
    use network, only: nspec
    use eos_module, only: eos
    use eos_type_module, only: eos_input_th, eos_t

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: h_lo(3), h_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: enthalpy(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3))
    real(rt), intent(inout) :: Linf_norm
    real(rt), intent(in   ), value :: actual_rho_max, actual_h_max, target_h_max

    integer  :: i, j, k
    real(rt) :: old_rho, drho

    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! We only want to call the EOS for zones with enthalpy > 0.
             ! For distances far enough from the center, the rotation
             ! term can overcome the other terms and make the enthalpy
             ! spuriously negative. If the enthalpy is negative, we just
             ! leave the zone alone -- this should be ambient material.

             if (enthalpy(i,j,k) > ZERO .and. state(i,j,k,URHO) > ZERO) then

                old_rho = state(i,j,k,URHO)

                ! Rescale the enthalpy by the maximum allowed value.

                enthalpy(i,j,k) = target_h_max * (enthalpy(i,j,k) / actual_h_max)

                eos_state%rho = state(i,j,k,URHO) ! Initial guess for the EOS
                eos_state%T   = state(i,j,k,UTEMP)
                eos_state%xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)
                eos_state%aux = state(i,j,k,UFX:UFX+naux-1) / state(i,j,k,URHO)
                eos_state%h   = enthalpy(i,j,k)

                call eos(eos_input_th, eos_state)

                state(i,j,k,URHO)  = eos_state % rho
                state(i,j,k,UTEMP) = eos_state % T
                state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
                state(i,j,k,UFS:UFS+nspec-1) = state(i,j,k,URHO) * eos_state % xn

                state(i,j,k,UMX:UMZ) = ZERO

                state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * sum(state(i,j,k,UMX:UMZ)**2) / state(i,j,k,URHO)

                ! Convergence test

                ! Zones only participate in this test if they have a density
                ! that is above a certain fraction of the peak, to avoid
                ! oscillations in low density zones stalling convergence.

                drho = abs( state(i,j,k,URHO) - old_rho ) / old_rho

                if (state(i,j,k,URHO) / actual_rho_max > 1.0e-3_rt) then
                   Linf_norm = max(Linf_norm, drho)
                end if

             end if

          enddo
       enddo
    enddo

  end subroutine scf_update_density



  subroutine scf_diagnostics(lo, hi, &
                             state, s_lo, s_hi, &
                             phi, phi_lo, phi_hi, &
                             phi_rot, phr_lo, phr_hi, &
                             dx, &
                             kin_eng, pot_eng, int_eng, mass) bind(C, name='scf_diagnostics')

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: NVAR, URHO, UTEMP, UFS
    use network, only: nspec
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: phr_lo(3), phr_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ) :: phi_rot(phr_lo(1):phr_hi(1),phr_lo(2):phr_hi(2),phr_lo(3):phr_hi(3))
    real(rt), intent(inout) :: kin_eng, pot_eng, int_eng, mass

    integer  :: i, j, k
    real(rt) :: dV, dm

    type (eos_t) :: eos_state

    dV = dx(1) * dx(2) * dx(3)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,URHO) > ZERO) then

                dm = state(i,j,k,URHO) * dV

                mass = mass + dm

                kin_eng = kin_eng + phi_rot(i,j,k) * dm

                pot_eng = pot_eng + HALF * phi(i,j,k) * dm

                eos_state%rho = state(i,j,k,URHO)
                eos_state%T   = state(i,j,k,UTEMP)
                eos_state%xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)
                eos_state%aux = state(i,j,k,UFX:UFX+naux-1) / state(i,j,k,URHO)

                call eos(eos_input_rt, eos_state)

                int_eng = int_eng + eos_state%p * dV

             end if

          enddo
       enddo
    enddo

  end subroutine scf_diagnostics



  subroutine scf_get_solar_mass(M_solar_out) bind(C, name="scf_get_solar_mass")

    use fundamental_constants_module, only: M_solar

    implicit none

    real(rt), intent(inout) :: M_solar_out

    M_solar_out = M_solar

  end subroutine scf_get_solar_mass

end module scf_relaxation_module
