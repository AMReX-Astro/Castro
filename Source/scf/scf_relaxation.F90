module scf_relaxation_module

  use amrex_fort_module, only: rt => amrex_real
  use eos_type_module, only: eos_t

  implicit none

  ! Internal data for SCF relaxation
  
  real(rt), save :: scf_h_max
  real(rt), save :: scf_enthalpy_min
  real(rt), save :: scf_r_A(3), scf_r_B(3)       ! Position of points A and B relative to system center

  type (eos_t), save :: ambient_state

contains

  subroutine scf_setup_relaxation() bind(C, name='scf_setup_relaxation')

    use amrex_constants_module, only: HALF, ONE, TWO
    use prob_params_module, only: problo, center, probhi
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t
    use meth_params_module, only: scf_maximum_density, scf_temperature, &
                                  scf_equatorial_radius, scf_polar_radius
    use network, only: nspec

    implicit none

    type (eos_t) :: eos_state

    ! We need to fix two points to uniquely determine an equilibrium 
    ! configuration for a rotating star. We can do this by specifying
    ! scf_equatorial_radius, the equatorial radius of the star, and
    ! scf_polar_radius, the polar radius of the star. In this configuration,
    ! we assume that the rotation axis is the z-axis, so the equatorial
    ! radius is along the xy-plane.

    scf_r_A(:) = center
    scf_r_B(:) = center

    ! Note that the choice of direction is arbitrary here: for the polar
    ! radius we could fix the point at positive z or negative z, and for
    ! the equatorial radius, we could use +/- x or y. Since the problem
    ! is symmetric, it should not matter what we actually choose.

    scf_r_A(1) = scf_r_A(1) + scf_equatorial_radius
    scf_r_B(3) = scf_r_B(3) + scf_polar_radius

    ! Convert the maximum density into a maximum enthalpy.

    eos_state % rho = scf_maximum_density
    eos_state % T   = scf_temperature
    eos_state % xn  = 1.0d0 / nspec

    call eos(eos_input_rt, eos_state)

    scf_h_max = eos_state % h

    ! Determine the lowest possible enthalpy that can be 
    ! obtained by the EOS; this is useful for protecting
    ! against trying to compute a corresponding temperature
    ! in zones where the enthalpy is just too low for convergence.

    ambient_state % T   = scf_temperature
    ambient_state % rho = 1.0d-4
    ambient_state % xn  = 1.0d0 / nspec

    call eos(eos_input_rt, ambient_state)

    eos_state = ambient_state

    scf_enthalpy_min = 1.0d100

    do while (eos_state % rho < 1.d11)

       eos_state % rho = eos_state % rho * 1.1

       call eos(eos_input_rt, eos_state)

       if (eos_state % h < scf_enthalpy_min) scf_enthalpy_min = eos_state % h

    enddo

  end subroutine scf_setup_relaxation



  ! Construct psi, the term that determines the rotation law.

  subroutine scf_construct_psi(lo, hi, &
                               psi, psi_lo, psi_hi, &
                               dx) bind(C, name='scf_construct_psi')

    use amrex_constants_module, only: HALF
    use meth_params_module, only: NVAR
    use prob_params_module, only: problo, center

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: psi_lo(3), psi_hi(3)
    real(rt), intent(inout) :: psi(psi_lo(1):psi_hi(1),psi_lo(2):psi_hi(2),psi_lo(3):psi_hi(3))
    real(rt), intent(in   ) :: dx(3)

    integer  :: i, j, k
    real(rt) :: r(3)

    ! The below assumes we are rotating on the z-axis.

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             psi(i,j,k) = -HALF * (r(1)**2 + r(2)**2)

          enddo
       enddo
    enddo

  end subroutine scf_construct_psi



  ! Calculate the phi and psi factors that go into
  ! updating the rotation frequency.

  subroutine scf_update_for_omegasq(lo, hi, &
                                    state, s_lo, s_hi, &
                                    phi, phi_lo, phi_hi, &
                                    psi, psi_lo, psi_hi, &
                                    dx, &
                                    phi_A, psi_A, phi_B, psi_B) bind(C, name='scf_update_for_omegasq')

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: NVAR
    use prob_params_module, only: problo, center

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: psi_lo(3), psi_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ) :: psi(psi_lo(1):psi_hi(1),psi_lo(2):psi_hi(2),psi_lo(3):psi_hi(3))
    real(rt), intent(inout) :: phi_A, phi_B, psi_A, psi_B

    integer  :: i, j, k
    integer  :: loc(3)
    real(rt) :: r(3), rr(3), scale

    ! The below assumes we are rotating on the z-axis.

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             ! Do a trilinear interpolation to find the contribution from
             ! this grid point. Limit so that only the nearest zone centers
             ! can participate. This implies that the maximum allowable
             ! distance from the target location is 0.5 * dx.

             rr = abs(r - scf_r_A) / dx
             rr = merge(rr, ZERO, rr <= HALF)

             scale = rr(1) * rr(2) * rr(3)

             phi_A = phi_A + scale * phi(i,j,k)
             psi_A = psi_A + scale * psi(i,j,k)

          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             rr = abs(r - scf_r_B) / dx
             rr = merge(rr, ZERO, rr <= HALF)

             scale = rr(1) * rr(2) * rr(3)

             phi_B = phi_B + scale * phi(i,j,k)
             psi_B = psi_B + scale * psi(i,j,k)

          enddo
       enddo
    enddo

  end subroutine scf_update_for_omegasq



  subroutine scf_get_bernoulli_const(lo, hi, &
                                     state, s_lo, s_hi, &
                                     phi, phi_lo, phi_hi, &
                                     psi, psi_lo, psi_hi, &
                                     dx, omega, bernoulli) bind(C, name='scf_get_bernoulli_const')

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: NVAR
    use prob_params_module, only: problo, center

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: psi_lo(3), psi_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ) :: psi(psi_lo(1):psi_hi(1),psi_lo(2):psi_hi(2),psi_lo(3):psi_hi(3))
    real(rt), intent(inout) :: bernoulli
    real(rt), intent(in   ), value :: omega

    integer  :: i, j, k
    integer  :: loc(3)
    real(rt) :: r(3), rr(3), scale

    ! The below assumes we are rotating on the z-axis.

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             rr = abs(r - scf_r_A) / dx
             rr = merge(rr, ZERO, rr <= HALF)

             scale = rr(1) * rr(2) * rr(3)

             bernoulli = bernoulli + scale * (phi(i,j,k) - omega**2 * psi(i,j,k))

          enddo
       enddo
    enddo

  end subroutine scf_get_bernoulli_const



  subroutine scf_construct_enthalpy(lo, hi, &
                                    state, s_lo, s_hi, &
                                    phi, phi_lo, phi_hi, &
                                    psi, psi_lo, psi_hi, &
                                    enthalpy, h_lo, h_hi, &
                                    dx, omega, bernoulli) bind(C, name='scf_construct_enthalpy')

    use amrex_constants_module, only: ZERO, HALF, ONE, TWO, M_PI
    use meth_params_module, only: NVAR
    use prob_params_module, only: problo, center, probhi

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: psi_lo(3), psi_hi(3)
    integer,  intent(in   ) :: h_lo(3), h_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ) :: psi(psi_lo(1):psi_hi(1),psi_lo(2):psi_hi(2),psi_lo(3):psi_hi(3))
    real(rt), intent(inout) :: enthalpy(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3))
    real(rt), intent(in   ), value :: omega, bernoulli

    integer  :: i, j, k
    real(rt) :: r(3)

    ! The Bernoulli equation says that energy is conserved:
    ! enthalpy + gravitational potential + rotational potential = const
    ! The rotational potential is equal to -1/2 | omega x r |^2.
    ! We already have the constant, so our goal is to construct the enthalpy field.

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             enthalpy(i,j,k) = bernoulli - phi(i,j,k) - omega**2 * psi(i,j,k)

          enddo
       enddo
    enddo

  end subroutine scf_construct_enthalpy



  subroutine scf_update_density(lo, hi, &
                                state, s_lo, s_hi, &
                                enthalpy, h_lo, h_hi, &
                                dx, omega, h_max, &
                                Linf_norm) bind(C, name='scf_update_density')

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: NVAR, URHO, UTEMP, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                  scf_temperature
    use network, only: nspec
    use prob_params_module, only: problo, center
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
    real(rt), intent(in   ), value :: omega, h_max

    integer  :: i, j, k
    real(rt) :: r(3)
    real(rt) :: old_rho, drho

    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             old_rho = state(i,j,k,URHO)

             ! Rescale the enthalpy by the maximum value.

             enthalpy(i,j,k) = scf_h_max * (enthalpy(i,j,k) / h_max)

             ! We only want to call the EOS for zones with enthalpy > 0,
             ! but for distances far enough from the center, the rotation
             ! term can overcome the other terms and make the enthalpy 
             ! spuriously negative. Avoid this by checking against an
             ! enthalpy floor.

             if (enthalpy(i,j,k) > scf_enthalpy_min) then

                eos_state%rho = state(i,j,k,URHO) ! Initial guess for the EOS
                eos_state%T   = scf_temperature
                eos_state%h   = enthalpy(i,j,k)
                eos_state%xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)

                call eos(eos_input_th, eos_state)

             else

                eos_state = ambient_state

             endif

             state(i,j,k,URHO)  = eos_state % rho
             state(i,j,k,UTEMP) = eos_state % T
             state(i,j,k,UEINT) = eos_state % rho * eos_state % e
             state(i,j,k,UFS:UFS+nspec-1) = eos_state % rho * eos_state % xn

             state(i,j,k,UMX:UMZ) = ZERO

             state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * sum(state(i,j,k,UMX:UMZ)**2) / state(i,j,k,URHO)

             ! Convergence test

             drho = abs( state(i,j,k,URHO) - old_rho ) / old_rho
             Linf_norm = max(Linf_norm, drho)

          enddo
       enddo
    enddo

  end subroutine scf_update_density



  subroutine scf_diagnostics(lo, hi, &
                             state, s_lo, s_hi, &
                             phi, phi_lo, phi_hi, &
                             psi, psi_lo, psi_hi, &
                             dx, omega, &
                             kin_eng, pot_eng, int_eng, mass) bind(C, name='scf_diagnostics')

    use amrex_constants_module, only: HALF
    use meth_params_module, only: NVAR, URHO, UTEMP, UFS
    use network, only: nspec
    use prob_params_module, only: problo, center
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: psi_lo(3), psi_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ) :: psi(psi_lo(1):psi_hi(1),psi_lo(2):psi_hi(2),psi_lo(3):psi_hi(3))
    real(rt), intent(inout) :: kin_eng, pot_eng, int_eng, mass
    real(rt), intent(in   ), value :: omega

    integer  :: i, j, k
    real(rt) :: r(3), dV, dm

    type (eos_t) :: eos_state

    dV = dx(1) * dx(2) * dx(3)

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             dm = state(i,j,k,URHO) * dV

             mass = mass + dm

             kin_eng = kin_eng - omega**2 * psi(i,j,k) * dm

             pot_eng = pot_eng + HALF * phi(i,j,k) * dm

             eos_state%rho = state(i,j,k,URHO)
             eos_state%T   = state(i,j,k,UTEMP)
             eos_state%xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)

             call eos(eos_input_rt, eos_state)

             int_eng = int_eng + eos_state%p * dV

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
