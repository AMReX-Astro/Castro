module scf_relaxation_module

  use amrex_fort_module, only: rt => amrex_real
  use eos_type_module, only: eos_t

  implicit none

  ! Internal data for SCF relaxation
  
  real(rt), save :: scf_h_max
  real(rt), save :: scf_enthalpy_min
  integer,  save :: scf_rloc_A(3), scf_rloc_B(3) ! Indices of zones nearby to these points
  real(rt), save :: scf_rpos_A(3), scf_rpos_B(3) ! Position of points A and B relative to lower left zone corner
  real(rt), save :: scf_r_A(3), scf_r_B(3)       ! Position of points A and B relative to system center
  real(rt), save :: scf_c_A(0:1,0:1,0:1), scf_c_B(0:1,0:1,0:1)  ! Interpolation coefficients for points

  type (eos_t), save :: ambient_state

contains

  subroutine scf_setup_relaxation(dx) bind(C)

    use amrex_constants_module, only: HALF, ONE, TWO
    use prob_params_module, only: problo, center, probhi
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t
    use meth_params_module, only: scf_d_A, scf_d_B, scf_relax_tol
    use network, only: nspec

    implicit none

    real(rt) :: dx(3)

    real(rt) :: x, y, z
    integer  :: n
    integer  :: ncell(3)

    real(rt) :: pos_l(3)
    type (eos_t) :: eos_state

    ! We need to fix two points to uniquely determine an equilibrium 
    ! configuration for a rotating star. We can do this by specifying
    ! scf_d_A, the equatorial radius of the star, and scf_d_B, the
    ! polar radius of the star. We give these widths in physical units,
    ! so in general these will be locations on the grid that are not
    ! coincident with a corner.  If a location is at point (x, y, z),
    ! then we find the eight  zone centers that surround this point,
    ! and do a tri-linear reconstruction to estimate the relative
    ! weight of each of the eight zone centers in determining the
    ! value at that point.

    ncell = NINT( (probhi - problo) / dx )

    scf_r_A(:) = center
    scf_r_B(:) = center

    scf_r_A(1) = scf_r_A(1) + scf_d_A
    scf_r_B(3) = scf_r_B(3) + scf_d_B

    ! Locate the zone centers that bracket each point at the 
    ! lower left corner. Note that the INT function rounds down,
    ! which is what we want here since we want the lower left.

    scf_rloc_A(:) = INT( (scf_r_A(:) + dx(:) / TWO) / dx(:)) + ncell / 2 - 1
    scf_rloc_B(:) = INT( (scf_r_B(:) + dx(:) / TWO) / dx(:)) + ncell / 2 - 1

    ! Obtain the location of these points relative to the cube surrounding them.
    ! The lower left corner is at (0,0,0) and the upper right corner is at (1,1,1).

    pos_l = (scf_rloc_A(:) - ncell / 2 + HALF) * dx(:)
    scf_rpos_A(:) = (scf_r_A(:) - pos_l) / dx(:)

    pos_l = (scf_rloc_B(:) - ncell / 2 + HALF) * dx(:)
    scf_rpos_B(:) = (scf_r_B(:) - pos_l) / dx(:)

    ! Determine the tri-linear coefficients

    x = scf_rpos_A(1)
    y = scf_rpos_A(2)
    z = scf_rpos_A(3)
    scf_c_A(0,0,0) = (ONE - x) * (ONE - y) * (ONE - z)
    scf_c_A(1,0,0) = x         * (ONE - y) * (ONE - z)
    scf_c_A(0,1,0) = (ONE - x) * y         * (ONE - z)
    scf_c_A(1,1,0) = x         * y         * (ONE - z)
    scf_c_A(0,0,1) = (ONE - x) * (ONE - y) * z
    scf_c_A(1,0,1) = x         * (ONE - y) * z
    scf_c_A(0,1,1) = (ONE - x) * y         * z
    scf_c_A(1,1,1) = x         * y         * z

    x = scf_rpos_B(1)
    y = scf_rpos_B(2)
    z = scf_rpos_B(3)
    scf_c_B(0,0,0) = (ONE - x) * (ONE - y) * (ONE - z)
    scf_c_B(1,0,0) = x         * (ONE - y) * (ONE - z)
    scf_c_B(0,1,0) = (ONE - x) * y         * (ONE - z)
    scf_c_B(1,1,0) = x         * y         * (ONE - z)
    scf_c_B(0,0,1) = (ONE - x) * (ONE - y) * z
    scf_c_B(1,0,1) = x         * (ONE - y) * z
    scf_c_B(0,1,1) = (ONE - x) * y         * z
    scf_c_B(1,1,1) = x         * y         * z

    ! Convert the maximum densities into maximum enthalpies.

    eos_state % T   = 1.0d7
    eos_state % rho = 5.0d6
    eos_state % xn  = 1.0d0 / nspec

    call eos(eos_input_rt, eos_state)

    scf_h_max = eos_state % h

    ! Determine the lowest possible enthalpy that can be 
    ! obtained by the EOS; this is useful for protecting
    ! against trying to compute a corresponding temperature
    ! in zones where the enthalpy is just too low for convergence.

    ambient_state % T   = 1.0d7
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



  ! Calculate the phi and psi factors that go into
  ! updating the rotation frequency.
  
  subroutine scf_update_for_omegasq(lo, hi, &
                                    domlo, domhi, &
                                    state, s_lo, s_hi, &
                                    phi, p_lo, p_hi, &
                                    dx, &
                                    phi_A, psi_A, phi_B, psi_B) bind(C, name='scf_update_for_omegasq')

    use amrex_constants_module, only: HALF
    use meth_params_module, only: NVAR
    use castro_util_module, only: position
    use prob_params_module, only: problo, center

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: p_lo(3), p_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    real(rt), intent(inout) :: phi_A, phi_B, psi_A, psi_B

    integer  :: i, j, k
    integer  :: loc(3)
    real(rt) :: omega(3), c(0:1,0:1,0:1), r(3), scale

    ! The below assumes we are rotating on the z-axis.

    loc = scf_rloc_A
    c = scf_c_A

    do k = loc(3), loc(3) + 1
       do j = loc(2), loc(2) + 1
          do i = loc(1), loc(1) + 1

             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                 i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then

                r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)
                r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
                r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)

                scale = c(i-loc(1),j-loc(2),k-loc(3))

                phi_A = phi_A + scale * phi(i,j,k)
                psi_A = psi_A + scale * (-HALF * (r(1)**2 + r(2)**2))

             endif

          enddo
       enddo
    enddo

    loc = scf_rloc_B
    c = scf_c_B

    do k = loc(3), loc(3) + 1
       do j = loc(2), loc(2) + 1
          do i = loc(1), loc(1) + 1

             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                 i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then

                r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)
                r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
                r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)

                scale = c(i-loc(1),j-loc(2),k-loc(3))

                phi_B = phi_B + scale * phi(i,j,k)
                psi_B = psi_B + scale * (-HALF * (r(1)**2 + r(2)**2))

             endif

          enddo
       enddo
    enddo

  end subroutine scf_update_for_omegasq



  subroutine scf_get_bernoulli_const(lo, hi, &
                                     domlo, domhi, &
                                     state, s_lo, s_hi, &
                                     phi, p_lo, p_hi, &
                                     dx, time, bernoulli) bind(C, name='scf_get_bernoulli_const')

    use amrex_constants_module, only: HALF, ONE, TWO, M_PI
    use meth_params_module, only: NVAR
    use rotation_frequency_module, only: get_omega
    use castro_util_module, only: position

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: p_lo(3), p_hi(3)
    real(rt), intent(in   ) :: dx(3), time
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    real(rt), intent(inout) :: bernoulli

    integer  :: i, j, k
    integer  :: loc(3)
    real(rt) :: omega(3), c(0:1,0:1,0:1), r(3)

    omega = get_omega(time)

    ! The below assumes we are rotating on the z-axis.

    loc = scf_rloc_A
    c = scf_c_A

    do k = loc(3), loc(3) + 1
       do j = loc(2), loc(2) + 1
          do i = loc(1), loc(1) + 1

             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                 i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then

                r(:) = position(i,j,k)

                bernoulli = bernoulli + c(i-loc(1),j-loc(2),k-loc(3)) * (phi(i,j,k) + HALF * (r(1)**2 + r(2)**2) * omega(3)**2)

             endif

          enddo
       enddo
    enddo

  end subroutine scf_get_bernoulli_const



  subroutine scf_construct_enthalpy(lo, hi, &
                                    domlo, domhi, &
                                    state, s_lo, s_hi, &
                                    phi, p_lo, p_hi, &
                                    enthalpy, h_lo, h_hi, &
                                    dx, time, &
                                    bernoulli, h_max) bind(C, name='scf_construct_enthalpy')

    use amrex_constants_module, only: ZERO, HALF, ONE, TWO, M_PI
    use meth_params_module, only: NVAR
    use prob_params_module, only: problo, center, probhi
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: p_lo(3), p_hi(3)
    integer,  intent(in   ) :: h_lo(3), h_hi(3)
    real(rt), intent(in   ) :: dx(3), time
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    real(rt), intent(inout) :: enthalpy(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3))
    real(rt), intent(in   ) :: bernoulli
    real(rt), intent(inout) :: h_max

    integer  :: i, j, k
    real(rt) :: r(3), omega(3)

    omega = get_omega(time)

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

             enthalpy(i,j,k) = bernoulli - phi(i,j,k) + HALF * sum(cross_product(omega, r)**2)

             if (enthalpy(i,j,k) > h_max) then
                h_max = enthalpy(i,j,k)
             end if

          enddo
       enddo
    enddo

  end subroutine scf_construct_enthalpy



  subroutine scf_update_density(lo, hi, &
                                domlo, domhi, &
                                state, s_lo, s_hi, &
                                phi, p_lo, p_hi, &
                                enthalpy, h_lo, h_hi, &
                                dx, time, &
                                h_max, &
                                kin_eng, pot_eng, int_eng, &
                                mass, &
                                delta_rho, l2_norm_resid, l2_norm_source) bind(C, name='scf_update_density')

    use amrex_constants_module, only: ZERO, HALF, ONE, TWO, M_PI
    use meth_params_module, only: NVAR, URHO, UTEMP, UMX, UMY, UMZ, UEDEN, UEINT, UFS
    use network, only: nspec
    use prob_params_module, only: problo, center, probhi
    use eos_module, only: eos
    use eos_type_module, only: eos_input_th, eos_t
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: p_lo(3), p_hi(3)
    integer,  intent(in   ) :: h_lo(3), h_hi(3)
    real(rt), intent(in   ) :: dx(3), time
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    real(rt), intent(inout) :: enthalpy(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3))
    real(rt), intent(in   ) :: h_max
    real(rt), intent(inout) :: kin_eng, pot_eng, int_eng
    real(rt), intent(inout) :: mass
    real(rt), intent(inout) :: delta_rho, l2_norm_resid, l2_norm_source

    integer  :: i, j, k
    real(rt) :: r(3)
    real(rt) :: old_rho, drho
    real(rt) :: dV
    real(rt) :: omega(3)

    type (eos_t) :: eos_state

    dV = dx(1) * dx(2) * dx(3)

    omega = get_omega(time)

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
             ! spuriously positive. So we'll only consider zones within 75%
             ! of the distance from the center.

             if (enthalpy(i,j,k) > scf_enthalpy_min) then

                eos_state % T   = state(i,j,k,UTEMP)
                eos_state % h   = enthalpy(i,j,k)
                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)
                eos_state % rho = state(i,j,k,URHO) ! Initial guess for the EOS

                call eos(eos_input_th, eos_state)

             else

                eos_state = ambient_state

             endif

             state(i,j,k,URHO) = eos_state % rho
             state(i,j,k,UEINT) = eos_state % rho * eos_state % e
             state(i,j,k,UFS:UFS+nspec-1) = eos_state % rho * eos_state % xn

             state(i,j,k,UMX:UMZ) = ZERO

             state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * (state(i,j,k,UMX)**2 + &
                  state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) / state(i,j,k,URHO)

             ! Convergence tests and diagnostic quantities

             drho = abs( state(i,j,k,URHO) - old_rho ) / old_rho
             if (drho > delta_rho) delta_rho = drho

             l2_norm_resid = l2_norm_resid + dV * (state(i,j,k,URHO) - old_rho)**2
             l2_norm_source = l2_norm_source + dV * old_rho**2

             kin_eng = kin_eng + HALF * sum(cross_product(omega, r)**2) * state(i,j,k,URHO) * dV

             pot_eng = pot_eng + HALF * state(i,j,k,URHO) * phi(i,j,k) * dV

             int_eng = int_eng + eos_state % p * dV

             mass = mass + state(i,j,k,URHO) * dV

          enddo
       enddo
    enddo

  end subroutine scf_update_density

  

  subroutine scf_check_convergence(kin_eng, pot_eng, int_eng, &
                                   mass, delta_rho, l2_norm, &
                                   is_relaxed, num_iterations) bind(C, name='scf_check_convergence')

    use meth_params_module, only: rot_period
    use amrex_constants_module, only: TWO, THREE
    use fundamental_constants_module, only: M_solar
    use meth_params_module, only: scf_relax_tol
    use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor

    implicit none

    integer,  intent(inout) :: is_relaxed
    integer,  intent(in   ) :: num_iterations

    real(rt), intent(in) :: kin_eng, pot_eng, int_eng
    real(rt), intent(in) :: mass, delta_rho, l2_norm

    real(rt) :: virial_error

    virial_error = abs(TWO * kin_eng + pot_eng + THREE * int_eng) / abs(pot_eng)

    if (l2_norm .lt. scf_relax_tol) then
       is_relaxed = 1
    endif

    if (amrex_pd_ioprocessor()) then

       write(*,*) ""
       write(*,*) ""
       write(*,'(A,I2)')      "   Relaxation iterations completed: ", num_iterations
       write(*,'(A,ES8.2)')   "   Maximum change in rho (g cm**-3): ", delta_rho
       write(*,'(A,ES8.2)')   "   L2 Norm of Residual (relative to old state): ", l2_norm
       write(*,'(A,f6.2)')    "   Rotational period (s): ", rot_period
       write(*,'(A,ES8.2)')   "   Kinetic energy: ", kin_eng
       write(*,'(A,ES9.2)')   "   Potential energy: ", pot_eng
       write(*,'(A,ES8.2)')   "   Internal energy: ", int_eng
       write(*,'(A,ES9.3)')   "   Virial error: ", virial_error
       write(*,'(A,f5.3,A)')  "   Mass: ", mass / M_solar, " solar masses"
       if (is_relaxed .eq. 1) write(*,*) "  Relaxation completed!"
       write(*,*) ""
       write(*,*) ""

    end if

  end subroutine scf_check_convergence

end module scf_relaxation_module
