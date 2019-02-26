module scf_relaxation_module

  implicit none

  ! Input parameters for SCF relaxation
  
  logical         , save :: do_scf_initial_models = .false.
  double precision, save :: scf_d_A = 1.0d9
  double precision, save :: scf_d_B = 1.0d9
  double precision, save :: scf_d_C = 1.8d9
  double precision, save :: scf_relax_tol = 1.d-3
  
  ! Internal data for SCF relaxation
  
  double precision, save :: scf_h_max_P, scf_h_max_S
  double precision, save :: scf_enthalpy_min
  double precision, save :: scf_rpos(3,3)         ! Relative position of points A, B, and C
  double precision, save :: scf_d_vector(3,3)     ! Positions of points relative to system center
  double precision, save :: scf_c(0:1,0:1,0:1,3)  ! Interpolation coefficients for points
  integer         , save :: scf_rloc(3,3)         ! Indices of zones nearby to these points

contains

  subroutine initialize_scf()

    implicit none

    call scf_read_namelist()

  end subroutine initialize_scf


  
  subroutine scf_read_namelist()

    use problem_io_module, only: probin
    
    implicit none

    namelist /scf/ do_scf_initial_models, &
                   scf_d_A, scf_d_B, scf_d_C, &
                   scf_relax_tol

    integer :: untin
    
    ! Read namelist to override the module defaults.

    untin = 9 
    open(untin,file=probin,form='formatted',status='old')
    read(untin,scf)
    close(unit=untin)    

  end subroutine scf_read_namelist

  
  
  subroutine scf_setup_relaxation(dx) bind(C)

    use bl_constants_module, only: HALF, ONE, TWO
    use prob_params_module, only: problo, center, probhi
    use eos_module
    use probdata_module

    implicit none

    double precision :: dx(3)

    double precision :: x, y, z
    integer          :: n
    integer          :: ncell(3)

    double precision :: pos_l(3)
    type (eos_t) :: eos_state, ambient_state

    call get_ambient(ambient_state)

    ! Following Swesty, Wang and Calder (2000), and Motl et al. (2002), 
    ! we need to fix three points to uniquely determine an equilibrium 
    ! configuration for two unequal mass stars. We can do this by 
    ! specifying scf_d_A, the distance from the center of mass to the 
    ! inner point of the primary; scf_d_B, the distance from the center 
    ! of mass to the inner point of the secondary; and scf_d_C, the 
    ! distance from the center of mass to the outer point of the secondary. 
    ! We give these widths in physical units, so in general these will 
    ! be locations on the grid that are not coincident with a corner. 
    ! If a location is at point (x, y, z), then we find the eight 
    ! zone centers that surround this point, and do a tri-linear
    ! reconstruction to estimate the relative weight of each of the eight 
    ! zone centers in determining the value at that point.

    ncell = NINT( (probhi - problo) / dx )

    do n = 1, 3
       scf_d_vector(:,n) = center
       if (n .eq. 1) then
          scf_d_vector(axis_1,n) = scf_d_vector(axis_1,n) - scf_d_A
       else if (n .eq. 2) then
          scf_d_vector(axis_1,n) = scf_d_vector(axis_1,n) + scf_d_B
       else if (n .eq. 3) then
          scf_d_vector(axis_1,n) = scf_d_vector(axis_1,n) + scf_d_C
       endif
    enddo

    ! Locate the zone centers that bracket each point at the 
    ! lower left corner. Note that the INT function rounds down,
    ! which is what we want here since we want the lower left.

    do n = 1, 3
       scf_rloc(:,n) = INT( (scf_d_vector(:,n) + dx(:)/TWO) / dx(:)) + ncell / 2 - 1
    enddo

    ! Obtain the location of these points relative to the cube surrounding them.
    ! The lower left corner is at (0,0,0) and the upper right corner is at (1,1,1).

    do n = 1, 3
       pos_l = (scf_rloc(:,n) - ncell / 2 + HALF) * dx
       scf_rpos(:,n) = (scf_d_vector(:,n) - pos_l) / dx(:)
    enddo

    ! Determine the tri-linear coefficients

    do n = 1, 3
       x = scf_rpos(1,n)
       y = scf_rpos(2,n)
       z = scf_rpos(3,n)
       scf_c(0,0,0,n) = (ONE - x) * (ONE - y) * (ONE - z)
       scf_c(1,0,0,n) = x         * (ONE - y) * (ONE - z)
       scf_c(0,1,0,n) = (ONE - x) * y         * (ONE - z)
       scf_c(1,1,0,n) = x         * y         * (ONE - z)
       scf_c(0,0,1,n) = (ONE - x) * (ONE - y) * z
       scf_c(1,0,1,n) = x         * (ONE - y) * z
       scf_c(0,1,1,n) = (ONE - x) * y         * z
       scf_c(1,1,1,n) = x         * y         * z
    enddo

    ! Convert the maximum densities into maximum enthalpies.

    eos_state % T   = model_P % central_temp
    eos_state % rho = model_P % central_density
    eos_state % xn  = model_P % core_comp

    call eos(eos_input_rt, eos_state)

    scf_h_max_P = eos_state % h

    eos_state % T   = model_S % central_temp
    eos_state % rho = model_S % central_density
    eos_state % xn  = model_S % core_comp

    call eos(eos_input_rt, eos_state)

    scf_h_max_S = eos_state % h

    ! Determine the lowest possible enthalpy that can be 
    ! obtained by the EOS; this is useful for protecting
    ! against trying to compute a corresponding temperature
    ! in zones where the enthalpy is just too low for convergence.

    eos_state = ambient_state

    scf_enthalpy_min = 1.0d100

    do while (eos_state % rho < 1.d11)

       eos_state % rho = eos_state % rho * 1.1

       call eos(eos_input_rt, eos_state)

       if (eos_state % h < scf_enthalpy_min) scf_enthalpy_min = eos_state % h

    enddo

  end subroutine scf_setup_relaxation



  subroutine scf_get_coeff_info(loc_A,loc_B,loc_C,c_A,c_B,c_C) bind(C)

    implicit none

    integer,          intent(inout) :: loc_A(3), loc_B(3), loc_C(3)
    double precision, intent(inout) :: c_A(2,2,2), c_B(2,2,2), c_C(2,2,2)

    loc_A(:) = scf_rloc(:,1)
    loc_B(:) = scf_rloc(:,2)
    loc_C(:) = scf_rloc(:,3)

    c_A(:,:,:) = scf_c(:,:,:,1)
    c_B(:,:,:) = scf_c(:,:,:,2)
    c_C(:,:,:) = scf_c(:,:,:,3)

  end subroutine scf_get_coeff_info



  subroutine scf_get_omegasq(lo,hi,domlo,domhi, &
                             state,s_lo,s_hi, &
                             phi,p_lo,p_hi, &
                             dx,time,omegasq) bind(C)

    use bl_constants_module, only: ONE, TWO
    use meth_params_module, only: NVAR
    use rotation_frequency_module, only: get_omega
    use prob_params_module, only: problo, probhi

    implicit none

    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: s_lo(3), s_hi(3)
    integer          :: p_lo(3), p_hi(3)
    double precision :: dx(3), time
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision :: omegasq

    integer          :: i, j, k
    double precision :: theta2, theta3, omega(3)

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))    

    omega = get_omega(time)

    ! To ensure that this process is fully general for any orientation of
    ! the rotation axis with respect to the stellar configuration, 
    ! we define the angle between the rotation axis and the position vector,
    ! using the standard formula:
    ! cos(theta) = a . b / ( |a| |b| )
    ! The rotational potential is then given by:
    ! (1/2) * | omega x r |^2 = (1/2) |omega|^2 |r|^2 sin^2(theta)

    theta2 = acos( dot_product(omega,scf_d_vector(:,2)) / (sqrt(sum(omega**2)) * sqrt(sum(scf_d_vector(:,2)**2))) )
    theta3 = acos( dot_product(omega,scf_d_vector(:,3)) / (sqrt(sum(omega**2)) * sqrt(sum(scf_d_vector(:,3)**2))) )

    do k = scf_rloc(3,2), scf_rloc(3,2) + 1
       do j = scf_rloc(2,2), scf_rloc(2,2) + 1
          do i = scf_rloc(1,2), scf_rloc(1,2) + 1
             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                  i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then

                omegasq = omegasq - scf_c(i-scf_rloc(1,2),j-scf_rloc(2,2),k-scf_rloc(3,2),2) &
                     * TWO * phi(i,j,k) / &
                     (sin(theta3)**2 * sum(scf_d_vector(:,3)**2) - sin(theta2)**2 * sum(scf_d_vector(:,2)**2))
             endif
          enddo
       enddo
    enddo

    do k = scf_rloc(3,3), scf_rloc(3,3) + 1
       do j = scf_rloc(2,3), scf_rloc(2,3) + 1
          do i = scf_rloc(1,3), scf_rloc(1,3) + 1
             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                  i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then

                omegasq = omegasq + scf_c(i-scf_rloc(1,3),j-scf_rloc(2,3),k-scf_rloc(3,3),3) &
                     * TWO * phi(i,j,k) / &
                     (sin(theta3)**2 * sum(scf_d_vector(:,3)**2) - sin(theta2)**2 * sum(scf_d_vector(:,2)**2))

             endif
          enddo
       enddo
    enddo

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

  end subroutine scf_get_omegasq



  subroutine scf_get_bernoulli_const(lo,hi,domlo,domhi, &
                                     state,s_lo,s_hi, &
                                     phi,p_lo,p_hi, &
                                     dx,time,bernoulli_1,bernoulli_2) bind(C)

    use bl_constants_module, only: HALF, ONE, TWO, M_PI
    use meth_params_module, only: NVAR
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product

    implicit none

    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: s_lo(3), s_hi(3)
    integer          :: p_lo(3), p_hi(3)
    double precision :: dx(3), time
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision :: bernoulli_1, bernoulli_2

    integer          :: i, j, k
    double precision :: omega(3)

    omega = get_omega(time)

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    do k = scf_rloc(3,1), scf_rloc(3,1)+1
       do j = scf_rloc(2,1), scf_rloc(2,1)+1
          do i = scf_rloc(1,1), scf_rloc(1,1)+1
             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                  i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then
                bernoulli_1 = bernoulli_1 + scf_c(i-scf_rloc(1,1),j-scf_rloc(2,1),k-scf_rloc(3,1),1) &
                     * (phi(i,j,k) - HALF * sum(cross_product(omega, scf_d_vector(:,1))**2))
             endif
          enddo
       enddo
    enddo

    do k = scf_rloc(3,2), scf_rloc(3,2) + 1
       do j = scf_rloc(2,2), scf_rloc(2,2) + 1
          do i = scf_rloc(1,2), scf_rloc(1,2) + 1
             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                  i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then
                bernoulli_2 = bernoulli_2 + scf_c(i-scf_rloc(1,2),j-scf_rloc(2,2),k-scf_rloc(3,2),2) &
                     * (phi(i,j,k) - HALF * sum(cross_product(omega, scf_d_vector(:,2))**2))
             endif
          enddo
       enddo
    enddo

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

  end subroutine scf_get_bernoulli_const



  subroutine scf_construct_enthalpy(lo,hi,domlo,domhi, &
                                    state,s_lo,s_hi, &
                                    phi,p_lo,p_hi, &
                                    enthalpy,h_lo,h_hi, &
                                    dx,time, &
                                    bernoulli_1,bernoulli_2,h_max_1,h_max_2) bind(C)

    use bl_constants_module, only: ZERO, HALF, ONE, TWO, M_PI
    use meth_params_module, only: NVAR
    use prob_params_module, only: problo, center, probhi
    use probdata_module
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product

    implicit none

    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: s_lo(3), s_hi(3)
    integer          :: p_lo(3), p_hi(3)
    integer          :: h_lo(3), h_hi(3)
    double precision :: dx(3), time
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision :: enthalpy(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3))
    double precision :: bernoulli_1, bernoulli_2, h_max_1,h_max_2

    integer          :: i, j, k
    double precision :: r(3), omega(3), max_dist

    omega = get_omega(time)

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    ! The Bernoulli equation says that energy is conserved:
    ! enthalpy + gravitational potential + rotational potential = const
    ! The rotational potential is equal to -1/2 | omega x r |^2.
    ! We already have the constant, so our goal is to construct the enthalpy field.

    ! We don't want to check regions that aren't part of the stars.

    max_dist = 0.75 * max(maxval(abs(probhi-center)), maxval(abs(problo-center)))

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             if (r(axis_1) < ZERO .and. sum(r**2) .lt. max_dist**2) then

                enthalpy(i,j,k) = bernoulli_1 - phi(i,j,k) + HALF * sum(cross_product(omega, r)**2)

                if (enthalpy(i,j,k) > h_max_1) then
                   h_max_1 = enthalpy(i,j,k)
                endif

             else if (r(axis_1) > ZERO .and. sum(r**2) .lt. max_dist**2) then

                enthalpy(i,j,k) = bernoulli_2 - phi(i,j,k) + HALF * sum(cross_product(omega, r)**2)

                if (enthalpy(i,j,k) > h_max_2) then
                   h_max_2 = enthalpy(i,j,k)
                endif

             endif

          enddo
       enddo
    enddo

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

  end subroutine scf_construct_enthalpy



  subroutine scf_update_density(lo,hi,domlo,domhi, &
                                state,s_lo,s_hi, &
                                phi,p_lo,p_hi, &
                                enthalpy,h_lo,h_hi, &
                                dx,time, &
                                h_max_1,h_max_2, &
                                kin_eng, pot_eng, int_eng, &
                                left_mass, right_mass, &
                                delta_rho, l2_norm_resid, l2_norm_source) bind(C)

    use bl_constants_module, only: ZERO, ONE, TWO, M_PI
    use meth_params_module, only: NVAR, URHO, UTEMP, UMX, UMY, UMZ, UEDEN, UEINT, UFS
    use network, only: nspec
    use prob_params_module, only: problo, center, probhi
    use probdata_module, only: axis_1, get_ambient
    use eos_module
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product

    implicit none

    integer :: lo(3), hi(3), domlo(3), domhi(3)
    integer :: s_lo(3), s_hi(3)
    integer :: p_lo(3), p_hi(3)
    integer :: h_lo(3), h_hi(3)
    double precision :: dx(3), time
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    double precision :: enthalpy(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3))
    double precision :: h_max_1,h_max_2
    double precision :: kin_eng, pot_eng, int_eng
    double precision :: left_mass, right_mass
    double precision :: delta_rho, l2_norm_resid, l2_norm_source

    integer :: i, j, k
    double precision :: r(3)
    double precision :: old_rho, drho
    double precision :: dV
    double precision :: omega(3)
    double precision :: max_dist

    type (eos_t) :: eos_state, ambient_state

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    dV = dx(1) * dx(2) * dx(3)

    max_dist = 0.75 * max(maxval(abs(probhi-center)), maxval(abs(problo-center)))

    omega = get_omega(time)

    call get_ambient(ambient_state)

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             old_rho = state(i,j,k,URHO)

             ! Rescale the enthalpy by the maximum value.

             if (r(axis_1) .lt. ZERO) then
                enthalpy(i,j,k) = scf_h_max_P * (enthalpy(i,j,k) / h_max_1)
             else
                enthalpy(i,j,k) = scf_h_max_S * (enthalpy(i,j,k) / h_max_2)
             endif

             ! We only want to call the EOS for zones with enthalpy > 0,
             ! but for distances far enough from the center, the rotation
             ! term can overcome the other terms and make the enthalpy 
             ! spuriously positive. So we'll only consider zones within 75%
             ! of the distance from the center.

             if (enthalpy(i,j,k) > scf_enthalpy_min .and. sum(r**2) .lt. max_dist**2) then

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

             if (r(axis_1) < ZERO) then
                left_mass = left_mass + state(i,j,k,URHO) * dV
             else
                right_mass = right_mass + state(i,j,k,URHO) * dV
             endif

          enddo
       enddo
    enddo

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

  end subroutine scf_update_density

  

  subroutine scf_check_convergence(kin_eng, pot_eng, int_eng, &
                                   left_mass, right_mass, &
                                   delta_rho, l2_norm, &
                                   is_relaxed, num_iterations) bind(C)

    use problem_io_module, only: ioproc
    use meth_params_module, only: rot_period
    use multifab_module
    use bl_constants_module, only: THREE
    use fundamental_constants_module, only: M_solar
    use eos_module

    implicit none

    integer :: is_relaxed
    integer :: num_iterations

    double precision :: kin_eng, pot_eng, int_eng
    double precision :: left_mass, right_mass
    double precision :: delta_rho, l2_norm

    double precision :: virial_error

    virial_error = abs(TWO * kin_eng + pot_eng + THREE * int_eng) / abs(pot_eng)

    if (l2_norm .lt. scf_relax_tol) then
       is_relaxed = 1
    endif

    if (ioproc) then
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
       write(*,'(A,f5.3,A)')  "   Primary mass: ", left_mass / M_solar, " solar masses"
       write(*,'(A,f5.3,A)')  "   Secondary mass: ", right_mass / M_solar, " solar masses"
       if (is_relaxed .eq. 1) write(*,*) "  Relaxation completed!"
       write(*,*) ""
       write(*,*) ""
    endif

  end subroutine scf_check_convergence



  ! Returns whether we are doing a relaxation step.
  
  subroutine get_do_scf_initial_models(do_scf_initial_models_out) bind(C)

    implicit none

    integer :: do_scf_initial_models_out

    if (do_scf_initial_models) then
       do_scf_initial_models_out = 1
    else
       do_scf_initial_models_out = 0
    endif

  end subroutine get_do_scf_initial_models  
  
end module scf_relaxation_module
