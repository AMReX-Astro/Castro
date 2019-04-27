! Problem-specific Fortran routines that are designed to interact with C++

subroutine problem_checkpoint(int_dir_name, len) bind(C, name="problem_checkpoint")

  ! called by the IO processor during checkpoint

  use amrex_IO_module
  use amrex_error_module, only: amrex_error
  use probdata_module, only: com_P, com_S, vel_P, vel_S, mass_P, mass_S, t_ff_P, t_ff_S, &
                             T_global_max, rho_global_max, ts_te_global_max
  use prob_params_module, only: center
  use meth_params_module, only: rot_period
  use probdata_module, only: jobIsDone, signalJobIsNotDone, num_previous_ener_timesteps, total_ener_array, &
                             relaxation_is_done

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i, un

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

  un = unit_new()
  open (unit=un, file=trim(dir)//"/COM", status="unknown")

100 format(1x, g30.20)
200 format(1x, g30.20, 1x, g30.20)
300 format(1x, g30.20, 1x, g30.20, 1x, g30.20)

  write (un,300) center(1), center(2), center(3)
  write (un,200) mass_P,   mass_S
  write (un,200) com_P(1), com_S(1)
  write (un,200) com_P(2), com_S(2)
  write (un,200) com_P(3), com_S(3)
  write (un,200) vel_P(1), vel_S(1)
  write (un,200) vel_P(2), vel_S(2)
  write (un,200) vel_P(3), vel_S(3)
  write (un,200) t_ff_P,   t_ff_S

  close (un)



  open (unit=un, file=trim(dir)//"/Rotation", status="unknown")

  write (un,100) rot_period

  close (un)



  open (unit=un, file=trim(dir)//"/Extrema", status="unknown")

  write (un,100) T_global_max
  write (un,100) rho_global_max
  write (un,100) ts_te_global_max

  close (un)



  open (unit=un, file=trim(dir)//"/Energy", status="unknown")

  do i = 1, num_previous_ener_timesteps

     write (un,100) total_ener_array(i)

  enddo

  close (un)



  open (unit=un, file=trim(dir)//"/Relaxation", status="unknown")

  write (un,100) relaxation_is_done

  close (un)



  ! If the job is done, write a file in the checkpoint indicating this.

  if (jobIsDone) then

     open  (unit=un, file=trim(dir)//"/jobIsDone", status="unknown")
     close (un)

  else

     if (signalJobIsNotDone .and. .not. jobIsDone) then

        open  (unit=un, file=trim(dir)//"/jobIsNotDone", status="unknown")
        close (un)

     endif

  endif

end subroutine problem_checkpoint



subroutine problem_restart(int_dir_name, len) bind(C, name="problem_restart")

  ! called by ALL processors during restart 

  use amrex_IO_module
  use probdata_module, only: com_P, com_S, vel_P, vel_S, mass_P, mass_S, t_ff_P, t_ff_S, problem, &
                             T_global_max, rho_global_max, ts_te_global_max, &
                             jobIsDone, num_previous_ener_timesteps, total_ener_array, &
                             problem, relaxation_is_done
  use wdmerger_util_module, only: turn_off_relaxation
  use problem_io_module, only: ioproc
  use prob_params_module, only: center
  use meth_params_module, only: rot_period
  use amrex_error_module, only: amrex_error

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i, un, stat
  logical :: fileExists

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

  un = unit_new()
  open (unit=un, file=trim(dir)//"/COM", status="old")

100 format(1x, g30.20)
200 format(1x, g30.20, 1x, g30.20)
300 format(1x, g30.20, 1x, g30.20, 1x, g30.20)

  read (un,300) center(1), center(2), center(3)
  read (un,200) mass_P,   mass_S
  read (un,200) com_P(1), com_S(1)
  read (un,200) com_P(2), com_S(2)
  read (un,200) com_P(3), com_S(3)
  read (un,200) vel_P(1), vel_S(1)
  read (un,200) vel_P(2), vel_S(2)
  read (un,200) vel_P(3), vel_S(3)
  read (un,200) t_ff_P,   t_ff_S

  close (un)



  if (problem == 1) then

     open (unit=un, file=trim(dir)//"/Rotation", status="old", IOSTAT = stat)

     if (stat .eq. 0) then

        read (un,100) rot_period

        if (ioproc) then
           print *, ""
           write(*,1001) rot_period
           print *, ""

1001       format ("  Based on the checkpoint, updating the rotational period to ", F7.3 " s.")
        endif

        close (un)

     else

        rot_period = -1.0d0

     endif

  endif



  open (unit=un, file=trim(dir)//"/Extrema", status="old", IOSTAT = stat)

  if (stat .eq. 0) then

     read (un,100) T_global_max
     read (un,100) rho_global_max
     read (un,100) ts_te_global_max

     close(un)

  else

     T_global_max = 0.0d0
     rho_global_max = 0.0d0
     ts_te_global_max = 0.0d0

  endif



  open (unit=un, file=trim(dir)//"/Energy", status="old", IOSTAT = stat)

  if (stat .eq. 0) then

     do i = 1, num_previous_ener_timesteps

        read (un,100) total_ener_array(i)

     enddo

     close (un)

  else

     total_ener_array(:) = -1.d200

  endif



  open (unit=un, file=trim(dir)//"/Relaxation", status="old", IOSTAT = stat)

  if (stat .eq. 0) then

     read (un,100) relaxation_is_done

     close (un)

  else

     if (problem == 1) then
        call amrex_error("Error: no Relaxation file found in the checkpoint.")
     endif

  endif

  if (relaxation_is_done == 1) then

     call turn_off_relaxation(-1.0d0)

  endif



  ! Pick up whether the job has been completed.

  inquire(file=trim(dir)//"/jobIsDone", EXIST = fileExists)

  if (fileExists) then

     jobIsDone = .true.

  else

     inquire(file=trim(dir)//"/jobIsNotDone", EXIST = fileExists)

     if (fileExists) then

        jobIsDone = .false.

     endif

  endif

end subroutine problem_restart



! Determine the critical Roche potential at the Lagrange point L1.
! We will use a tri-linear interpolation that gets a contribution
! from all the zone centers that bracket the Lagrange point.

subroutine get_critical_roche_potential(phiEff,p_lo,p_hi,lo,hi,L1,potential) &
                                        bind(C,name='get_critical_roche_potential')

  use amrex_constants_module, only: ZERO, HALF, ONE
  use castro_util_module, only: position
  use prob_params_module, only: dim, dx_level
  use amrinfo_module, only: amr_level

  implicit none

  integer          :: lo(3), hi(3)
  integer          :: p_lo(3), p_hi(3)
  double precision :: phiEff(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
  double precision :: L1(3), potential

  double precision :: r(3), dx(3)
  integer          :: i, j, k

  dx = dx_level(:,amr_level)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           r = position(i,j,k) - L1

           ! Scale r by dx (in dimensions we're actually simulating).

           r(1:dim) = r(1:dim) / dx(1:dim)
           r(dim+1:3) = ZERO

           ! We want a contribution from this zone if it is
           ! less than one zone width away from the Lagrange point.

           if (sum(r**2) < ONE) then

              potential = potential + product(ONE - abs(r)) * phiEff(i,j,k)

           endif

        enddo
     enddo
  enddo

end subroutine get_critical_roche_potential



! Given state data in the rotating frame, transform it to the inertial frame.

subroutine transform_to_inertial_frame(state, s_lo, s_hi, lo, hi, time) &
                                       bind(C,name='transform_to_inertial_frame')

  use meth_params_module, only: NVAR, URHO, UMX, UMZ
  use wdmerger_util_module, only: inertial_velocity
  use castro_util_module, only: position

  implicit none

  integer          :: lo(3), hi(3)
  integer          :: s_lo(3), s_hi(3)
  double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
  double precision :: time

  double precision :: loc(3), vel(3)
  integer          :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           loc = position(i,j,k)
           vel = state(i,j,k,UMX:UMZ) / state(i,j,k,URHO)

           state(i,j,k,UMX:UMZ) = state(i,j,k,URHO) * inertial_velocity(loc, vel, time)

        enddo
     enddo
  enddo

end subroutine transform_to_inertial_frame



! Given the above quadrupole tensor, calculate the strain tensor.

subroutine gw_strain_tensor(h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3, Qtt, time) &
                            bind(C,name='gw_strain_tensor')

  use amrex_constants_module, only: ZERO, HALF, ONE, TWO
  use fundamental_constants_module, only: Gconst, c_light, parsec
  use probdata_module, only: gw_dist, axis_1, axis_2, axis_3
  use prob_params_module, only: dim

  implicit none

  double precision, intent(inout) :: h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3
  double precision, intent(in   ) :: Qtt(3,3)
  double precision, intent(in   ) :: time

  integer          :: i, j, k, l, dir
  double precision :: h(3,3), proj(3,3,3,3), delta(3,3), n(3), r
  double precision :: dist(3)

  ! Standard Kronecker delta.

  delta(:,:) = ZERO

  do i = 1, 3
     delta(i,i) = ONE
  enddo

  ! Unit vector for the wave is simply the distance 
  ! vector to the observer normalized by the total distance.
  ! We are going to repeat this process by looking along 
  ! all three coordinate axes.

  do dir = 1, 3

     dist(:) = ZERO
     dist(dir) = gw_dist

     r = sqrt(sum(dist**2))

     n(:) = dist(:) / r

     h = ZERO

     ! Projection operator onto the unit vector n.

     do l = 1, 3
        do k = 1, 3
           do j = 1, 3
              do i = 1, 3
                 proj(i,j,k,l) = (delta(i,k) - n(i) * n(k)) * (delta(j,l) - n(j) * n(l)) &
                               - HALF * (delta(i,j) - n(i) * n(j)) * (delta(k,l) - n(k) * n(l))
              enddo
           enddo
        enddo
     enddo

     ! Now we can calculate the strain tensor.

     do l = 1, 3
        do k = 1, 3
           do j = 1, 3
              do i = 1, 3
                 h(i,j) = h(i,j) + proj(i,j,k,l) * Qtt(k,l)
              enddo
           enddo
        enddo
     enddo

     ! Finally multiply by the coefficients.

     r = r * parsec * 1d3 ! Convert from kpc to cm.

     h(:,:) = h(:,:) * TWO * Gconst / (c_light**4 * r)

     if (dim .eq. 3) then

        ! If rot_axis == 3, then h_+ = h_{11} = -h_{22} and h_x = h_{12} = h_{21}.
        ! Analogous statements hold along the other axes.

        ! We are adding here so that this calculation makes sense on multiple levels.

        if (dir .eq. axis_1) then

           h_plus_1  = h_plus_1  + h(axis_2,axis_2)
           h_cross_1 = h_cross_1 + h(axis_2,axis_3)

        else if (dir .eq. axis_2) then

           h_plus_2  = h_plus_2  + h(axis_3,axis_3)
           h_cross_2 = h_cross_2 + h(axis_3,axis_1)

        else if (dir .eq. axis_3) then

           h_plus_3  = h_plus_3  + h(axis_1,axis_1)
           h_cross_3 = h_cross_3 + h(axis_1,axis_2)

        endif

     else

        ! In 2D axisymmetric coordinates, enforce that axis_1 is the x-axis,
        ! axis_2 is the y-axis, and axis_3 is the z-axis.

        if (dir .eq. 1) then

           h_plus_1  = h_plus_1  + h(2,2)
           h_cross_1 = h_cross_1 + h(2,3)

        else if (dir .eq. 2) then

           h_plus_2  = h_plus_2  + h(3,3)
           h_cross_2 = h_cross_2 + h(3,1)

        else if (dir .eq. 3) then

           h_plus_3  = h_plus_3  + h(1,1)
           h_cross_3 = h_cross_3 + h(1,2)

        endif

     endif

  enddo

end subroutine gw_strain_tensor



subroutine update_center(time) bind(C,name='update_center')

  use amrex_constants_module, only: ZERO
  use amrex_error_module, only: amrex_error
  use probdata_module, only: bulk_velx, bulk_vely, bulk_velz, &
                             center_fracx, center_fracy, center_fracz
  use prob_params_module, only: center, problo, probhi, dim

  implicit none

  double precision, intent(in) :: time

  ! Determine the original location of the center.

  if (dim .eq. 3) then

     center(1) = problo(1) + center_fracx * (probhi(1) - problo(1))
     center(2) = problo(2) + center_fracy * (probhi(2) - problo(2))
     center(3) = problo(3) + center_fracz * (probhi(3) - problo(3))

  else if (dim .eq. 2) then

     center(1) = problo(1)
     center(2) = problo(2) + center_fracz * (probhi(2) - problo(2))
     center(3) = ZERO

  else if (dim .eq. 1) then

     center(1) = problo(1) + center_fracx * (probhi(1) - problo(1))
     center(2) = ZERO
     center(3) = ZERO

  else

     call amrex_error("Error: unknown dim in subroutine update_center.")

  endif

  ! Now update using the time passed since the beginning of the simulation.

  center(1) = center(1) + bulk_velx * time
  center(2) = center(2) + bulk_vely * time
  center(3) = center(3) + bulk_velz * time

end subroutine update_center



! Updates the CASTRO rotational period.

subroutine set_period(period) bind(C,name='set_period')

  use meth_params_module, only: rot_period

  implicit none

  double precision :: period

  rot_period = period

end subroutine set_period



! Returns the CASTRO rotational period.

subroutine get_period(period) bind(C,name='get_period')

  use meth_params_module, only: rot_period

  implicit none

  double precision :: period

  period = rot_period

end subroutine get_period



! Returns the CASTRO rotation frequency vector.

subroutine get_omega_vec(omega_in, time) bind(C,name='get_omega_vec')

  use rotation_frequency_module, only: get_omega

  implicit none

  double precision, intent(inout) :: omega_in(3)
  double precision, intent(in   ), value :: time

  omega_in = get_omega(time)

end subroutine get_omega_vec



! Updates the global extrema.

subroutine set_extrema(T_max, rho_max, ts_te_max) bind(C,name='set_extrema')

  use probdata_module, only: T_global_max, rho_global_max, ts_te_global_max

  implicit none

  double precision, intent(in) :: T_max, rho_max, ts_te_max

  T_global_max     = T_max
  rho_global_max   = rho_max
  ts_te_global_max = ts_te_max

end subroutine set_extrema



! Retrieves the global extrema.

subroutine get_extrema(T_max, rho_max, ts_te_max) bind(C,name='get_extrema')

  use probdata_module, only: T_global_max, rho_global_max, ts_te_global_max

  implicit none

  double precision, intent(inout) :: T_max, rho_max, ts_te_max

  T_max     = T_global_max
  rho_max   = rho_global_max
  ts_te_max = ts_te_global_max

end subroutine get_extrema



! Returns whether the simulation is done.

subroutine get_job_status(jobDoneStatus) bind(C,name='get_job_status')

  use probdata_module, only: jobIsDone

  implicit none

  integer, intent(inout) :: jobDoneStatus

  if (jobIsDone) then
     jobDoneStatus = 1
  else
     jobDoneStatus = 0
  endif

end subroutine get_job_status



! Sets whether the simulation is done.

subroutine set_job_status(jobDoneStatus) bind(C,name='set_job_status')

  use probdata_module, only: jobIsDone

  implicit none

  integer, intent(in) :: jobDoneStatus

  if (jobDoneStatus == 1) then
     jobIsDone = .true.
  else
     jobIsDone = .false.
  endif

end subroutine set_job_status



! Get the relaxation_cutoff_time parameter.

subroutine get_relaxation_cutoff_time(relaxation_cutoff_time_in) bind(C,name='get_relaxation_cutoff_time')

  use amrex_fort_module, only: rt => amrex_real
  use probdata_module, only: relaxation_cutoff_time

  implicit none

  real(rt), intent(inout) :: relaxation_cutoff_time_in

  relaxation_cutoff_time_in = relaxation_cutoff_time

end subroutine get_relaxation_cutoff_time



! Gets whether the relaxation is done.

subroutine get_relaxation_status(relaxation_status) bind(C,name='get_relaxation_status')

  use probdata_module, only: relaxation_is_done

  implicit none

  integer, intent(inout) :: relaxation_status

  relaxation_status = relaxation_is_done

end subroutine get_relaxation_status



! Sets whether the relaxation is done.

subroutine set_relaxation_status(relaxation_status) bind(C,name='set_relaxation_status')

  use probdata_module, only: relaxation_is_done

  implicit none

  integer, intent(in) :: relaxation_status

  relaxation_is_done = relaxation_status

end subroutine set_relaxation_status



! Retrieve the total energy array.

subroutine get_total_ener_array(ener_array_in) bind(C,name='get_total_ener_array')

  use probdata_module, only: num_previous_ener_timesteps, total_ener_array

  implicit none

  double precision, intent(inout) :: ener_array_in(num_previous_ener_timesteps)

  ener_array_in(:) = total_ener_array(:)

end subroutine get_total_ener_array



! Set the total energy array.

subroutine set_total_ener_array(ener_array_in) bind(C,name='set_total_ener_array')

  use probdata_module, only: num_previous_ener_timesteps, total_ener_array

  implicit none

  double precision, intent(in) :: ener_array_in(num_previous_ener_timesteps)

  total_ener_array(:) = ener_array_in(:)

end subroutine set_total_ener_array
