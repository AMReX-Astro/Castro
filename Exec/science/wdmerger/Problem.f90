! Problem-specific Fortran routines that are designed to interact with C++

subroutine problem_checkpoint(int_dir_name, len) bind(C, name="problem_checkpoint")

  ! called by the IO processor during checkpoint

  use amrex_IO_module
  use castro_error_module, only: castro_error
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
  use wdmerger_util_module, only: set_relaxation_damping_factor
  use problem_io_module, only: ioproc
  use prob_params_module, only: center
  use meth_params_module, only: rot_period
  use castro_error_module, only: castro_error

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
        call castro_error("Error: no Relaxation file found in the checkpoint.")
     endif

  endif

  if (relaxation_is_done == 1) then

     call set_relaxation_damping_factor(-1.0_rt)

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
