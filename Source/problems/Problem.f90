! problem-specific Fortran stuff goes here


!> @brief called by the IO processor during checkpoint
!!
!! @note Binds to C function ``problem_checkpoint``
!!
subroutine problem_checkpoint(int_dir_name, len) bind(C, name="problem_checkpoint")


  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo



end subroutine problem_checkpoint



!> @brief called by ALL processors during restart
!!
!! @note Binds to C function ``problem_restart``
!!
subroutine problem_restart(int_dir_name, len) bind(C, name="problem_restart")


  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

end subroutine problem_restart
