! problem-specific Fortran stuff goes here

subroutine problem_checkpoint(int_dir_name, len) bind(C, name="problem_checkpoint")

  ! called by the IO processor during checkpoint

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


subroutine problem_restart(int_dir_name, len) bind(C, name="problem_restart")

  ! called by ALL processors during restart 

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



subroutine check_stopping_criteria(lo, hi, &
                                   xvel, v_lo, v_hi, &
                                   temp, t_lo, t_hi, &
                                   ts_te, s_lo, s_hi, &
                                   T_criterion, ts_te_criterion, to_stop) bind(c, name="check_stopping_criteria")

  use amrex_constants_module, only: ZERO
  use amrex_fort_Module, only: rt => amrex_real
  use meth_params_module, only: NVAR, UMX, UTEMP
  use probdata_module, only: vel

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: v_lo(3), v_hi(3)
  integer,  intent(in   ) :: t_lo(3), t_hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)

  real(rt), intent(in   ) :: xvel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
  real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
  real(rt), intent(in   ) :: ts_te(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

  real(rt), intent(in   ), value :: T_criterion, ts_te_criterion
  integer,  intent(inout) :: to_stop

  integer :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           ! Note that this stopping criterion only makes sense for the collision problem.

           if (vel > ZERO .and. xvel(i,j,k) > ZERO) then
              if (temp(i,j,k) >= T_criterion .or. ts_te(i,j,k) >= ts_te_criterion) then
                 to_stop = 1
              end if
           end if

        end do
     end do
  end do

end subroutine check_stopping_criteria
