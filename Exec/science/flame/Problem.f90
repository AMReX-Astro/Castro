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


! Calculate temperature properties, used for determining flame width.

subroutine flame_width_temp(temp, t_lo, t_hi, &
                            lo, hi, dx, time, &
                            T_max, T_min, grad_T_max) bind(C)

  use bl_constants_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer         , intent(in   ) :: t_lo(3), t_hi(3)

  real(rt)        , intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))

  integer         , intent(in   ) :: lo(3), hi(3)
  real(rt)        , intent(in   ) :: dx(3), time

  real(rt)        , intent(inout) :: T_max, T_min, grad_T_max

  ! Local variables
  
  integer :: i, j, k
  real(rt)         :: T, grad_T
  
  ! Assumes 1D simulation, right now. Also assumes that
  ! we have at least one ghost cell in the x dimension.
  
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           T = temp(i,j,k)
           grad_T = abs(temp(i+1,j,k) - temp(i-1,j,k)) / (TWO * dx(1))

           ! Ignore problem zones where we have a negative temperature

           if (T .lt. ZERO) continue

           if (T > T_max) then
              T_max = T
           endif

           if (T < T_min) then
              T_min = T
           endif

           if (grad_T > grad_T_max) then
              grad_T_max = grad_T
           endif
           
        enddo
     enddo
  enddo

end subroutine flame_width_temp



subroutine flame_speed_data(omegadot, od_lo, od_hi, &
                            lo, hi, dx, &
                            rho_X_dot) bind(C)

  use bl_constants_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer         , intent(in   ) :: od_lo(3), od_hi(3)

  real(rt)        , intent(in   ) :: omegadot(od_lo(1):od_hi(1),od_lo(2):od_hi(2),od_lo(3):od_hi(3))

  integer         , intent(in   ) :: lo(3), hi(3)
  real(rt)        , intent(in   ) :: dx(3)

  real(rt)        , intent(inout) :: rho_X_dot

  ! Local variables
  
  integer :: i, j, k
  
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           rho_X_dot = rho_X_dot + omegadot(i,j,k) * dx(1)
           
        enddo
     enddo
  enddo

end subroutine flame_speed_data
