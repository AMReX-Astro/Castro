module mhd_state_module

  use meth_params_module, only: NVAR
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: UMAGX = NVAR+1
  integer, parameter :: UMAGY = NVAR+2
  integer, parameter :: UMAGZ = NVAR+3

contains

  pure function epsilon_ijk(i, j, k) result (val)

    integer, intent(in) :: i, j, k
    real(rt) :: val

    if (i == j .or. j == k .or. i == k) then
       val = 0.0_rt
    else if ((i == 1 .and. j == 2 .and. k == 3) .or. &
             (i == 2 .and. j == 3 .and. k == 1) .or. &
             (i == 3 .and. j == 1 .and. k == 2)) then
       val = 1.0_rt
    else
       val = -1.0_rt
    end if

  end function epsilon_ijk

end module mhd_state_module

