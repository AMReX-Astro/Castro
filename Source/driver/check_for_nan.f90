module nan_check

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

  subroutine check_for_nan(s, s_lo, s_hi, lo, hi, ncomp, comp)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: ncomp
    real(rt), intent(in) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)

    integer :: i, j, k, n, comp

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (isnan(s(i,j,k,comp))) then
                print *, "NaN: ", i, j, k, comp
                call castro_error("NaN")
             endif
          end do
       end do
    end do

  end subroutine check_for_nan

end module nan_check
