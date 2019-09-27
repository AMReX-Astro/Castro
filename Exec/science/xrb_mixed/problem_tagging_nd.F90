module problem_tagging_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine set_problem_tags(lo, hi, &
                              tag, tag_lo, tag_hi, &
                              state, state_lo, state_hi, &
                              dx, problo, &
                              set, clear, time, level) &
                              bind(C, name="set_problem_tags")

    use meth_params_module, only: URHO, NVAR, UFS
    use probdata_module, only: H_min, cutoff_density

    use amrex_fort_module, only : rt => amrex_real
    use iso_c_binding, only : c_int8_t

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: tag_lo(3), tag_hi(3)
    integer,    intent(in   ) :: state_lo(3), state_hi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(tag_lo(1):tag_hi(1), tag_lo(2):tag_hi(2), tag_lo(3):tag_hi(3))
    real(rt),   intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3), NVAR)
    real(rt),   intent(in   ) :: dx(3), problo(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: level
    real(rt),   intent(in   ), value :: time

    integer :: i, j, k

    ! Tag on regions of with H > H_min and rho < cutoff_density.
    ! Note that H is the first species variable and so is in index UFS of the state array.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if ( state(i,j,k,URHO) > cutoff_density .and. state(i,j,k,UFS) > H_min) then
                tag(i,j,k) = set
             end if
          end do
       end do
    end do

  end subroutine set_problem_tags

end module problem_tagging_module

