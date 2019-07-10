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

    use amrex_constants_module
    use meth_params_module, only: URHO, NVAR, UFS
    use probdata_module, only: refine_cutoff_height
    use prob_params_module, only : dim
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: tag_lo(3), tag_hi(3)
    integer,    intent(in   ) :: state_lo(3), state_hi(3)
    integer(1), intent(inout) :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    real(rt),   intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    real(rt),   intent(in   ) :: dx(3), problo(3)
    integer(1), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: level
    real(rt),   intent(in   ), value :: time

    integer  :: i, j, k
    real(rt) :: y, z, height

    do k = lo(3), hi(3)
       z = problo(3) + (dble(k) + HALF)*dx(3)
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j) + HALF)*dx(2)
          do i = lo(1), hi(1)
             if (dim == 2) then
                height = y
             else
                height = z
             endif
             if ( height > refine_cutoff_height) then
                tag(i,j,k) = clear
             end if
          end do
       end do
    end do

  end subroutine set_problem_tags

end module problem_tagging_module

