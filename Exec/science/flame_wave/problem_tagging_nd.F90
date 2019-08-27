module problem_tagging_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains

  subroutine set_problem_tags(lo, hi, &
                              tag, tag_lo, tag_hi, &
                              state, state_lo, state_hi, &
                              dx, problo, &
                              set, clear, time, level) &
                              bind(C, name="set_problem_tags")

    use amrex_constants_module, only: HALF
    use meth_params_module, only: URHO, NVAR, UFS
    use probdata_module, only: X_min, cutoff_density, &
                               max_hse_tagging_level, max_base_tagging_level, x_refine_distance
    use prob_params_module, only: center

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: tag_lo(3), tag_hi(3)
    integer,    intent(in   ) :: state_lo(3), state_hi(3)
    integer(1), intent(inout) :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    real(rt),   intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    real(rt),   intent(in   ) :: problo(3), dx(3)
    integer(1), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: level
    real(rt),   intent(in   ), value :: time

    integer  :: i, j, k
    real(rt) :: x, xdist

    !$gpu

    ! Tag on regions of with X > X_min and rho < cutoff_density.  Note
    ! that X is the first species variable and so is in index UFS of
    ! the state array.  This should refine just the fueld layer

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i) + HALF)*dx(1)

             if ( state(i,j,k,URHO) > cutoff_density .and. state(i,j,k,UFS)/state(i,j,k,URHO) > X_min) then
                xdist = abs(x - center(1))
                if (level < max_hse_tagging_level .and. xdist < x_refine_distance) then
                   tag(i,j,k) = set
                end if
             end if
          end do
       end do
    end do

    ! Always tag on the atmosphere up to some level, where we assume that
    ! max_base_tagging_level <= max_hse_tagging_level

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if ( state(i,j,k,URHO) > cutoff_density) then
                if (level < max_base_tagging_level) then
                   tag(i,j,k) = set
                end if
             end if
          end do
       end do
    end do

  end subroutine set_problem_tags

end module problem_tagging_module
