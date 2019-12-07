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
    use probdata_module, only: max_base_tagging_level, tag_density, tag_max_density_fraction
    use eos_type_module, only: maxdens
    use prob_params_module, only: center
    use iso_c_binding, only : c_int8_t

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: tag_lo(3), tag_hi(3)
    integer,    intent(in   ) :: state_lo(3), state_hi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    real(rt),   intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    real(rt),   intent(in   ) :: problo(3), dx(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: level
    real(rt),   intent(in   ), value :: time

    integer  :: i, j, k
    real(rt) :: x, y, dist

    !$gpu

    ! Tag cells for refinement based on tag_density and the EOS maximum density

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (state(i,j,k,URHO) > tag_density) then
                if (level < max_base_tagging_level) then
                   tag(i,j,k) = set
                end if
             else if (state(i,j,k,URHO) >= maxdens * tag_max_density_fraction) then
                 tag(i,j,k) = set
             end if
          end do
       end do
    end do

  end subroutine set_problem_tags

end module problem_tagging_module
