module problem_tagging_module

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding, only : c_int8_t

  implicit none

  public

contains

  subroutine set_problem_tags(lo, hi, &
                              tag, tag_lo, tag_hi, &
                              state, state_lo, state_hi, &
                              dx, problo, &
                              set, clear, time, level) &
                              bind(C, name="set_problem_tags")
    ! This is a template routine for users to set their own tags based on the state.
    ! It will be overwritten by having a copy of this file in the user's problem setup.

    use meth_params_module, only: NVAR

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: tag_lo(3), tag_hi(3)
    integer,    intent(in   ) :: state_lo(3), state_hi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(tag_lo(1):tag_hi(1), tag_lo(2):tag_hi(2), tag_lo(3):tag_hi(3))
    real(rt),   intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3), NVAR)
    real(rt),   intent(in   ) :: problo(3), dx(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: level
    real(rt),   intent(in   ), value :: time

    !$gpu

  end subroutine set_problem_tags

end module problem_tagging_module
