module problem_tagging_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  ! This is a template routine for users to set their own tags based on the state.
  ! It will be overwritten by having a copy of this file in the user's problem setup.
  
  subroutine set_problem_tags(lo, hi, &
                              tag, tagl1, tagh1, &
                              state, state_l1, state_h1, &
                              dx, problo, &
                              set, clear, time, level) &
                              bind(C, name="set_problem_tags")

    use meth_params_module, only: NVAR
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,    intent(in   ) :: lo(1), hi(1)
    integer,    intent(in   ) :: tagl1, tagh1
    integer,    intent(in   ) :: state_l1, state_h1
    integer(1), intent(inout) :: tag(tagl1:tagh1)
    real(rt),   intent(in   ) :: state(state_l1:state_h1, NVAR)
    real(rt),   intent(in   ) :: dx(1), problo(1)
    integer(1), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: level
    real(rt),   intent(in   ), value :: time

  end subroutine set_problem_tags

end module problem_tagging_module
