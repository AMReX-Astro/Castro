module problem_tagging_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  ! This is a template routine for users to set their own tags based on the state.
  ! It will be overwritten by having a copy of this file in the user's problem setup.
  
  subroutine set_problem_tags(tag,tag_lo,tag_hi, &
                              state,state_lo,state_hi, &
                              set,clear,&
                              lo,hi,&
                              dx,problo,time,level) &
                              bind(C, name="set_problem_tags")

    use meth_params_module, only : NVAR
    
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3),hi(3)
    integer          :: state_lo(3),state_hi(3)
    integer          :: tag_lo(3),tag_hi(3)
    real(rt)         :: state(state_lo(1):state_hi(1), &
                        state_lo(2):state_hi(2), &
                        state_lo(3):state_hi(3),NVAR)
    integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    real(rt)         :: problo(3),dx(3),time
    integer          :: level,set,clear

  end subroutine set_problem_tags

end module problem_tagging_module
