module problem_tagging_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  ! This is a template routine for users to set their own tags based on the state.
  ! It will be overwritten by having a copy of this file in the user's problem setup.
  
  subroutine set_problem_tags(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                              state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,&
                              set,clear,&
                              lo,hi,&
                              dx,problo,time,level) &
                              bind(C, name="set_problem_tags")

    use meth_params_module, only: NVAR

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer         ,intent(in   ) :: lo(3),hi(3)
    integer         ,intent(in   ) :: state_l1,state_l2,state_l3, &
                                      state_h1,state_h2,state_h3
    integer         ,intent(in   ) :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
    real(rt)        ,intent(in   ) :: state(state_l1:state_h1, &
                                      state_l2:state_h2, &
                                      state_l3:state_h3,NVAR)
    integer         ,intent(inout) :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
    real(rt)        ,intent(in   ) :: problo(3),dx(3),time
    integer         ,intent(in   ) :: level,set,clear

  end subroutine set_problem_tags

end module problem_tagging_module
