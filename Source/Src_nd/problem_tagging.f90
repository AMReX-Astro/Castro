! This is a template routine for users to set their own tags based on the state.
! It will be overwritten by having a copy of this file in the user's problem setup.

subroutine set_problem_tags(tag,t_lo,t_hi, &
                            state,s_lo,s_hi,,&
                            set,clear,&
                            lo,hi,&
                            dx,time,level)

  use prob_params_module, only: problo
 
  implicit none
  
  integer         ,intent(in   ) :: lo(3),hi(3)
  integer         ,intent(in   ) :: s_lo(3), s_hi(3)
  integer         ,intent(in   ) :: t_lo(3), t_hi(3)
  double precision,intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
  integer         ,intent(inout) :: tag(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
  double precision,intent(in   ) :: dx(3),time
  integer         ,intent(in   ) :: level,set,clear
  
end subroutine set_problem_tags

