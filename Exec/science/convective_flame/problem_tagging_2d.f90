module problem_tagging_module

  implicit none

  public

contains

  subroutine set_problem_tags(tag,tagl1,tagl2,tagh1,tagh2, &
                              state,state_l1,state_l2,state_h1,state_h2,&
                              set,clear,&
                              lo,hi,&
                              dx,problo,time,level) bind(C, name="set_problem_tags")


    use probdata_module, only: refine_cutoff_height
    use meth_params_module, only: NVAR
    
    implicit none

    integer         ,intent(in   ) :: lo(2),hi(2)
    integer         ,intent(in   ) :: state_l1,state_l2,state_h1,state_h2
    integer         ,intent(in   ) :: tagl1,tagl2,tagh1,tagh2
    double precision,intent(in   ) :: state(state_l1:state_h1,state_l2:state_h2, NVAR)
    integer         ,intent(inout) :: tag(tagl1:tagh1,tagl2:tagh2)
    double precision,intent(in   ) :: problo(2),dx(2),time
    integer         ,intent(in   ) :: level,set,clear

    integer :: i, j

    double precision :: y
    
    do j = lo(2), hi(2)
       y = problo(2) + dble(j + 0.5d0)*dx(2)
       do i = lo(1), hi(1)
          if (y > refine_cutoff_height) then
             tag(i,j) = clear
          endif
       enddo
    enddo

  end subroutine set_problem_tags

end module problem_tagging_module

