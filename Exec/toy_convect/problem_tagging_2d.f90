module problem_tagging_module

  implicit none

  public

contains

  subroutine set_problem_tags(tag,tagl1,tagl2,tagh1,tagh2, &
                              state,state_l1,state_l2,state_h1,state_h2,&
                              set,clear,&
                              lo,hi,&
                              dx,problo,time,level) bind(C)

    use meth_params_module, only: URHO, NVAR, UFS
    use probdata_module, only: H_min, cutoff_density

    implicit none

    integer         ,intent(in   ) :: lo(2),hi(2)
    integer         ,intent(in   ) :: state_l1,state_l2,state_h1,state_h2
    integer         ,intent(in   ) :: tagl1,tagl2,tagh1,tagh2
    double precision,intent(in   ) :: state(state_l1:state_h1,state_l2:state_h2, NVAR)
    integer         ,intent(inout) :: tag(tagl1:tagh1,tagl2:tagh2)
    double precision,intent(in   ) :: problo(2),dx(2),time
    integer         ,intent(in   ) :: level,set,clear

    integer :: i, j

    ! Tag on regions of with H > H_min and rho < cutoff_density.
    ! Note that H is the first species variable and so is in index UFS of the state array.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if ( state(i,j,URHO) > cutoff_density .and. state(i,j,UFS) > H_min) then
             tag(i,j) = set
          endif
       enddo
    enddo

  end subroutine set_problem_tags

end module problem_tagging_module

