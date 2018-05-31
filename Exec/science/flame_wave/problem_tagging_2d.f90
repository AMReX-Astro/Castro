module problem_tagging_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine set_problem_tags(lo, hi, tag,tagl1,tagl2,tagh1,tagh2, &
                              state,state_l1,state_l2,state_h1,state_h2,&
                              set,clear,&
                              dx,problo,time,level) bind(C, name="set_problem_tags")

    use meth_params_module, only: URHO, NVAR, UFS
    use probdata_module, only: X_min, cutoff_density, &
         burn_tagging_min, burn_tagging_max, max_hse_tagging_level

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer         ,intent(in   ) :: lo(2),hi(2)
    integer         ,intent(in   ) :: state_l1,state_l2,state_h1,state_h2
    integer         ,intent(in   ) :: tagl1,tagl2,tagh1,tagh2
    real(rt)        ,intent(in   ) :: state(state_l1:state_h1,state_l2:state_h2, NVAR)
    integer         ,intent(inout) :: tag(tagl1:tagh1,tagl2:tagh2)
    real(rt)        ,intent(in   ) :: problo(2),dx(2),time
    integer         ,intent(in   ) :: level,set,clear

    integer :: i, j

    ! Tag on regions of with X > X_min (note: this is the first
    ! species variable) and rho < cutoff_density.  Note that X is the
    ! first species variable and so is in index UFS of the state
    ! array.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if ( state(i,j,URHO) > cutoff_density .and. state(i,j,UFS) > X_min) then
             if (level < max_hse_tagging_level) then
                tag(i,j) = set
             endif
          endif
       enddo
    enddo

    ! additional tagging for just a layer where we are burning
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          if ( state(i,j,URHO) >= burn_tagging_min .and. &
               state(i,j,URHO) <= burn_tagging_max) then
             tag(i,j) = set

          else if (level > max_hse_tagging_level) then
             tag(i,j) = clear
          endif

       enddo
    enddo

  end subroutine set_problem_tags

end module problem_tagging_module

