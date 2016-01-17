subroutine set_problem_tags(tag,tag_lo,tag_hi, &
                            state,state_lo,state_hi, &
                            set,clear,&
                            lo,hi,&
                            dx,problo,time,level) bind(C)

  use meth_params_module, only : NVAR, URHO
  use probdata_module, only : cutoff_density, H_min

  implicit none
  
  integer          :: lo(3),hi(3)
  integer          :: state_lo(3),state_hi(3)
  integer          :: tag_lo(3),tag_hi(3)
  double precision :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3),NVAR)
  integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
  double precision :: problo(3),dx(3),time
  integer          :: level,set,clear

  ! Tag on regions of with H > H_min and rho < cutoff_density.
  ! Note that H is the first species variable and so is in index UFS of the state array.

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           if ( state(i,j,k,URHO) > cutoff_density .and. var(i,j,k,UFS) > H_min) then
              tag(i,j,k) = set
           endif
        enddo
     enddo
  enddo
  
end subroutine set_problem_tags

