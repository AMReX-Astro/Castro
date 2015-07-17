! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the density
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: den       => density array
! ::: nd        => number of components in den array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
subroutine ca_state_error(tag,tagl1,tagl2,tagh1,tagh2, &
                          set,clear, &
                          var,varl1,varl2,varh1,varh2, &
                          lo,hi,nd,domlo,domhi, &
                          delta,xlo,problo,time,level)
  use probdata_module
  implicit none
  
  integer          :: set, clear, nd, level
  integer          :: tagl1,tagl2,tagh1,tagh2
  integer          :: varl1,varl2,varh1,varh2
  integer          :: lo(2), hi(2), domlo(2), domhi(2)
  integer          ::tag(tagl1:tagh1,tagl2:tagh2)
  double precision :: var(varl1:varh1,varl2:varh2,nd)
  double precision :: delta(2), xlo(2), problo(2), time
  
  integer i, j
      
  ! DENSITY     is component 1                                                            
  ! TEMPERATURE is component 2                                                            
  ! H FRAC      is component 3                                                            
  
  ! Tag on regions of with H > H_min and rho < cutoff_density
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        if ( var(i,j,1) > cutoff_density .and. var(i,j,3) > H_min) then
           tag(i,j) = set
        endif
     enddo
  enddo
      
end subroutine ca_state_error

