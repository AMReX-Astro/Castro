module tagging_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , save ::    denerr,   dengrad
  real(rt)        , save ::    enterr,   entgrad
  real(rt)        , save ::    velerr,   velgrad
  real(rt)        , save ::   temperr,  tempgrad
  real(rt)        , save ::  presserr, pressgrad
  real(rt)        , save ::    raderr,   radgrad
  integer         , save ::  max_denerr_lev,   max_dengrad_lev
  integer         , save ::  max_enterr_lev,   max_entgrad_lev
  integer         , save ::  max_velerr_lev,   max_velgrad_lev
  integer         , save ::  max_temperr_lev,  max_tempgrad_lev
  integer         , save ::  max_presserr_lev, max_pressgrad_lev
  integer         , save ::  max_raderr_lev,   max_radgrad_lev

  public

contains

  ! All tagging subroutines in this file must be threadsafe because
  ! they are called inside OpenMP parallel regions.


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the temperature
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: temp      => temperature array
! ::: np        => number of components in temp array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
  
  subroutine ca_temperror(tag,taglo,taghi, &
       set,clear, &
       temp,templo,temphi, &
       lo,hi,np,domlo,domhi, &
       delta,xlo,problo,time,level) bind(C, name="ca_temperror")

    use prob_params_module, only: dim

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: set, clear, np, level
    integer          :: taglo(3), taghi(3)
    integer          :: templo(3), temphi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt)         :: temp(templo(1):temphi(1),templo(2):temphi(2),templo(3):temphi(3),np)
    real(rt)         :: delta(3), xlo(3), problo(3), time

    real(rt)         :: ax, ay, az, gradT
    integer          :: i, j, k

    !     Tag on regions of high temperature
    if (level .lt. max_temperr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (temp(i,j,k,1) .ge. temperr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high temperature gradient
    if (level .lt. max_tempgrad_lev) then
       if (dim .eq. 1) then
          j = lo(2)
          k = lo(3)
          do i = lo(1), hi(1)
             ax = ABS(temp(i+1,j,k,1) - temp(i,j,k,1))
             ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1,j,k,1)))
             if ( ax / temp(i,j,k,1) .ge. tempgrad) &
                  tag(i,j,k) = set
          enddo
       else if (dim .eq. 2) then
          k = lo(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(temp(i+1,j,k,1) - temp(i,j,k,1))
                ay = ABS(temp(i,j+1,k,1) - temp(i,j,k,1))
                ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1,j,k,1)))
                ay = MAX(ay,ABS(temp(i,j,k,1) - temp(i,j-1,k,1)))
                gradT = sqrt(ax**2 + ay**2)
                if (gradT / temp(i,j,k,1) .ge. tempgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       else
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   ax = ABS(temp(i+1,j,k,1) - temp(i,j,k,1))
                   ay = ABS(temp(i,j+1,k,1) - temp(i,j,k,1))
                   az = ABS(temp(i,j,k+1,1) - temp(i,j,k,1))
                   ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1,j,k,1)))
                   ay = MAX(ay,ABS(temp(i,j,k,1) - temp(i,j-1,k,1)))
                   az = MAX(az,ABS(temp(i,j,k,1) - temp(i,j,k-1,1)))
                   gradT = sqrt(ax**2 + ay**2 + az**2)
                   if (gradT / temp(i,j,k,1) .ge. tempgrad .and. &
                        temp(i,j,k,1) .ge. temperr) then
                      tag(i,j,k) = set
                   endif
                enddo
             enddo
          enddo
       end if
    end if

  end subroutine ca_temperror

end module tagging_module
