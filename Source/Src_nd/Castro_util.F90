module castro_util_module

  implicit none

contains

  ! Given 3D indices (i,j,k), return the cell-centered spatial position.
  ! Optionally we can also be edge-centered in any of the directions.
  
  function position(i, j, k, ccx, ccy, ccz) result(loc)

    use amrinfo_module, only: amr_level
    use prob_params_module, only: problo, probhi, physbc_lo, physbc_hi, dx_level, &
                                  domlo_level, domhi_level, Interior
    use bl_constants_module, only: ZERO, HALF

    ! Input arguments
    integer :: i, j, k
    logical, optional :: ccx, ccy, ccz

    ! Local variables
    double precision :: loc(3), dx(3), offset(3)
    integer :: idx(3)
    logical :: cc(3)
    integer :: domlo(3), domhi(3)
    integer :: dir
    
    idx = (/ i, j, k /)
    
    dx(:) = dx_level(:,amr_level)
    domlo = domlo_level(:,amr_level)
    domhi = domhi_level(:,amr_level)

    offset(:) = problo(:)       
    
    cc(:) = .true.
    
    if (present(ccx)) then
       cc(1) = ccx
    endif

    if (present(ccy)) then
       cc(2) = ccy
    endif

    if (present(ccz)) then
       cc(3) = ccz
    endif

    do dir = 1, 3
       if (cc(i)) then
          ! If we're cell-centered, we want to be in the middle of the zone.
          
          offset(dir) = offset(dir) + HALF * dx(dir)
       else
          ! Take care of the fact that for edge-centered indexing,
          ! we actually range from (domlo, domhi+1).
          
          domhi(dir) = domhi(dir) + 1
       endif
    enddo

    ! Be careful when using periodic boundary conditions. In that case,
    ! we need to loop around to the other side of the domain.

    do dir = 1, 3
       if      (physbc_lo(dir) .eq. Interior .and. idx(dir) .lt. domlo(dir)) then
          offset(dir) = offset(dir) + (probhi(dir) - problo(dir))
       else if (physbc_hi(dir) .eq. Interior .and. idx(dir) .gt. domhi(dir)) then
          offset(dir) = offset(dir) + (problo(dir) - probhi(dir))
       endif
    enddo
       
    loc(:) = offset(:) + dble(idx(:)) * dx(:)
    
  end function position
 
end module castro_util_module
