module castro_util_module

  implicit none

contains

  ! Given 3D indices (i,j,k), return the cell-centered spatial position.
  ! Optionally we can also be edge-centered in any of the directions.
  
  function position(i, j, k, ccx, ccy, ccz)

    use amrinfo_module, only: amr_level
    use prob_params_module, only: problo, probhi, physbc_lo, physbc_hi, dx_level, &
                                  domlo_level, domhi_level, Interior
    use bl_constants_module, only: ZERO, HALF

    ! Input arguments
    integer :: i, j, k
    logical, optional :: ccx, ccy, ccz

    ! Local variables
    double precision :: position(3), dx(3), offset(3)
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
       if (cc(dir)) then
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

    position(:) = offset(:) + dble(idx(:)) * dx(:)

  end function position



  ! Given 3D indices (i,j,k) and a direction dir, return the
  ! area of the face perpendicular to direction d. We assume
  ! the coordinates perpendicular to the dir axies are edge-centered.
  ! Note that Castro has no support for angular coordinates, so 
  ! this function only provides Cartesian in 1D/2D/3D, Cylindrical (R-Z)
  ! in 2D, and Spherical in 1D.

  function area(i, j, k, dir)

    use amrinfo_module, only: amr_level
    use bl_constants_module, only: ZERO, ONE, TWO, M_PI, FOUR
    use prob_params_module, only: dim, coord_type, dx_level

    implicit none

    integer, intent(in) :: i, j, k, dir

    double precision :: area

    logical :: cc(3) = .true.
    double precision :: dx(3), loc(3)

    ! Force edge-centering along the direction of interest

    cc(dir) = .false.

    dx = dx_level(:,amr_level)

    if (coord_type .eq. 0) then

       ! Cartesian (1D/2D/3D)

       if (dim .eq. 1) then

          select case (dir)

          case (1)
             area = ONE
          case default
             area = ZERO

          end select

       else if (dim .eq. 2) then

          select case (dir)

          case (1)
             area = dx(2)
          case (2)
             area = dx(1)
          case default
             area = ZERO

          end select

       else if (dim .eq. 3) then

          select case (dir)

          case (1)
             area = dx(2) * dx(3)
          case (2)
             area = dx(1) * dx(3)
          case (3)
             area = dx(1) * dx(2)
          case default
             area = ZERO

          end select

       endif

    else if (coord_type .eq. 1) then

       ! Cylindrical (2D only)

       ! Get edge-centered position

       loc = position(i,j,k,cc(1),cc(2),cc(3))

       if (dim .eq. 2) then

          select case (dir)

          case (1)
             area = TWO * M_PI * loc(1) * dx(2)
          case (2)
             area = TWO * M_PI * loc(1) * dx(1)
          case default
             area = ZERO

          end select

       else

          call bl_error("Cylindrical coordinates only supported in 2D.")

       endif

    else if (coord_type .eq. 2) then

       ! Spherical (1D only)

       ! Get edge-centered position

       loc = position(i,j,k,cc(1),cc(2),cc(3))

       if (dim .eq. 1) then

          select case (dir)

          case (1)
             area = FOUR * M_PI * loc(1)**2
          case default
             area = ZERO

          end select

       else

          call bl_error("Spherical coordinates only supported in 1D.")

       endif

    endif

  end function area



  ! Given 3D cell-centered indices (i,j,k), return the volume of the zone.
  ! Note that Castro has no support for angular coordinates, so 
  ! this function only provides Cartesian in 1D/2D/3D, Cylindrical (R-Z)
  ! in 2D, and Spherical in 1D.

  function volume(i, j, k)

    use amrinfo_module, only: amr_level
    use bl_constants_module, only: ZERO, HALF, FOUR3RD, TWO, M_PI
    use prob_params_module, only: dim, coord_type, dx_level

    implicit none

    integer, intent(in) :: i, j, k

    double precision :: volume

    double precision :: dx(3), loc_l(3), loc_r(3)

    dx = dx_level(:,amr_level)

    if (coord_type .eq. 0) then

       ! Cartesian (1D/2D/3D)

       if (dim .eq. 1) then

          volume = dx(1)

       else if (dim .eq. 2) then

          volume = dx(1) * dx(2)

       else if (dim .eq. 3) then

          volume = dx(1) * dx(2) * dx(3)

       endif

    else if (coord_type .eq. 1) then

       ! Cylindrical (2D only)

       ! Get inner and outer radii

       loc_l = position(i  ,j,k,ccx=.true.)
       loc_r = position(i+1,j,k,ccx=.true.)

       if (dim .eq. 2) then

          volume = TWO * M_PI * (HALF * (loc_l(1) + loc_r(1))) * dx(1) * dx(2)

       else

          call bl_error("Cylindrical coordinates only supported in 2D.")

       endif

    else if (coord_type .eq. 2) then

       ! Spherical (1D only)

       ! Get inner and outer radii

       loc_l = position(i  ,j,k,ccx=.true.)
       loc_r = position(i+1,j,k,ccx=.true.)

       if (dim .eq. 1) then

          volume = FOUR3RD * M_PI * (loc_r(1)**3 - loc_l(1)**3)

       else

          call bl_error("Spherical coordinates only supported in 1D.")

       endif

    endif

  end function volume

end module castro_util_module
