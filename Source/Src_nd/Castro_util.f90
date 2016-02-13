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

  
  
  subroutine ca_enforce_consistent_e(lo,hi,state,s_lo,s_hi) &
       bind(C, name="ca_enforce_consistent_e")

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    ! Local variables
    integer          :: i,j,k
    double precision :: u, v, w, rhoInv

    ! 
    ! Enforces (rho E) = (rho e) + 1/2 rho (u^2 +_ v^2 + w^2)
    !
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE / state(i,j,k,URHO)
             u = state(i,j,k,UMX) * rhoInv
             v = state(i,j,k,UMY) * rhoInv
             w = state(i,j,k,UMZ) * rhoInv

             state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                  HALF * state(i,j,k,URHO) * (u*u + v*v + w*w)

          end do
       end do
    end do

  end subroutine ca_enforce_consistent_e



  subroutine reset_internal_e(lo,hi,u,u_lo,u_hi,verbose) &
       bind(C, name="reset_internal_e")

    use eos_module 
    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFX, &
         small_temp, allow_negative_energy, &
         dual_energy_eta2, dual_energy_update_E_from_e
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3), verbose
    integer          :: u_lo(3), u_hi(3)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    ! Local variables
    integer          :: i,j,k
    double precision :: Up, Vp, Wp, ke, rho_eint, eint_new, rhoInv

    type (eos_t) :: eos_state

    ! Reset internal energy
    if (allow_negative_energy .eq. 0) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                rhoInv = ONE/u(i,j,k,URHO)
                Up = u(i,j,k,UMX) * rhoInv
                Vp = u(i,j,k,UMY) * rhoInv
                Wp = u(i,j,k,UMZ) * rhoInv
                ke = HALF * (Up**2 + Vp**2 + Wp**2)

                rho_eint = u(i,j,k,UEDEN) - u(i,j,k,URHO) * ke

                ! Reset (e from e) if it's greater than eta * E.
                if (rho_eint .gt. ZERO .and. rho_eint / u(i,j,k,UEDEN) .gt. dual_energy_eta2) then

                   u(i,j,k,UEINT) = rho_eint

                   ! If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
                else if (u(i,j,k,UEINT) .gt. ZERO .and. dual_energy_update_E_from_e) then

                   u(i,j,k,UEDEN) = u(i,j,k,UEINT) + u(i,j,k,URHO) * ke

                   ! If not resetting and little e is negative ...
                else if (u(i,j,k,UEINT) .le. ZERO) then

                   eos_state % rho   = u(i,j,k,URHO)
                   eos_state % T     = small_temp
                   eos_state % xn(:) = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
                   eos_state % aux(1:naux) = u(i,j,k,UFX:UFX+naux-1) * rhoInv

                   call eos(eos_input_rt, eos_state)

                   eint_new = eos_state % e

                   if (verbose .gt. 0) then
                      print *,'   '
                      print *,'>>> Warning: Castro_3d::reset_internal_energy  ',i,j,k
                      print *,'>>> ... resetting neg. e from EOS using small_temp'
                      print *,'>>> ... from ',u(i,j,k,UEINT)/u(i,j,k,URHO),' to ', eint_new
                      print *,'    '
                   end if

                   u(i,j,k,UEDEN) = u(i,j,k,UEDEN) + (u(i,j,k,URHO) * eint_new - u(i,j,k,UEINT))
                   u(i,j,k,UEINT) = u(i,j,k,URHO) * eint_new

                end if
             enddo
          enddo
       enddo

       ! If (allow_negative_energy .eq. 1) then just reset (rho e) from (rho E)
    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                rhoInv = ONE/u(i,j,k,URHO)
                Up = u(i,j,k,UMX) * rhoInv
                Vp = u(i,j,k,UMY) * rhoInv
                Wp = u(i,j,k,UMZ) * rhoInv
                ke = HALF * (Up**2 + Vp**2 + Wp**2)

                u(i,j,k,UEINT) = u(i,j,k,UEDEN) - u(i,j,k,URHO) * ke

             enddo
          enddo
       enddo

    endif

  end subroutine reset_internal_e  



  subroutine compute_temp(lo,hi,state,s_lo,s_hi) &
       bind(C, name="compute_temp")

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, UTEMP, &
         UFS, UFX, allow_negative_energy
    use bl_constants_module

    implicit none

    integer         , intent(in   ) :: lo(3),hi(3)
    integer         , intent(in   ) :: s_lo(3),s_hi(3)
    double precision, intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer          :: i,j,k
    double precision :: rhoInv

    type (eos_t) :: eos_state

    ! First check the inputs for validity.

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             if (state(i,j,k,URHO) <= ZERO) then
                print *,'   '
                print *,'>>> Error: Castro_3d::compute_temp ',i,j,k
                print *,'>>> ... negative density ',state(i,j,k,URHO)
                print *,'    '
                call bl_error("Error:: compute_temp_nd.f90")
             end if

             if (allow_negative_energy .eq. 0 .and. state(i,j,k,UEINT) <= ZERO) then
                print *,'   '
                print *,'>>> Warning: Castro_3d::compute_temp ',i,j,k
                print *,'>>> ... negative (rho e) ',state(i,j,k,UEINT)
                print *,'   '
                call bl_error("Error:: compute_temp_nd.f90")
             end if

          enddo
       enddo
    enddo

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE / state(i,j,k,URHO)

             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   = state(i,j,k,UTEMP) ! Initial guess for the EOS
             eos_state % e   = state(i,j,k,UEINT) * rhoInv
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_re, eos_state)

             state(i,j,k,UTEMP) = eos_state % T

             ! In case we've floored, or otherwise allowed the energy to change, update the energy accordingly.

             state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + (state(i,j,k,URHO) * eos_state % e - state(i,j,k,UEINT))
             state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e

          enddo
       enddo
    enddo

  end subroutine compute_temp 



  subroutine ca_normalize_species(u,u_lo,u_hi,lo,hi) &
       bind(C, name="ca_normalize_species")

    use network, only : nspec, smallx
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module, only: ONE

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    ! Local variables
    integer          :: i, j, k
    double precision :: xn(nspec)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             xn = u(i,j,k,UFS:UFS+nspec-1)

             xn = max(smallx * u(i,j,k,URHO), min(u(i,j,k,URHO), xn))

             xn = u(i,j,k,URHO) * (xn / sum(xn))

             u(i,j,k,UFS:UFS+nspec-1) = xn

          enddo
       enddo
    enddo

  end subroutine ca_normalize_species



  ! Given 3D spatial coordinates, return the cell-centered zone indices closest to it.
  ! Optionally we can also be edge-centered in any of the directions.
  
  function position_to_index(loc) result(index)

    use amrinfo_module, only: amr_level
    use prob_params_module, only: dx_level, dim
    
    double precision :: loc(3)

    integer :: index(3)

    index(1:dim)   = NINT(loc(1:dim) / dx_level(1:dim,amr_level))
    index(dim+1:3) = 0
    
  end function position_to_index  



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
