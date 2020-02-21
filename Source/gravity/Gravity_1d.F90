module gravity_1D_module

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_test_residual(lo, hi, &
       rhs, rhl1, rhh1,  &
       ecx, ecxl1, ecxh1, &
       dx, problo, coord_type) bind(C, name="ca_test_residual")

    use amrex_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(1),hi(1)
    integer , intent(in   ) :: rhl1, rhh1
    integer , intent(in   ) :: ecxl1, ecxh1
    integer , intent(in   ) :: coord_type
    real(rt), intent(inout) :: rhs(rhl1:rhh1)
    real(rt), intent(in   ) :: ecx(ecxl1:ecxh1)
    real(rt), intent(in   ) :: dx(1),problo(1)

    real(rt)         :: lapphi
    real(rt)         :: rlo,rhi,rcen
    integer          :: i

    ! Cartesian
    if (coord_type .eq. 0) then

       do i=lo(1),hi(1)
          lapphi = (ecx(i+1)-ecx(i)) / dx(1)
          rhs(i) = rhs(i) - lapphi
       enddo

       ! r-z
    else if (coord_type .eq. 1) then

       do i=lo(1),hi(1)
          rlo  = problo(1) + dble(i)*dx(1)
          rhi  = rlo + dx(1)
          rcen = HALF * (rlo+rhi)
          lapphi = (rhi*ecx(i+1)-rlo*ecx(i)) / (rcen*dx(1))
          rhs(i) = rhs(i) - lapphi
       enddo

       ! spherical
    else if (coord_type .eq. 2) then

       do i=lo(1),hi(1)
          rlo  = problo(1) + dble(i)*dx(1)
          rhi  = rlo + dx(1)
          rcen = HALF * (rlo+rhi)
          lapphi = (rhi**2 * ecx(i+1)-rlo**2 * ecx(i)) / (rcen**2 * dx(1))
          rhs(i) = rhs(i) - lapphi
       enddo

    else
       print *,'Bogus coord_type in test_residual ' ,coord_type
       call castro_error("Error:: Gravity_1d.f90 :: ca_test_residual")
    end if

  end subroutine ca_test_residual



  subroutine ca_put_radial_grav(lo,hi,dx,dr,&
                                grav,g_l1,g_h1, &
                                radial_grav,problo,n1d,level) bind(C, name="ca_put_radial_grav")

    use prob_params_module, only: center
    use amrex_constants_module, only: HALF, TWO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(1), hi(1)
    real(rt), intent(in   ) :: dx(1), dr
    real(rt), intent(in   ) :: problo(1)

    integer , intent(in   ) :: n1d, level
    real(rt), intent(in   ) :: radial_grav(0:n1d-1)

    integer , intent(in   ) :: g_l1, g_h1
    real(rt), intent(inout) :: grav(g_l1:g_h1)

    integer          :: i, index
    real(rt)         :: r, mag_grav
    real(rt)         :: cen, xi, slope, glo, gmd, ghi, minvar, maxvar

    ! Note that we are interpolating onto the entire range of grav,
    ! including the ghost cells. Taking the absolute value of r ensures
    ! that we will get the correct behavior even for the ghost zones with
    ! negative indices, which have a reflecting boundary condition.

    do i = lo(1), hi(1)

       r = abs(problo(1) + (dble(i) + HALF) * dx(1) - center(1))

       index = int(r / dr)

       cen = (dble(index) + HALF) * dr
       xi = r - cen

       if (index == 0) then

          ! Linear interpolation or extrapolation
          slope = ( radial_grav(index+1) - radial_grav(index) ) / dr
          mag_grav = radial_grav(index) + slope * xi

       else if (index == n1d-1) then

          ! Linear interpolation or extrapolation
          slope = ( radial_grav(index) - radial_grav(index-1) ) / dr
          mag_grav = radial_grav(index) + slope * xi

       else if (index .gt. n1d-1) then

          if (level .eq. 0) then
             print *,'PUT_RADIAL_GRAV: INDEX TOO BIG ', index, ' > ', n1d-1
             print *,'AT i ', i
             print *,'R / DR IS ', r, dr
             call castro_error("Error:: Gravity_1d.f90 :: ca_put_radial_grav")
          else
             ! NOTE: we don't do anything to this point if it's outside the
             !       radial grid and level > 0
          end if

       else

          ! Quadratic interpolation
          ghi = radial_grav(index+1)
          gmd = radial_grav(index  )
          glo = radial_grav(index-1)
          mag_grav = ( ghi -  TWO*gmd + glo)*xi**2/(TWO*dr**2) + &
                     ( ghi       - glo     )*xi   /(TWO*dr   ) + &
                     (-ghi + 26.e0_rt*gmd - glo)/24.e0_rt

          minvar = min(gmd, min(glo,ghi))
          maxvar = max(gmd, max(glo,ghi))
          mag_grav = max(mag_grav,minvar)
          mag_grav = min(mag_grav,maxvar)

       end if

       if (index .le. n1d-1) then

          grav(i) = mag_grav

       end if

    enddo

  end subroutine ca_put_radial_grav

end module gravity_1D_module
