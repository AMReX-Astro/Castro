module gravity_3D_module

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_test_residual(lo, hi, &
       rhs, rhl1, rhl2, rhl3, rhh1, rhh2, rhh3, &
       ecx, ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3, &
       ecy, ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3, &
       ecz, eczl1, eczl2, eczl3, eczh1, eczh2, eczh3, &
       dx,problo,coord_type) bind(C, name="ca_test_residual")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    integer , intent(in   ) :: coord_type
    integer , intent(in   ) :: rhl1, rhl2, rhl3, rhh1, rhh2, rhh3
    integer , intent(in   ) :: ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3
    integer , intent(in   ) :: ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3
    integer , intent(in   ) :: eczl1, eczl2, eczl3, eczh1, eczh2, eczh3
    real(rt), intent(inout) :: rhs(rhl1:rhh1,rhl2:rhh2,rhl3:rhh3)
    real(rt), intent(in   ) :: ecx(ecxl1:ecxh1,ecxl2:ecxh2, ecxl3:ecxh3)
    real(rt), intent(in   ) :: ecy(ecyl1:ecyh1,ecyl2:ecyh2, ecyl3:ecyh3)
    real(rt), intent(in   ) :: ecz(eczl1:eczh1,eczl2:eczh2, eczl3:eczh3)
    real(rt), intent(in   ) :: dx(3),problo(3)

    ! Local variables
    integer          :: i,j,k
    real(rt)         :: lapphi

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             lapphi = (ecx(i+1,j,k)-ecx(i,j,k)) / dx(1) + &
                  (ecy(i,j+1,k)-ecy(i,j,k)) / dx(2) + &
                  (ecz(i,j,k+1)-ecz(i,j,k)) / dx(3)
             rhs(i,j,k) = rhs(i,j,k) - lapphi
          enddo
       enddo
    enddo

  end subroutine ca_test_residual



  subroutine ca_put_radial_grav (lo,hi,dx,dr,&
       grav,g_l1,g_l2,g_l3,g_h1,g_h2,g_h3, &
       radial_grav,problo,n1d,level) bind(C, name="ca_put_radial_grav")

    use amrex_constants_module
    use prob_params_module, only: center

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    real(rt), intent(in   ) :: dx(3),dr
    real(rt), intent(in   ) :: problo(3)

    integer , intent(in   ) :: n1d,level
    real(rt), intent(in   ) :: radial_grav(0:n1d-1)

    integer , intent(in   ) :: g_l1,g_l2,g_l3,g_h1,g_h2,g_h3
    real(rt), intent(inout) :: grav(g_l1:g_h1,g_l2:g_h2,g_l3:g_h3,3)

    integer          :: i,j,k,index
    real(rt)         :: x,y,z,r,mag_grav
    real(rt)         :: cen,xi,slope,glo,gmd,ghi,minvar,maxvar
    !
    ! Note that we are interpolating onto the entire range of grav,
    ! including the ghost cells.
    !
    do k = lo(3), hi(3)
       z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)

       do j = lo(2), hi(2)
          y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)

          do i = lo(1), hi(1)
             x     = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
             r     = sqrt(x**2 + y**2 + z**2)
             index = int(r/dr)
             cen   = (dble(index)+HALF)*dr
             xi    = r - cen

             if (index == 0) then
                !
                ! Linear interpolation or extrapolation
                !
                slope = ( radial_grav(index+1) - radial_grav(index) ) / dr
                mag_grav = radial_grav(index) + slope * xi
             else if (index == n1d-1) then
                !
                ! Linear interpolation or extrapolation
                !
                slope = ( radial_grav(index) - radial_grav(index-1) ) / dr
                mag_grav = radial_grav(index) + slope * xi
             else if (index .gt. n1d-1) then

                if (level .eq. 0) then
                   print *,'PUT_RADIAL_GRAV: INDEX TOO BIG ',index,' > ',n1d-1
                   print *,'AT (i,j,k) ',i,j,k
                   print *,'X Y Z ',x,y,z
                   print *,'R / DR ',r,dr
                   call castro_error("Error:: Gravity_3d.f90 :: ca_put_radial_grav")
                else
                   ! NOTE: we don't do anything to this point if it's outside the
                   !       radial grid and level > 0
                end if

             else
                !
                ! Quadratic interpolation
                !
                ghi = radial_grav(index+1)
                gmd = radial_grav(index  )
                glo = radial_grav(index-1)
                mag_grav = &
                     ( ghi -   TWO*gmd + glo)*xi**2/(TWO*dr**2) + &
                     ( ghi       - glo      )*xi   /(TWO*dr   ) + &
                     (-ghi + 26.e0_rt*gmd - glo)/24.e0_rt

                minvar = min(gmd, min(glo,ghi))
                maxvar = max(gmd, max(glo,ghi))
                mag_grav = max(mag_grav,minvar)
                mag_grav = min(mag_grav,maxvar)
             end if

             if (index .le. n1d-1) then
                grav(i,j,k,1) = mag_grav * (x/r)
                grav(i,j,k,2) = mag_grav * (y/r)
                grav(i,j,k,3) = mag_grav * (z/r)
             end if
          enddo
       enddo
    enddo

  end subroutine ca_put_radial_grav

end module gravity_3D_module
