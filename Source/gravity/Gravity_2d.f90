module gravity_2D_module

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_test_residual(lo, hi, &
       rhs, rhl1, rhl2, rhh1, rhh2, &
       ecx, ecxl1, ecxl2, ecxh1, ecxh2, &
       ecy, ecyl1, ecyl2, ecyh1, ecyh2, &
       dx,problo,coord_type) bind(C, name="ca_test_residual")

    use amrex_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(2),hi(2)
    integer , intent(in   ) :: coord_type
    integer , intent(in   ) :: rhl1, rhl2, rhh1, rhh2
    integer , intent(in   ) :: ecxl1, ecxl2, ecxh1, ecxh2
    integer , intent(in   ) :: ecyl1, ecyl2, ecyh1, ecyh2
    real(rt), intent(inout) :: rhs(rhl1:rhh1,rhl2:rhh2)
    real(rt), intent(in   ) :: ecx(ecxl1:ecxh1,ecxl2:ecxh2)
    real(rt), intent(in   ) :: ecy(ecyl1:ecyh1,ecyl2:ecyh2)
    real(rt), intent(in   ) :: dx(2), problo(2)

    ! Local variables
    real(rt)         :: lapphi
    real(rt)         :: ri,ro,rc
    integer          :: i,j

    ! Cartesian
    if (coord_type .eq. 0) then

       do i=lo(1),hi(1)
          do j=lo(2),hi(2)
             lapphi = (ecx(i+1,j)-ecx(i,j)) / dx(1) + &
                  (ecy(i,j+1)-ecy(i,j)) / dx(2)
             rhs(i,j) = rhs(i,j) - lapphi
          enddo
       enddo

       ! r-z
    else if (coord_type .eq. 1) then

       do i=lo(1),hi(1)

          ri = problo(1) + dble(i  )*dx(1)
          ro = problo(1) + dble(i+1)*dx(1)
          rc = HALF * (ri + ro)

          do j=lo(2),hi(2)

             lapphi = (ro*ecx(i+1,j)-ri*ecx(i,j)) / (rc*dx(1)) + &
                  (   ecy(i,j+1)-   ecy(i,j)) / dx(2)

             rhs(i,j) = rhs(i,j) - lapphi
          enddo
       enddo

    else

       print *,'Bogus coord_type in test_residual ' ,coord_type
       call castro_error("Error:: Gravity_2d.f90 :: ca_test_residual")

    end if

  end subroutine ca_test_residual



  subroutine ca_compute_radial_mass (lo,hi,dx,dr,&
       state,r_l1,r_l2,r_h1,r_h2, &
       radial_mass,radial_vol,problo, &
       n1d,drdxfac,level) bind(C, name="ca_compute_radial_mass")
    
    use amrex_constants_module
    use prob_params_module, only: center
    use meth_params_module, only: NVAR, URHO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(2),hi(2)
    real(rt), intent(in   ) :: dx(2),dr
    real(rt), intent(in   ) :: problo(2)

    integer , intent(in   ) :: n1d,drdxfac,level
    real(rt), intent(inout) :: radial_mass(0:n1d-1)
    real(rt), intent(inout) :: radial_vol (0:n1d-1)

    integer , intent(in   ) :: r_l1,r_l2,r_h1,r_h2
    real(rt), intent(in   ) :: state(r_l1:r_h1,r_l2:r_h2,NVAR)

    integer          :: i,j,index
    integer          :: ii,jj
    real(rt)         :: xc,yc,r,octant_factor
    real(rt)         :: fac,xx,yy,dx_frac,dy_frac,vol_frac,vol_frac_fac
    real(rt)         :: lo_i,lo_j,rlo,rhi

    if (abs(center(2) - problo(2)) .lt. 1.e-2 * dx(2)) then
       octant_factor = TWO
    else
       octant_factor = ONE
    end if

    fac  = dble(drdxfac)
    dx_frac = dx(1) / fac
    dy_frac = dx(2) / fac

    do j = lo(2), hi(2)
       yc = problo(2) + (dble(j)+HALF) * dx(2) - center(2)

       do i = lo(1), hi(1)
          xc = problo(1) + (dble(i)+HALF) * dx(1) - center(1)

          vol_frac_fac = TWO * M_PI * dx_frac * dy_frac * octant_factor

          r = sqrt(xc**2  + yc**2)
          index = int(r/dr)

          if (index .gt. n1d-1) then

             if (level .eq. 0) then
                print *,'   '
                print *,'>>> Error: Gravity_2d::ca_compute_radial_mass ',i,j
                print *,'>>> ... index too big: ', index,' > ',n1d-1
                print *,'>>> ... at (i,j)     : ',i,j
                print *,'    ' 
                call castro_error("Error:: Gravity_2d.f90 :: ca_compute_radial_mass")
             end if

          else

             ! Note that we assume we are in r-z coordinates in 2d or we wouldn't be 
             !      doing monopole gravity
             lo_i = problo(1) + dble(i)*dx(1) - center(1)
             lo_j = problo(2) + dble(j)*dx(2) - center(2)
             do ii = 0,drdxfac-1
                xx  = lo_i + (dble(ii  )+HALF)*dx_frac
                rlo = lo_i +  dble(ii  )       *dx_frac
                rhi = lo_i +  dble(ii+1)       *dx_frac
                vol_frac = vol_frac_fac * xx 
                do jj = 0,drdxfac-1
                   yy = lo_j + (dble(jj)+HALF)*dy_frac
                   r = sqrt(xx**2  + yy**2)
                   index = int(r/dr)
                   if (index .le. n1d-1) then
                      radial_mass(index) = radial_mass(index) + vol_frac*state(i,j,URHO)
                      radial_vol (index) = radial_vol (index) + vol_frac
                   end if
                end do
             end do

          end if
       enddo
    enddo

  end subroutine ca_compute_radial_mass



  subroutine ca_put_radial_grav (lo,hi,dx,dr,&
       grav,g_l1,g_l2,g_h1,g_h2, &
       radial_grav,problo,n1d,level) bind(C, name="ca_put_radial_grav")

    use prob_params_module, only: center
    use amrex_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    integer , intent(in   ) :: lo(2),hi(2)
    real(rt), intent(in   ) :: dx(2),dr
    real(rt), intent(in   ) :: problo(2)

    integer , intent(in   ) :: n1d,level
    real(rt), intent(in   ) :: radial_grav(0:n1d-1)

    integer , intent(in   ) :: g_l1,g_l2,g_h1,g_h2
    real(rt), intent(inout) :: grav(g_l1:g_h1,g_l2:g_h2,2)

    integer          :: i,j,index
    real(rt)         :: x,y,r,mag_grav
    real(rt)         :: cen,xi,slope,glo,gmd,ghi,minvar,maxvar

    ! Note that we are interpolating onto the entire range of grav,
    ! including the ghost cells

    do j = lo(2), hi(2)
       y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
       do i = lo(1), hi(1)
          x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
          r = sqrt( x**2 + y**2 )
          index = int(r/dr)
          cen = (dble(index)+HALF)*dr
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
                print *,'PUT_RADIAL_GRAV: INDEX TOO BIG ',index,' > ',n1d-1
                print *,'AT (i,j) ',i,j
                print *,'R / DR IS ',r,dr
                call castro_error("Error:: Gravity_2d.f90 :: ca_put_radial_grav")
             else 
                ! NOTE: we don't do anything to this point if it's outside the
                !       radial grid and level > 0
             end if
          else 
             ! Quadratic interpolation
             ghi = radial_grav(index+1)
             gmd = radial_grav(index  )
             glo = radial_grav(index-1)
             mag_grav = &
                  ( ghi -  TWO*gmd + glo)*xi**2/(TWO*dr**2) + &
                  ( ghi       - glo     )*xi   /(TWO*dr   ) + &
                  (-ghi + 26.e0_rt*gmd - glo)/24.e0_rt

             minvar = min(gmd, min(glo,ghi))
             maxvar = max(gmd, max(glo,ghi))
             mag_grav = max(mag_grav,minvar)
             mag_grav = min(mag_grav,maxvar)
          end if

          if (index .le. n1d-1) then
             grav(i,j,1) = mag_grav * (x/r)
             grav(i,j,2) = mag_grav * (y/r)
          end if
       enddo
    enddo

  end subroutine ca_put_radial_grav



  subroutine ca_put_radial_phi (lo,hi,domlo,domhi,dx,dr,&
       phi,p_l1,p_l2,p_h1,p_h2, &
       radial_phi,problo, &
       numpts_1d,fill_interior) bind(C, name="ca_put_radial_phi")

    use prob_params_module, only: center
    use amrex_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    integer , intent(in   ) :: lo(2),hi(2)
    integer , intent(in   ) :: domlo(2),domhi(2)
    real(rt), intent(in   ) :: dx(2),dr
    real(rt), intent(in   ) :: problo(2)

    integer , intent(in   ) :: numpts_1d
    real(rt), intent(in   ) :: radial_phi(0:numpts_1d-1)
    integer , intent(in   ) :: fill_interior

    integer , intent(in   ) :: p_l1,p_l2,p_h1,p_h2
    real(rt), intent(inout) :: phi(p_l1:p_h1,p_l2:p_h2)

    integer          :: i,j,index
    real(rt)         :: x,y,r
    real(rt)         :: cen,xi,slope,p_lo,p_md,p_hi,minvar,maxvar

    ! Note that when we interpolate into the ghost cells we use the
    ! location of the edge, not the cell center

    do j = lo(2), hi(2)
       if (j .gt. domhi(2)) then
          y = problo(2) + (dble(j  )     ) * dx(2) - center(2)
       else if (j .lt. domlo(2)) then
          y = problo(2) + (dble(j+1)     ) * dx(2) - center(2)
       else 
          y = problo(2) + (dble(j  )+HALF) * dx(2) - center(2)
       end if
       do i = lo(1), hi(1)
          if (i .gt. domhi(1)) then
             x = problo(1) + (dble(i  )     ) * dx(1) - center(1)
          else if (i .lt. domlo(1)) then
             x = problo(1) + (dble(i+1)     ) * dx(1) - center(1)
          else 
             x = problo(1) + (dble(i  )+HALF) * dx(1) - center(1)
          end if
          r = sqrt( x**2 + y**2)
          index = int(r/dr)
          if (index .gt. numpts_1d-1) then
             print *,'PUT_RADIAL_PHI: INDEX TOO BIG ',index,' > ',numpts_1d-1
             print *,'AT (i,j) ',i,j
             print *,'R / DR IS ',r,dr
             call castro_error("Error:: Gravity_2d.f90 :: ca_put_radial_phi")
          end if

          if (  (fill_interior .eq. 1) .or. &
               ( j.lt.domlo(2).or.j.gt.domhi(2)  .or. &
               i.lt.domlo(1).or.i.gt.domhi(1)  ) ) then  

             cen = (dble(index)+HALF)*dr
             xi = r - cen
             if (index == 0) then
                ! Linear interpolation or extrapolation
                slope = ( radial_phi(index+1) - radial_phi(index) ) / dr
                phi(i,j) = radial_phi(index) + slope * xi
             else if (index == numpts_1d-1) then
                ! Linear interpolation or extrapolation 
                slope = ( radial_phi(index) - radial_phi(index-1) ) / dr
                phi(i,j) = radial_phi(index) + slope * xi
             else 
                ! Quadratic interpolation
                p_hi = radial_phi(index+1)
                p_md = radial_phi(index  )
                p_lo = radial_phi(index-1)
                phi(i,j) = &
                     ( p_hi -  TWO*p_md + p_lo)*xi**2/(TWO*dr**2) + &
                     ( p_hi       - p_lo      )*xi   /(TWO*dr   ) + &
                     (-p_hi + 26.e0_rt*p_md - p_lo)/24.e0_rt
                minvar = min(p_md, min(p_lo,p_hi))
                maxvar = max(p_md, max(p_lo,p_hi))
                phi(i,j) = max(phi(i,j),minvar)
                phi(i,j) = min(phi(i,j),maxvar)
             end if

          end if
       enddo
    enddo

  end subroutine ca_put_radial_phi

end module gravity_2D_module
