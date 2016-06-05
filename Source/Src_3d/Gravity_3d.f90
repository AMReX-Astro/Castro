module gravity_3D_module

  implicit none

  public
  
contains
  
  subroutine ca_test_residual(lo, hi, &
       rhs, rhl1, rhl2, rhl3, rhh1, rhh2, rhh3, &
       ecx, ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3, &
       ecy, ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3, &
       ecz, eczl1, eczl2, eczl3, eczh1, eczh2, eczh3, &
       dx,problo,coord_type) bind(C, name="ca_test_residual")

    implicit none

    integer          :: lo(3),hi(3)
    integer          :: coord_type
    integer          :: rhl1, rhl2, rhl3, rhh1, rhh2, rhh3
    integer          :: ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3
    integer          :: ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3
    integer          :: eczl1, eczl2, eczl3, eczh1, eczh2, eczh3
    double precision :: rhs(rhl1:rhh1,rhl2:rhh2,rhl3:rhh3)
    double precision :: ecx(ecxl1:ecxh1,ecxl2:ecxh2, ecxl3:ecxh3)
    double precision :: ecy(ecyl1:ecyh1,ecyl2:ecyh2, ecyl3:ecyh3)
    double precision :: ecz(eczl1:eczh1,eczl2:eczh2, eczl3:eczh3)
    double precision :: dx(3),problo(3)

    ! Local variables
    integer          :: i,j,k
    double precision :: lapphi

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



  subroutine ca_compute_radial_mass (lo,hi,dx,dr,&
       rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,&
       radial_mass,radial_vol,problo,&
       n1d,drdxfac,level) bind(C, name="ca_compute_radial_mass")
    
    use bl_constants_module
    use prob_params_module, only: center

    implicit none

    integer          :: lo(3),hi(3)
    double precision :: dx(3),dr
    double precision :: problo(3)

    integer          :: n1d,drdxfac,level
    double precision :: radial_mass(0:n1d-1)
    double precision :: radial_vol (0:n1d-1)

    integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
    double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)

    integer          :: i,j,k,index
    integer          :: ii,jj,kk
    double precision :: xc,yc,zc,r,xxsq,yysq,zzsq,octant_factor
    double precision :: fac,xx,yy,zz,dx_frac,dy_frac,dz_frac
    double precision :: vol_frac, drinv
    double precision :: lo_i,lo_j,lo_k

    if (( abs(center(1) - problo(1)) .lt. 1.d-2 * dx(1) ) .and. &
         ( abs(center(2) - problo(2)) .lt. 1.d-2 * dx(2) ) .and. &
         ( abs(center(3) - problo(3)) .lt. 1.d-2 * dx(3) ) ) then
       octant_factor = EIGHT
    else
       octant_factor = ONE
    end if

    drinv = ONE/dr

    fac     = dble(drdxfac)
    dx_frac = dx(1) / fac
    dy_frac = dx(2) / fac
    dz_frac = dx(3) / fac

    vol_frac = octant_factor * dx_frac * dy_frac * dz_frac

    do k = lo(3), hi(3)
       zc = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
       lo_k =  problo(3) + dble(k)*dx(3) - center(3)

       do j = lo(2), hi(2)
          yc = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
          lo_j =  problo(2) + dble(j)*dx(2) - center(2)

          do i = lo(1), hi(1)
             xc  = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
             lo_i =  problo(1) + dble(i)*dx(1) - center(1)

             r = sqrt(xc**2 + yc**2 + zc**2)
             index = int(r*drinv)

             if (index .gt. n1d-1) then

                if (level .eq. 0) then
                   print *,'   '  
                   print *,'>>> Error: Gravity_3d::ca_compute_radial_mass ',i,j,k
                   print *,'>>> ... index too big: ', index,' > ',n1d-1
                   print *,'>>> ... at (i,j,k)   : ',i,j,k
                   call bl_error("Error:: Gravity_3d.f90 :: ca_compute_radial_mass")
                end if

             else

                do kk = 0,drdxfac-1
                   zz   = lo_k + (dble(kk)+HALF)*dz_frac
                   zzsq = zz*zz
                   do jj = 0,drdxfac-1
                      yy   = lo_j + (dble(jj)+HALF)*dy_frac
                      yysq = yy*yy
                      do ii = 0,drdxfac-1

                         xx    = lo_i + (dble(ii)+HALF)*dx_frac
                         xxsq  = xx*xx
                         r     = sqrt(xxsq  + yysq + zzsq)
                         index = int(r*drinv)

                         if (index .le. n1d-1) then
                            radial_mass(index) = radial_mass(index) + vol_frac * rho(i,j,k)
                            radial_vol (index) = radial_vol (index) + vol_frac
                         end if
                      end do
                   end do
                end do

             end if
          enddo
       enddo
    enddo

  end subroutine ca_compute_radial_mass



  subroutine ca_put_radial_grav (lo,hi,dx,dr,&
       grav,g_l1,g_l2,g_l3,g_h1,g_h2,g_h3, &
       radial_grav,problo,n1d,level) bind(C, name="ca_put_radial_grav")

    use bl_constants_module
    use prob_params_module, only: center

    implicit none

    integer          :: lo(3),hi(3)
    double precision :: dx(3),dr
    double precision :: problo(3)

    integer          :: n1d,level
    double precision :: radial_grav(0:n1d-1)

    integer          :: g_l1,g_l2,g_l3,g_h1,g_h2,g_h3
    double precision :: grav(g_l1:g_h1,g_l2:g_h2,g_l3:g_h3,3)

    integer          :: i,j,k,index
    double precision :: x,y,z,r,mag_grav
    double precision :: cen,xi,slope,glo,gmd,ghi,minvar,maxvar
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
                   call bl_error("Error:: Gravity_3d.f90 :: ca_put_radial_grav")
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
                     (-ghi + 26.d0*gmd - glo)/24.d0

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

  

  subroutine ca_put_radial_phi (lo,hi,domlo,domhi,dx,dr,&
       phi,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
       radial_phi,problo,&
       numpts_1d,fill_interior) bind(C, name="ca_put_radial_phi")
    
    use bl_constants_module
    use prob_params_module, only: center

    implicit none

    integer          :: lo(3),hi(3)
    integer          :: domlo(3),domhi(3)
    double precision :: dx(3),dr
    double precision :: problo(3)

    integer          :: numpts_1d
    double precision :: radial_phi(0:numpts_1d-1)
    integer          :: fill_interior

    integer          :: p_l1,p_l2,p_l3,p_h1,p_h2,p_h3
    double precision :: phi(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3)

    integer          :: i,j,k,index
    double precision :: x,y,z,r
    double precision :: cen,xi,slope,p_lo,p_md,p_hi,minvar,maxvar
    !
    ! Note that when we interpolate into the ghost cells we use the
    ! location of the edge, not the cell center
    !
    do k = lo(3), hi(3)
       if (k .gt. domhi(3)) then
          z = problo(3) + (dble(k  )       ) * dx(3) - center(3)
       else if (k .lt. domlo(3)) then
          z = problo(3) + (dble(k+1)       ) * dx(3) - center(3)
       else 
          z = problo(3) + (dble(k  )+HALF) * dx(3) - center(3)
       end if

       do j = lo(2), hi(2)
          if (j .gt. domhi(2)) then
             y = problo(2) + (dble(j  )       ) * dx(2) - center(2)
          else if (j .lt. domlo(2)) then
             y = problo(2) + (dble(j+1)       ) * dx(2) - center(2)
          else 
             y = problo(2) + (dble(j  )+HALF) * dx(2) - center(2)
          end if

          do i = lo(1), hi(1)
             if (i .gt. domhi(1)) then
                x = problo(1) + (dble(i  )       ) * dx(1) - center(1)
             else if (i .lt. domlo(1)) then
                x = problo(1) + (dble(i+1)       ) * dx(1) - center(1)
             else 
                x = problo(1) + (dble(i  )+HALF) * dx(1) - center(1)
             end if

             r     = sqrt( x**2 + y**2 + z**2 )
             index = int(r/dr)

             if (index .gt. numpts_1d-1) then
                print *,'PUT_RADIAL_PHI: INDEX TOO BIG ',index,' > ',numpts_1d-1
                print *,'AT (i,j) ',i,j,k
                print *,'R / DR IS ',r,dr
                call bl_error("Error:: Gravity_3d.f90 :: ca_put_radial_phi")
             end if

             if ( (fill_interior .eq. 1) .or. &
                  ( i.lt.domlo(1).or.i.gt.domhi(1)  .or. &
                  j.lt.domlo(2).or.j.gt.domhi(2)  .or. &
                  k.lt.domlo(3).or.k.gt.domhi(3)  ) ) then  
                cen = (dble(index)+HALF)*dr
                xi  = r - cen
                if (index == 0) then
                   !
                   ! Linear interpolation or extrapolation
                   !
                   slope      = ( radial_phi(index+1) - radial_phi(index) ) / dr
                   phi(i,j,k) = radial_phi(index) + slope * xi
                else if (index == numpts_1d-1) then
                   !
                   ! Linear interpolation or extrapolation
                   !
                   slope      = ( radial_phi(index) - radial_phi(index-1) ) / dr
                   phi(i,j,k) = radial_phi(index) + slope * xi
                else
                   !
                   ! Quadratic interpolation
                   !
                   p_hi = radial_phi(index+1)
                   p_md = radial_phi(index  )
                   p_lo = radial_phi(index-1)
                   phi(i,j,k) = &
                        ( p_hi -   TWO*p_md + p_lo)*xi**2/(TWO*dr**2) + &
                        ( p_hi       - p_lo      )*xi    /(TWO*dr   ) + &
                        (-p_hi + 26.d0*p_md - p_lo)/24.d0
                   minvar     = min(p_md, min(p_lo,p_hi))
                   maxvar     = max(p_md, max(p_lo,p_hi))
                   phi(i,j,k) = max(phi(i,j,k),minvar)
                   phi(i,j,k) = min(phi(i,j,k),maxvar)
                end if
             end if
          enddo
       enddo
    enddo

  end subroutine ca_put_radial_phi



  subroutine ca_compute_direct_sum_bc (lo,hi,domlo,domhi,&
       symmetry_type,lo_bc,hi_bc, &
       dx,rho,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
       problo, probhi, &
       bcXYLo,bcXYHi,bcXZLo,bcXZHi,bcYZLo,bcYZHi) &
       bind(C, name="ca_compute_direct_sum_bc")
    
    use fundamental_constants_module, only: Gconst
    use bl_constants_module

    implicit none

    integer          :: lo(3),hi(3)
    integer          :: domlo(3),domhi(3)
    double precision :: dx(3),dV
    double precision :: problo(3),probhi(3)

    integer          :: symmetry_type
    integer          :: lo_bc(3), hi_bc(3)

    double precision :: bcXYLo(domlo(1)-1:domhi(1)+1,domlo(2)-1:domhi(2)+1)
    double precision :: bcXYHi(domlo(1)-1:domhi(1)+1,domlo(2)-1:domhi(2)+1)
    double precision :: bcXZLo(domlo(1)-1:domhi(1)+1,domlo(3)-1:domhi(3)+1)
    double precision :: bcXZHi(domlo(1)-1:domhi(1)+1,domlo(3)-1:domhi(3)+1)
    double precision :: bcYZLo(domlo(2)-1:domhi(2)+1,domlo(3)-1:domhi(3)+1)
    double precision :: bcYZHi(domlo(2)-1:domhi(2)+1,domlo(3)-1:domhi(3)+1)

    integer          :: p_l1,p_l2,p_l3,p_h1,p_h2,p_h3
    double precision :: rho(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3)

    integer          :: i,j,k,l,m,n,b
    double precision :: r
    double precision :: loc(3), locb(3), dx2, dy2, dz2

    logical          :: doSymmetricAddLo(3), doSymmetricAddHi(3), doSymmetricAdd

    dV = dx(1) * dx(2) * dx(3)

    ! Determine if we need to add contributions from any symmetric boundaries

    doSymmetricAddLo(:) = .false.
    doSymmetricAddHi(:) = .false.
    doSymmetricAdd      = .false.

    do b = 1, 3
       if ( lo_bc(b) .eq. symmetry_type ) then
          doSymmetricAddLo(b) = .true.
          doSymmetricAdd           = .true.
       endif

       if ( hi_bc(b) .eq. symmetry_type ) then
          doSymmetricAddHi(b) = .true.
          doSymmetricAdd           = .true.
       endif
    enddo

    do k = lo(3), hi(3)
       loc(3) = problo(3) + (dble(k)+HALF) * dx(3)

       do j = lo(2), hi(2)
          loc(2) = problo(2) + (dble(j)+HALF) * dx(2)

          do i = lo(1), hi(1)
             loc(1) = problo(1) + (dble(i)+HALF) * dx(1)

             ! Do xy interfaces first.

             do l = domlo(1) - 1, domhi(1) + 1
                if     ( l .lt. domlo(1) ) then
                   locb(1) = problo(1) + (dble(l+1)     ) * dx(1)
                elseif ( l .gt. domhi(1) ) then
                   locb(1) = problo(1) + (dble(l  )     ) * dx(1)
                else
                   locb(1) = problo(1) + (dble(l  )+HALF) * dx(1)
                endif

                dx2 = (loc(1) - locb(1))**2

                locb(3) = problo(3)
                dz2 = (loc(3) - locb(3))**2

                do m = domlo(2) - 1, domhi(2) + 1
                   if     ( m .lt. domlo(1) ) then
                      locb(2) = problo(2) + (dble(m+1)     ) * dx(2)
                   elseif ( m .gt. domhi(1) ) then
                      locb(2) = problo(2) + (dble(m  )     ) * dx(2)
                   else
                      locb(2) = problo(2) + (dble(m  )+HALF) * dx(2)
                   endif

                   r = ( dx2 + (loc(2) - locb(2))**2 + dz2 )**HALF

                   bcXYLo(l,m) = bcXYLo(l,m) + Gconst * rho(i,j,k) * dV / r

                   ! Now, add any contributions from mass that is hidden behind
                   ! a symmetric boundary.

                   if ( doSymmetricAdd ) then

                      bcXYLo(l,m) = bcXYLo(l,m) + &
                           direct_sum_symmetric_add(loc,locb,problo,probhi, &
                           rho(i,j,k),dV,doSymmetricAddLo,doSymmetricAddHi)

                   endif

                enddo



                locb(3) = probhi(3)
                dz2 = (loc(3) - locb(3))**2

                do m = domlo(2) - 1, domhi(2) + 1
                   if     ( m .lt. domlo(1) ) then
                      locb(2) = problo(2) + (dble(m+1)     ) * dx(2)
                   elseif ( m .gt. domhi(1) ) then
                      locb(2) = problo(2) + (dble(m  )     ) * dx(2)
                   else
                      locb(2) = problo(2) + (dble(m  )+HALF) * dx(2)
                   endif

                   r = ( dx2 + (loc(2) - locb(2))**2 + dz2 )**HALF

                   bcXYHi(l,m) = bcXYHi(l,m) + Gconst * rho(i,j,k) * dV / r

                   if ( doSymmetricAdd ) then

                      bcXYHi(l,m) = bcXYHi(l,m) + &
                           direct_sum_symmetric_add(loc,locb,problo,probhi, &
                           rho(i,j,k),dV,doSymmetricAddLo,doSymmetricAddHi)

                   endif

                enddo



                ! Now do xz interfaces.

                locb(2) = problo(2)
                dy2 = (loc(2) - locb(2))**2

                do n = domlo(3) - 1, domhi(3) + 1
                   if     ( n .lt. domlo(3) ) then
                      locb(3) = problo(3) + (dble(n+1)     ) * dx(3)
                   elseif ( n .gt. domhi(3) ) then
                      locb(3) = problo(3) + (dble(n  )     ) * dx(3)
                   else
                      locb(3) = problo(3) + (dble(n  )+HALF) * dx(3)
                   endif

                   r = ( dx2 + dy2 + (loc(3) - locb(3))**2 )**HALF

                   bcXZLo(l,n) = bcXZLo(l,n) + Gconst * rho(i,j,k) * dV / r

                   if ( doSymmetricAdd ) then

                      bcXZLo(l,n) = bcXZLo(l,n) + &
                           direct_sum_symmetric_add(loc,locb,problo,probhi, &
                           rho(i,j,k),dV,doSymmetricAddLo,doSymmetricAddHi)

                   endif

                enddo

                locb(2) = probhi(2)
                dy2 = (loc(2) - locb(2))**2

                do n = domlo(3) - 1, domhi(3) + 1
                   if     ( n .lt. domlo(3) ) then
                      locb(3) = problo(3) + (dble(n+1)     ) * dx(3)
                   elseif ( n .gt. domhi(3) ) then
                      locb(3) = problo(3) + (dble(n  )     ) * dx(3)
                   else
                      locb(3) = problo(3) + (dble(n  )+HALF) * dx(3)
                   endif

                   r = ( dx2 + dy2 + (loc(3) - locb(3))**2 )**HALF

                   bcXZHi(l,n) = bcXZHi(l,n) + Gconst * rho(i,j,k) * dV / r

                   if ( doSymmetricAdd ) then

                      bcXZHi(l,n) = bcXZHi(l,n) + &
                           direct_sum_symmetric_add(loc,locb,problo,probhi, &
                           rho(i,j,k),dV,doSymmetricAddLo,doSymmetricAddHi)

                   endif

                enddo

             enddo

             ! Finally, do yz interfaces.

             do m = domlo(2) - 1, domhi(2) + 1
                if     ( m .lt. domlo(2) ) then
                   locb(2) = problo(2) + (dble(m+1)     ) * dx(2)
                elseif ( m .gt. domhi(2) ) then
                   locb(2) = problo(2) + (dble(m  )     ) * dx(2)
                else
                   locb(2) = problo(2) + (dble(m  )+HALF) * dx(2)
                endif

                dy2 = (loc(2) - locb(2))**2

                locb(1) = problo(1)
                dx2 = (loc(1) - locb(1))**2

                do n = domlo(3) - 1, domhi(3) + 1
                   if     ( n .lt. domlo(3) ) then
                      locb(3) = problo(3) + (dble(n+1)     ) * dx(3)
                   elseif ( n .gt. domhi(3) ) then
                      locb(3) = problo(3) + (dble(n  )     ) * dx(3)
                   else
                      locb(3) = problo(3) + (dble(n  )+HALF) * dx(3)
                   endif

                   r = ( dx2 + dy2 + (loc(3) - locb(3))**2 )**HALF

                   bcYZLo(m,n) = bcYZLo(m,n) + Gconst * rho(i,j,k) * dV / r

                   if ( doSymmetricAdd ) then

                      bcYZLo(m,n) = bcYZLo(m,n) + &
                           direct_sum_symmetric_add(loc,locb,problo,probhi, &
                           rho(i,j,k),dV,doSymmetricAddLo,doSymmetricAddHi)
                   endif

                enddo

                locb(1) = probhi(1)
                dx2 = (loc(1) - locb(1))**2

                do n = domlo(3) - 1, domhi(3) + 1
                   if     ( n .lt. domlo(3) ) then
                      locb(3) = problo(3) + (dble(n+1)     ) * dx(3)
                   elseif ( n .gt. domhi(3) ) then
                      locb(3) = problo(3) + (dble(n  )     ) * dx(3)
                   else
                      locb(3) = problo(3) + (dble(n  )+HALF) * dx(3)
                   endif

                   r = ( dx2 + dy2 + (loc(3) - locb(3))**2 )**HALF

                   bcYZHi(m,n) = bcYZHi(m,n) + Gconst * rho(i,j,k) * dV / r

                   if ( doSymmetricAdd ) then

                      bcYZHi(m,n) = bcYZHi(m,n) + &
                           direct_sum_symmetric_add(loc,locb,problo,probhi, &
                           rho(i,j,k),dV,doSymmetricAddLo,doSymmetricAddHi)

                   endif

                enddo

             enddo

          enddo
       enddo
    enddo

  end subroutine ca_compute_direct_sum_bc



  subroutine ca_put_direct_sum_bc (lo,hi,domlo,domhi,phi,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
       bcXYLo,bcXYHi,bcXZLo,bcXZHi,bcYZLo,bcYZHi) &
       bind(C, name="ca_put_direct_sum_bc")
    
    implicit none

    integer          :: lo(3),hi(3)
    integer          :: domlo(3),domhi(3)

    double precision :: bcXYLo(domlo(1)-1:domhi(1)+1,domlo(2)-1:domhi(2)+1)
    double precision :: bcXYHi(domlo(1)-1:domhi(1)+1,domlo(2)-1:domhi(2)+1)
    double precision :: bcXZLo(domlo(1)-1:domhi(1)+1,domlo(3)-1:domhi(3)+1)
    double precision :: bcXZHi(domlo(1)-1:domhi(1)+1,domlo(3)-1:domhi(3)+1)
    double precision :: bcYZLo(domlo(2)-1:domhi(2)+1,domlo(3)-1:domhi(3)+1)
    double precision :: bcYZHi(domlo(2)-1:domhi(2)+1,domlo(3)-1:domhi(3)+1)

    integer          :: p_l1,p_l2,p_l3,p_h1,p_h2,p_h3
    double precision :: phi(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3)

    integer          :: i,j,k

    i = lo(1)
    if (i .eq. domlo(1)-1) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             phi(i,j,k) = bcYZLo(j,k)
          end do
       end do
    end if

    i = hi(1)
    if (i .eq. domhi(1)+1) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             phi(i,j,k) = bcYZHi(j,k)
          end do
       end do
    end if

    j = lo(2)
    if (j .eq. domlo(2)-1) then
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             phi(i,j,k) = bcXZLo(i,k)
          end do
       end do
    end if

    j = hi(2)
    if (j .eq. domhi(2)+1) then
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             phi(i,j,k) = bcXZHi(i,k)
          end do
       end do
    end if

    k = lo(3)
    if (k .eq. domlo(3)-1) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             phi(i,j,k) = bcXYLo(i,j)
          end do
       end do
    end if

    k = hi(3)
    if (k .eq. domhi(3)+1) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             phi(i,j,k) = bcXYHi(i,j)
          end do
       end do
    end if

  end subroutine ca_put_direct_sum_bc



  function direct_sum_symmetric_add (loc,locb,problo,probhi, &
       rho,dV,doSymmetricAddLo,doSymmetricAddHi) result(bcTerm)

    use fundamental_constants_module, only: Gconst
    use bl_constants_module

    implicit none

    double precision :: loc(3), locb(3)
    double precision :: problo(3), probhi(3)
    logical          :: doSymmetricAddLo(3), doSymmetricAddHi(3)

    double precision :: x, y, z, r
    double precision :: rho, dV
    double precision :: bcTerm

    ! Add contributions from any symmetric boundaries.

    bcTerm = ZERO

    if ( doSymmetricAddLo(1) ) then

       x = TWO * problo(1) - loc(1)
       y = loc(2)
       z = loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm + Gconst * rho * dV / r

       if ( doSymmetricAddLo(2) ) then

          x = TWO * problo(1) - loc(1)
          y = TWO * problo(2) - loc(2)
          z = loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm + Gconst * rho * dV / r

       endif

       if ( doSymmetricAddLo(3) ) then

          x = TWO * problo(1) - loc(1)
          y = loc(2)
          z = TWO * problo(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm + Gconst * rho * dV / r

       endif

       if ( doSymmetricAddLo(2) .and. doSymmetricAddLo(3) ) then

          x = TWO * problo(1) - loc(1)
          y = TWO * problo(2) - loc(2)
          z = TWO * problo(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm + Gconst * rho * dV / r

       endif

    endif

    if ( doSymmetricAddLo(2) ) then

       x = loc(1)
       y = TWO * problo(2) - loc(2)
       z = loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm + Gconst * rho * dV / r

       if ( doSymmetricAddLo(3) ) then

          x = loc(1)
          y = TWO * problo(2) - loc(2)
          z = TWO * problo(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm + Gconst * rho * dV / r

       endif

    endif

    if ( doSymmetricAddLo(3) ) then

       x = loc(1)
       y = loc(2)
       z = TWO * problo(3) - loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm + Gconst * rho * dV / r

    endif



    if ( doSymmetricAddHi(1) ) then

       x = TWO * probhi(1) - loc(1)
       y = loc(2)
       z = loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm + Gconst * rho * dV / r

       if ( doSymmetricAddHi(2) ) then

          x = TWO * probhi(1) - loc(1)
          y = TWO * probhi(2) - loc(2)
          z = loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm + Gconst * rho * dV / r

       endif

       if ( doSymmetricAddHi(3) ) then

          x = TWO * probhi(1) - loc(1)
          y = loc(2)
          z = TWO * probhi(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm + Gconst * rho * dV / r

       endif

       if ( doSymmetricAddHi(2) .and. doSymmetricAddHi(3) ) then

          x = TWO * probhi(1) - loc(1)
          y = TWO * probhi(2) - loc(2)
          z = TWO * probhi(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm + Gconst * rho * dV / r

       endif

    endif

    if ( doSymmetricAddHi(2) ) then

       x = loc(1)
       y = TWO * probhi(2) - loc(2)
       z = loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm + Gconst * rho * dV / r

       if ( doSymmetricAddHi(3) ) then

          x = loc(1)
          y = TWO * probhi(2) - loc(2)
          z = TWO * probhi(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm + Gconst * rho * dV / r

       endif

    endif

    if ( doSymmetricAddHi(3) ) then

       x = loc(1)
       y = loc(2)
       z = TWO * probhi(3) - loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm + Gconst * rho * dV / r

    endif

  end function direct_sum_symmetric_add

end module gravity_3D_module
