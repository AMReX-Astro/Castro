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



  subroutine ca_compute_radial_mass (lo,hi,dx,dr,&
       state,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,&
       radial_mass,radial_vol,problo,&
       n1d,drdxfac,level) bind(C, name="ca_compute_radial_mass")

    use amrex_constants_module
    use prob_params_module, only: center
    use meth_params_module, only: NVAR, URHO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    real(rt), intent(in   ) :: dx(3),dr
    real(rt), intent(in   ) :: problo(3)

    integer , intent(in   ) :: n1d,drdxfac,level
    real(rt), intent(inout) :: radial_mass(0:n1d-1)
    real(rt), intent(inout) :: radial_vol (0:n1d-1)

    integer , intent(in   ) :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
    real(rt), intent(in   ) :: state(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3,NVAR)

    integer          :: i,j,k,index
    integer          :: ii,jj,kk
    real(rt)         :: xc,yc,zc,r,xxsq,yysq,zzsq,octant_factor
    real(rt)         :: fac,xx,yy,zz,dx_frac,dy_frac,dz_frac
    real(rt)         :: vol_frac, drinv
    real(rt)         :: lo_i,lo_j,lo_k

    if (( abs(center(1) - problo(1)) .lt. 1.e-2_rt * dx(1) ) .and. &
         ( abs(center(2) - problo(2)) .lt. 1.e-2_rt * dx(2) ) .and. &
         ( abs(center(3) - problo(3)) .lt. 1.e-2_rt * dx(3) ) ) then
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
                   call castro_error("Error:: Gravity_3d.f90 :: ca_compute_radial_mass")
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
                            radial_mass(index) = radial_mass(index) + vol_frac * state(i,j,k,URHO)
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



  subroutine ca_put_radial_phi (lo,hi,domlo,domhi,dx,dr,&
       phi,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
       radial_phi,problo,&
       numpts_1d,fill_interior) bind(C, name="ca_put_radial_phi")

    use amrex_constants_module
    use prob_params_module, only: center

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(3),hi(3)
    integer , intent(in   ) :: domlo(3),domhi(3)
    real(rt), intent(in   ) :: dx(3),dr
    real(rt), intent(in   ) :: problo(3)

    integer , intent(in   ) :: numpts_1d
    real(rt), intent(in   ) :: radial_phi(0:numpts_1d-1)
    integer , intent(in   ) :: fill_interior

    integer , intent(in   ) :: p_l1,p_l2,p_l3,p_h1,p_h2,p_h3
    real(rt), intent(inout) :: phi(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3)

    integer          :: i,j,k,index
    real(rt)         :: x,y,z,r
    real(rt)         :: cen,xi,slope,p_lo,p_md,p_hi,minvar,maxvar
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
                call castro_error("Error:: Gravity_3d.f90 :: ca_put_radial_phi")
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
                        (-p_hi + 26.e0_rt*p_md - p_lo)/24.e0_rt
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



  subroutine ca_compute_direct_sum_bc (lo, hi, dx, &
                                       symmetry_type, lo_bc, hi_bc, &
                                       rho, r_lo, r_hi, &
                                       vol, v_lo, v_hi, &
                                       problo, probhi, &
                                       bcXYLo, bcXYHi, &
                                       bcXZLo, bcXZHi, &
                                       bcYZLo, bcYZHi, &
                                       bclo, bchi, bcdx) bind(C, name="ca_compute_direct_sum_bc")

    use fundamental_constants_module, only: Gconst
    use amrex_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer , intent(in   ) :: lo(3), hi(3)
    integer , intent(in   ) :: bclo(3), bchi(3)
    integer , intent(in   ) :: r_lo(3), r_hi(3)
    integer , intent(in   ) :: v_lo(3), v_hi(3)
    real(rt), intent(in   ) :: dx(3), bcdx(3)
    real(rt), intent(in   ) :: problo(3), probhi(3)

    integer , intent(in   ) :: symmetry_type
    integer , intent(in   ) :: lo_bc(3), hi_bc(3)

    real(rt), intent(out) :: bcXYLo(bclo(1):bchi(1),bclo(2):bchi(2))
    real(rt), intent(out) :: bcXYHi(bclo(1):bchi(1),bclo(2):bchi(2))
    real(rt), intent(out) :: bcXZLo(bclo(1):bchi(1),bclo(3):bchi(3))
    real(rt), intent(out) :: bcXZHi(bclo(1):bchi(1),bclo(3):bchi(3))
    real(rt), intent(out) :: bcYZLo(bclo(2):bchi(2),bclo(3):bchi(3))
    real(rt), intent(out) :: bcYZHi(bclo(2):bchi(2),bclo(3):bchi(3))

    real(rt), intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer          :: i, j, k, l, m, n, b
    real(rt)         :: r
    real(rt)         :: loc(3), locb(3), dx2, dy2, dz2

    logical          :: doSymmetricAddLo(3), doSymmetricAddHi(3), doSymmetricAdd

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

             ! Do xy interfaces first. Note that the boundary conditions
             ! on phi are expected to live directly on the interface.
             ! We also have to handle the domain corners correctly. We are
             ! assuming that bclo = domlo - 1 and bchi = domhi + 1, where
             ! domlo and domhi are the coarse domain extent.

             do m = bclo(2), bchi(2)
                if (m .eq. bclo(2)) then
                   locb(2) = problo(2)
                else if (m .eq. bchi(2)) then
                   locb(2) = probhi(2)
                else
                   locb(2) = problo(2) + (dble(m)+HALF) * bcdx(2)
                endif
                dy2 = (loc(2) - locb(2))**2

                do l = bclo(1), bchi(1)
                   if (l .eq. bclo(1)) then
                      locb(1) = problo(1)
                   else if (l .eq. bchi(1)) then
                      locb(1) = probhi(2)
                   else
                      locb(1) = problo(1) + (dble(l)+HALF) * bcdx(1)
                   endif
                   dx2 = (loc(1) - locb(1))**2

                   locb(3) = problo(3)
                   dz2 = (loc(3) - locb(3))**2
                   
                   r = ( dx2 + dy2 + dz2 )**HALF

                   bcXYLo(l,m) = bcXYLo(l,m) - Gconst * rho(i,j,k) * vol(i,j,k) / r

                   ! Now, add any contributions from mass that is hidden behind
                   ! a symmetric boundary.

                   if ( doSymmetricAdd ) then

                      bcXYLo(l,m) = bcXYLo(l,m) + &
                                    direct_sum_symmetric_add(loc,locb,problo,probhi, &
                                                             rho(i,j,k),vol(i,j,k), &
                                                             doSymmetricAddLo,doSymmetricAddHi)

                   endif

                   locb(3) = probhi(3)
                   dz2 = (loc(3) - locb(3))**2

                   r = ( dx2 + dy2 + dz2 )**HALF

                   bcXYHi(l,m) = bcXYHi(l,m) - Gconst * rho(i,j,k) * vol(i,j,k) / r

                   if ( doSymmetricAdd ) then

                      bcXYHi(l,m) = bcXYHi(l,m) + &
                                    direct_sum_symmetric_add(loc,locb,problo,probhi, &
                                                             rho(i,j,k),vol(i,j,k), &
                                                             doSymmetricAddLo,doSymmetricAddHi)

                   endif

                enddo

             enddo


             ! Now do xz interfaces.

             do n = bclo(3), bchi(3)
                if (n .eq. bclo(3)) then
                   locb(3) = problo(3)
                else if (n .eq. bchi(3)) then
                   locb(3) = probhi(3)
                else
                   locb(3) = problo(3) + (dble(n)+HALF) * bcdx(3)
                endif
                dz2 = (loc(3) - locb(3))**2
                
                do l = bclo(1), bchi(1)
                   if (l .eq. bclo(1)) then
                      locb(1) = problo(1)
                   else if (l .eq. bchi(1)) then
                      locb(1) = probhi(1)
                   else
                      locb(1) = problo(1) + (dble(l)+HALF) * bcdx(1)
                   endif
                   dx2 = (loc(1) - locb(1))**2

                   locb(2) = problo(2)
                   dy2 = (loc(2) - locb(2))**2

                   r = ( dx2 + dy2 + dz2 )**HALF

                   bcXZLo(l,n) = bcXZLo(l,n) - Gconst * rho(i,j,k) * vol(i,j,k) / r

                   if ( doSymmetricAdd ) then

                      bcXZLo(l,n) = bcXZLo(l,n) + &
                                    direct_sum_symmetric_add(loc,locb,problo,probhi, &
                                                             rho(i,j,k),vol(i,j,k), &
                                                             doSymmetricAddLo,doSymmetricAddHi)

                   endif

                   locb(2) = probhi(2)
                   dy2 = (loc(2) - locb(2))**2

                   r = ( dx2 + dy2 + dz2 )**HALF

                   bcXZHi(l,n) = bcXZHi(l,n) - Gconst * rho(i,j,k) * vol(i,j,k) / r

                   if ( doSymmetricAdd ) then

                      bcXZHi(l,n) = bcXZHi(l,n) + &
                                    direct_sum_symmetric_add(loc,locb,problo,probhi, &
                                                             rho(i,j,k),vol(i,j,k), &
                                                             doSymmetricAddLo,doSymmetricAddHi)

                   endif

                enddo

             enddo

             ! Finally, do yz interfaces.

             do n = bclo(3), bchi(3)
                if (n .eq. bclo(3)) then
                   locb(3) = problo(3)
                else if (n .eq. bchi(3)) then
                   locb(3) = probhi(3)
                else
                   locb(3) = problo(3) + (dble(n)+HALF) * bcdx(3)
                endif
                dz2 = (loc(3) - locb(3))**2

                do m = bclo(2), bchi(2)
                   if (m .eq. bclo(2)) then
                      locb(2) = problo(2)
                   else if (m .eq. bchi(2)) then
                      locb(2) = probhi(2)
                   else
                      locb(2) = problo(2) + (dble(m)+HALF) * bcdx(2)
                   endif
                   dy2 = (loc(2) - locb(2))**2

                   locb(1) = problo(1)
                   dx2 = (loc(1) - locb(1))**2

                   r = ( dx2 + dy2 + dz2 )**HALF

                   bcYZLo(m,n) = bcYZLo(m,n) - Gconst * rho(i,j,k) * vol(i,j,k) / r

                   if ( doSymmetricAdd ) then

                      bcYZLo(m,n) = bcYZLo(m,n) + &
                                    direct_sum_symmetric_add(loc,locb,problo,probhi, &
                                                             rho(i,j,k),vol(i,j,k), &
                                                             doSymmetricAddLo,doSymmetricAddHi)
                   endif

                   locb(1) = probhi(1)
                   dx2 = (loc(1) - locb(1))**2

                   r = ( dx2 + dy2 + dz2 )**HALF

                   bcYZHi(m,n) = bcYZHi(m,n) - Gconst * rho(i,j,k) * vol(i,j,k) / r

                   if ( doSymmetricAdd ) then

                      bcYZHi(m,n) = bcYZHi(m,n) + &
                                    direct_sum_symmetric_add(loc,locb,problo,probhi, &
                                                             rho(i,j,k),vol(i,j,k), &
                                                             doSymmetricAddLo,doSymmetricAddHi)

                   endif

                enddo

             enddo

          enddo
       enddo
    enddo

  end subroutine ca_compute_direct_sum_bc



  subroutine ca_put_direct_sum_bc (lo, hi, &
                                   phi, p_lo, p_hi, &
                                   bcXYLo, bcXYHi, &
                                   bcXZLo, bcXZHi, &
                                   bcYZLo, bcYZHi, &
                                   bclo, bchi) bind(C, name="ca_put_direct_sum_bc")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in   ) :: lo(3), hi(3)
    integer, intent(in   ) :: bclo(3), bchi(3)
    integer, intent(in   ) :: p_lo(3), p_hi(3)

    real(rt), intent(in   ) :: bcXYLo(bclo(1):bchi(1),bclo(2):bchi(2))
    real(rt), intent(in   ) :: bcXYHi(bclo(1):bchi(1),bclo(2):bchi(2))
    real(rt), intent(in   ) :: bcXZLo(bclo(1):bchi(1),bclo(3):bchi(3))
    real(rt), intent(in   ) :: bcXZHi(bclo(1):bchi(1),bclo(3):bchi(3))
    real(rt), intent(in   ) :: bcYZLo(bclo(2):bchi(2),bclo(3):bchi(3))
    real(rt), intent(in   ) :: bcYZHi(bclo(2):bchi(2),bclo(3):bchi(3))

    real(rt), intent(inout) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))

    integer          :: i, j, k

    ! Note that we assume phi has one ghost zone relative to the domain,
    ! and since we use a grown box for this loop, so the boxes that lie
    ! on the perimeter of the domain will indeed have lo = domlo - 1.

    i = lo(1)
    if (i .eq. bclo(3)) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             phi(i,j,k) = bcYZLo(j,k)
          end do
       end do
    end if

    i = hi(1)
    if (i .eq. bchi(3)) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             phi(i,j,k) = bcYZHi(j,k)
          end do
       end do
    end if

    j = lo(2)
    if (j .eq. bclo(2)) then
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             phi(i,j,k) = bcXZLo(i,k)
          end do
       end do
    end if

    j = hi(2)
    if (j .eq. bchi(2)) then
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             phi(i,j,k) = bcXZHi(i,k)
          end do
       end do
    end if

    k = lo(3)
    if (k .eq. bclo(3)) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             phi(i,j,k) = bcXYLo(i,j)
          end do
       end do
    end if

    k = hi(3)
    if (k .eq. bchi(3)) then
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
    use amrex_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: loc(3), locb(3)
    real(rt)         :: problo(3), probhi(3)
    logical          :: doSymmetricAddLo(3), doSymmetricAddHi(3)

    real(rt)         :: x, y, z, r
    real(rt)         :: rho, dV
    real(rt)         :: bcTerm

    ! Add contributions from any symmetric boundaries.

    bcTerm = ZERO

    if ( doSymmetricAddLo(1) ) then

       x = TWO * problo(1) - loc(1)
       y = loc(2)
       z = loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm - Gconst * rho * dV / r

       if ( doSymmetricAddLo(2) ) then

          x = TWO * problo(1) - loc(1)
          y = TWO * problo(2) - loc(2)
          z = loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm - Gconst * rho * dV / r

       endif

       if ( doSymmetricAddLo(3) ) then

          x = TWO * problo(1) - loc(1)
          y = loc(2)
          z = TWO * problo(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm - Gconst * rho * dV / r

       endif

       if ( doSymmetricAddLo(2) .and. doSymmetricAddLo(3) ) then

          x = TWO * problo(1) - loc(1)
          y = TWO * problo(2) - loc(2)
          z = TWO * problo(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm - Gconst * rho * dV / r

       endif

    endif

    if ( doSymmetricAddLo(2) ) then

       x = loc(1)
       y = TWO * problo(2) - loc(2)
       z = loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm - Gconst * rho * dV / r

       if ( doSymmetricAddLo(3) ) then

          x = loc(1)
          y = TWO * problo(2) - loc(2)
          z = TWO * problo(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm - Gconst * rho * dV / r

       endif

    endif

    if ( doSymmetricAddLo(3) ) then

       x = loc(1)
       y = loc(2)
       z = TWO * problo(3) - loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm - Gconst * rho * dV / r

    endif



    if ( doSymmetricAddHi(1) ) then

       x = TWO * probhi(1) - loc(1)
       y = loc(2)
       z = loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm - Gconst * rho * dV / r

       if ( doSymmetricAddHi(2) ) then

          x = TWO * probhi(1) - loc(1)
          y = TWO * probhi(2) - loc(2)
          z = loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm - Gconst * rho * dV / r

       endif

       if ( doSymmetricAddHi(3) ) then

          x = TWO * probhi(1) - loc(1)
          y = loc(2)
          z = TWO * probhi(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm - Gconst * rho * dV / r

       endif

       if ( doSymmetricAddHi(2) .and. doSymmetricAddHi(3) ) then

          x = TWO * probhi(1) - loc(1)
          y = TWO * probhi(2) - loc(2)
          z = TWO * probhi(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm - Gconst * rho * dV / r

       endif

    endif

    if ( doSymmetricAddHi(2) ) then

       x = loc(1)
       y = TWO * probhi(2) - loc(2)
       z = loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm - Gconst * rho * dV / r

       if ( doSymmetricAddHi(3) ) then

          x = loc(1)
          y = TWO * probhi(2) - loc(2)
          z = TWO * probhi(3) - loc(3)

          r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

          bcTerm = bcTerm - Gconst * rho * dV / r

       endif

    endif

    if ( doSymmetricAddHi(3) ) then

       x = loc(1)
       y = loc(2)
       z = TWO * probhi(3) - loc(3)

       r = ( (x - locb(1))**2 + (y - locb(2))**2 + (z - locb(3))**2 )**HALF

       bcTerm = bcTerm - Gconst * rho * dV / r

    endif

  end function direct_sum_symmetric_add

end module gravity_3D_module
