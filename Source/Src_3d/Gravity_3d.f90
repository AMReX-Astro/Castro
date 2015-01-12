! :: ----------------------------------------------------------
! :: Average the fine grid phi onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! :: Note this differs from ca_avgdown in that there is no volume weighting.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  rfine      => (ignore) used in 2-D RZ calc
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
      subroutine ca_avgdown_phi (crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3, &
                                 fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                                 lo,hi,lrat)

      use bl_constants_module

      implicit none

      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer lo(3), hi(3)
      integer lrat(3)
      double precision crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3)
      double precision fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)

      integer i, j, k, ic, jc, kc, ioff, joff, koff
      integer lratx, lraty, lratz
      double precision volfrac

      lratx   = lrat(1)
      lraty   = lrat(2)
      lratz   = lrat(3)
      volfrac = ONE/float(lrat(1)*lrat(2)*lrat(3))
      !
      ! ::::: set coarse grid to zero on overlap
      !
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,kc) = ZERO
            enddo
         enddo
      enddo
      !
      ! ::::: sum fine data
      !
      do koff = 0, lratz-1
        do kc = lo(3), hi(3)
          k = kc*lratz + koff
          do joff = 0, lraty-1
            do jc = lo(2), hi(2)
              j = jc*lraty + joff
              do ioff = 0, lratx-1
                do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  crse(ic,jc,kc) = crse(ic,jc,kc) + fine(i,j,k)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,kc) = volfrac*crse(ic,jc,kc)
            enddo
         enddo
      enddo

      end subroutine ca_avgdown_phi

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_edge_interp(flo, fhi, nc, ratio, dir, &
           fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)

      implicit none
      integer flo(0:3-1), fhi(0:3-1), nc, ratio(0:3-1), dir
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      double precision &
           fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,n,P,M,L
      double precision val, df
      !
      ! Do linear in dir, pc transverse to dir, leave alone the fine values
      ! lining up with coarse edges--assume these have been set to hold the 
      ! values you want to interpolate to the rest.
      !
      if (dir.eq.0) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                     df = fine(i+ratio(dir),j,k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) &
                             + df*dble(M)/dble(ratio(dir))
                        do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(i+M,P,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1)-ratio(dir),ratio(1)
                  do i=flo(0),fhi(0)
                     df = fine(i,j+ratio(dir),k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) &
                             + df*dble(M)/dble(ratio(dir))
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(P,j+M,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=flo(2),fhi(2)-ratio(dir),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0),ratio(0)
                     df = fine(i,j,k+ratio(dir),n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) &
                             + df*dble(M)/dble(ratio(dir))
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                              fine(P,L,k+M,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      endif

      end subroutine ca_edge_interp

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_pc_edge_interp(lo, hi, nc, ratio, dir, &
           crse, crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2, &
           fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)
      implicit none
      integer lo(3),hi(3), nc, ratio(0:3-1), dir
      integer crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      double precision &
           crse(crse_l0:crse_h0,crse_l1:crse_h1,crse_l2:crse_h2,nc)
      double precision &
           fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,ii,jj,kk,n,L, P
      !
      ! For edge-based data, fill fine values with piecewise-constant interp of coarse data.
      ! Operate only on faces that overlap--ie, only fill the fine faces that make up each
      ! coarse face, leave the in-between faces alone.
      !
      if (dir.eq.0) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(1)-1
                           fine(ii,jj+L,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(1)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj+P,kk,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      end subroutine ca_pc_edge_interp

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_avg_ec_to_cc(lo, hi, bc_lo, bc_hi, &
           symmetry_type, &
           cc, ccl1, ccl2, ccl3, cch1, cch2, cch3, &
           ecx, ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3, &
           ecy, ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3, &
           ecz, eczl1, eczl2, eczl3, eczh1, eczh2, eczh3, &
           dx, problo, coord_type)

      use bl_constants_module

      implicit none
      integer          :: lo(3),hi(3)
      integer          :: bc_lo(3),bc_hi(3)
      integer          :: symmetry_type, coord_type
      integer          :: ccl1, ccl2, ccl3, cch1, cch2, cch3
      integer          :: ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3
      integer          :: ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3
      integer          :: eczl1, eczl2, eczl3, eczh1, eczh2, eczh3
      double precision :: cc(ccl1:cch1,ccl2:cch2,ccl3:cch3, 3)
      double precision :: ecx(ecxl1:ecxh1,ecxl2:ecxh2, ecxl3:ecxh3)
      double precision :: ecy(ecyl1:ecyh1,ecyl2:ecyh2, ecyl3:ecyh3)
      double precision :: ecz(eczl1:eczh1,eczl2:eczh2, eczl3:eczh3)
      double precision :: dx(3), problo(3)

      ! Local variables
      integer          :: i,j,k

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,k,1) = HALF * ( ecx(i+1,j,k) + ecx(i,j,k) )
               cc(i,j,k,2) = HALF * ( ecy(i,j+1,k) + ecy(i,j,k) )
               cc(i,j,k,3) = HALF * ( ecz(i,j,k+1) + ecz(i,j,k) )
            enddo
         enddo
      enddo
         
      end subroutine ca_avg_ec_to_cc

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
 
      subroutine ca_test_residual(lo, hi, &
           rhs, rhl1, rhl2, rhl3, rhh1, rhh2, rhh3, &
           ecx, ecxl1, ecxl2, ecxl3, ecxh1, ecxh2, ecxh3, &
           ecy, ecyl1, ecyl2, ecyl3, ecyh1, ecyh2, ecyh3, &
           ecz, eczl1, eczl2, eczl3, eczh1, eczh2, eczh3, &
           dx,problo,coord_type)

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

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_average_ec ( &
           f, fl1, fl2, fl3, fh1, fh2, fh3, &
           c, cl1, cl2, cl3, ch1, ch2, ch3, &
           lo, hi, rr, idir)

      use bl_constants_module

      implicit none

      integer lo(3),hi(3)
      integer fl1, fl2, fl3, fh1, fh2, fh3
      integer cl1, cl2, cl3, ch1, ch2, ch3
      double precision f(fl1:fh1,fl2:fh2,fl3:fh3)
      double precision c(cl1:ch1,cl2:ch2,cl3:ch3)
      integer rr(3), idir

      integer i,j,k,n,m,facx,facy,facz

      facx = rr(1)
      facy = rr(2)
      facz = rr(3)

      if (idir .eq. 0) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j,k) = ZERO
                  do n = 0,facy-1
                     do m = 0,facz-1
                        c(i,j,k) = c(i,j,k) + f(facx*i,facy*j+n,facz*k+m)
                     end do
                  end do
                  c(i,j,k) = c(i,j,k) / (facy*facz)
               end do
            end do
         end do
      else if (idir .eq. 1) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j,k) = ZERO
                  do n = 0,facx-1
                     do m = 0,facz-1
                        c(i,j,k) = c(i,j,k) + f(facx*i+n,facy*j,facz*k+m)
                     end do
                  end do
                  c(i,j,k) = c(i,j,k) / (facx*facz)
               end do
            end do
         end do
      else
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j,k) = ZERO
                  do n = 0,facx-1
                     do m = 0,facy-1
                        c(i,j,k) = c(i,j,k) + f(facx*i+n,facy*j+m,facz*k)
                     end do
                  end do
                  c(i,j,k) = c(i,j,k) / (facx*facy)
               end do
            end do
         end do
      end if

      end subroutine ca_average_ec

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_compute_radial_mass (lo,hi,dx,dr,&
                                         rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,&
                                         radial_mass,radial_vol,problo,&
                                         n1d,drdxfac,level)
      use bl_constants_module
      use probdata_module

      implicit none

      integer          :: lo(3),hi(3)
      double precision :: dx(3),dr
      double precision :: problo(3)

      integer          :: n1d,drdxfac,level
      double precision :: radial_mass(0:n1d-1)
      double precision :: radial_vol (0:n1d-1)

      integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
      double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)

      integer          :: i,j,k,index,tid,nthreads
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

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_put_radial_grav (lo,hi,dx,dr,&
                                     grav,g_l1,g_l2,g_l3,g_h1,g_h2,g_h3, &
                                     radial_grav,problo,n1d,level)

      use bl_constants_module
      use probdata_module

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

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_put_radial_phi (lo,hi,domlo,domhi,dx,dr,&
                                    phi,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
                                    radial_phi,problo,&
                                    numpts_1d,fill_interior)
        use bl_constants_module
        use probdata_module

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

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_put_multipole_bc (lo,hi,domlo,domhi,dx,&
                                      phi,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
                                      problo,probhi,lnum,q0,qC,qS)
        use probdata_module
        use fundamental_constants_module, only: Gconst
        use bl_constants_module

        implicit none

        integer          :: lo(3),hi(3)
        integer          :: domlo(3),domhi(3)
        double precision :: dx(3), dV
        double precision :: problo(3),probhi(3)

        integer          :: lnum
        double precision :: q0(0:lnum), qC(0:lnum,0:lnum), qS(0:lnum,0:lnum)

        integer          :: p_l1,p_l2,p_l3,p_h1,p_h2,p_h3
        double precision :: phi(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3)

        integer          :: i,j,k
        integer          :: l,m
        double precision :: x,y,z,r,cosTheta,phiAngle
        double precision :: rmax
        double precision :: legPolyArr(0:lnum), assocLegPolyArr(0:lnum,0:lnum)
        double precision :: r_to_mlm1

        dV = dx(1) * dx(2) * dx(3)

        rmax = probhi(1)
        if ( probhi(2) > rmax ) then
          rmax = probhi(2)
        endif
        if ( probhi(3) > rmax ) then
          rmax = probhi(3)
        endif

        rmax = rmax * sqrt(THREE) / TWO
        
        !$OMP PARALLEL DO PRIVATE(i,j,k,l,m,x,y,z,r,cosTheta,phiAngle,legPolyArr,assocLegPolyArr,r_to_mlm1)
        do k = p_l3,p_h3
           if (k .gt. domhi(3)) then
              z = problo(3) + (dble(k  )     ) * dx(3) - center(3)
           else if (k .lt. domlo(3)) then
              z = problo(3) + (dble(k+1)     ) * dx(3) - center(3)
           else 
              z = problo(3) + (dble(k  )+HALF) * dx(3) - center(3)
           end if

           z = z / rmax

           do j = p_l2,p_h2
              if (j .gt. domhi(2)) then
                 y = problo(2) + (dble(j  )     ) * dx(2) - center(2)
              else if (j .lt. domlo(2)) then
                 y = problo(2) + (dble(j+1)     ) * dx(2) - center(2)
              else 
                 y = problo(2) + (dble(j  )+HALF) * dx(2) - center(2)
              end if

              y = y / rmax

              do i = p_l1,p_h1
                 if (i .gt. domhi(1)) then
                    x = problo(1) + (dble(i  )     ) * dx(1) - center(1)
                 else if (i .lt. domlo(1)) then
                    x = problo(1) + (dble(i+1)     ) * dx(1) - center(1)
                 else 
                    x = problo(1) + (dble(i  )+HALF) * dx(1) - center(1)
                 end if

                 x = x / rmax
        
                 ! Only adjust ghost zones here

                 if ( i .lt. domlo(1) .or. i .gt. domhi(1) .or. &
                      j .lt. domlo(2) .or. j .gt. domhi(2) .or. &
                      k .lt. domlo(3) .or. k .gt. domhi(3) ) then

                   ! There are some cases where r == 0. This might occur, for example,
                   ! when we have symmetric BCs and our corner is at one edge.
                   ! In this case, we'll set phi to zero for safety, to avoid NaN issues.
                   ! These cells should not be called anyway in the gravity BCs.

                   r = sqrt( x**2 + y**2 + z**2 )

                   if ( r < 1.0d-12 ) then
                     phi(i,j,k) = ZERO
                     cycle
                   endif

                   cosTheta = z / r
                   phiAngle = atan2(y,x)

                   phi(i,j,k) = ZERO

                   ! First, calculate the Legendre polynomials.

                   legPolyArr(:) = ZERO
                   assocLegPolyArr(:,:) = ZERO

                   call fill_legendre_arrays(legPolyArr, assocLegPolyArr, cosTheta, lnum)

                   ! Now compute the potentials on the ghost cells.

                   do l = 0, lnum

                     r_to_mlm1 = r**(-l-1)
 
                     phi(i,j,k) = phi(i,j,k) + q0(l) * legPolyArr(l) * r_to_mlm1

                     do m = 1, l
       
                       phi(i,j,k) = phi(i,j,k) + (qC(l,m) * cos(m * phiAngle) + qS(l,m) * sin(m * phiAngle)) * &
                                                 assocLegPolyArr(l,m) * r_to_mlm1

                     enddo

                   enddo

                   phi(i,j,k) = Gconst * phi(i,j,k) * dV / rmax

                 endif

              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO

      end subroutine ca_put_multipole_bc

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_compute_multipole_moments (lo,hi,domlo,domhi,symmetry_type,lo_bc,hi_bc,&
                                               dx,rho,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3,&
                                               problo,probhi,lnum,q0,qC,qS)
        use probdata_module
        use bl_constants_module
        use meth_params_module, only : deterministic

        implicit none

        integer          :: lo(3),hi(3)
        integer          :: lo_bc(3),hi_bc(3)
        integer          :: symmetry_type
        integer          :: domlo(3),domhi(3)
        double precision :: dx(3)
        double precision :: problo(3), probhi(3)

        integer          :: lnum
        double precision :: q0(0:lnum), qC(0:lnum,0:lnum), qS(0:lnum,0:lnum)
        double precision :: factArray(0:lnum,0:lnum)

        integer          :: p_l1,p_l2,p_l3,p_h1,p_h2,p_h3
        double precision :: rho(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3)

        integer          :: i,j,k
        integer          :: ilo, ihi, jlo, jhi, klo, khi
        integer          :: l,m,b

        double precision :: factorial

        double precision :: x,y,z,r,cosTheta,phiAngle

        double precision :: volumeFactor, parityFactor
        double precision, parameter :: edgeTolerance = 1.0d-2
        double precision :: rmax
        double precision :: legPolyArr(0:lnum), assocLegPolyArr(0:lnum,0:lnum)
        double precision :: rho_r_to_l
        double precision :: parity_q0(0:lnum), parity_qC_qS(0:lnum,0:lnum)

        logical          :: doSymmetricAddLo(3), doSymmetricAddHi(3), doSymmetricAdd
        logical          :: doReflectionLo(3), doReflectionHi(3)

        ! If any of the boundaries are symmetric, we need to account for the mass that is assumed
        ! to lie on the opposite side of the symmetric axis. If the center in any direction 
        ! coincides with the boundary, then we can simply double the mass as a result of that reflection.
        ! Otherwise, we need to do a more general solve. We include a logical that is set to true
        ! if any boundary is symmetric, so that we can avoid unnecessary function calls.

        volumeFactor = ONE
        parityFactor = ONE

        doSymmetricAddLo(:) = .false.
        doSymmetricAddHi(:) = .false.

        doSymmetricAdd      = .false.

        doReflectionLo(:)   = .false.
        doReflectionHi(:)   = .false.

        do b = 1, 3

          if ( (lo_bc(b) .eq. symmetry_type) ) then
            if ( abs(center(b) - problo(b)) < edgeTolerance ) then
              volumeFactor = volumeFactor * TWO
              doReflectionLo(b) = .true.
            else
              doSymmetricAddLo(b) = .true.
              doSymmetricAdd      = .true.
            endif
          endif

          if ( (hi_bc(b) .eq. symmetry_type) ) then
            if ( abs(center(b) - probhi(b)) < edgeTolerance ) then
              volumeFactor = volumeFactor * TWO
              doReflectionHi(b) = .true.
            else
              doSymmetricAddHi(b) = .true.
              doSymmetricAdd      = .true.
            endif
          endif

        enddo 


        ! Compute pre-factors now to save computation time, for qC and qS

        do l = 0, lnum

          ! The odd l Legendre polynomials are odd in their argument, so
          ! a symmetric reflection about the z axis leads to a total cancellation.

          parity_q0(l) = ONE

          if ( MODULO(l,2) /= 0 .and. ( doReflectionLo(3) .or. doReflectionHi(3) ) ) then
            parity_q0(l) = ZERO
          endif

          do m = 1, l

            ! The parity properties of the associated Legendre polynomials are:
            ! P_l^m (-x) = (-1)^(l+m) P_l^m (x)
            ! Therefore, a complete cancellation occurs if l+m is odd and
            ! we are reflecting about the z axis.

            ! Additionally, the cosine and sine terms flip sign when reflected
            ! about the x or y axis, so if we have a reflection about x or y
            ! then the terms have a complete cancellation.

            parity_qC_qS(l,m) = ONE

            if ( MODULO(l+m,2) /= 0 .and. ( doReflectionLo(3) .or. doReflectionHi(3) ) ) then
              parity_qC_qS(l,m) = ZERO
            endif

            if ( doReflectionLo(1) .or. doReflectionLo(2) .or. doReflectionHi(1) .or. doReflectionHi(2) ) then
              parity_qC_qS(l,m) = ZERO
            endif

            factArray(l,m) = TWO * factorial(l-m) / factorial(l+m) * volumeFactor

          enddo

        enddo


        ! Don't add to multipole moments if we're outside the physical domain

        if ( p_l3 .lt. domlo(3) ) then
          klo = domlo(3)
        else
          klo = p_l3
        endif

        if ( p_h3 .gt. domhi(3) ) then
          khi = domhi(3)
        else
          khi = p_h3
        endif

        if ( p_l2 .lt. domlo(2) ) then
          jlo = domlo(2)
        else
          jlo = p_l2
        endif

        if ( p_h2 .gt. domhi(2) ) then
          jhi = domhi(2)
        else
          jhi = p_h2
        endif

        if ( p_l1 .lt. domlo(1) ) then
          ilo = domlo(1)
        else
          ilo = p_l1
        endif

        if ( p_h1 .gt. domhi(1) ) then
          ihi = domhi(1)
        else
          ihi = p_h1
        endif

        ! Now let's take care of a safety issue. The multipole calculation involves taking powers of r^l, 
        ! which can overflow the double precision exponent limit if lnum is very large. Therefore,
        ! we will normalize all distances to the maximum possible physical distance from the center,
        ! which is the diagonal from the center to the edge of the box. Then r^l will always be
        ! less than or equal to one. For large enough lnum, this will still result in roundoff
        ! errors that don't make your answer any more precise, but at least it avoids
        ! possible NaN issues from having numbers that are too large for double precision.
        ! We will put the rmax factor back in at the end of ca_put_multipole_bc.

        ! Another constraint is that the method used to calculate the Legendre polynomials
        ! is numerically unstable for large l, just because of how big they get, so
        ! for lnum > 50 we'll throw an error (this number was determined by experiment).

        rmax = probhi(1)
        if ( probhi(2) > rmax ) then
          rmax = probhi(2)
        endif
        if ( probhi(3) > rmax ) then
          rmax = probhi(3)
        endif

        rmax = rmax * sqrt(THREE) / TWO ! This is the distance from the center to the corner of a cube. 

        if ( lnum > 50 ) then
          print *, ">>> CA_COMPUTE_MULTIPOLE_MOMENTS: The value of l you have chosen is too large."
          print *, ">>> Try again with max_multipole_order <= 50."
          call bl_error("Error: Gravity_3d.f90: ca_compute_multipole_moments")
        endif

        !$OMP PARALLEL DO PRIVATE(i,j,k,l,m,legPolyArr,assocLegPolyArr) &
        !$OMP PRIVATE(x,y,z,r,cosTheta,phiAngle,parityFactor,rho_r_to_l) &
        !$OMP REDUCTION(+:q0,qC,qS) if (.not.deterministic)
        do k = klo, khi
           z = ( problo(3) + (dble(k)+HALF) * dx(3) - center(3) ) / rmax

           do j = jlo, jhi
              y = ( problo(2) + (dble(j)+HALF) * dx(2) - center(2) ) / rmax
 
              do i = ilo, ihi
                 x = ( problo(1) + (dble(i)+HALF) * dx(1) - center(1) ) / rmax

                 r = sqrt( x**2 + y**2 + z**2 )
                 cosTheta = z / r
                 phiAngle = atan2(y,x)

                 ! Compute contribution from this cell to all multipole moments.
                 ! First, calculate the Legendre polynomials.

                 call fill_legendre_arrays(legPolyArr, assocLegPolyArr, cosTheta, lnum)

                 ! Absorb the factorial terms into the Legendre polynomials to save on multiplications later.

                 assocLegPolyArr = assocLegPolyArr * factArray

                 ! Now, compute the multipole moments using the tabulated polynomials.

                 do l = 0, lnum 
                 
                   rho_r_to_l = rho(i,j,k) * (r**dble(l))

                   q0(l) = q0(l) + legPolyArr(l) * rho_r_to_l * volumeFactor * parity_q0(l)

                   do m = 1, l

                     qC(l,m) = qC(l,m) + assocLegPolyArr(l,m) * cos(m * phiAngle) * &
                                         rho_r_to_l * parity_qC_qS(l,m)

                     qS(l,m) = qS(l,m) + assocLegPolyArr(l,m) * sin(m * phiAngle) * &
                                         rho_r_to_l * parity_qC_qS(l,m)

                   enddo

                 enddo

                 ! Now add in contributions if we have any symmetric boundaries

                 if ( doSymmetricAdd ) then

                   call multipole_symmetric_add(doSymmetricAddLo, doSymmetricAddHi, &
                        x, y, z, problo, probhi, rmax, &
                        rho(i,j,k), factArray, &
                        q0, qC, qS, lnum)

                 endif

              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO

      end subroutine ca_compute_multipole_moments


! ::
! :: ----------------------------------------------------------
! ::

      double precision function factorial(n)

        use bl_constants_module      

        implicit none

        integer :: n, i

        factorial = ONE
 
        do i = 2, n
          factorial = factorial * dble(i)
        enddo

      end function factorial

! ::
! :: ----------------------------------------------------------
! ::

      subroutine fill_legendre_arrays(legPolyArr, assocLegPolyArr, x, lnum)

        use bl_constants_module
      
        implicit none

        integer :: lnum
        integer :: l, m, n
        double precision :: x
        double precision :: legPolyArr(0:lnum), assocLegPolyArr(0:lnum,0:lnum)

        legPolyArr(:)        = ZERO
        assocLegPolyArr(:,:) = ZERO

        ! First we'll do the associated Legendre polynomials. There are a number of
        ! recurrence relations, but many are unstable. We'll use one that is known
        ! to be stable for the reasonably low values of l we care about in a simulation:
        ! (l-m)P_l^m(x) = x(2l-1)P_{l-1}^m(x) - (l+m-1)P_{l-2}^m(x).
        ! This uses the following two expressions as initial conditions:
        ! P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2)
        ! P_{m+1}^m(x) = x (2m+1) P_m^m (x)

        do m = 1, lnum

          ! P_m^m

          assocLegPolyArr(m,m) = (-1)**m * ( (ONE - x) * (ONE + x) )**(dble(m)/TWO)

          ! Multiply by the double factorial term

          do n = (2*m-1), 3, -2

            assocLegPolyArr(m,m) = assocLegPolyArr(m,m) * n

          enddo

          ! P_{m+1}^m
          ! We need to be careful at m = lnum, because we could reference a non-existent array subscript.

          if ( m < lnum ) then

            assocLegPolyArr(m+1,m) = x * (2*m + 1) * assocLegPolyArr(m,m)

          endif

          ! All other l
          ! The loop will automatically be avoided if we're close to lnum here.

          do l = m+2, lnum

            assocLegPolyArr(l,m) = (x * (2*l - 1) * assocLegPolyArr(l-1,m) - (l + m - 1) * assocLegPolyArr(l-2,m) ) / (l-m)

          enddo

        enddo


        ! Now we'll do the normal Legendre polynomials. We use a stable recurrence relation:
        ! (l+1) P_{l+1}(x) = (2l+1) x P_l(x) - l P_{l-1}(x). This uses initial conditions:
        ! P_0(x) = 1
        ! P_1(x) = x

        do l = 0, lnum

          if ( l == 0 ) then

            legPolyArr(0) = ONE

          elseif ( l == 1 ) then

            legPolyArr(1) = x

          else

            legPolyArr(l) = ( (2*l - 1) * x * legPolyArr(l-1) - (l-1) * legPolyArr(l-2) ) / l

          endif

        enddo

      end subroutine fill_legendre_arrays

! ::
! :: ----------------------------------------------------------
! ::

      subroutine multipole_symmetric_add(doSymmetricAddLo, doSymmetricAddHi, &
                      x, y, z, problo, probhi, rmax, &
                      rho, factArray, &
                      q0, qC, qS, lnum)

        use probdata_module
        use bl_constants_module

        implicit none

        integer,          intent(in) :: lnum
        double precision, intent(in) :: factArray(0:lnum,0:lnum)
        double precision, intent(in) :: x, y, z
        double precision, intent(in) :: problo(3), probhi(3), rmax
        double precision, intent(in) :: rho

        logical,          intent(in) :: doSymmetricAddLo(3), doSymmetricAddHi(3)

        double precision, intent(inout) :: q0(0:lnum)
        double precision, intent(inout) :: qC(0:lnum,0:lnum), qS(0:lnum,0:lnum)

        double precision :: cosTheta, phiAngle, r
        double precision :: xLo, yLo, zLo, xHi, yHi, zHi

        xLo = ( TWO * (problo(1) - center(1)) ) / rmax - x
        xHi = ( TWO * (probhi(1) - center(1)) ) / rmax - x

        yLo = ( TWO * (problo(2) - center(2)) ) / rmax - y
        yHi = ( TWO * (probhi(2) - center(2)) ) / rmax - y

        zLo = ( TWO * (problo(3) - center(3)) ) / rmax - z
        zHi = ( TWO * (probhi(3) - center(3)) ) / rmax - z

        if ( doSymmetricAddLo(1) ) then

          r        = sqrt( xLo**2 + y**2 + z**2 )
          phiAngle = atan2(y, xLo)
          cosTheta = z / r

          call multipole_add(cosTheta, phiAngle, r, rho, factArray, q0, qC, qS, lnum)

          if ( doSymmetricAddLo(2) ) then

            r        = sqrt( xLo**2 + yLo**2 + z**2 )
            phiAngle = atan2(yLo, xLo)
            cosTheta = z / r

            call multipole_add(cosTheta, phiAngle, r, rho, factArray, q0, qC, qS, lnum)

          endif

          if ( doSymmetricAddLo(3) ) then

            r        = sqrt( xLo**2 + y**2 + zLo**2 )
            phiAngle = atan2(y, xLo)
            cosTheta = zLo / r

            call multipole_add(cosTheta, phiAngle, r, rho, factArray, q0, qC, qS, lnum)

          endif

          if ( doSymmetricAddLo(2) .and. doSymmetricAddLo(3) ) then

            r        = sqrt( xLo**2 + yLo**2 + zLo**2 )
            phiAngle = atan2(yLo, xLo)
            cosTheta = zLo / r

            call multipole_add(cosTheta, phiAngle, r, rho, factArray, q0, qC, qS, lnum)

          endif

        endif

        if ( doSymmetricAddLo(2) ) then

          r        = sqrt( x**2 + yLo**2 + z**2 )
          phiAngle = atan2(yLo, x)
          cosTheta = z / r

          call multipole_add(cosTheta, phiAngle, r, rho, factArray, q0, qC, qS, lnum)

          if ( doSymmetricAddLo(3) ) then

            r        = sqrt( x**2 + yLo**2 + zLo**2 )
            phiAngle = atan2(yLo, x)
            cosTheta = zLo / r

            call multipole_add(cosTheta, phiAngle, r, rho, factArray, q0, qC, qS, lnum)

          endif

        endif

        if ( doSymmetricAddLo(3) ) then
 
          r        = sqrt( x**2 + y**2 + zLo**2 )
          phiAngle = atan2(y, x)
          cosTheta = zLo / r

          call multipole_add(cosTheta, phiAngle, r, rho, factArray, q0, qC, qS, lnum)

        endif

      end subroutine multipole_symmetric_add

! ::
! :: ----------------------------------------------------------
! ::

      subroutine multipole_add(cosTheta, phiAngle, r, rho, factArray, q0, qC, qS, lnum)

        implicit none

        integer,          intent(in) :: lnum
        double precision, intent(in) :: cosTheta, phiAngle, r, rho, factArray(0:lnum,0:lnum)

        double precision, intent(inout) :: q0(0:lnum), qC(0:lnum,0:lnum), qS(0:lnum,0:lnum)

        integer :: l, m

        double precision :: legPolyArr(0:lnum), assocLegPolyArr(0:lnum,0:lnum)

        double precision :: rho_r_to_l

        call fill_legendre_arrays(legPolyArr, assocLegPolyArr, cosTheta, lnum)

        ! Absorb factorial terms into associated Legendre polynomials

        assocLegPolyArr = assocLegPolyArr * factArray

        do l = 0, lnum

          rho_r_to_l = rho * (r ** l)

          q0(l) = q0(l) + legPolyArr(l) * rho_r_to_l

          do m = 1, l
            
            qC(l,m) = qC(l,m) + assocLegPolyArr(l,m) * cos(m * phiAngle) * rho_r_to_l
            qS(l,m) = qS(l,m) + assocLegPolyArr(l,m) * sin(m * phiAngle) * rho_r_to_l

          enddo

        enddo

      end subroutine multipole_add

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_compute_direct_sum_bc (lo,hi,domlo,domhi,&
                                           symmetry_type,lo_bc,hi_bc, &
                                           dx,rho,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
                                           problo, probhi, &
                                           bcXYLo,bcXYHi,bcXZLo,bcXZHi,bcYZLo,bcYZHi)
        use probdata_module
        use fundamental_constants_module, only: Gconst
        use bl_constants_module
        use meth_params_module, only : deterministic

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
        integer          :: ilo,ihi,jlo,jhi,klo,khi
        double precision :: r
        double precision :: loc(3), locb(3), dx2, dy2, dz2

        logical          :: doSymmetricAddLo(3), doSymmetricAddHi(3), doSymmetricAdd
        double precision :: direct_sum_symmetric_add
        
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

        ! We only contribute to the BCs if we're inside the physical domain.
        
        if ( p_l3 .lt. domlo(3) ) then
          klo = domlo(3)
        else
          klo = p_l3
        endif

        if ( p_h3 .gt. domhi(3) ) then
          khi = domhi(3)
        else
          khi = p_h3
        endif

        if ( p_l2 .lt. domlo(2) ) then
          jlo = domlo(2)
        else
          jlo = p_l2
        endif

        if ( p_h2 .gt. domhi(2) ) then
          jhi = domhi(2)
        else
          jhi = p_h2
        endif

        if ( p_l1 .lt. domlo(1) ) then
          ilo = domlo(1)
        else
          ilo = p_l1
        endif

        if ( p_h1 .gt. domhi(1) ) then
          ihi = domhi(1)
        else
          ihi = p_h1
        endif

 
        !$OMP PARALLEL DO PRIVATE(i,j,k,loc,locb,dx2,dy2,dz2,r,l,m,n) &
        !$OMP REDUCTION(+:bcXYLo,bcXYHi,bcXZLo,bcXZHi,bcYZLo,bcYZHi) if(.not.deterministic)
        do k = klo, khi
           loc(3) = problo(3) + (dble(k)+HALF) * dx(3)

           do j = jlo, jhi
              loc(2) = problo(2) + (dble(j)+HALF) * dx(2)

              do i = ilo, ihi
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
        !$OMP END PARALLEL DO

      end subroutine ca_compute_direct_sum_bc

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_put_direct_sum_bc (lo,hi,domlo,domhi,phi,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
                                       bcXYLo,bcXYHi,bcXZLo,bcXZHi,bcYZLo,bcYZHi)
        use probdata_module

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

    
        !$OMP PARALLEL DO PRIVATE(i,j,k)
        do k = p_l3,p_h3

           do j = p_l2,p_h2

              do i = p_l1,p_h1

                 if     ( k .lt. domlo(3) ) then

                   phi(i,j,k) = bcXYLo(i,j)

                 elseif ( k .gt. domhi(3) ) then
 
                   phi(i,j,k) = bcXYHi(i,j)

                 elseif ( j .lt. domlo(2) ) then

                   phi(i,j,k) = bcXZLo(i,k)

                 elseif ( j .gt. domhi(2) ) then

                   phi(i,j,k) = bcXZHi(i,k)

                 elseif ( i .lt. domlo(1) ) then

                   phi(i,j,k) = bcYZLo(j,k)

                 elseif ( i .gt. domhi(1) ) then

                   phi(i,j,k) = bcYZHi(j,k)

                endif

              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO

      end subroutine ca_put_direct_sum_bc

! ::
! :: ----------------------------------------------------------
! ::

      double precision function direct_sum_symmetric_add (loc,locb,problo,probhi, &
                       rho,dV,doSymmetricAddLo,doSymmetricAddHi) result(bcTerm)

        use probdata_module
        use fundamental_constants_module, only: Gconst
        use bl_constants_module

        implicit none

        double precision :: loc(3), locb(3)
        double precision :: problo(3), probhi(3)
        logical          :: doSymmetricAddLo(3), doSymmetricAddHi(3)

        double precision :: x, y, z, r
        double precision :: rho, dV

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

      
