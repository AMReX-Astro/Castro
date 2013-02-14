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
      volfrac = 1.d0/float(lrat(1)*lrat(2)*lrat(3))
      !
      ! ::::: set coarse grid to zero on overlap
      !
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,kc) = 0.d0
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

      !$OMP PARALLEL PRIVATE(i,j,k)
      !$OMP DO
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,k,1) = 0.5d0 * ( ecx(i+1,j,k) + ecx(i,j,k) )
            enddo
         enddo
      enddo
      !$OMP END DO NOWAIT

      !$OMP DO
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,k,2) = 0.5d0 * ( ecy(i,j+1,k) + ecy(i,j,k) )
            enddo
         enddo
      enddo
      !$OMP END DO NOWAIT

      !$OMP DO
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,k,3) = 0.5d0 * ( ecz(i,j,k+1) + ecz(i,j,k) )
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      ! Note this assumes the lo end of the domain is 0.
      !
      if (ccl1 .lt. 0 .and. bc_lo(1) .eq. symmetry_type) then
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),-1
                  cc(i,j,k,1) = -cc(-i-1,j,k,1)
                  cc(i,j,k,2) =  cc(-i-1,j,k,2)
                  cc(i,j,k,3) =  cc(-i-1,j,k,3)
               enddo
            enddo
         enddo
      endif
      !
      ! Note this assumes the lo end of the domain is 0.
      !
      if (ccl2 .lt. 0 .and. bc_lo(2) .eq. symmetry_type) then
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               do j=lo(2),-1
                  cc(i,j,k,1) =  cc(i,-j-1,k,1)
                  cc(i,j,k,2) = -cc(i,-j-1,k,2)
                  cc(i,j,k,3) =  cc(i,-j-1,k,3)
               enddo
            enddo
         enddo
      endif
      !
      ! Note this assumes the lo end of the domain is 0.
      !
      if (ccl3 .lt. 0 .and. bc_lo(3) .eq. symmetry_type) then
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do k=lo(3),-1
                  cc(i,j,k,1) =  cc(i,j,-k-1,1)
                  cc(i,j,k,2) =  cc(i,j,-k-1,2)
                  cc(i,j,k,3) = -cc(i,j,-k-1,3)
               enddo
            enddo
         enddo
      endif
         
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

      !$OMP PARALLEL DO PRIVATE(i,j,k,lapphi) 
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
      !$OMP END PARALLEL DO
 
      end subroutine ca_test_residual

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_average_ec ( &
           fx, fxl1, fxl2, fxl3, fxh1, fxh2, fxh3, &
           fy, fyl1, fyl2, fyl3, fyh1, fyh2, fyh3, &
           fz, fzl1, fzl2, fzl3, fzh1, fzh2, fzh3, &
           cx, cxl1, cxl2, cxl3, cxh1, cxh2, cxh3, &
           cy, cyl1, cyl2, cyl3, cyh1, cyh2, cyh3, &
           cz, czl1, czl2, czl3, czh1, czh2, czh3, &
           lo, hi, rr)

      integer lo(3),hi(3)
      integer fxl1, fxl2, fxl3, fxh1, fxh2, fxh3
      integer fyl1, fyl2, fyl3, fyh1, fyh2, fyh3
      integer fzl1, fzl2, fzl3, fzh1, fzh2, fzh3
      integer cxl1, cxl2, cxl3, cxh1, cxh2, cxh3
      integer cyl1, cyl2, cyl3, cyh1, cyh2, cyh3
      integer czl1, czl2, czl3, czh1, czh2, czh3
      double precision fx(fxl1:fxh1,fxl2:fxh2,fxl3:fxh3)
      double precision fy(fyl1:fyh1,fyl2:fyh2,fyl3:fyh3)
      double precision fz(fzl1:fzh1,fzl2:fzh2,fzl3:fzh3)
      double precision cx(cxl1:cxh1,cxl2:cxh2,cxl3:cxh3)
      double precision cy(cyl1:cyh1,cyl2:cyh2,cyl3:cyh3)
      double precision cz(czl1:czh1,czl2:czh2,czl3:czh3)
      integer rr(3)

      integer i,j,k,n,m,facx,facy,facz

      facx = rr(1)
      facy = rr(2)
      facz = rr(3)

      !$OMP PARALLEL PRIVATE(i,j,k,m,n)
      !$OMP DO
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               cx(i,j,k) = 0.d0
               do n = 0,facy-1
                  do m = 0,facz-1
                     cx(i,j,k) = cx(i,j,k) + fx(facx*i,facy*j+n,facz*k+m)
                  end do
               end do
               cx(i,j,k) = cx(i,j,k) / (facy*facz)
            end do
         end do
      end do
      !$OMP END DO NOWAIT

      !$OMP DO
      do k = lo(3), hi(3)
         do i = lo(1), hi(1)
            do j = lo(2), hi(2)+1
               cy(i,j,k) = 0.d0
               do n = 0,facx-1
                  do m = 0,facz-1
                     cy(i,j,k) = cy(i,j,k) + fy(facx*i+n,facy*j,facz*k+m)
                  end do
               end do
               cy(i,j,k) = cy(i,j,k) / (facx*facz)
            end do
         end do
      end do
      !$OMP END DO NOWAIT

      !$OMP DO
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do k = lo(3), hi(3)+1
               cz(i,j,k) = 0.d0
               do n = 0,facx-1
                  do m = 0,facy-1
                     cz(i,j,k) = cz(i,j,k) + fz(facx*i+n,facy*j+m,facz*k)
                  end do
               end do
               cz(i,j,k) = cz(i,j,k) / (facx*facy)
            end do
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      end subroutine ca_average_ec

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_compute_radial_mass (lo,hi,dx,dr,&
                                         rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,&
                                         radial_mass,radial_vol,problo,&
                                         n1d,drdxfac,level)
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

      integer          :: i,j,k,index
      integer          :: ii,jj,kk
      double precision :: xc,yc,zc,r,xxsq,yysq,zzsq,octant_factor
      double precision :: fac,xx,yy,zz,dx_frac,dy_frac,dz_frac
      double precision :: vol_frac
      double precision :: lo_i,lo_j,lo_k

      if (( abs(center(1) - problo(1)) .lt. 1.e-2 * dx(1) ) .and. &
          ( abs(center(2) - problo(2)) .lt. 1.e-2 * dx(2) ) .and. &
          ( abs(center(3) - problo(3)) .lt. 1.e-2 * dx(3) ) ) then
         octant_factor = 8.d0
      else
         octant_factor = 1.d0
      end if


      fac     = dble(drdxfac)
      dx_frac = dx(1) / fac
      dy_frac = dx(2) / fac
      dz_frac = dx(3) / fac

      vol_frac = octant_factor * dx_frac * dy_frac * dz_frac
      !
      ! Don't OMP this.
      !
      do k = lo(3), hi(3)
         zc = problo(3) + (dble(k)+0.50d0) * dx(3) - center(3)

         do j = lo(2), hi(2)
            yc = problo(2) + (dble(j)+0.50d0) * dx(2) - center(2)

            do i = lo(1), hi(1)
               xc  = problo(1) + (dble(i)+0.50d0) * dx(1) - center(1)

               r = sqrt(xc**2 + yc**2 + zc**2)
               index = int(r/dr)

               if (index .gt. n1d-1) then

                  if (level .eq. 0) then
                     print *,'   '  
                     print *,'>>> Error: Gravity_3d::ca_compute_radial_mass ',i,j,k
                     print *,'>>> ... index too big: ', index,' > ',n1d-1
                     print *,'>>> ... at (i,j,k)   : ',i,j,k
                     call bl_error("Error:: Gravity_3d.f90 :: ca_compute_radial_mass")
                  end if

               else

                  lo_i =  problo(1) + dble(i)*dx(1) - center(1)
                  lo_j =  problo(2) + dble(j)*dx(2) - center(2)
                  lo_k =  problo(3) + dble(k)*dx(3) - center(3)

                  do kk = 0,drdxfac-1
                     zz   = lo_k + (dble(kk)+0.5d0)*dz_frac
                     zzsq = zz*zz
                     do jj = 0,drdxfac-1
                        yy   = lo_j + (dble(jj)+0.5d0)*dy_frac
                        yysq = yy*yy
                        do ii = 0,drdxfac-1

                           xx    = lo_i + (dble(ii)+0.5d0)*dx_frac
                           xxsq  = xx*xx
                           r     = sqrt(xxsq  + yysq + zzsq)
                           index = int(r/dr)

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
      !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,r,index,cen,xi,slope,mag_grav,ghi,gmd,glo,minvar,maxvar)
      do k = g_l3,g_h3
         z = problo(3) + (dble(k)+0.50d0) * dx(3) - center(3)

         do j = g_l2,g_h2
            y = problo(2) + (dble(j)+0.50d0) * dx(2) - center(2)

            do i = g_l1,g_h1
               x     = problo(1) + (dble(i)+0.50d0) * dx(1) - center(1)
               r     = sqrt(x**2 + y**2 + z**2)
               index = int(r/dr)
               cen   = (dble(index)+0.5d0)*dr
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
                   ( ghi -  2.d0*gmd + glo)*xi**2/(2.d0*dr**2) + &
                   ( ghi       - glo      )*xi   /(2.d0*dr   ) + &
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
      !$OMP END PARALLEL DO

      end subroutine ca_put_radial_grav

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_put_radial_phi (lo,hi,domlo,domhi,dx,dr,&
                                    phi,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
                                    radial_phi,problo,&
                                    numpts_1d,fill_interior)
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
        do k = p_l3,p_h3
           if (k .gt. domhi(3)) then
              z = problo(3) + (dble(k  )       ) * dx(3) - center(3)
           else if (k .lt. domlo(3)) then
              z = problo(3) + (dble(k+1)       ) * dx(3) - center(3)
           else 
              z = problo(3) + (dble(k  )+0.50d0) * dx(3) - center(3)
           end if

           do j = p_l2,p_h2
              if (j .gt. domhi(2)) then
                 y = problo(2) + (dble(j  )       ) * dx(2) - center(2)
              else if (j .lt. domlo(2)) then
                 y = problo(2) + (dble(j+1)       ) * dx(2) - center(2)
              else 
                 y = problo(2) + (dble(j  )+0.50d0) * dx(2) - center(2)
              end if

              do i = p_l1,p_h1
                 if (i .gt. domhi(1)) then
                    x = problo(1) + (dble(i  )       ) * dx(1) - center(1)
                 else if (i .lt. domlo(1)) then
                    x = problo(1) + (dble(i+1)       ) * dx(1) - center(1)
                 else 
                    x = problo(1) + (dble(i  )+0.50d0) * dx(1) - center(1)
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
                    cen = (dble(index)+0.5d0)*dr
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
                            ( p_hi -  2.d0*p_md + p_lo)*xi**2/(2.d0*dr**2) + &
                            ( p_hi       - p_lo      )*xi   /(2.d0*dr   ) + &
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
