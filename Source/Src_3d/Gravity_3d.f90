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

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_put_multipole_bc (lo,hi,domlo,domhi,dx,&
                                      phi,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3, &
                                      problo,lnum,lmax,q0,qC,qS)
        use probdata_module
        use fundamental_constants_module, only: Gconst

        implicit none

        integer          :: lo(3),hi(3)
        integer          :: domlo(3),domhi(3)
        double precision :: dx(3)
        double precision :: problo(3)

        integer          :: lnum, lmax
        double precision :: q0(0:lmax), qC(0:lmax,0:lmax), qS(0:lmax,0:lmax)

        integer          :: p_l1,p_l2,p_l3,p_h1,p_h2,p_h3
        double precision :: phi(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3)

        integer          :: i,j,k
        integer          :: l,m
        double precision :: x,y,z,r,cosTheta,phiAngle
        double precision :: legPolyArr(0:lmax), assocLegPolyArr(0:lmax,0:lmax)
        
        !$OMP PARALLEL DO PRIVATE(i,j,k,l,m,x,y,z,r,cosTheta,phiAngle,legPolyArr,assocLegPolyArr)
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
        
                 ! Only adjust ghost zones here

                 if ( i .lt. domlo(1) .or. i .gt. domhi(1) .or. &
                      j .lt. domlo(2) .or. j .gt. domhi(2) .or. &
                      k .lt. domlo(3) .or. k .gt. domhi(3) ) then

                   r = sqrt( x**2 + y**2 + z**2 )
                   cosTheta = z / r
                   phiAngle = atan2(y,x)

                   phi(i,j,k) = 0.0d0

                   ! First, calculate the Legendre polynomials.

                   legPolyArr(:) = 0.0d0
                   assocLegPolyArr(:,:) = 0.0d0

                   call fill_legendre_arrays(legPolyArr, assocLegPolyArr, cosTheta, lnum, lmax)

                   ! Now compute the potentials on the ghost cells.

                   do l = 0, lnum
 
                     phi(i,j,k) = phi(i,j,k) + q0(l) * legPolyArr(l) / r**(l+1)

                     do m = 1, l
       
                       phi(i,j,k) = phi(i,j,k) + (qC(l,m) * cos(m * phiAngle) + qS(l,m) * sin(m * phiAngle)) * &
                                                 assocLegPolyArr(l,m) / r**(l+1)

                     enddo

                   enddo

                   phi(i,j,k) = Gconst * phi(i,j,k)

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
                                               problo,probhi,lnum,lmax,q0,qC,qS)
        use probdata_module

        implicit none

        integer          :: lo(3),hi(3)
        integer          :: lo_bc(3),hi_bc(3)
        integer          :: symmetry_type
        integer          :: domlo(3),domhi(3)
        double precision :: dx(3)
        double precision :: problo(3), probhi(3)

        ! lmax only controls the size of the arrays;
        ! lnum is the actual maximum value of l we calculate.

        integer          :: lnum, lmax
        double precision :: q0(0:lmax), qC(0:lmax,0:lmax), qS(0:lmax,0:lmax)
        double precision :: factArray(0:lmax,0:lmax)

        integer          :: p_l1,p_l2,p_l3,p_h1,p_h2,p_h3
        double precision :: rho(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3)

        integer          :: i,j,k
        integer          :: ilo, ihi, jlo, jhi, klo, khi
        integer          :: l,m,b

        double precision :: factorial

        double precision :: x,y,z,r,cosTheta,phiAngle,dV

        double precision :: volumeFactor, parityFactor
        double precision :: edgeTolerance = 1.0d-2
        double precision :: legPolyArr(0:lmax), assocLegPolyArr(0:lmax,0:lmax)

        ! Variables for symmetric BCs

        double precision :: xLo, yLo, zLo, xHi, yHi, zHi
        double precision :: cosThetaLoX(4), phiAngleLoX(4)
        double precision :: cosThetaLoY(2), phiAngleLoY(2)
        double precision :: cosThetaLoZ(1), phiAngleLoZ(1)
        double precision :: cosThetaHiX(4), phiAngleHiX(4)
        double precision :: cosThetaHiY(2), phiAngleHiY(2)
        double precision :: cosThetaHiZ(1), phiAngleHiZ(1)
        double precision :: legPolyArrLoX(0:lmax,4), assocLegPolyArrLoX(0:lmax,0:lmax,4)
        double precision :: legPolyArrLoY(0:lmax,2), assocLegPolyArrLoY(0:lmax,0:lmax,2)
        double precision :: legPolyArrLoZ(0:lmax,1), assocLegPolyArrLoZ(0:lmax,0:lmax,1)
        double precision :: legPolyArrHiX(0:lmax,4), assocLegPolyArrHiX(0:lmax,0:lmax,4)
        double precision :: legPolyArrHiY(0:lmax,2), assocLegPolyArrHiY(0:lmax,0:lmax,2)
        double precision :: legPolyArrHiZ(0:lmax,1), assocLegPolyArrHiZ(0:lmax,0:lmax,1)

        logical          :: doSymmetricAddLo(3), doSymmetricAddHi(3)
        logical          :: doReflectionLo(3), doReflectionHi(3)

        dV = dx(1) * dx(2) * dx(3)

        ! Safety check: if lnum > lmax, set it equal to lmax to avoid array overflows.
        ! lmax is set so that you really could never need that many terms.

        if ( lnum > lmax ) then
          lnum = lmax
        endif

        ! If any of the boundaries are symmetric, we need to account for the mass that is assumed
        ! to lie on the opposite side of the symmetric axis. If the center in any direction 
        ! coincides with the boundary, then we can simply double the mass as a result of that reflection.
        ! Otherwise, we need to do a more general solve.

        volumeFactor = 1.0d0
        parityFactor = 1.0d0

        doSymmetricAddLo(:) = .false.
        doSymmetricAddHi(:) = .false.

        doReflectionLo(:)   = .false.
        doReflectionHi(:)   = .false.

        do b = 1, 3

          if ( (lo_bc(b) .eq. symmetry_type) ) then
            if ( abs(center(b) - problo(b)) < edgeTolerance ) then
              volumeFactor = volumeFactor * 2.0d0
              doReflectionLo(b) = .true.
            else
              doSymmetricAddLo(b) = .true.
            endif
          endif

          if ( (hi_bc(b) .eq. symmetry_type) ) then
            if ( abs(center(b) - probhi(b)) < edgeTolerance ) then
              volumeFactor = volumeFactor * 2.0d0
              doReflectionHi(b) = .true.
            else
              doSymmetricAddHi(b) = .true.
            endif
          endif

        enddo 


        ! Compute pre-factors now to save computation time, for qC and qS

        do l = 0, lnum
          do m = 1, l
            factArray(l,m) = 2.0 * factorial(l-m) / factorial(l+m) * dV * volumeFactor
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

        !$OMP PARALLEL DO PRIVATE(i,j,k,legPolyArr,assocLegPolyArr,x,y,z,r,cosTheta,phiAngle ) &
        !$OMP PRIVATE(l,m,parityFactor,xLo,yLo,zLo,xHi,yHi,zHi                               ) &
        !$OMP PRIVATE(cosThetaLoX,phiAngleLoX,cosThetaLoY,phiAngleLoY,cosThetaLoZ,phiAngleLoZ) &
        !$OMP PRIVATE(legPolyArrLoX,legPolyArrLoY,legPolyArrLoZ                              ) &
        !$OMP PRIVATE(assocLegPolyArrLoX,assocLegPolyArrLoY,assocLegPolyArrLoZ               ) &
        !$OMP PRIVATE(cosThetaHiX,phiAngleHiX,cosThetaHiY,phiAngleHiY,cosThetaHiZ,phiAngleHiZ) &
        !$OMP PRIVATE(legPolyArrHiX,legPolyArrHiY,legPolyArrHiZ                              ) &
        !$OMP PRIVATE(assocLegPolyArrHiX,assocLegPolyArrHiY,assocLegPolyArrHiZ               ) &
        !$OMP REDUCTION(+:q0,qC,qS)
        do k = klo, khi
           z = problo(3) + (dble(k)+0.50d0) * dx(3) - center(3)

           ! Compute reflections of x, y and z about their respective axes,
           ! modulo any difference of center from zero. These are used
           ! if we have symmetric boundary conditions, though in those cases
           ! the center will most likely coincide with the symmetric axis.
 

           do j = jlo, jhi
              y = problo(2) + (dble(j)+0.50d0) * dx(2) - center(2)
 
              do i = ilo, ihi
                 x = problo(1) + (dble(i)+0.50d0) * dx(1) - center(1)

                 r = sqrt( x**2 + y**2 + z**2 )
                 cosTheta = z / r
                 phiAngle = atan2(y,x)

                 ! Compute contribution from this cell to all multipole moments.
                 ! First, calculate the Legendre polynomials.

                 legPolyArr(:) = 0.0d0
                 assocLegPolyArr(:,:) = 0.0d0

                 call fill_legendre_arrays(legPolyArr, assocLegPolyArr, cosTheta, lnum, lmax)

                 ! Now, determine the contributions from any symmetric axes. First we do the 
                 ! lower boundaries. We pre-populate the arrays now instead of doing it
                 ! in the loop over multipole orders. This minimizes function calls at
                 ! the expense of longer code.

                 if ( doSymmetricAddLo(1) ) then

                   legPolyArrLoX(:,:) = 0.0d0
                   assocLegPolyArrLoX(:,:,:) = 0.0d0

                   cosThetaLoX(1) = z / r
                   phiAngleLoX(1) = atan2(y,xLo)
                   call fill_legendre_arrays(legPolyArrLoX(:,1), assocLegPolyArrLoX(:,:,1), cosThetaLoX(1), lnum, lmax) 

                   if ( doSymmetricAddLo(2) ) then

                     cosThetaLoX(2) = z / r
                     phiAngleLoX(2) = atan2(yLo,xLo)
                     call fill_legendre_arrays(legPolyArrLoX(:,2), assocLegPolyArrLoX(:,:,2), cosThetaLoX(2), lnum, lmax)

                   endif

                   if ( doSymmetricAddLo(3) ) then

                     cosThetaLoX(3) = zLo / r
                     phiAngleLoX(3) = atan2(y,xLo)
                     call fill_legendre_arrays(legPolyArrLoX(:,3), assocLegPolyArrLoX(:,:,3), cosThetaLoX(3), lnum, lmax)

                   endif

                   if ( doSymmetricAddLo(2) .and. doSymmetricAddLo(3) ) then
 
                     cosThetaLoX(4) = zLo / r
                     phiAngleLoX(4) = atan2(yLo,xLo)
                     call fill_legendre_arrays(legPolyArrLoX(:,4), assocLegPolyArrLoX(:,:,4), cosThetaLoX(4), lnum, lmax)
                
                   endif

                 endif

                 if ( doSymmetricAddLo(2) ) then

                   legPolyArrLoY(:,:) = 0.0d0
                   assocLegPolyArrLoY(:,:,:) = 0.0d0

                   cosThetaLoY(1) = z / r
                   phiAngleLoY(1) = atan2(yLo,x)
                   call fill_legendre_arrays(legPolyArrLoY(:,1), assocLegPolyArrLoY(:,:,1), cosThetaLoY(1), lnum, lmax)

                   if ( doSymmetricAddLo(3) ) then

                     cosThetaLoY(2) = zLo / r
                     phiAngleLoY(2) = atan2(yLo,x)
                     call fill_legendre_arrays(legPolyArrLoY(:,2), assocLegPolyArrLoY(:,:,2), cosThetaLoY(2), lnum, lmax)

                   endif

                 endif

                 if ( doSymmetricAddLo(3) ) then

                   legPolyArrLoZ(:,:) = 0.0d0
                   assocLegPolyArrLoZ(:,:,:) = 0.0d0

                   cosThetaLoZ(1) = zLo / r
                   phiAngleLoZ(1) = atan2(y,x)
                   call fill_legendre_arrays(legPolyArrLoZ(:,1), assocLegPolyArrLoZ(:,:,1), cosThetaLoZ(1), lnum, lmax)

                 endif


                 ! Do the same for the upper boundaries.

                 if ( doSymmetricAddHi(1) ) then

                   legPolyArrHiX(:,:) = 0.0d0
                   assocLegPolyArrHiX(:,:,:) = 0.0d0

                   cosThetaHiX(1) = z / r
                   phiAngleHiX(1) = atan2(y,xHi)
                   call fill_legendre_arrays(legPolyArrHiX(:,1), assocLegPolyArrHiX(:,:,1), cosThetaHiX(1), lnum, lmax) 

                   if ( doSymmetricAddHi(2) ) then

                     cosThetaHiX(2) = z / r
                     phiAngleHiX(2) = atan2(yHi,xHi)
                     call fill_legendre_arrays(legPolyArrHiX(:,2), assocLegPolyArrHiX(:,:,2), cosThetaHiX(2), lnum, lmax)

                   endif

                   if ( doSymmetricAddHi(3) ) then

                     cosThetaHiX(3) = zHi / r
                     phiAngleHiX(3) = atan2(y,xHi)
                     call fill_legendre_arrays(legPolyArrHiX(:,3), assocLegPolyArrHiX(:,:,3), cosThetaHiX(3), lnum, lmax)

                   endif

                   if ( doSymmetricAddHi(2) .and. doSymmetricAddHi(3) ) then
 
                     cosThetaHiX(4) = zHi / r
                     phiAngleHiX(4) = atan2(yHi,xHi)
                     call fill_legendre_arrays(legPolyArrHiX(:,4), assocLegPolyArrHiX(:,:,4), cosThetaHiX(4), lnum, lmax)
                
                   endif

                 endif

                 if ( doSymmetricAddHi(2) ) then

                   legPolyArrHiY(:,:) = 0.0d0
                   assocLegPolyArrHiY(:,:,:) = 0.0d0

                   cosThetaHiY(1) = z / r
                   phiAngleHiY(1) = atan2(yHi,x)
                   call fill_legendre_arrays(legPolyArrHiY(:,1), assocLegPolyArrHiY(:,:,1), cosThetaHiY(1), lnum, lmax)

                   if ( doSymmetricAddHi(3) ) then

                     cosThetaHiY(2) = zHi / r
                     phiAngleHiY(2) = atan2(yHi,x)
                     call fill_legendre_arrays(legPolyArrHiY(:,2), assocLegPolyArrHiY(:,:,2), cosThetaHiY(2), lnum, lmax)

                   endif

                 endif

                 if ( doSymmetricAddHi(3) ) then

                   legPolyArrHiZ(:,:) = 0.0d0
                   assocLegPolyArrHiZ(:,:,:) = 0.0d0

                   cosThetaHiZ(1) = zHi / r
                   phiAngleHiZ(1) = atan2(y,x)
                   call fill_legendre_arrays(legPolyArrHiZ(:,1), assocLegPolyArrHiZ(:,:,1), cosThetaHiZ(1), lnum, lmax)

                 endif



                 ! Now, compute the multipole moments using the tabulated polynomials.

                 do l = 0, lnum 

                   ! The odd l Legendre polynomials are odd in their argument, so
                   ! a symmetric reflection about the z axis leads to a total cancellation.

                   parityFactor = 1.0d0

                   if ( MODULO(l,2) /= 0 .and. ( doReflectionLo(3) .or. doReflectionHi(3) ) ) then
                     parityFactor = 0.0d0
                   endif

                   q0(l) = q0(l) + legPolyArr(l) * (r ** l) * rho(i,j,k) * dV * volumeFactor * parityFactor



                   ! Add contributions from any symmetric boundaries that were not handled
                   ! by the reflections.       

                   if ( doSymmetricAddLo(1) ) then

                     q0(l) = q0(l) + legPolyArrLoX(l,1) * (r ** l) * rho(i,j,k) * dV

                     if ( doSymmetricAddLo(2) ) then
                       q0(l) = q0(l) + legPolyArrLoX(l,2) * (r ** l) * rho(i,j,k) * dV
                     endif

                     if ( doSymmetricAddLo(3) ) then
                       q0(l) = q0(l) + legPolyArrLoX(l,3) * (r ** l) * rho(i,j,k) * dV
                     endif

                     if ( doSymmetricAddLo(2) .and. doSymmetricAddLo(3) ) then
                       q0(l) = q0(l) + legPolyArrLoX(l,4) * (r ** l) * rho(i,j,k) * dV
                     endif

                   endif

                   if ( doSymmetricAddLo(2) ) then
 
                     q0(l) = q0(l) + legPolyArrLoY(l,1) * (r ** l) * rho(i,j,k) * dV

                     if ( doSymmetricAddLo(3) ) then
                       q0(l) = q0(l) + legPolyArrLoY(l,2) * (r ** l) * rho(i,j,k) * dV
                     endif

                   endif

                   if ( doSymmetricAddLo(3) ) then

                     q0(l) = q0(l) + legPolyArrLoZ(l,1) * (r ** l) * rho(i,j,k) * dV

                   endif





                   if ( doSymmetricAddHi(1) ) then

                     q0(l) = q0(l) + legPolyArrHiX(l,1) * (r ** l) * rho(i,j,k) * dV

                     if ( doSymmetricAddHi(2) ) then
                       q0(l) = q0(l) + legPolyArrHiX(l,2) * (r ** l) * rho(i,j,k) * dV
                     endif

                     if ( doSymmetricAddHi(3) ) then
                       q0(l) = q0(l) + legPolyArrHiX(l,3) * (r ** l) * rho(i,j,k) * dV
                     endif

                     if ( doSymmetricAddHi(2) .and. doSymmetricAddHi(3) ) then
                       q0(l) = q0(l) + legPolyArrHiX(l,4) * (r ** l) * rho(i,j,k) * dV
                     endif

                   endif

                   if ( doSymmetricAddHi(2) ) then
 
                     q0(l) = q0(l) + legPolyArrHiY(l,1) * (r ** l) * rho(i,j,k) * dV

                     if ( doSymmetricAddHi(3) ) then
                       q0(l) = q0(l) + legPolyArrHiY(l,2) * (r ** l) * rho(i,j,k) * dV
                     endif

                   endif

                   if ( doSymmetricAddHi(3) ) then

                     q0(l) = q0(l) + legPolyArrHiZ(l,1) * (r ** l) * rho(i,j,k) * dV

                   endif


                   do m = 1, l

                     ! The parity properties of the associated Legendre polynomials are:
                     ! P_l^m (-x) = (-1)^(l+m) P_l^m (x)
                     ! Therefore, a complete cancellation occurs if l+m is odd and
                     ! we are reflecting about the z axis.

                     ! Additionally, the cosine and sine terms flip sign when reflected
                     ! about the x or y axis, so if we have a reflection about x or y
                     ! then the terms have a complete cancellation.

                     parityFactor = 1.0d0

                     if ( MODULO(l+m,2) /= 0 .and. ( doReflectionLo(3) .or. doReflectionHi(3) ) ) then
                       parityFactor = 0.0d0
                     endif

                     if ( doReflectionLo(1) .or. doReflectionLo(2) .or. doReflectionHi(1) .or. doReflectionHi(2) ) then
                       parityFactor = 0.0d0
                     endif

                     qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArr(l,m) * &
                                         cos(m * phiAngle) * (r ** l) * rho(i,j,k) * parityFactor

                     qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArr(l,m) * &
                                         sin(m * phiAngle) * (r ** l) * rho(i,j,k) * parityFactor


                     ! Now add in contributions if we have any symmetric boundaries

                     if ( doSymmetricAddLo(1) ) then

                       qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrLoX(l,m,1) * &
                                           cos(m * phiAngleLoX(1)) * (r ** l) * rho(i,j,k)
                       qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrLoX(l,m,1) * &
                                           sin(m * phiAngleLoX(1)) * (r ** l) * rho(i,j,k)

                       if ( doSymmetricAddLo(2) ) then

                         qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrLoX(l,m,2) * &
                                             cos(m * phiAngleLoX(2)) * (r ** l) * rho(i,j,k)
                         qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrLoX(l,m,2) * &
                                             sin(m * phiAngleLoX(2)) * (r ** l) * rho(i,j,k)

                       endif

                       if ( doSymmetricAddLo(3) ) then

                         qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrLoX(l,m,3) * &
                                             cos(m * phiAngleLoX(3)) * (r ** l) * rho(i,j,k)
                         qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrLoX(l,m,3) * &
                                             sin(m * phiAngleLoX(3)) * (r ** l) * rho(i,j,k)

                       endif

                       if ( doSymmetricAddLo(2) .and. doSymmetricAddLo(3) ) then

                         qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrLoX(l,m,4) * &
                                             cos(m * phiAngleLoX(4)) * (r ** l) * rho(i,j,k)
                         qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrLoX(l,m,4) * &
                                             sin(m * phiAngleLoX(4)) * (r ** l) * rho(i,j,k)

                       endif

                     endif

                     if ( doSymmetricAddLo(2) ) then

                       qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrLoY(l,m,1) * &
                                           cos(m * phiAngleLoY(1)) * (r ** l) * rho(i,j,k)
                       qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrLoY(l,m,1) * &
                                           sin(m * phiAngleLoY(1)) * (r ** l) * rho(i,j,k)

                       if ( doSymmetricAddLo(3) ) then

                         qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrLoY(l,m,2) * &
                                             cos(m * phiAngleLoY(2)) * (r ** l) * rho(i,j,k)
                         qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrLoY(l,m,2) * &
                                             sin(m * phiAngleLoY(2)) * (r ** l) * rho(i,j,k)

                       endif

                     endif

                     if ( doSymmetricAddLo(3) ) then

                       qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrLoZ(l,m,1) * &
                                           cos(m * phiAngleLoZ(1)) * (r ** l) * rho(i,j,k)
                       qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrLoX(l,m,1) * &
                                           sin(m * phiAngleLoZ(1)) * (r ** l) * rho(i,j,k)

                     endif


                     ! Repeat for the upper boundaries.

                     if ( doSymmetricAddHi(1) ) then

                       qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrHiX(l,m,1) * &
                                           cos(m * phiAngleHiX(1)) * (r ** l) * rho(i,j,k)
                       qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrHiX(l,m,1) * &
                                           sin(m * phiAngleHiX(1)) * (r ** l) * rho(i,j,k)

                       if ( doSymmetricAddHi(2) ) then

                         qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrHiX(l,m,2) * &
                                             cos(m * phiAngleHiX(2)) * (r ** l) * rho(i,j,k)
                         qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrHiX(l,m,2) * &
                                             sin(m * phiAngleHiX(2)) * (r ** l) * rho(i,j,k)

                       endif

                       if ( doSymmetricAddHi(3) ) then

                         qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrHiX(l,m,3) * &
                                             cos(m * phiAngleHiX(3)) * (r ** l) * rho(i,j,k)
                         qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrHiX(l,m,3) * &
                                             sin(m * phiAngleHiX(3)) * (r ** l) * rho(i,j,k)

                       endif

                       if ( doSymmetricAddHi(2) .and. doSymmetricAddHi(3) ) then

                         qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrHiX(l,m,4) * &
                                             cos(m * phiAngleHiX(4)) * (r ** l) * rho(i,j,k)
                         qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrHiX(l,m,4) * &
                                             sin(m * phiAngleHiX(4)) * (r ** l) * rho(i,j,k)

                       endif

                     endif

                     if ( doSymmetricAddHi(2) ) then

                       qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrHiY(l,m,1) * &
                                           cos(m * phiAngleHiY(1)) * (r ** l) * rho(i,j,k)
                       qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrHiY(l,m,1) * &
                                           sin(m * phiAngleHiY(1)) * (r ** l) * rho(i,j,k)

                       if ( doSymmetricAddHi(3) ) then

                         qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrHiY(l,m,2) * &
                                             cos(m * phiAngleHiY(2)) * (r ** l) * rho(i,j,k)
                         qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrHiY(l,m,2) * &
                                             sin(m * phiAngleHiY(2)) * (r ** l) * rho(i,j,k)

                       endif

                     endif

                     if ( doSymmetricAddHi(3) ) then

                       qC(l,m) = qC(l,m) + factArray(l,m) * assocLegPolyArrHiZ(l,m,1) * &
                                           cos(m * phiAngleHiZ(1)) * (r ** l) * rho(i,j,k)
                       qS(l,m) = qS(l,m) + factArray(l,m) * assocLegPolyArrHiX(l,m,1) * &
                                           sin(m * phiAngleHiZ(1)) * (r ** l) * rho(i,j,k)

                     endif

                   enddo

                 enddo

              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO

      end subroutine ca_compute_multipole_moments


! ::
! :: ----------------------------------------------------------
! ::

      double precision function factorial(n)
      
        implicit none

        integer :: n, i

        factorial = 1
 
        do i = 1, n
          factorial = factorial * i
        enddo

      end function factorial

! ::
! :: ----------------------------------------------------------
! ::

      subroutine fill_legendre_arrays(legPolyArr, assocLegPolyArr, x, lnum, lmax)
      
        implicit none

        integer :: lnum, lmax
        integer :: l, m, ll
        double precision :: x
        double precision :: legPolyArr(0:lmax), assocLegPolyArr(0:lmax,0:lmax)

        do l = 0, lnum

          ! If l = 0 or 1, use explicit expressions for the polynomials.
          ! Otherwise, use the recurrence relations.

          if ( l == 0 ) then

            legPolyArr(0) = 1.0d0
            assocLegPolyArr(0,0) = 1.0d0

          elseif ( l == 1 ) then

            legPolyArr(1) = x

            assocLegPolyArr(1,0) = x
            assocLegPolyArr(1,1) = -(1 - x**2)**(0.5d0)

          else

            legPolyArr(l) = ( (2*l - 1) * x * legPolyArr(l-1) - (l-1) * legPolyArr(l-2) ) / l

            ! The m = l term has an explicit expression using a double factorial:
            ! P_l^l(x) = (-1)^l * (2l-1)!! * (1-x^2)^(l/2)

            assocLegPolyArr(l,l) = (-1)**l * (1 - x**2)**(l / 2.0d0)

            do ll = (2*l-1), 3, -2

              assocLegPolyArr(l,l) = assocLegPolyArr(l,l) * (2*l - 1)

            enddo

            ! For 1 < m < l, use the following recurrence relation (Wikipedia):
            ! P_l^m(x) = -(1 - x^2)^(1/2) / (2*m) * ( P_{l-1}^{m+1}(x) + (l+m-1)(l+m)*P_{l-1}^{m-1}(x) )

            do m = 1, l-1

              assocLegPolyArr(l,m) = -(1.0d0 - x**2)**(0.5d0) / (2*m) * &
                                      ( assocLegPolyArr(l-1,m+1) + &
                                      (l+m-1)*(l+m)*assocLegPolyArr(l-1,m-1) )

            enddo

            ! That one fails for m = 0, which we need to populate for the next l
            ! even though it's not used in the expansion directly. Wikipedia:
            ! l * P_l^0 (x) = x * (2l - 1)*P_{l-1}^0 (x) - (l - 1) P_{l-2}^0(x)

            assocLegPolyArr(l,0) = 1.0d0 / l * ( x * (2*l - 1) * assocLegPolyArr(l-1,0) &
                                   - (l-1) * assocLegPolyArr(l-2,0) )

          endif

        enddo

      end subroutine fill_legendre_arrays
