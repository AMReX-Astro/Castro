! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine ca_edge_interp(flo, fhi, nc, ratio, dir, &
           fine, fine_l0, fine_h0)
      implicit none
      integer flo(0:2-1), fhi(0:2-1), nc, ratio(0:2-1), dir
      integer fine_l0, fine_h0
      double precision fine(fine_l0:fine_h0,nc)
      integer i,n,M
      double precision val, df

!     Do linear in dir, pc transverse to dir, leave alone the fine values
!     lining up with coarse edges--assume these have been set to hold the 
!     values you want to interpolate to the rest.

      do n=1,nc
         do i=flo(0),fhi(0)-ratio(dir),ratio(0)
            df = fine(i+ratio(dir),n)-fine(i,n)
            do M=1,ratio(dir)-1
               val = fine(i,n) + df*dble(M)/dble(ratio(dir))
               fine(i+M,n) = val
            enddo                     
         enddo
      enddo

      end subroutine ca_edge_interp

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine ca_pc_edge_interp(lo, hi, nc, ratio, dir, &
           crse, crse_l0, crse_h0,  &
           fine, fine_l0, fine_h0)
      implicit none
      integer lo(1),hi(1), nc, ratio(0:2-1), dir
      integer crse_l0, crse_h0
      integer fine_l0, fine_h0
      double precision crse(crse_l0:crse_h0,nc)
      double precision fine(fine_l0:fine_h0,nc)
      integer i,ii,n

!     For edge-based data, fill fine values with piecewise-constant interp of coarse data.
!     Operate only on faces that overlap--ie, only fill the fine faces that make up each
!     coarse face, leave the in-between faces alone.
      do n=1,nc
         do i=lo(1),hi(1)
            ii = ratio(0)*i
            fine(ii,n) = crse(i,n)
         enddo
      enddo

      end subroutine ca_pc_edge_interp

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine ca_test_residual(lo, hi, &
           rhs, rhl1, rhh1,  &
           ecx, ecxl1, ecxh1, &
           dx, problo, coord_type)

      use bl_constants_module

      implicit none

      integer          :: lo(1),hi(1)
      integer          :: rhl1, rhh1
      integer          :: ecxl1, ecxh1
      integer          :: coord_type
      double precision :: rhs(rhl1:rhh1)
      double precision :: ecx(ecxl1:ecxh1)
      double precision :: dx(1),problo(1)

      double precision :: lapphi
      double precision :: rlo,rhi,rcen
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
         call bl_error("Error:: Gravity_1d.f90 :: ca_test_residual")
      end if

      end subroutine ca_test_residual

! :::
! ::: ----------------------------------------------------------------
! :::

      ! Note that we come into this routine with the full 1D density
      ! array, so we can compute the gravity in one pass.

      subroutine ca_compute_1d_grav(rho, r_l1, r_h1, lo, hi, grav, phi, dx, problo)

      use fundamental_constants_module, only : Gconst
      use bl_constants_module
      use meth_params_module, only : get_g_from_phi

      implicit none

      integer         , intent(in   ) :: r_l1, r_h1
      integer         , intent(in   ) :: lo,   hi
      double precision, intent(in   ) ::  rho(r_l1:r_h1)
      double precision, intent(  out) :: grav(r_l1:r_h1)
      double precision, intent(  out) ::  phi(r_l1:r_h1)
      double precision, intent(in   ) :: dx, problo(1)

      double precision :: phi_temp(r_l1-1:r_h1+1)
      
      double precision, parameter ::  fourthirdspi = FOUR3RD * M_PI
      double precision :: rc,rlo,mass_encl,halfdx,dm,rloj,rcj,rhij
      integer          :: i,j

      halfdx = HALF * dx

      if (get_g_from_phi) then

         phi_temp = ZERO
         grav = ZERO

         ! First do all the zones in the physical domain and in the
         ! upper ghost cells, using the standard approach of integrating
         ! the Green's function for the potential.
         
         do i = lo, r_h1+1

            rlo = problo(1) + dble(i) * dx
            rc = rlo + halfdx

            mass_encl = ZERO                  
            
            do j = lo, hi

               rloj = problo(1) + dble(j) * dx
               rcj = rloj + halfdx
               rhij = rcj + halfdx

               dm = fourthirdspi * (rhij**3 - rloj**3) * rho(j)

               mass_encl = mass_encl + dm
               
               ! If the mass shell is interior to us, the shell theorem
               ! (or an expansion in in the potential) tells us
               ! that its contribution to the potential is given by
               ! a point mass located at the origin.
                  
               if (j .lt. i) then

                  phi_temp(i) = phi_temp(i) + Gconst * dm / rc

               ! If the mass shell is exterior, the potential is G * M / R where
               ! R is the radius of the shell.
                     
               else if (j .gt. i) then

                  phi_temp(i) = phi_temp(i) + Gconst * dm / rcj

               endif

            enddo

         enddo

         ! We want to do even reflection of phi for the lower ghost cells on a
         ! symmetry axis, to ensure that the gradient at r == 0 vanishes.
         
         if (problo(1) .eq. ZERO) then
            do i = r_l1-1, lo-1
               phi_temp(i) = phi_temp(-i-1)
            enddo
         endif

         ! For the outermost zones, use phi = G * M / r again.

         do i = hi+1, r_h1
            rc = problo(1) + dble(i) * dx + halfdx
            phi_temp(i) = Gconst * mass_encl / rc
         enddo

         ! Now that we have phi, construct g by taking the gradient.
         ! We use simple second-order centered differencing.

         do i = r_l1, r_h1

            grav(i) = (phi_temp(i+1) - phi_temp(i-1)) / (TWO * dx)

         enddo

         phi = phi_temp(r_l1:r_h1)         
         
      else

         mass_encl = ZERO      
         
         do i = 0, r_h1
            rlo = problo(1) + dble(i) * dx
            rc  = rlo + halfdx
            if (i .gt. 0) then
               dm = fourthirdspi * halfdx * (rlo**2 + rlo*(rlo-halfdx) + (rlo-halfdx)**2) * rho(i-1) + &
                    fourthirdspi * halfdx * ( rc**2 +  rc* rlo         +  rlo**2        ) * rho(i  )
            else
               dm = fourthirdspi * halfdx * (rc**2 + rc*rlo +  rlo**2) * rho(i)
            endif
            mass_encl = mass_encl + dm
            grav(i) = -Gconst * mass_encl / rc**2

         enddo

         if (problo(1) .eq. ZERO) then
            do i = r_l1,-1
               grav(i) = -grav(-i-1)
            enddo
         endif

      endif
         
      end subroutine ca_compute_1d_grav
