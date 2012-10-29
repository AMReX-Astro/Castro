
! ::
! :: ----------------------------------------------------------
! ::

      ! Returns the gravitational constant, G

      subroutine get_grav_const(Gconst_out)

         use fundamental_constants_module, only: Gconst

         double precision :: Gconst_out

         Gconst_out = Gconst

      end subroutine get_grav_const

! ::
! :: ----------------------------------------------------------
! ::

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Given a radial mass distribution, this computes the gravitational
      ! acceleration as a function of radius by computing the mass enclosed
      ! in successive spherical shells.
      ! Inputs: mass(r), dr, numpts_1d
      ! Outputs: grav(r)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ca_integrate_grav (mass,grav,dr,numpts_1d)

      use fundamental_constants_module, only : Gconst
      use bl_constants_module

      implicit none
      integer          :: numpts_1d
      double precision :: mass(0:numpts_1d-1)
      double precision :: grav(0:numpts_1d-1)
      double precision :: dr

      integer          :: i
      double precision :: rc,rlo,rhi,halfdr
      double precision :: mass_encl
      double precision :: vol_inner_shell, vol_outer_shell
      double precision :: vol_lower_shell, vol_upper_shell
      double precision :: vol_total_im1, vol_total_i

      double precision, parameter ::  fourthirdspi = 4.d0 * M_PI / 3.d0

      halfdr = 0.5d0 * dr

      mass_encl = 0.d0
      do i = 0,numpts_1d-1
         rlo = (dble(i)      ) * dr
         rc  = (dble(i)+0.5d0) * dr
         rhi = (dble(i)+1.0d0) * dr

         if (i.eq.0) then

            ! The mass at (i) is distributed into these two regions
            vol_outer_shell = fourthirdspi * rc**3 
            vol_upper_shell = fourthirdspi * (rhi**3 - rc**3)
            vol_total_i     = vol_outer_shell + vol_upper_shell

            mass_encl = vol_outer_shell * mass(i)  / vol_total_i

         else

            ! The mass at (i-1) is distributed into these two shells
            vol_lower_shell = vol_outer_shell   ! This copies from the previous i
            vol_inner_shell = vol_upper_shell   ! This copies from the previous i
            vol_total_im1   = vol_total_i       ! This copies from the previous i

            ! The mass at (i)   is distributed into these two shells
            vol_outer_shell = fourthirdspi * halfdr * ( rc**2 + rlo*rc + rlo**2)
            vol_upper_shell = fourthirdspi * halfdr * ( rc**2 + rhi*rc + rhi**2)
            vol_total_i     = vol_outer_shell + vol_upper_shell

            mass_encl = mass_encl + (vol_inner_shell/vol_total_im1) * mass(i-1) + & 
                                    (vol_outer_shell/vol_total_i  ) * mass(i  ) 
         end if
         grav(i) = -Gconst * mass_encl / rc**2
      enddo

      end subroutine ca_integrate_grav

! ::
! :: ----------------------------------------------------------
! ::

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Integrates radial mass elements of a spherically symmetric          
      ! mass distribution to calculate both the gravitational acceleration,  
      ! grav, and the gravitational potential, phi. Here the mass variable   
      ! gives the mass contained in each radial shell.                      
      !                                                                     
      ! The convention in Castro for Poisson's equation is                  
      !                                                                     
      !     laplacian(phi) = -4*pi*G*rho                                    
      !
      ! The gravitational acceleration is then
      !
      !     g(r) = -G*M(r) / r**2
      !
      ! where M(r) is the mass interior to radius r.
      !
      ! The strategy for calculating the potential is to calculate the potential
      ! at the boundary assuming all the mass is enclosed:
      !
      !     phi(R) = G * M / R 
      !
      ! Then, the potential in all other zones can be found using
      !
      !     d(phi)/dr = g    ==>    phi(r < R) = phi(R) - int(g * dr)
      !
      ! Inputs: mass, grav, dr, numpts_1d
      ! Outputs: phi
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ca_integrate_phi (mass,grav,phi,dr,numpts_1d)

      use fundamental_constants_module, only : Gconst

      implicit none
      integer          :: numpts_1d
      double precision :: mass(0:numpts_1d-1)
      double precision :: grav(0:numpts_1d-1)
      double precision :: phi(0:numpts_1d-1)
      double precision :: dr
      double precision :: gravBC, phiBC

      integer          :: i
      double precision :: mass_encl,rhi

      mass_encl = 0.d0
      grav(0)   = 0.d0
      do i = 1,numpts_1d-1
         rhi = dble(i) * dr
         mass_encl = mass_encl + mass(i-1)
         grav(i) = -Gconst * mass_encl / rhi**2
      enddo

      mass_encl = mass_encl + mass(numpts_1d-1)
      phiBC = Gconst * mass_encl / (numpts_1d*dr)
      gravBC = -Gconst * mass_encl / (numpts_1d*dr)**2
      phi(numpts_1d-1) = phiBC - gravBC * dr
       
      do i = numpts_1d-2,0,-1
        phi(i) = phi(i+1) - grav(i+1) * dr
      enddo
      
      end subroutine ca_integrate_phi

! ::
! :: ----------------------------------------------------------
! ::

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Same as ca_integrate_grav above, but includes general relativistic effects.
      ! Inputs: rho, mass, pressure, dr, numpts_1d
      ! Outputs: grav
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ca_integrate_gr_grav (rho,mass,pres,grav,dr,numpts_1d)

      use fundamental_constants_module, only : Gconst
      use bl_constants_module

      implicit none
      integer          :: numpts_1d
      double precision ::  rho(0:numpts_1d-1)
      double precision :: mass(0:numpts_1d-1)
      double precision :: pres(0:numpts_1d-1)
      double precision :: grav(0:numpts_1d-1)
      double precision :: dr

      integer          :: i
      double precision :: mass_encl,rc,rlo,rhi,halfdr
      double precision :: ga, gb, gc, P,R
      double precision :: vol_inner_shell, vol_outer_shell
      double precision :: vol_lower_shell, vol_upper_shell
      double precision :: vol_total_im1, vol_total_i

      double precision, parameter ::  fourpi       = 4.d0 * M_PI
      double precision, parameter ::  fourthirdspi = 4.d0 * M_PI / 3.d0
      double precision, parameter ::  sqvc         = 29979245800.d0**2

      halfdr = 0.5d0 * dr

      mass_encl = 0.d0
      do i = 0,numpts_1d-1
         rlo = (dble(i)      ) * dr
         rc  = (dble(i)+0.5d0) * dr
         rhi = (dble(i)+1.0d0) * dr

         if (i.eq.0) then

            ! The mass at (i) is distributed into these two regions
            vol_outer_shell = fourthirdspi * rc**3 
            vol_upper_shell = fourthirdspi * (rhi**3 - rc**3)
            vol_total_i     = vol_outer_shell + vol_upper_shell

            mass_encl = vol_outer_shell * mass(i)  / vol_total_i

         else

            ! The mass at (i-1) is distributed into these two shells
            vol_lower_shell = vol_outer_shell   ! This copies from the previous i
            vol_inner_shell = vol_upper_shell   ! This copies from the previous i
            vol_total_im1   = vol_total_i       ! This copies from the previous i

            ! The mass at (i)   is distributed into these two shells
            vol_outer_shell = fourthirdspi * halfdr * ( rc**2 + rlo*rc + rlo**2)
            vol_upper_shell = fourthirdspi * halfdr * ( rc**2 + rhi*rc + rhi**2)
            vol_total_i     = vol_outer_shell + vol_upper_shell

            mass_encl = mass_encl + vol_inner_shell / vol_total_im1 * mass(i-1) + & 
                                    vol_outer_shell / vol_total_i   * mass(i  ) 
         end if
         grav(i) = -Gconst * mass_encl / rc**2

!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       This adds the post-Newtonian correction
!!       Added by Ken Chen, 6/9 2010
!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!       Tolman-Oppenheimer-Volkoff(TOV) case

         if (rho(i) .gt. 0.d0) then
            P =  pres(i)
            R =  rho(i)
            ga = (1.d0 + P/(R*sqvc))
            gb = (1.d0 + fourpi * rc**3 * P / (mass_encl*sqvc))
            gc = 1.d0 / (1.d0 - 2.d0 * Gconst * mass_encl / (rc*sqvc))

            grav(i) = grav(i)*ga*gb*gc
         end if

!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       This ends the post-Newtonian correction
!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo

      end subroutine ca_integrate_gr_grav
