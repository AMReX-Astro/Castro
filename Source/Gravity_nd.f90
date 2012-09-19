
! ::
! :: ----------------------------------------------------------
! ::

      subroutine get_grav_const(Gconst_out)

         use fundamental_constants_module, only: Gconst

         double precision :: Gconst_out

         Gconst_out = Gconst

      end subroutine get_grav_const

! ::
! :: ----------------------------------------------------------
! ::

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

      subroutine ca_integrate_phi (mass,grav,phi,dr,numpts_1d)

      use fundamental_constants_module, only : Gconst

      implicit none
      integer          :: numpts_1d
      double precision :: mass(0:numpts_1d-1)
      double precision :: grav(0:numpts_1d-1)
      double precision ::  phi(0:numpts_1d-1)
      double precision :: dr

      integer          :: i
      double precision :: mass_encl,rhi

      mass_encl = 0.d0
      grav(0)   = 0.d0
      do i = 1,numpts_1d-1
         rhi = dble(i) * dr
         mass_encl = mass_encl + mass(i-1)
         grav(i) = -Gconst * mass_encl / rhi**2
      enddo

      phi(0) = 0.d0
      do i = 1,numpts_1d-1
        phi(i) = phi(i-1) + grav(i) * dr
      enddo

      end subroutine ca_integrate_phi

! ::
! :: ----------------------------------------------------------
! ::

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
