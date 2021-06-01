module gravity_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains

  subroutine ca_integrate_grav (mass,den,grav,max_radius,dr,numpts_1d) &
       bind(C, name="ca_integrate_grav")
    ! Given a radial mass distribution, this computes the gravitational
    ! acceleration as a function of radius by computing the mass enclosed
    ! in successive spherical shells.

    use fundamental_constants_module, only : Gconst
    use amrex_constants_module

    implicit none
    integer , intent(in   ) :: numpts_1d
    real(rt), intent(in   ) :: mass(0:numpts_1d-1)   ! radial mass distribution
    real(rt), intent(in   ) ::  den(0:numpts_1d-1)
    real(rt), intent(inout) :: grav(0:numpts_1d-1)   ! gravitational acceleration as a function of radius
    real(rt), intent(in   ) :: max_radius,dr

    integer          :: i
    real(rt)         :: rc,rlo,rhi,halfdr
    real(rt)         :: mass_encl
    real(rt)         :: vol_inner_shell, vol_outer_shell
    real(rt)         :: vol_lower_shell, vol_upper_shell
    real(rt)         :: vol_total_im1, vol_total_i

    real(rt)        , parameter ::  fourthirdspi = 4.e0_rt * M_PI / 3.e0_rt

    halfdr = 0.5e0_rt * dr

    mass_encl = 0.e0_rt
    do i = 0,numpts_1d-1
       rlo = (dble(i)      ) * dr
       rc  = (dble(i)+0.5e0_rt) * dr
       rhi = (dble(i)+1.0e0_rt) * dr

       if (i.eq.0) then

          ! The mass at (i) is distributed into these two regions
          vol_outer_shell = fourthirdspi * rc**3
          vol_upper_shell = fourthirdspi * (rhi**3 - rc**3)
          vol_total_i     = vol_outer_shell + vol_upper_shell

          mass_encl = vol_outer_shell * mass(i)  / vol_total_i

       else if (rc .lt. max_radius) then

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

       else

          ! The mass at (i-1) is distributed into these two shells
          vol_lower_shell = vol_outer_shell   ! This copies from the previous i
          vol_inner_shell = vol_upper_shell   ! This copies from the previous i
          vol_total_im1   = vol_total_i       ! This copies from the previous i

          ! The mass at (i)   is distributed into these two shells
          vol_outer_shell = fourthirdspi * halfdr * ( rc**2 + rlo*rc + rlo**2)
          vol_upper_shell = fourthirdspi * halfdr * ( rc**2 + rhi*rc + rhi**2)
          vol_total_i     = vol_outer_shell + vol_upper_shell

          mass_encl = mass_encl + vol_inner_shell*den(i-1) + vol_outer_shell*den(i  )
       end if

       grav(i) = -Gconst * mass_encl / rc**2
       !        print *,'GRAV MASS ',rc, mass_encl

    enddo

  end subroutine ca_integrate_grav



  subroutine ca_integrate_gr_grav (rho,mass,pres,grav,dr,numpts_1d) &
       bind(C, name="ca_integrate_gr_grav")
    ! Same as ca_integrate_grav above, but includes general relativistic effects.

    use fundamental_constants_module, only : Gconst, c_light
    use amrex_constants_module

    implicit none
    integer , intent(in   ) :: numpts_1d
    real(rt), intent(in   ) ::  rho(0:numpts_1d-1)
    real(rt), intent(in   ) :: mass(0:numpts_1d-1)
    real(rt), intent(in   ) :: pres(0:numpts_1d-1)
    real(rt), intent(inout) :: grav(0:numpts_1d-1)
    real(rt), intent(in   ):: dr

    integer          :: i
    real(rt)         :: mass_encl,rc,rlo,rhi,halfdr
    real(rt)         :: ga, gb, gc, P,R
    real(rt)         :: vol_inner_shell, vol_outer_shell
    real(rt)         :: vol_lower_shell, vol_upper_shell
    real(rt)         :: vol_total_im1, vol_total_i

    real(rt)        , parameter ::  fourpi       = 4.e0_rt * M_PI
    real(rt)        , parameter ::  fourthirdspi = 4.e0_rt * M_PI / 3.e0_rt
    real(rt)        , parameter ::  sqvc         = c_light**2

    halfdr = 0.5e0_rt * dr

    mass_encl = 0.e0_rt
    do i = 0,numpts_1d-1
       rlo = (dble(i)      ) * dr
       rc  = (dble(i)+0.5e0_rt) * dr
       rhi = (dble(i)+1.0e0_rt) * dr

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

       !       !!!!!!!!!!!!!!!!!!!!!!!!!
       !       This adds the post-Newtonian correction
       !       Added by Ken Chen, 6/9 2010
       !       !!!!!!!!!!!!!!!!!!!!!!!!!

       !       Tolman-Oppenheimer-Volkoff(TOV) case

       if (rho(i) .gt. 0.e0_rt) then
          P =  pres(i)
          R =  rho(i)
          ga = (1.e0_rt + P/(R*sqvc))
          gb = (1.e0_rt + fourpi * rc**3 * P / (mass_encl*sqvc))
          gc = 1.e0_rt / (1.e0_rt - 2.e0_rt * Gconst * mass_encl / (rc*sqvc))

          grav(i) = grav(i)*ga*gb*gc
       end if

       !       !!!!!!!!!!!!!!!!!!!!!!!!!
       !       This ends the post-Newtonian correction
       !       !!!!!!!!!!!!!!!!!!!!!!!!!

    enddo

  end subroutine ca_integrate_gr_grav

end module gravity_module
