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

  ! ::
  ! :: ----------------------------------------------------------
  ! ::

  subroutine ca_integrate_phi (mass,grav,phi,dr,numpts_1d) &
       bind(C, name="ca_integrate_phi")
    ! Integrates a radial gravitational acceleration, grav, to get the
    ! gravitational potential, phi. Here the mass variable
    ! gives the mass contained in each radial shell.
    !
    ! The convention in Castro for Poisson's equation is
    !
    !     laplacian(phi) = 4*pi*G*rho
    !
    ! The gravitational acceleration is then
    ! \f[
    !     g(r) = -G M(r) / r^2
    ! \f]
    ! where \f$M(r)\f$ is the mass interior to radius r.
    !
    ! The strategy for calculating the potential is to calculate the potential
    ! at the boundary assuming all the mass is enclosed:
    ! \f[
    !     \phi(R) = -G M / R
    ! \f]
    ! Then, the potential in all other zones can be found using
    ! \f[
    !     \frac{d(\phi)}{dr} = g   \; \Rightarrow \;
    !     \phi(r < R) = \phi(R) - \int(g \dr)
    ! \f]

    use fundamental_constants_module, only: Gconst

    implicit none

    integer , intent(in   ) :: numpts_1d   ! number of points in radial direction
    real(rt), intent(in   ) :: mass(0:numpts_1d-1)   ! radial mass distribution
    real(rt), intent(inout) :: grav(0:numpts_1d-1)   ! radial gravitational acceleration
    real(rt), intent(inout) :: phi(0:numpts_1d-1)   ! radial gravitational potential
    real(rt), intent(in   ) :: dr   ! radial cell spacing

    real(rt) :: gravBC, phiBC
    integer  :: i
    real(rt) :: mass_encl, rhi

    mass_encl = sum(mass)
    phiBC = -Gconst * mass_encl / (numpts_1d * dr)
    gravBC = -Gconst * mass_encl / (numpts_1d * dr)**2
    phi(numpts_1d-1) = phiBC + gravBC * dr

    do i = numpts_1d-2, 0, -1
       phi(i) = phi(i+1) + grav(i+1) * dr
    enddo

  end subroutine ca_integrate_phi

  ! ::
  ! :: ----------------------------------------------------------
  ! ::

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

          print *, i, grav(i), ga, gb, gc

          grav(i) = grav(i)*ga*gb*gc
       end if

       !       !!!!!!!!!!!!!!!!!!!!!!!!!
       !       This ends the post-Newtonian correction
       !       !!!!!!!!!!!!!!!!!!!!!!!!!

    enddo

  end subroutine ca_integrate_gr_grav



  subroutine ca_compute_avgpres (lo, hi, dx, dr, &
                                 var, r_lo, r_hi, &
                                 radial_pres, problo, &
                                 n1d, drdxfac, level) &
                                 bind(C,name='ca_compute_avgpres')

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: HALF, ONE, FOUR3RD, M_PI
    use meth_params_module, only: NVAR, URHO, UEINT, UTEMP, UFS, UFX
    use prob_params_module, only: center
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos
    use network, only: nspec, naux
    use castro_error_module, only: castro_error

    implicit none

    integer , intent(in   ) :: lo(3), hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: problo(3)
    real(rt), intent(in   ), value :: dr
 
    real(rt), intent(inout) :: radial_pres(0:n1d-1)
    integer , intent(in   ), value :: n1d, drdxfac, level

    integer , intent(in   ) :: r_lo(3), r_hi(3)
    real(rt), intent(in   ) :: var(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),NVAR)

    integer  :: i, j, k, n, index
    integer  :: ii, jj, kk
    real(rt) :: xc, yc, zc, r, rlo, rhi
    real(rt) :: fac, xx, yy, zz, dx_frac, dy_frac, dz_frac, vol
    real(rt) :: lo_i, lo_j, lo_k

    type (eos_t) :: eos_state
    real(rt)     :: rhoInv

    fac     = dble(drdxfac)
    dx_frac = dx(1) / fac
    dy_frac = dx(2) / fac
    dz_frac = dx(3) / fac
    !
    ! Don't OMP this.
    !
    do k = lo(3), hi(3)
       zc = problo(3) + (dble(k) + HALF) * dx(3) - center(3)

       do j = lo(2), hi(2)
          yc = problo(2) + (dble(j) + HALF) * dx(2) - center(2)

          do i = lo(1), hi(1)
             xc = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             r  = sqrt(xc**2 + yc**2 + zc**2)

             index = int(r / dr)

             if (index .gt. n1d-1) then

                if (level .eq. 0) then
                   print *, '   '  
                   print *, '>>> Error: Gravity_nd.F90::ca_compute_avgpres ', i, j, k
                   print *, '>>> ... index too big: ', index, ' > ', n1d-1
                   print *, '>>> ... at (i,j,k)   : ', i, j, k
                   call castro_error("Error:: Gravity_nd.F90 :: ca_compute_avgpres")
                end if

             else

                ! We may be coming in here with a masked out zone (in a zone on a coarse
                ! level underlying a fine level). We don't want to be calling the EOS in
                ! this case, so we'll skip these masked out zones (which will have rho
                ! exactly equal to zero).

                if (var(i,j,k,URHO) == 0.0_rt) cycle

                rhoInv = ONE / var(i,j,k,URHO)

                eos_state % rho = var(i,j,k,URHO)
                eos_state % e   = var(i,j,k,UEINT) * rhoInv
                eos_state % T   = var(i,j,k,UTEMP)
                eos_state % xn  = var(i,j,k,UFS:UFS+nspec-1) * rhoInv
                eos_state % aux = var(i,j,k,UFX:UFX+naux-1) * rhoInv

                ! Compute pressure from the EOS
                call eos(eos_input_re, eos_state)

                lo_i = problo(1) + dble(i) * dx(1) - center(1)
                lo_j = problo(2) + dble(j) * dx(2) - center(2)
                lo_k = problo(3) + dble(k) * dx(3) - center(3)

                do kk = 0, drdxfac-1
                   zz = lo_k + (dble(kk) + HALF) * dz_frac
                   do jj = 0, drdxfac-1
                      yy = lo_j + (dble(jj) + HALF) * dy_frac
                      do ii = 0, drdxfac-1
                         xx = lo_i + (dble(ii) + HALF) * dx_frac

                         r  = sqrt(xx**2 + yy**2 + zz**2)
                         index = int(r / dr)

                         if (index .le. n1d-1) then

                            if (AMREX_SPACEDIM == 1) then

                               ! Spherical
                               rlo = abs(lo_i +  dble(ii  ) * dx_frac)
                               rhi = abs(lo_i +  dble(ii+1) * dx_frac)

                               vol = FOUR3RD * M_PI * (rhi**3 - rlo**3)

                            else if (AMREX_SPACEDIM == 2) then

                               ! Cylindrical
                               rlo = abs(lo_i + dble(ii  ) * dx_frac)
                               rhi = abs(lo_i + dble(ii+1) * dx_frac)

                               vol = (rhi**2 - rlo**2) * dy_frac

                            else

                               ! Cartesian
                               vol = ONE

                            end if

                            radial_pres(index) = radial_pres(index) + vol * eos_state % p

                         end if

                      end do
                   end do
                end do

             end if

          end do
       end do
    end do

  end subroutine ca_compute_avgpres



  subroutine ca_put_radial_phi(lo, hi, &
                               domlo, domhi, &
                               dx, dr, &
                               phi, p_lo, p_hi, &
                               radial_phi, problo, &
                               numpts_1d, fill_interior) &
                               bind(C, name="ca_put_radial_phi")

    use amrex_constants_module, only: HALF, TWO
    use prob_params_module, only: center
    use amrex_fort_module, only: rt => amrex_real
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    implicit none

    integer , intent(in   ) :: lo(3), hi(3)
    integer , intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: p_lo(3), p_hi(3)
    real(rt), intent(in   ) :: radial_phi(0:numpts_1d-1)
    real(rt), intent(in   ) :: dx(3), problo(3)
    real(rt), intent(inout) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    real(rt), intent(in   ), value :: dr
    integer , intent(in   ), value :: numpts_1d, fill_interior

    integer  :: i, j, k, index
    real(rt) :: x, y, z, r
    real(rt) :: cen, xi, slope, phi_lo, phi_md, phi_hi, minvar, maxvar

    !$gpu

    ! Note that when we interpolate into the ghost cells we use the
    ! location of the edge, not the cell center

    do k = lo(3), hi(3)
       if (k .gt. domhi(3)) then
          z = problo(3) + (dble(k  )       ) * dx(3) - center(3)
       else if (k .lt. domlo(3)) then
          z = problo(3) + (dble(k+1)       ) * dx(3) - center(3)
       else
          z = problo(3) + (dble(k  ) + HALF) * dx(3) - center(3)
       end if

       do j = lo(2), hi(2)
          if (j .gt. domhi(2)) then
             y = problo(2) + (dble(j  )       ) * dx(2) - center(2)
          else if (j .lt. domlo(2)) then
             y = problo(2) + (dble(j+1)       ) * dx(2) - center(2)
          else
             y = problo(2) + (dble(j  ) + HALF) * dx(2) - center(2)
          end if

          do i = lo(1), hi(1)
             if (i .gt. domhi(1)) then
                x = problo(1) + (dble(i  )       ) * dx(1) - center(1)
             else if (i .lt. domlo(1)) then
                x = problo(1) + (dble(i+1)       ) * dx(1) - center(1)
             else
                x = problo(1) + (dble(i  ) + HALF) * dx(1) - center(1)
             end if

             r     = sqrt(x**2 + y**2 + z**2)
             index = int(r / dr)

#ifndef AMREX_USE_GPU
             if (index .gt. numpts_1d - 1) then
                print *, 'PUT_RADIAL_PHI: INDEX TOO BIG ', index, ' > ', numpts_1d - 1
                print *, 'AT (i,j) ', i, j, k
                print *, 'R, DR IS ', r, dr
                call castro_error("Error:: Gravity_nd.F90 :: ca_put_radial_phi")
             end if
#endif

             if ((fill_interior .eq. 1) .or. &
                 (i .lt. domlo(1) .or. i .gt. domhi(1) .or. &
                  j .lt. domlo(2) .or. j .gt. domhi(2) .or. &
                  k .lt. domlo(3) .or. k .gt. domhi(3))) then

                cen = (dble(index) + HALF) * dr
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
                   phi_hi = radial_phi(index+1)
                   phi_md = radial_phi(index  )
                   phi_lo = radial_phi(index-1)
                   phi(i,j,k) = &
                        ( phi_hi -      TWO * phi_md + phi_lo) * xi**2 / (TWO * dr**2) + &
                        ( phi_hi                     - phi_lo) * xi    / (TWO * dr   ) + &
                        (-phi_hi + 26.e0_rt * phi_md - phi_lo) / 24.e0_rt
                   minvar     = min(phi_md, min(phi_lo, phi_hi))
                   maxvar     = max(phi_md, max(phi_lo, phi_hi))
                   phi(i,j,k) = max(phi(i,j,k), minvar)
                   phi(i,j,k) = min(phi(i,j,k), maxvar)

                end if

             end if

          end do
       end do
    end do

  end subroutine ca_put_radial_phi

end module gravity_module
