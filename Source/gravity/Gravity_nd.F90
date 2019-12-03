module gravity_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

  ! Data for the multipole gravity

  real(rt), allocatable :: volumeFactor, parityFactor
  real(rt), parameter :: edgeTolerance = 1.0e-2_rt
  real(rt), allocatable :: rmax
  logical,  allocatable :: doSymmetricAddLo(:), doSymmetricAddHi(:), doSymmetricAdd
  logical,  allocatable :: doReflectionLo(:), doReflectionHi(:)
  integer,  allocatable :: lnum_max
  real(rt), allocatable :: factArray(:,:)
  real(rt), allocatable :: parity_q0(:), parity_qC_qS(:,:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: volumeFactor, parityFactor
  attributes(managed) :: rmax
  attributes(managed) :: doSymmetricAddLo, doSymmetricAddHi, doSymmetricAdd
  attributes(managed) :: doReflectionLo, doReflectionHi
  attributes(managed) :: lnum_max
  attributes(managed) :: factArray
  attributes(managed) :: parity_q0, parity_qC_qS
#endif

contains

  ! ::
  ! :: ----------------------------------------------------------
  ! ::


  subroutine get_grav_const(Gconst_out) bind(C, name="get_grav_const")
    ! Returns the gravitational constant, G
    !
    ! .. note::
    !    Binds to C function ``get_grav_const``

    use fundamental_constants_module, only: Gconst

    real(rt), intent(inout) :: Gconst_out

    Gconst_out = Gconst

  end subroutine get_grav_const

  ! ::
  ! :: ----------------------------------------------------------
  ! ::

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
    ! Integrates radial mass elements of a spherically symmetric
    ! mass distribution to calculate both the gravitational acceleration,
    ! grav, and the gravitational potential, phi. Here the mass variable
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

    use fundamental_constants_module, only : Gconst

    implicit none
    integer , intent(in   ) :: numpts_1d   ! number of points in radial direction
    real(rt), intent(in   ) :: mass(0:numpts_1d-1)   ! radial mass distribution
    real(rt), intent(inout) :: grav(0:numpts_1d-1)   ! radial gravitational acceleration
    real(rt), intent(inout) :: phi(0:numpts_1d-1)   ! radial gravitational potential
    real(rt), intent(in   ) :: dr   ! radial cell spacing

    real(rt)         :: gravBC, phiBC
    integer          :: i
    real(rt)         :: mass_encl,rhi

    mass_encl = 0.e0_rt
    grav(0)   = 0.e0_rt
    do i = 1,numpts_1d-1
       rhi = dble(i) * dr
       mass_encl = mass_encl + mass(i-1)
       grav(i) = -Gconst * mass_encl / rhi**2
    enddo

    mass_encl = mass_encl + mass(numpts_1d-1)
    phiBC = -Gconst * mass_encl / (numpts_1d*dr)
    gravBC = -Gconst * mass_encl / (numpts_1d*dr)**2
    phi(numpts_1d-1) = phiBC + gravBC * dr

    do i = numpts_1d-2,0,-1
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

  ! ::
  ! :: ----------------------------------------------------------
  ! ::

  subroutine init_multipole_gravity(lnum, lo_bc, hi_bc) bind(C,name="init_multipole_gravity")
    ! If any of the boundaries are symmetric, we need to account for the mass that is assumed
    ! to lie on the opposite side of the symmetric axis. If the center in any direction
    ! coincides with the boundary, then we can simply double the mass as a result of that reflection.
    ! Otherwise, we need to do a more general solve. We include a logical that is set to true
    ! if any boundary is symmetric, so that we can avoid unnecessary function calls.
    !
    ! .. note::
    !    Binds to C function ``init_multipole_gravity``

    use amrex_constants_module
    use prob_params_module, only: coord_type, Symmetry, problo, probhi, center, dim

    implicit none

    integer, intent(in) :: lnum, lo_bc(3), hi_bc(3)

    integer :: b, l, m

    allocate(volumeFactor)
    allocate(parityFactor)

    volumeFactor = ONE
    parityFactor = ONE

    allocate(doSymmetricAddLo(3))
    allocate(doSymmetricAddHi(3))

    doSymmetricAddLo(:) = .false.
    doSymmetricAddHi(:) = .false.

    allocate(doSymmetricAdd)

    doSymmetricAdd      = .false.

    allocate(doReflectionLo(3))
    allocate(doReflectionHi(3))
    
    doReflectionLo(:)   = .false.
    doReflectionHi(:)   = .false.

    do b = 1, 3

       if ( (lo_bc(b) .eq. Symmetry) .and. (coord_type .eq. 0) ) then
          if ( abs(center(b) - problo(b)) < edgeTolerance ) then
             volumeFactor = volumeFactor * TWO
             doReflectionLo(b) = .true.
          else
             doSymmetricAddLo(b) = .true.
             doSymmetricAdd      = .true.
          endif
       endif

       if ( (hi_bc(b) .eq. Symmetry) .and. (coord_type .eq. 0) ) then
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

    allocate(lnum_max)

    lnum_max = lnum

    allocate(factArray(0:lnum_max, 0:lnum_max))
    allocate(parity_q0(0:lnum_max))
    allocate(parity_qC_qS(0:lnum_max, 0:lnum_max))

    factArray(:,:) = ZERO
    parity_q0 = ONE
    parity_qC_qS = ONE

    do l = 0, lnum_max

       ! The odd l Legendre polynomials are odd in their argument, so
       ! a symmetric reflection about the z axis leads to a total cancellation.

       parity_q0(l) = ONE

       if ( MODULO(l,2) /= 0 ) then
          if ( dim .eq. 3 .and. ( doReflectionLo(3) .or. doReflectionHi(3) ) ) then
             parity_q0(l) = ZERO
          else if ( dim .eq. 2 .and. coord_type .eq. 1 ) then
             parity_q0(l) = ZERO
          endif
       endif

       do m = 1, l

          ! The parity properties of the associated Legendre polynomials are:
          ! P_l^m (-x) = (-1)^(l+m) P_l^m (x)
          ! Therefore, a complete cancellation occurs if l+m is odd and
          ! we are reflecting about the z axis.

          ! Additionally, the cosine and sine terms flip sign when reflected
          ! about the x or y axis, so if we have a reflection about x or y
          ! then the terms have a complete cancellation.

          if (dim .eq. 3) then
             parity_qC_qS(l,m) = ONE
          else if (dim .eq. 2 .and. coord_type .eq. 1) then
             parity_qC_qS(l,m) = ZERO
          endif

          if ( MODULO(l+m,2) /= 0 .and. ( doReflectionLo(3) .or. doReflectionHi(3) ) ) then
             parity_qC_qS(l,m) = ZERO
          endif

          if ( doReflectionLo(1) .or. doReflectionLo(2) .or. doReflectionHi(1) .or. doReflectionHi(2) ) then
             parity_qC_qS(l,m) = ZERO
          endif

          factArray(l,m) = TWO * factorial(l-m) / factorial(l+m) * volumeFactor

       enddo

    enddo

    ! Now let's take care of a safety issue. The multipole calculation involves taking powers of r^l,
    ! which can overflow the real(rt)         exponent limit if lnum is very large. Therefore,
    ! we will normalize all distances to the maximum possible physical distance from the center,
    ! which is the diagonal from the center to the edge of the box. Then r^l will always be
    ! less than or equal to one. For large enough lnum, this may still result in roundoff
    ! errors that don't make your answer any more precise, but at least it avoids
    ! possible NaN issues from having numbers that are too large for real(rt)        .
    ! We will put the rmax factor back in at the end of ca_put_multipole_phi.

    allocate(rmax)

    rmax = HALF * maxval(probhi(1:dim) - problo(1:dim)) * sqrt(dble(dim))

  end subroutine init_multipole_gravity



  subroutine ca_put_multipole_phi (lo, hi, domlo, domhi, dx, &
                                   phi, p_lo, p_hi, &
                                   lnum, &
                                   qL0, qLC, qLS, qU0, qUC, qUS, &
                                   npts, boundary_only) &
                                   bind(C, name="ca_put_multipole_phi")

#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use prob_params_module, only: problo, center, dim, coord_type
    use fundamental_constants_module, only: Gconst
    use amrex_constants_module

    implicit none

    integer , intent(in   ) :: lo(3), hi(3)
    integer , intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3)

    integer , intent(in   ), value :: lnum, npts, boundary_only
    real(rt), intent(in   ) :: qL0(0:lnum,0:npts-1), qLC(0:lnum,0:lnum,0:npts-1), qLS(0:lnum,0:lnum,0:npts-1)
    real(rt), intent(in   ) :: qU0(0:lnum,0:npts-1), qUC(0:lnum,0:lnum,0:npts-1), qUS(0:lnum,0:lnum,0:npts-1)

    integer , intent(in   ) :: p_lo(3), p_hi(3)
    real(rt), intent(inout) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))

    integer  :: i, j, k
    integer  :: l, m, n, nlo
    real(rt) :: x, y, z, r, cosTheta, phiAngle
    real(rt) :: legPolyL, legPolyL1, legPolyL2
    real(rt) :: assocLegPolyLM, assocLegPolyLM1, assocLegPolyLM2
    real(rt) :: r_U
    real(rt) :: rmax_cubed

    !$gpu

    ! If we're using this to construct boundary values, then only use
    ! the outermost bin.

    if (boundary_only .eq. 1) then
       nlo = npts-1
    else
       nlo = 0
    endif

#ifndef AMREX_USE_CUDA
    if (lnum > lnum_max) then
       call castro_error("Error: ca_put_multipole_phi: requested more multipole moments than we allocated data for.")
    endif
#endif

    rmax_cubed = rmax**3

    do k = lo(3), hi(3)
       if (k .gt. domhi(3)) then
          z = problo(3) + (dble(k  )     ) * dx(3) - center(3)
       else if (k .lt. domlo(3)) then
          z = problo(3) + (dble(k+1)     ) * dx(3) - center(3)
       else
          z = problo(3) + (dble(k  )+HALF) * dx(3) - center(3)
       end if

       z = z / rmax

       do j = lo(2), hi(2)
          if (j .gt. domhi(2)) then
             y = problo(2) + (dble(j  )     ) * dx(2) - center(2)
          else if (j .lt. domlo(2)) then
             y = problo(2) + (dble(j+1)     ) * dx(2) - center(2)
          else
             y = problo(2) + (dble(j  )+HALF) * dx(2) - center(2)
          end if

          y = y / rmax

          do i = lo(1), hi(1)
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
                ! These cells should not be accessed anyway during the gravity solve.

                r = sqrt( x**2 + y**2 + z**2 )

                if ( r < 1.0e-12_rt ) then
                   phi(i,j,k) = ZERO
                   cycle
                endif

                if (dim .eq. 3) then
                   cosTheta = z / r
                   phiAngle = atan2(y,x)
                else if (dim .eq. 2 .and. coord_type .eq. 1) then
                   cosTheta = y / r
                   phiAngle = ZERO
                endif

                phi(i,j,k) = ZERO

                ! Compute the potentials on the ghost cells.

                do n = nlo, npts-1

                   do l = 0, lnum

                      call calcLegPolyL(l, legPolyL, legPolyL1, legPolyL2, cosTheta)

                      r_U = r**dble(-l-1)

                      ! Make sure we undo the volume scaling here.

                      phi(i,j,k) = phi(i,j,k) + qL0(l,n) * legPolyL * r_U * rmax_cubed

                   end do

                   do m = 1, lnum
                      do l = 1, lnum

                         if (m > l) continue

                         call calcAssocLegPolyLM(l, m, assocLegPolyLM, assocLegPolyLM1, assocLegPolyLM2, cosTheta)

                         r_U = r**dble(-l-1)

                         ! Make sure we undo the volume scaling here.

                         phi(i,j,k) = phi(i,j,k) + (qLC(l,m,n) * cos(m * phiAngle) + qLS(l,m,n) * sin(m * phiAngle)) * &
                                                   assocLegPolyLM * r_U * rmax_cubed

                      enddo

                   enddo

                enddo

                phi(i,j,k) = -Gconst * phi(i,j,k) / rmax

             endif

          enddo
       enddo
    enddo

  end subroutine ca_put_multipole_phi



  subroutine ca_compute_multipole_moments(lo, hi, domlo, domhi, &
                                          dx, rho, r_lo, r_hi, &
                                          vol, v_lo, v_hi, &
                                          lnum, &
                                          qL0, qLC, qLS, qU0, qUC, qUS, &
                                          npts, boundary_only) &
                                          bind(C, name="ca_compute_multipole_moments")

#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use prob_params_module, only: problo, center, probhi, dim, coord_type
    use amrex_constants_module

    implicit none

    integer , intent(in   ) :: lo(3), hi(3)
    integer , intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3)
    integer , intent(in   ), value :: boundary_only, npts, lnum

    real(rt), intent(inout) :: qL0(0:lnum,0:npts-1), qLC(0:lnum,0:lnum,0:npts-1), qLS(0:lnum,0:lnum,0:npts-1)
    real(rt), intent(inout) :: qU0(0:lnum,0:npts-1), qUC(0:lnum,0:lnum,0:npts-1), qUS(0:lnum,0:lnum,0:npts-1)

    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    real(rt), intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer  :: i, j, k
    integer  :: nlo, index

    real(rt) :: x, y, z, r, drInv, cosTheta, phiAngle
    real(rt) :: rmax_cubed_inv

    !$gpu

    ! If we're using this to construct boundary values, then only fill
    ! the outermost bin.

    if (boundary_only .eq. 1) then
       nlo = npts-1
    else
       nlo = 0
    endif

    ! Note that we don't currently support dx != dy != dz, so this is acceptable.

    drInv = rmax / dx(1)

    ! Sanity check

#ifndef AMREX_USE_CUDA
    if (lnum > lnum_max) then
       call castro_error("Error: ca_compute_multipole_moments: requested more multipole moments than we allocated data for.")
    endif
#endif

    rmax_cubed_inv = ONE / rmax**3

    do k = lo(3), hi(3)
       z = ( problo(3) + (dble(k)+HALF) * dx(3) - center(3) ) / rmax

       do j = lo(2), hi(2)
          y = ( problo(2) + (dble(j)+HALF) * dx(2) - center(2) ) / rmax

          do i = lo(1), hi(1)
             x = ( problo(1) + (dble(i)+HALF) * dx(1) - center(1) ) / rmax

             r = sqrt( x**2 + y**2 + z**2 )

             if (dim .eq. 3) then
                index = int(r * drInv)
                cosTheta = z / r
                phiAngle = atan2(y, x)
             else if (dim .eq. 2 .and. coord_type .eq. 1) then
                index = nlo ! We only do the boundary potential in 2D.
                cosTheta = y / r
                phiAngle = z
             endif

             ! Now, compute the multipole moments.

             call multipole_add(cosTheta, phiAngle, r, rho(i,j,k), vol(i,j,k) * rmax_cubed_inv, &
                                qL0, qLC, qLS, qU0, qUC, qUS, lnum, npts, nlo, index, .true.)

             ! Now add in contributions if we have any symmetric boundaries in 3D.
             ! The symmetric boundary in 2D axisymmetric is handled separately.

             if ( doSymmetricAdd ) then

                call multipole_symmetric_add(doSymmetricAddLo, doSymmetricAddHi, &
                                             x, y, z, problo, probhi, &
                                             rho(i,j,k), vol(i,j,k) * rmax_cubed_inv, &
                                             qL0, qLC, qLS, qU0, qUC, qUS, &
                                             lnum, npts, nlo, index)

             endif

          enddo
       enddo
    enddo

  end subroutine ca_compute_multipole_moments



  function factorial(n) result(fact)

    use amrex_constants_module, only: ONE

    implicit none

    integer :: n, i
    real(rt) :: fact

    fact = ONE

    do i = 2, n
       fact = fact * dble(i)
    enddo

  end function factorial



  subroutine calcAssocLegPolyLM(l, m, assocLegPolyLM, assocLegPolyLM1, assocLegPolyLM2, x)

    use amrex_constants_module, only: ONE, TWO

    implicit none

    integer,  intent(in   ) :: l, m
    real(rt), intent(inout) :: assocLegPolyLM, assocLegPolyLM1, assocLegPolyLM2
    real(rt), intent(in   ) :: x

    integer :: n

    !$gpu

    ! Calculate the associated Legendre polynomials. There are a number of
    ! recurrence relations, but many are unstable. We'll use one that is known
    ! to be stable for the reasonably low values of l we care about in a simulation:
    ! (l-m)P_l^m(x) = x(2l-1)P_{l-1}^m(x) - (l+m-1)P_{l-2}^m(x).
    ! This uses the following two expressions as initial conditions:
    ! P_m^m(x) = (-1)^m (2m-1)! (1-x^2)^(m/2)
    ! P_{m+1}^m(x) = x (2m+1) P_m^m (x)

    if (l == m) then

       ! P_m^m

       assocLegPolyLM = (-1)**m * ((ONE - x) * (ONE + x))**(dble(m)/TWO)

       do n = (2*m-1), 3, -2

          assocLegPolyLM = assocLegPolyLM * n

       end do

    else if (l == m + 1) then

       ! P_{m+1}^m

       assocLegPolyLM1 = assocLegPolyLM
       assocLegPolyLM  = x * (2*m + 1) * assocLegPolyLM1

    else

       assocLegPolyLM2 = assocLegPolyLM1
       assocLegPolyLM1 = assocLegPolyLM
       assocLegPolyLM  = (x * (2*l - 1) * assocLegPolyLM1 - (l + m - 1) * assocLegPolyLM2) / (l-m)

    end if

  end subroutine calcAssocLegPolyLM



  subroutine calcLegPolyL(l, legPolyL, legPolyL1, legPolyL2, x)

    use amrex_constants_module, only: ONE

    implicit none

    integer,  intent(in   ) :: l
    real(rt), intent(inout) :: legPolyL, legPolyL1, legPolyL2
    real(rt), intent(in   ) :: x

    !$gpu

    ! Calculate the Legendre polynomials. We use a stable recurrence relation:
    ! (l+1) P_{l+1}(x) = (2l+1) x P_l(x) - l P_{l-1}(x).
    ! This uses initial conditions:
    ! P_0(x) = 1
    ! P_1(x) = x

    if (l == 0) then

       legPolyL = ONE

    elseif (l == 1) then

       legPolyL1 = legPolyL
       legPolyL  = x

    else

       legPolyL2 = legPolyL1
       legPolyL1 = legPolyL
       legPolyL  = ((2*l - 1) * x * legPolyL1 - (l-1) * legPolyL2) / l

    end if

  end subroutine calcLegPolyL


  subroutine multipole_symmetric_add(doSymmetricAddLo, doSymmetricAddHi, &
                                     x, y, z, problo, probhi, &
                                     rho, vol, &
                                     qU0, qUC, qUS, qL0, qLC, qLS, &
                                     lnum, npts, nlo, index)

    use prob_params_module, only: center
    use amrex_constants_module

    implicit none

    integer,  intent(in) :: lnum, npts, nlo, index
    real(rt), intent(in) :: x, y, z
    real(rt), intent(in) :: problo(3), probhi(3)
    real(rt), intent(in) :: rho, vol

    logical,  intent(in) :: doSymmetricAddLo(3), doSymmetricAddHi(3)

    real(rt), intent(inout) :: qL0(0:lnum,0:npts-1)
    real(rt), intent(inout) :: qLC(0:lnum,0:lnum,0:npts-1), qLS(0:lnum,0:lnum,0:npts-1)

    real(rt), intent(inout) :: qU0(0:lnum,0:npts-1)
    real(rt), intent(inout) :: qUC(0:lnum,0:lnum,0:npts-1), qUS(0:lnum,0:lnum,0:npts-1)

    real(rt) :: cosTheta, phiAngle, r
    real(rt) :: xLo, yLo, zLo, xHi, yHi, zHi

    !$gpu

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

       call multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, lnum, npts, nlo, index)

       if ( doSymmetricAddLo(2) ) then

          r        = sqrt( xLo**2 + yLo**2 + z**2 )
          phiAngle = atan2(yLo, xLo)
          cosTheta = z / r

          call multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, lnum, npts, nlo, index)

       endif

       if ( doSymmetricAddLo(3) ) then

          r        = sqrt( xLo**2 + y**2 + zLo**2 )
          phiAngle = atan2(y, xLo)
          cosTheta = zLo / r

          call multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, lnum, npts, nlo, index)

       endif

       if ( doSymmetricAddLo(2) .and. doSymmetricAddLo(3) ) then

          r        = sqrt( xLo**2 + yLo**2 + zLo**2 )
          phiAngle = atan2(yLo, xLo)
          cosTheta = zLo / r

          call multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, lnum, npts, nlo, index)

       endif

    endif

    if ( doSymmetricAddLo(2) ) then

       r        = sqrt( x**2 + yLo**2 + z**2 )
       phiAngle = atan2(yLo, x)
       cosTheta = z / r

       call multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, lnum, npts, nlo, index)

       if ( doSymmetricAddLo(3) ) then

          r        = sqrt( x**2 + yLo**2 + zLo**2 )
          phiAngle = atan2(yLo, x)
          cosTheta = zLo / r

          call multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, lnum, npts, nlo, index)

       endif

    endif

    if ( doSymmetricAddLo(3) ) then

       r        = sqrt( x**2 + y**2 + zLo**2 )
       phiAngle = atan2(y, x)
       cosTheta = zLo / r

       call multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, lnum, npts, nlo, index)

    endif

  end subroutine multipole_symmetric_add


  subroutine multipole_add(cosTheta, phiAngle, r, rho, vol, &
                           qL0, qLC, qLS, qU0, qUC, qUS, &
                           lnum, npts, nlo, index, do_parity)

    use amrex_constants_module, only: ONE
    use reduction_module, only: reduce_add

    implicit none

    integer,  intent(in)    :: lnum, npts, nlo, index
    real(rt), intent(in)    :: cosTheta, phiAngle, r, rho, vol

    real(rt), intent(inout) :: qL0(0:lnum,0:npts-1), qLC(0:lnum,0:lnum,0:npts-1), qLS(0:lnum,0:lnum,0:npts-1)
    real(rt), intent(inout) :: qU0(0:lnum,0:npts-1), qUC(0:lnum,0:lnum,0:npts-1), qUS(0:lnum,0:lnum,0:npts-1)

    logical, optional, intent(in)   :: do_parity

    integer :: l, m, n

    real(rt) :: legPolyL, legPolyL1, legPolyL2
    real(rt) :: assocLegPolyLM, assocLegPolyLM1, assocLegPolyLM2

    real(rt) :: rho_r_L, rho_r_U

    real(rt) :: dQ

    logical :: parity

    !$gpu

    parity = .false.

    if (present(do_parity)) then
       parity = do_parity
    endif

    do n = nlo, npts-1

       do l = 0, lnum

          call calcLegPolyL(l, legPolyL, legPolyL1, legPolyL2, cosTheta)

          if (index .le. n) then

             rho_r_L = rho * (r ** dble( l  ))

             dQ = legPolyL * rho_r_L * vol * volumeFactor
             if (parity) then
                dQ = dQ * parity_q0(l)
             end if

             call reduce_add(qL0(l,n), dQ)

          else

             rho_r_U = rho * (r ** dble(-l-1))

             dQ = legPolyL * rho_r_U * vol * volumeFactor
             if (parity) then
                dQ = dQ * parity_q0(l)
             end if
             call reduce_add(qU0(l,n), dQ)

          end if

       end do

       ! For the associated Legendre polynomial loop, we loop over m and then l.
       ! It means that we have to recompute rho_r_L or rho_r_U again, but the
       ! recursion relation we use for the polynomials depends on l being the
       ! innermost loop index.

       do m = 1, lnum
          do l = 1, lnum

             if (m > l) continue

             call calcAssocLegPolyLM(l, m, assocLegPolyLM, assocLegPolyLM1, assocLegPolyLM2, cosTheta)

             if (index .le. n) then

                rho_r_L = rho * (r ** dble( l  ))

                dQ = assocLegPolyLM * cos(m * phiAngle) * rho_r_L * vol * factArray(l,m)
                if (parity) then
                   dQ = dQ * parity_qC_qS(l,m)
                end if
                call reduce_add(qLC(l,m,n), dQ)

                dQ = assocLegPolyLM * sin(m * phiAngle) * rho_r_L * vol * factArray(l,m)
                if (parity) then
                   dQ = dQ * parity_qC_qS(l,m)
                end if
                call reduce_add(qLS(l,m,n), dQ)

             else

                rho_r_U = rho * (r ** dble(-l-1))

                dQ = assocLegPolyLM * cos(m * phiAngle) * rho_r_U * vol * factArray(l,m)
                if (parity) then
                   dQ = dQ * parity_qC_qS(l,m)
                end if
                call reduce_add(qUC(l,m,n), dQ)

                dQ = assocLegPolyLM * sin(m * phiAngle) * rho_r_U * vol * factArray(l,m)
                if (parity) then
                   dQ = dQ * parity_qC_qS(l,m)
                end if
                call reduce_add(qUS(l,m,n), dQ)

             end if

          end do
       end do

    enddo

  end subroutine multipole_add

end module gravity_module
