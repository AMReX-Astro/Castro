module advection_util_3d_module

  implicit none

  private

  public normalize_species_fluxes, divu, &
         limit_hydro_fluxes_on_small_dens

contains

  subroutine normalize_species_fluxes(flux1,flux1_lo,flux1_hi, &
                                      flux2,flux2_lo,flux2_hi, &
                                      flux3,flux3_lo,flux3_hi, &
                                      lo, hi)

    ! here we normalize the fluxes of the mass fractions so that
    ! they sum to 0.  This is essentially the CMA procedure that is
    ! defined in Plewa & Muller, 1999, A&A, 342, 179

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    double precision, intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    double precision, intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    double precision, intent(inout) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)

    ! Local variables
    integer          :: i, j, k, n
    double precision :: sum, fac

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux1(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux1(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux1(i,j,k,n) = flux1(i,j,k,n) * fac
             end do
          end do
       end do
    end do

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux2(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux2(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux2(i,j,k,n) = flux2(i,j,k,n) * fac
             end do
          end do
       end do
    end do

    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux3(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux3(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux3(i,j,k,n) = flux3(i,j,k,n) * fac
             end do
          end do
       end do
    end do

  end subroutine normalize_species_fluxes

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine limit_hydro_fluxes_on_small_dens(u,u_lo,u_hi, &
                                              q,q_lo,q_hi, &
                                              vol,vol_lo,vol_hi,   &
                                              lo,hi,dt,dx, &
                                              flux1,flux1_lo,flux1_hi, &
                                              area1,area1_lo,area1_hi, &
                                              flux2,flux2_lo,flux2_hi, &
                                              area2,area2_lo,area2_hi, &
                                              flux3,flux3_lo,flux3_hi, &
                                              area3,area3_lo,area3_hi)

    use bl_constants_module, only: ZERO, HALF, ONE, TWO
    use meth_params_module, only: NVAR, QVAR, URHO, small_dens, cfl
    use prob_params_module, only: dim, coord_type
    use mempool_module, only: bl_allocate, bl_deallocate
    use advection_util_module, only: dflux

    implicit none

    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)

    double precision, intent(in   ) :: dt, dx(3)

    double precision, intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    double precision, intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    double precision, intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    double precision, intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    double precision, intent(in   ) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    double precision, intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    double precision, intent(in   ) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    double precision, intent(inout) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    double precision, intent(in   ) :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))

    double precision, pointer :: thetap(:,:,:), thetam(:,:,:)

    integer          :: i, j, k

    double precision :: alpha_x, alpha_y, alpha_z
    double precision :: rho, drho, fluxLF(NVAR), fluxL(NVAR), fluxR(NVAR), rhoLF, drhoLF, dtdx, theta
    integer          :: dir
    logical          :: include_pressure

    ! q is a temporary array that is designed to hold primitive data. The cell-centered flux
    ! routine expects a QVAR-sized array,
    ! The following algorithm comes from Hu, Adams, and Shu (2013), JCP, 242, 169,
    ! "Positivity-preserving method for high-order conservative schemes solving
    ! compressible Euler equations." It has been modified to enforce not only positivity
    ! but also the stronger requirement that rho > small_dens.

    ! We implement the flux limiter on a dimension-by-dimension basis, starting with the x-direction.

    call bl_allocate(thetap,flux1_lo(1)-1,flux1_hi(1)+1,flux1_lo(2)-1,flux1_hi(2)+1,flux1_lo(3)-1,flux1_hi(3)+1)
    call bl_allocate(thetam,flux1_lo(1)-1,flux1_hi(1)+1,flux1_lo(2)-1,flux1_hi(2)+1,flux1_lo(3)-1,flux1_hi(3)+1)

    thetap(:,:,:) = ONE
    thetam(:,:,:) = ONE

    dir = 1

    dtdx = dt / dx(1)

    ! Whether or not to include pressure in the cell-centered fluxes we calculate
    ! will depend on which dimensionality and coordinate system we are in.

    if (dim .eq. 1 .or. (dim .eq. 2 .and. coord_type .eq. 1)) then
       include_pressure = .false.
    else
       include_pressure = .true.
    endif

    alpha_x = (ONE / dim)
    alpha_y = (ONE / dim)
    alpha_z = (ONE / dim)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1) - 1, hi(1) + 1

             ! Note that this loop includes one ghost zone on either side of the current
             ! bounds, but we only need a one-sided limiter for lo(1)-1 and hi(1)+1.

             ! First we'll do the plus state, which is on the left edge of the zone.

             if (i .ge. lo(1)) then

                ! Obtain the one-sided update to the density, based on Hu et al., Eq. 11.
                ! Note that the sign convention for the notation is opposite to our convention
                ! for the edge states for the flux limiter, that is, the "plus" limiter is on
                ! the left edge of the zone and so is the "minus" rho. The flux limiter convention
                ! is analogous to the convention for the hydro reconstruction edge states.

                rho = u(i,j,k,URHO) + TWO * (dt / alpha_x) * (area1(i,j,k) / vol(i,j,k)) * flux1(i,j,k,URHO)

                if (rho < small_dens) then

                   ! Construct the Lax-Friedrichs flux on the interface (Equation 12).
                   ! Note that we are using the information from Equation 9 to obtain the
                   ! effective maximum wave speed, (|u| + c)_max = CFL / lambda where
                   ! lambda = dt/(dx * alpha); alpha = 1 in 1D and may be chosen somewhat
                   ! freely in multi-D as long as alpha_x + alpha_y + alpha_z = 1.

                   fluxL = dflux(u(i-1,j,k,:), q(i-1,j,k,:), dir, [i-1, j, k], include_pressure)
                   fluxR = dflux(u(i  ,j,k,:), q(i  ,j,k,:), dir, [i  , j, k], include_pressure)
                   fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha_x) * (u(i-1,j,k,:) - u(i,j,k,:)))

                   ! Limit the Lax-Friedrichs flux so that it doesn't cause a density < small_dens.
                   ! To do this, first, construct the density change corresponding to the LF density flux.
                   ! Then, if this update would create a density that is less than small_dens, scale all
                   ! fluxes linearly such that the density flux gives small_dens when applied.

                   drhoLF = TWO * (dt / alpha_x) * (area1(i,j,k) / vol(i,j,k)) * fluxLF(URHO)

                   if (u(i,j,k,URHO) + drhoLF < small_dens) then
                      fluxLF = fluxLF * abs((small_dens - u(i,j,k,URHO)) / drhoLF)
                   endif

                   ! Obtain the final density corresponding to the LF flux.

                   rhoLF = u(i,j,k,URHO) + TWO * (dt / alpha_x) * (area1(i,j,k) / vol(i,j,k)) * fluxLF(URHO)

                   ! Solve for theta from (1 - theta) * rhoLF + theta * rho = small_dens.

                   thetap(i,j,k) = (small_dens - rhoLF) / (rho - rhoLF)

                endif

             endif

             ! Now do the minus state, which is on the right edge of the zone.
             ! This uses the same logic as the above, so we don't replicate the comments.

             if (i .le. hi(1)) then

                rho = u(i,j,k,URHO) - TWO * (dt / alpha_x) * (area1(i+1,j,k) / vol(i,j,k)) * flux1(i+1,j,k,URHO)

                if (rho < small_dens) then

                   fluxL = dflux(u(i  ,j,k,:), q(i  ,j,k,:), dir, [i  , j, k], include_pressure)
                   fluxR = dflux(u(i+1,j,k,:), q(i+1,j,k,:), dir, [i+1, j, k], include_pressure)
                   fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha_x) * (u(i,j,k,:) - u(i+1,j,k,:)))

                   drhoLF = -TWO * (dt / alpha_x) * (area1(i+1,j,k) / vol(i,j,k)) * fluxLF(URHO)

                   if (u(i,j,k,URHO) + drhoLF < small_dens) then
                      fluxLF = fluxLF * abs((small_dens - u(i,j,k,URHO)) / drhoLF)
                   endif

                   rhoLF = u(i,j,k,URHO) - TWO * (dt / alpha_x) * (area1(i+1,j,k) / vol(i,j,k)) * fluxLF(URHO)

                   thetam(i,j,k) = (small_dens - rhoLF) / (rho - rhoLF)

                endif

             endif

          enddo
       enddo
    enddo

    ! Now figure out the limiting values of theta. Each zone center has a thetap and thetam,
    ! but we want a nodal value of theta that is the strongest of the two limiters in each case.
    ! Then, limit the flux accordingly.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1) + 1

             theta = min(thetam(i-1,j,k), thetap(i,j,k))

             fluxL = dflux(u(i-1,j,k,:), q(i-1,j,k,:), dir, [i-1, j, k], include_pressure)
             fluxR = dflux(u(i  ,j,k,:), q(i  ,j,k,:), dir, [i  , j, k], include_pressure)
             fluxLF(:) = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha_x) * (u(i-1,j,k,:) - u(i,j,k,:)))

             drhoLF = TWO * (dt / alpha_x) * (area1(i,j,k) / vol(i,j,k)) * fluxLF(URHO)

             if (u(i,j,k,URHO) + drhoLF < small_dens) then
                fluxLF(:) = fluxLF(:) * abs((small_dens - u(i,j,k,URHO)) / drhoLF)
             else if (u(i-1,j,k,URHO) - drhoLF < small_dens) then
                fluxLF(:) = fluxLF(:) * abs((small_dens - u(i-1,j,k,URHO)) / drhoLF)
             endif

             flux1(i,j,k,:) = (ONE - theta) * fluxLF(:) + theta * flux1(i,j,k,:)

             drho = TWO * (dt / alpha_x) * (area1(i,j,k) / vol(i,j,k)) * flux1(i,j,k,URHO)

             if (u(i,j,k,URHO) + drho < small_dens) then
                 flux1(i,j,k,:) = flux1(i,j,k,:) * abs((small_dens - u(i,j,k,URHO)) / drho)
             else if (u(i-1,j,k,URHO) - drho < small_dens) then
                 flux1(i,j,k,:) = flux1(i,j,k,:) * abs((small_dens - u(i-1,j,k,URHO)) / drho)
             endif

          enddo
       enddo
    enddo

    ! Now do the y-direction. The logic is all the same as for the x-direction,
    ! so the comments are skipped.

    thetap(:,:,:) = ONE
    thetam(:,:,:) = ONE

    dir = 2

    dtdx = dt / dx(2)

    if (dim .ge. 2) then

       do k = lo(3), hi(3)
          do j = lo(2) - 1, hi(2) + 1
             do i = lo(1), hi(1)

                if (j .ge. lo(2)) then

                   rho = u(i,j,k,URHO) + TWO * (dt / alpha_y) * (area2(i,j,k) / vol(i,j,k)) * flux2(i,j,k,URHO)

                   if (rho < small_dens) then

                      fluxL = dflux(u(i,j-1,k,:), q(i,j-1,k,:), dir, [i, j-1, k], include_pressure)
                      fluxR = dflux(u(i,j  ,k,:), q(i,j  ,k,:), dir, [i, j  , k], include_pressure)
                      fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha_y) * (u(i,j-1,k,:) - u(i,j,k,:)))

                      drhoLF = TWO * (dt / alpha_y) * (area2(i,j,k) / vol(i,j,k)) * fluxLF(URHO)

                      if (u(i,j,k,URHO) + drhoLF < small_dens) then
                         fluxLF = fluxLF * abs((small_dens - u(i,j,k,URHO)) / drhoLF)
                      endif

                      rhoLF = u(i,j,k,URHO) + TWO * (dt / alpha_y) * (area2(i,j,k) / vol(i,j,k)) * fluxLF(URHO)

                      thetap(i,j,k) = (small_dens - rhoLF) / (rho - rhoLF)

                   endif

                endif

                if (j .le. hi(2)) then

                   rho = u(i,j,k,URHO) - TWO * (dt / alpha_y) * (area2(i,j+1,k) / vol(i,j,k)) * flux2(i,j+1,k,URHO)

                   if (rho < small_dens) then

                      fluxL = dflux(u(i,j  ,k,:), q(i,j  ,k,:), dir, [i, j  , k], include_pressure)
                      fluxR = dflux(u(i,j+1,k,:), q(i,j+1,k,:), dir, [i, j+1, k], include_pressure)
                      fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha_y) * (u(i,j,k,:) - u(i,j+1,k,:)))

                      drhoLF = -TWO * (dt / alpha_y) * (area2(i,j+1,k) / vol(i,j,k)) * fluxLF(URHO)

                      if (u(i,j,k,URHO) + drhoLF < small_dens) then
                         fluxLF = fluxLF * abs((small_dens - u(i,j,k,URHO)) / drhoLF)
                      endif

                      rhoLF = u(i,j,k,URHO) - TWO * (dt / alpha_y) * (area2(i,j+1,k) / vol(i,j,k)) * fluxLF(URHO)

                      thetam(i,j,k) = (small_dens - rhoLF) / (rho - rhoLF)

                   endif

                endif

             enddo
          enddo
       enddo

       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)

                theta = min(thetam(i,j-1,k), thetap(i,j,k))

                fluxL = dflux(u(i,j-1,k,:), q(i,j-1,k,:), dir, [i, j-1, k], include_pressure)
                fluxR = dflux(u(i,j  ,k,:), q(i,j  ,k,:), dir, [i, j  , k], include_pressure)
                fluxLF(:) = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha_y) * (u(i,j-1,k,:) - u(i,j,k,:)))

                drhoLF = TWO * (dt / alpha_y) * (area2(i,j,k) / vol(i,j,k)) * fluxLF(URHO)

                if (u(i,j,k,URHO) + drhoLF < small_dens) then
                   fluxLF(:) = fluxLF(:) * abs((small_dens - u(i,j,k,URHO)) / drhoLF)
                else if (u(i,j-1,k,URHO) - drhoLF < small_dens) then
                   fluxLF(:) = fluxLF(:) * abs((small_dens - u(i,j-1,k,URHO)) / drhoLF)
                endif

                flux2(i,j,k,:) = (ONE - theta) * fluxLF(:) + theta * flux2(i,j,k,:)

                drho = TWO * (dt / alpha_y) * (area2(i,j,k) / vol(i,j,k)) * flux2(i,j,k,URHO)

                if (u(i,j,k,URHO) + drho < small_dens) then
                   flux2(i,j,k,:) = flux2(i,j,k,:) * abs((small_dens - u(i,j,k,URHO)) / drho)
                else if (u(i,j-1,k,URHO) - drho < small_dens) then
                   flux2(i,j,k,:) = flux2(i,j,k,:) * abs((small_dens - u(i,j-1,k,URHO)) / drho)
                endif

             enddo
          enddo
       enddo

    endif

    ! Now do the z-direction. The logic is all the same as for the x-direction,
    ! so the comments are skipped.

    thetap(:,:,:) = ONE
    thetam(:,:,:) = ONE

    dir = 3

    dtdx = dt / dx(3)

    if (dim .eq. 3) then

       do k = lo(3) - 1, hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (k .ge. lo(3)) then

                   rho = u(i,j,k,URHO) + TWO * (dt / alpha_z) * (area3(i,j,k) / vol(i,j,k)) * flux3(i,j,k,URHO)

                   if (rho < small_dens) then

                      fluxL = dflux(u(i,j,k-1,:), q(i,j,k-1,:), dir, [i, j, k-1], include_pressure)
                      fluxR = dflux(u(i,j,k  ,:), q(i,j,k  ,:), dir, [i, j, k  ], include_pressure)
                      fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha_z) * (u(i,j,k-1,:) - u(i,j,k,:)))

                      drhoLF = TWO * (dt / alpha_z) * (area3(i,j,k) / vol(i,j,k)) * fluxLF(URHO)

                      if (u(i,j,k,URHO) + drhoLF < small_dens) then
                         fluxLF = fluxLF * abs((small_dens - u(i,j,k,URHO)) / drhoLF)
                      endif

                      rhoLF = u(i,j,k,URHO) + TWO * (dt / alpha_z) * (area3(i,j,k) / vol(i,j,k)) * fluxLF(URHO)

                      thetap(i,j,k) = (small_dens - rhoLF) / (rho - rhoLF)

                   endif

                endif

                if (k .le. hi(3)) then

                   rho = u(i,j,k,URHO) - TWO * (dt / alpha_z) * (area3(i,j,k+1) / vol(i,j,k)) * flux3(i,j,k+1,URHO)

                   if (rho < small_dens) then

                      fluxL = dflux(u(i,j,k  ,:), q(i,j,k  ,:), dir, [i, j, k  ], include_pressure)
                      fluxR = dflux(u(i,j,k+1,:), q(i,j,k+1,:), dir, [i, j, k+1], include_pressure)
                      fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha_z) * (u(i,j,k,:) - u(i,j,k+1,:)))

                      drhoLF = -TWO * (dt / alpha_z) * (area3(i,j,k+1) / vol(i,j,k)) * fluxLF(URHO)

                      if (u(i,j,k,URHO) + drhoLF < small_dens) then
                         fluxLF = fluxLF * abs((small_dens - u(i,j,k,URHO)) / drhoLF)
                      endif

                      rhoLF = u(i,j,k,URHO) - TWO * (dt / alpha_z) * (area3(i,j,k+1) / vol(i,j,k)) * fluxLF(URHO)

                      thetam(i,j,k) = (small_dens - rhoLF) / (rho - rhoLF)

                   endif

                endif

             enddo
          enddo
       enddo

       do k = lo(3), hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                theta = min(thetam(i,j,k-1), thetap(i,j,k))

                fluxL = dflux(u(i,j,k-1,:), q(i,j,k-1,:), dir, [i, j, k-1], include_pressure)
                fluxR = dflux(u(i,j,k  ,:), q(i,j,k  ,:), dir, [i, j, k  ], include_pressure)
                fluxLF(:) = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha_z) * (u(i,j,k-1,:) - u(i,j,k,:)))

                drhoLF = TWO * (dt / alpha_z) * (area3(i,j,k) / vol(i,j,k)) * fluxLF(URHO)

                if (u(i,j,k,URHO) + drhoLF < small_dens) then
                   fluxLF(:) = fluxLF(:) * abs((small_dens - u(i,j,k,URHO)) / drhoLF)
                else if (u(i,j,k-1,URHO) - drhoLF < small_dens) then
                   fluxLF(:) = fluxLF(:) * abs((small_dens - u(i,j,k-1,URHO)) / drhoLF)
                endif

                flux3(i,j,k,:) = (ONE - theta) * fluxLF(:) + theta * flux3(i,j,k,:)

                drho = TWO * (dt / alpha_z) * (area3(i,j,k) / vol(i,j,k)) * flux3(i,j,k,URHO)

                if (u(i,j,k,URHO) + drho < small_dens) then
                   flux3(i,j,k,:) = flux3(i,j,k,:) * abs((small_dens - u(i,j,k,URHO)) / drho)
                else if (u(i,j,k-1,URHO) - drho < small_dens) then
                   flux3(i,j,k,:) = flux3(i,j,k,:) * abs((small_dens - u(i,j,k-1,URHO)) / drho)
                endif

             enddo
          enddo
       enddo

    endif

    call bl_deallocate(thetap)
    call bl_deallocate(thetam)

  end subroutine limit_hydro_fluxes_on_small_dens

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine divu(lo,hi,q,q_lo,q_hi,dx,div,div_lo,div_hi)

    use meth_params_module, only : QU, QV, QW, QVAR
    use bl_constants_module

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: div_lo(3), div_hi(3)
    double precision, intent(in) :: dx(3)
    double precision, intent(inout) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    double precision, intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    integer          :: i, j, k
    double precision :: ux, vy, wz, dxinv, dyinv, dzinv

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1

             ux = FOURTH*( &
                    + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                    + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                    + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                    + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) ) * dxinv

             vy = FOURTH*( &
                    + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                    + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                    + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                    + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) ) * dyinv

             wz = FOURTH*( &
                    + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                    + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                    + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                    + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) ) * dzinv

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo

  end subroutine divu

end module advection_util_3d_module
