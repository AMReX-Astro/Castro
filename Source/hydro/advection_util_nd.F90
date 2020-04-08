module advection_util_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_enforce_minimum_density(lo, hi, &
                                        state, s_lo, s_hi, &
                                        frac_change, verbose) bind(c,name='ca_enforce_minimum_density')

    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, small_dens, density_reset_method
    use amrex_constants_module, only: ZERO
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif
    use amrex_fort_module, only: rt => amrex_real
    use reduction_module, only: reduce_min

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: frac_change
    integer,  intent(in   ), value :: verbose

    ! Local variables
    integer  :: i, j, k
    integer  :: ii, jj, kk
    integer  :: i_set, j_set, k_set
    real(rt) :: max_dens, old_rho
    real(rt) :: uold(NVAR), unew(NVAR)
    integer  :: num_positive_zones
    real(rt) :: frac_change_tmp

    !$gpu

    max_dens = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             frac_change_tmp = 1.0_rt

             if (state(i,j,k,URHO) .eq. ZERO) then

#ifndef AMREX_USE_GPU
                print *,'DENSITY EXACTLY ZERO AT CELL ', i, j, k
                print *,'  in grid ',lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
                call castro_error("Error :: ca_enforce_minimum_density")
#endif

             else if (state(i,j,k,URHO) < small_dens) then

                old_rho = state(i,j,k,URHO)

                if (density_reset_method == 1) then

                   ! Reset to the characteristics of the adjacent state with the highest density.

                   max_dens = state(i,j,k,URHO)
                   i_set = i
                   j_set = j
                   k_set = k
                   do kk = -1, 1
                      do jj = -1, 1
                         do ii = -1, 1

                            if (i+ii >= s_lo(1) .and. j+jj >= s_lo(2) .and. k+kk >= s_lo(3) .and. &
                                i+ii <= s_hi(1) .and. j+jj <= s_hi(2) .and. k+kk <= s_hi(3)) then

                               if (state(i+ii,j+jj,k+kk,URHO) .gt. max_dens) then

                                  i_set = i+ii
                                  j_set = j+jj
                                  k_set = k+kk
                                  max_dens = state(i_set,j_set,k_set,URHO)

                               end if

                            end if

                         end do
                      end do
                   end do

                   if (max_dens < small_dens) then

                      ! We could not find any nearby zones with sufficient density.

                      uold = state(i,j,k,:)
                      call reset_to_small_state(uold, [i, j, k], s_lo, s_hi, verbose)
                      state(i,j,k,:) = uold

                   else

                      uold = state(i,j,k,:)
                      unew = state(i_set,j_set,k_set,:)

                      call reset_to_zone_state(uold, unew, [i, j, k], s_lo, s_hi, verbose)

                      state(i,j,k,:) = uold

                   endif

                else if (density_reset_method == 2) then

                   ! Reset to the average of adjacent zones. The median is independently calculated for each variable.

                   num_positive_zones = 0
                   unew(:) = ZERO

                   do kk = -1, 1
                      do jj = -1, 1
                         do ii = -1, 1

                            if (i+ii >= s_lo(1) .and. j+jj >= s_lo(2) .and. k+kk >= s_lo(3) .and. &
                                i+ii <= s_hi(1) .and. j+jj <= s_hi(2) .and. k+kk <= s_hi(3)) then

                               if (state(i+ii,j+jj,k+kk,URHO) .ge. small_dens) then

                                  unew(:) = unew(:) + state(i+ii,j+jj,k+kk,:)
                                  num_positive_zones = num_positive_zones + 1

                               end if

                            end if

                         end do
                      end do
                   end do

                   if (num_positive_zones == 0) then

                      ! We could not find any nearby zones with sufficient density.

                      uold = state(i,j,k,:)
                      call reset_to_small_state(uold, [i, j, k], s_lo, s_hi, verbose)
                      state(i,j,k,:) = uold

                   else

                      uold = state(i,j,k,:)
                      unew(:) = unew(:) / num_positive_zones

                      call reset_to_zone_state(uold, unew, [i, j, k], s_lo, s_hi, verbose)

                      state(i,j,k,:) = uold

                   endif

#ifndef AMREX_USE_CUDA
                else

                   call castro_error("Unknown density_reset_method in subroutine ca_enforce_minimum_density.")
#endif
                endif

                ! Store the maximum (negative) fractional change in the density from this reset.

                if (old_rho < ZERO) then
                   frac_change_tmp = 1.0_rt
                else
                   frac_change_tmp = (state(i,j,k,URHO) - old_rho) / old_rho
                end if

             end if

             call reduce_min(frac_change, frac_change_tmp)

          end do
       end do
    end do

  end subroutine ca_enforce_minimum_density


  subroutine reset_to_small_state(state, idx, lo, hi, verbose)
    ! If no neighboring zones are above small_dens, our only recourse
    ! is to set the density equal to small_dens, and the temperature
    ! equal to small_temp. We set the velocities to zero,
    ! though any choice here would be arbitrary.
    !

    use amrex_constants_module, only: ZERO
    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEINT, UEDEN, UFS, small_temp, small_dens, npassive, upass_map
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_module, only: eos
    use castro_util_module, only: position ! function
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only: linear_to_hybrid ! function
    use meth_params_module, only: UMR, UMP
#endif

    use amrex_fort_module, only: rt => amrex_real
    implicit none

    real(rt)         :: state(NVAR)
    integer          :: idx(3), lo(3), hi(3), verbose

    integer          :: n, ipassive
    type (eos_t)     :: eos_state

#ifdef HYBRID_MOMENTUM
    real(rt)         :: loc(3)
#endif

    !$gpu

#ifndef AMREX_USE_CUDA
    if (verbose .gt. 0) then
       print *,'   '
       if (state(URHO) < ZERO) then
          print *,'>>> RESETTING NEG.  DENSITY AT ', idx(1), idx(2), idx(3)
       else
          print *,'>>> RESETTING SMALL DENSITY AT ', idx(1), idx(2), idx(3)
       endif
       print *,'>>> FROM ', state(URHO), ' TO ', small_dens
       print *,'>>> IN GRID ', lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
       print *,'   '
    end if
#endif

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       state(n) = state(n) * (small_dens / state(URHO))
    end do

    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = state(UFS:UFS+nspec-1) / small_dens
    eos_state % aux = state(UFS:UFS+naux-1) / small_dens

    call eos(eos_input_rt, eos_state)

    state(URHO ) = eos_state % rho
    state(UTEMP) = eos_state % T

    state(UMX  ) = ZERO
    state(UMY  ) = ZERO
    state(UMZ  ) = ZERO

    state(UEINT) = eos_state % rho * eos_state % e
    state(UEDEN) = state(UEINT)

#ifdef HYBRID_MOMENTUM
    loc = position(idx(1),idx(2),idx(3))
    state(UMR:UMP) = linear_to_hybrid(loc, state(UMX:UMZ))
#endif

  end subroutine reset_to_small_state



  subroutine reset_to_zone_state(state, input_state, idx, lo, hi, verbose)

    use amrex_constants_module, only: ZERO
    use meth_params_module, only: NVAR, URHO
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt) :: state(NVAR), input_state(NVAR)
    integer  :: idx(3), lo(3), hi(3), verbose

    !$gpu

#ifndef AMREX_USE_CUDA
    if (verbose .gt. 0) then
       if (state(URHO) < ZERO) then
          print *,'   '
          print *,'>>> RESETTING NEG.  DENSITY AT ',idx(1), idx(2), idx(3)
          print *,'>>> FROM ', state(URHO) ,' TO ', input_state(URHO)
          print *,'>>> IN GRID ', lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
          print *,'   '
       else
          print *,'   '
          print *,'>>> RESETTING SMALL DENSITY AT ', idx(1), idx(2), idx(3)
          print *,'>>> FROM ', state(URHO), ' TO ', input_state(URHO)
          print *,'>>> IN GRID ', lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
          print *,'   '
       end if
    end if
#endif

    state(:) = input_state(:)

  end subroutine reset_to_zone_state


  function dflux(u, q, dir, idx) result(flux)
    ! Given a conservative state and its corresponding primitive state, calculate the
    ! corresponding flux in a given direction.
    !

    use amrex_constants_module, only: ZERO
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, &
         NQ, QU, QPRES, &
         npassive, upass_map
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only: compute_hybrid_flux
    use meth_params_module, only: NGDNV, GDRHO, GDU, GDW, GDPRES, QRHO, QW
#endif
    use prob_params_module, only: mom_flux_has_p
    use amrex_fort_module, only: rt => amrex_real
    implicit none

    integer :: dir, idx(3)
    real(rt)         :: u(NVAR), q(NQ), flux(NVAR)

    real(rt)         :: v_adv
    integer :: ipassive, n
#ifdef HYBRID_MOMENTUM
    real(rt)         :: qgdnv(NGDNV)
    logical :: cell_centered
#endif

    !$gpu

    ! Set everything to zero; this default matters because some
    ! quantities like temperature are not updated through fluxes.

    flux = ZERO

    ! Determine the advection speed based on the flux direction.

    v_adv = q(QU + dir - 1)

    ! Core quantities (density, momentum, energy).

    flux(URHO) = u(URHO) * v_adv
    flux(UMX:UMZ) = u(UMX:UMZ) * v_adv
    flux(UEDEN) = (u(UEDEN) + q(QPRES)) * v_adv
    flux(UEINT) = u(UEINT) * v_adv

    ! Optionally include the pressure term in the momentum flux.
    ! It is optional because for some geometries we cannot write
    ! the pressure term in a conservative form.

    if (mom_flux_has_p(dir)%comp(UMX+dir-1)) then
       flux(UMX + dir - 1) = flux(UMX + dir - 1) + q(QPRES)
    endif

    ! Hybrid flux.

#ifdef HYBRID_MOMENTUM
    ! Create a temporary edge-based q for this routine.
    qgdnv(:) = ZERO
    qgdnv(GDRHO) = q(QRHO)
    qgdnv(GDU:GDW) = q(QU:QW)
    qgdnv(GDPRES) = q(QPRES)
    cell_centered = .true.
    call compute_hybrid_flux(qgdnv, flux, dir, idx, cell_centered)
#endif

    ! Passively advected quantities.

    do ipassive = 1, npassive

       n = upass_map(ipassive)
       flux(n) = u(n) * v_adv

    enddo

  end function dflux


  subroutine limit_hydro_fluxes_on_small_dens(lo, hi, &
                                              idir, &
                                              u, u_lo, u_hi, &
                                              q, q_lo, q_hi, &
                                              vol, vol_lo, vol_hi, &
                                              flux, flux_lo, flux_hi, &
                                              area, area_lo, area_hi, &
                                              dt, dx) bind(c, name="limit_hydro_fluxes_on_small_dens")
    ! The following algorithm comes from Hu, Adams, and Shu (2013), JCP, 242, 169,
    ! "Positivity-preserving method for high-order conservative schemes solving
    ! compressible Euler equations." It has been modified to enforce not only positivity
    ! but also the stronger requirement that rho > small_dens. We do not limit on pressure
    ! (or, similarly, internal energy) because those cases are easily fixed by calls to
    ! reset_internal_energy that enforce a thermodynamic floor. The density limiter, by
    ! contrast, is very important because calls to enforce_minimum_density can yield
    ! hydrodynamic states that are inconsistent (there is no clear strategy for what to do
    ! when a density is negative).

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE, TWO
    use meth_params_module, only: NVAR, NQ, URHO, UTEMP, small_dens, cfl
#ifdef SHOCK_VAR
    use meth_params_module, only: USHK
#endif
    use prob_params_module, only: dim
    use amrex_mempool_module, only: bl_allocate, bl_deallocate

    implicit none

    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in), value :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: flux_lo(3), flux_hi(3)
    integer, intent(in) :: area_lo(3), area_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: flux(flux_lo(1):flux_hi(1),flux_lo(2):flux_hi(2),flux_lo(3):flux_hi(3),NVAR)
    real(rt), intent(in   ) :: area(area_lo(1):area_hi(1),area_lo(2):area_hi(2),area_lo(3):area_hi(3))

    integer  :: i, j, k

    real(rt) :: rhoL, rhoR, drhoL, drhoR, fluxLF(NVAR), fluxL(NVAR), fluxR(NVAR), rhoLF, drhoLF, dtdx, theta, alpha
    real(rt) :: uL(NVAR), uR(NVAR), qL(NQ), qR(NQ), volL, volR, flux_coefL, flux_coefR
    integer  :: idxL(3), idxR(3)

    real(rt), parameter :: density_floor_tolerance = 1.1_rt
    real(rt) :: density_floor

    !$gpu

    ! The density floor is the small density, modified by a small factor.
    ! In practice numerical error can cause the density that is created
    ! by this flux limiter to be slightly lower than the target density,
    ! so we set the target to be slightly larger than the real density floor
    ! to avoid density resets.

    density_floor = small_dens * density_floor_tolerance

    dtdx = dt / dx(idir)
    alpha = ONE / DIM

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! Grab the states on either side of the interface we are working with,
             ! depending on which dimension we're currently calling this with.

             uR = u(i,j,k,:)
             qR = q(i,j,k,:)
             volR = vol(i,j,k)
             idxR = [i,j,k]

             if (idir == 1) then
                uL = u(i-1,j,k,:)
                qL = q(i-1,j,k,:)
                volL = vol(i-1,j,k)
                idxL = [i-1,j,k]
             else if (idir == 2) then
                uL = u(i,j-1,k,:)
                qL = q(i,j-1,k,:)
                volL = vol(i,j-1,k)
                idxL = [i,j-1,k]
             else
                uL = u(i,j,k-1,:)
                qL = q(i,j,k-1,:)
                volL = vol(i,j,k-1)
                idxL = [i,j,k-1]
             end if

             ! If an adjacent zone has a floor-violating density, set the flux to zero and move on.
             ! At that point, the only thing to do is wait for a reset at a later point.

             if (uR(URHO) < density_floor .or. uL(URHO) < density_floor) then

                flux(i,j,k,:) = ZERO
                cycle

             endif

             ! Construct cell-centered fluxes.

             fluxL = dflux(uL, qL, idir, idxL)
             fluxR = dflux(uR, qR, idir, idxR)

             ! Construct the Lax-Friedrichs flux on the interface (Equation 12).
             ! Note that we are using the information from Equation 9 to obtain the
             ! effective maximum wave speed, (|u| + c)_max = CFL / lambda where
             ! lambda = dt/(dx * alpha); alpha = 1 in 1D and may be chosen somewhat
             ! freely in multi-D as long as alpha_x + alpha_y + alpha_z = 1.

             fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha) * (uL(:) - uR(:)))

             ! Coefficients of fluxes on either side of the interface.

             flux_coefR = TWO * (dt / alpha) * area(i,j,k) / volR
             flux_coefL = TWO * (dt / alpha) * area(i,j,k) / volL

             ! Obtain the one-sided update to the density, based on Hu et al., Eq. 11.
             ! If we would violate the floor, then we need to limit the flux. Since the
             ! flux adds to the density on one side and subtracts from the other, the floor
             ! can only be violated in at most one direction, so we'll do an if-else test
             ! below. This means that we can simplify the approach of Hu et al. -- whereas
             ! they constructed two thetas for each interface (corresponding to either side)
             ! we can complete the operation in one step with a single theta.

             drhoL = flux_coefL * flux(i,j,k,URHO)
             rhoL = uL(URHO) - drhoL

             drhoR = flux_coefR * flux(i,j,k,URHO)
             rhoR = uR(URHO) + drhoR

             theta = ONE

             if (rhoL < density_floor) then

                ! Obtain the final density corresponding to the LF flux.

                drhoLF = flux_coefL * fluxLF(URHO)
                rhoLF = uL(URHO) - drhoLF

                ! Solve for theta from (1 - theta) * rhoLF + theta * rho = density_floor.

                theta = (density_floor - rhoLF) / (rhoL - rhoLF)

                ! Limit theta to the valid range (this will deal with roundoff issues).

                theta = min(ONE, max(theta, ZERO))

             else if (rhoR < density_floor) then

                drhoLF = flux_coefR * fluxLF(URHO)
                rhoLF = uR(URHO) + drhoLF

                theta = (density_floor - rhoLF) / (rhoR - rhoLF)

                theta = min(ONE, max(theta, ZERO))

             endif

             ! Assemble the limited flux (Equation 16).

             flux(i,j,k,:) = (ONE - theta) * fluxLF(:) + theta * flux(i,j,k,:)

             ! Zero out fluxes for quantities that don't advect.

             flux(i,j,k,UTEMP) = ZERO
#ifdef SHOCK_VAR
             flux(i,j,k,USHK) = ZERO
#endif

             ! Now, apply our requirement that the final flux cannot violate the density floor.

             drhoR = flux_coefR * flux(i,j,k,URHO)
             drhoL = flux_coefL * flux(i,j,k,URHO)

             if (uR(URHO) + drhoR < density_floor) then
                flux(i,j,k,:) = flux(i,j,k,:) * abs((density_floor - uR(URHO)) / drhoR)
             else if (uL(URHO) - drhoL < density_floor) then
                flux(i,j,k,:) = flux(i,j,k,:) * abs((density_floor - uL(URHO)) / drhoL)
             endif

          enddo
       enddo
    enddo

  end subroutine limit_hydro_fluxes_on_small_dens



  subroutine limit_hydro_fluxes_on_large_vel(lo, hi, &
                                             idir, &
                                             u, u_lo, u_hi, &
                                             q, q_lo, q_hi, &
                                             vol, vol_lo, vol_hi, &
                                             flux, flux_lo, flux_hi, &
                                             area, area_lo, area_hi, &
                                             dt, dx) bind(c, name="limit_hydro_fluxes_on_large_vel")
    ! This limiter is similar to the density-based limiter above, but limits
    ! on velocities that are too large instead. The comments are minimal since
    ! the algorithm is effectively the same.

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE, TWO
    use meth_params_module, only: NVAR, NQ, URHO, UMX, UTEMP, cfl, speed_limit
#ifdef SHOCK_VAR
    use meth_params_module, only: USHK
#endif
    use prob_params_module, only: dim

    implicit none

    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in), value :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: flux_lo(3), flux_hi(3)
    integer, intent(in) :: area_lo(3), area_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: flux(flux_lo(1):flux_hi(1),flux_lo(2):flux_hi(2),flux_lo(3):flux_hi(3),NVAR)
    real(rt), intent(in   ) :: area(area_lo(1):area_hi(1),area_lo(2):area_hi(2),area_lo(3):area_hi(3))

    integer  :: i, j, k, n

    real(rt) :: fluxLF(NVAR), fluxL(NVAR), fluxR(NVAR), dtdx, theta, alpha
    real(rt) :: rhouL, rhouR, drhouL, drhouR, rhouLF, drhouLF
    real(rt) :: rhoL, rhoR, drhoL, drhoR
    real(rt) :: uL(NVAR), uR(NVAR), qL(NQ), qR(NQ), volL, volR, flux_coefL, flux_coefR
    integer  :: idxL(3), idxR(3), UMOM

    !$gpu

    dtdx = dt / dx(idir)
    alpha = ONE / DIM

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             uR = u(i,j,k,:)
             qR = q(i,j,k,:)
             volR = vol(i,j,k)
             idxR = [i,j,k]

             if (idir == 1) then
                uL = u(i-1,j,k,:)
                qL = q(i-1,j,k,:)
                volL = vol(i-1,j,k)
                idxL = [i-1,j,k]
             else if (idir == 2) then
                uL = u(i,j-1,k,:)
                qL = q(i,j-1,k,:)
                volL = vol(i,j-1,k)
                idxL = [i,j-1,k]
             else
                uL = u(i,j,k-1,:)
                qL = q(i,j,k-1,:)
                volL = vol(i,j,k-1)
                idxL = [i,j,k-1]
             end if

             ! Construct cell-centered fluxes.

             fluxL = dflux(uL, qL, idir, idxL)
             fluxR = dflux(uR, qR, idir, idxR)

             ! Construct the Lax-Friedrichs flux on the interface.

             fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha) * (uL(:) - uR(:)))

             ! Coefficients of fluxes on either side of the interface.

             flux_coefR = TWO * (dt / alpha) * area(i,j,k) / volR
             flux_coefL = TWO * (dt / alpha) * area(i,j,k) / volL

             theta = ONE

             ! Loop over all three momenta, and choose the strictest
             ! limiter among them.

             do n = 1, 3

                UMOM = UMX + n - 1

                ! Obtain the one-sided update to the momentum.

                drhouL = flux_coefL * flux(i,j,k,UMOM)
                rhouL = abs(uL(UMOM) - drhouL)

                drhoL = flux_coefL * flux(i,j,k,URHO)
                rhoL = uL(URHO) - drhoL

                drhouR = flux_coefR * flux(i,j,k,UMOM)
                rhouR = abs(uR(UMOM) + drhouR)

                drhoR = flux_coefR * flux(i,j,k,URHO)
                rhoR = uR(uRHO) + drhoR

                if (abs(rhouL) > rhoL * speed_limit) then

                   ! Obtain the final density corresponding to the LF flux.

                   drhouLF = flux_coefL * fluxLF(UMOM)
                   rhouLF = abs(uL(UMOM) - drhouLF)

                   ! Solve for theta from (1 - theta) * rhouLF + theta * rhou = rhoL * speed_limit.

                   theta = abs(rhoL * speed_limit - rhouLF) / abs(rhouL - rhouLF)

                   ! Limit theta to the valid range (this will deal with roundoff issues).

                   theta = min(ONE, max(theta, ZERO))

                else if (abs(rhouR) > rhoR * speed_limit) then

                   drhouLF = flux_coefR * fluxLF(UMOM)
                   rhouLF = abs(uR(UMOM) + drhouLF)

                   theta = abs(rhoR * speed_limit - rhouLF) / abs(rhouR - rhouLF)

                   theta = min(ONE, max(theta, ZERO))

                endif

             end do

             ! Assemble the limited flux (Equation 16).

             flux(i,j,k,:) = (ONE - theta) * fluxLF(:) + theta * flux(i,j,k,:)

             ! Zero out fluxes for quantities that don't advect.

             flux(i,j,k,UTEMP) = ZERO
#ifdef SHOCK_VAR
             flux(i,j,k,USHK) = ZERO
#endif

          enddo
       enddo
    enddo

  end subroutine limit_hydro_fluxes_on_large_vel

end module advection_util_module
