module hse_bc_module

  use bl_types
  use bl_error_module
  use prob_params_module
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use bl_constants_module, only : ZERO, HALF

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, parameter :: MAX_ITER = 100
  real(rt)        , parameter :: TOL = 1.e-8_rt
  
contains

  ! this expects a 3-d problem
  subroutine hse_bc_zlo(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc,density_only)

    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   UFS, UTEMP, const_grav
    use network, only: nspec
    use amrex_fort_module, only : rt => amrex_real
    implicit none
    integer :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    real(rt)         :: delta(3), xlo(3), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)
    logical, intent(in), optional :: density_only

    integer :: i, j, k, iter
    real(rt)         :: z, z0, slope
    real(rt)         :: pres_above, p_want, pres_zone, A
    real(rt)         :: drho,dpdr, temp_zone, eint, X_zone(nspec), dens_zone

    logical :: just_density

    logical :: converged_hse

    type (eos_t) :: eos_state

    if (present(density_only)) then
       just_density = density_only
    else
       just_density = .false.
    endif

    if (just_density) then
       ! we cannot rely on any of the other state variables being
       ! initialized, so the only thing we can do here is extrapolate
       ! the density from the domain interior
       do j=adv_l2,adv_h2
          do i=adv_l1,adv_h1
             slope = (adv(i,j,domlo(3)+1,URHO) - adv(i,j,domlo(3),URHO))/delta(3)
             z0 = problo(3) + dble(domlo(3) + HALF)*delta(3)

             do k=domlo(3)-1,adv_l3,-1
                z = problo(3) + dble(k + HALF)*delta(3)
                adv(i,j,k,URHO) = adv(i,j,domlo(3),URHO) + slope*(z - z0)
             enddo
          enddo
       enddo
       return
    endif

    ! HSE integration to get density, pressure
    do j = adv_l2, adv_h2
       do i = adv_l1, adv_h1

          ! temperature and species held constant in BCs
          temp_zone = adv(i,j,domlo(3),UTEMP)
          X_zone(:) = adv(i,j,domlo(3),UFS:UFS-1+nspec)/adv(i,j,domlo(3),URHO)

          ! get pressure in zone above
          eos_state%rho = adv(i,j,domlo(3),URHO)
          eos_state%T = temp_zone
          eos_state%xn(:) = X_zone(:)

          call eos(eos_input_rt, eos_state)

          eint = eos_state%e
          pres_above = eos_state%p

          ! this do loop counts backwards since we want to work downward
          do k = domlo(3)-1, adv_l3, -1

             ! initial guesses
             dens_zone = adv(i,j,k+1,URHO)

             converged_hse = .FALSE.

             do iter = 1, MAX_ITER

                ! pressure needed from HSE
                p_want = pres_above - &
                     delta(3)*HALF*(dens_zone + adv(i,j,k+1,URHO))*const_grav

                ! pressure from EOS
                eos_state%rho = dens_zone
                eos_state%T = temp_zone
                eos_state%xn(:) = X_zone(:)

                call eos(eos_input_rt, eos_state)

                pres_zone = eos_state%p
                dpdr = eos_state%dpdr
                eint = eos_state%e

                ! Newton-Raphson - we want to zero A = p_want - p(rho)
                A = p_want - pres_zone
                drho = A/(dpdr + HALF*delta(3)*const_grav)

                dens_zone = max(0.9_rt*dens_zone, &
                                min(dens_zone + drho, 1.1_rt*dens_zone))

                ! convergence?
                if (abs(drho) < TOL*dens_zone) then
                   converged_hse = .TRUE.
                   exit
                endif

             enddo

             if (.not. converged_hse) call bl_error("ERROR: failure to converge in -Z BC")

             ! velocity
             if (zero_vels) then

                ! zero normal momentum causes pi waves to pass through
                adv(i,j,k,UMZ) = ZERO

                ! zero transverse momentum
                adv(i,j,k,UMX) = ZERO
                adv(i,j,k,UMY) = ZERO
             else

                ! zero gradient velocity
                adv(i,j,k,UMX) = dens_zone*(adv(i,j,domlo(3),UMX)/adv(i,j,domlo(3),URHO))
                adv(i,j,k,UMY) = dens_zone*(adv(i,j,domlo(3),UMY)/adv(i,j,domlo(3),URHO))
                adv(i,j,k,UMZ) = dens_zone*(adv(i,j,domlo(3),UMZ)/adv(i,j,domlo(3),URHO))
             endif
          
             adv(i,j,k,URHO) = dens_zone

             eos_state%rho = dens_zone
             eos_state%T = temp_zone
             eos_state%xn(:) = X_zone

             call eos(eos_input_rt, eos_state)

             pres_zone = eos_state%p
             eint = eos_state%e

             adv(i,j,k,UEINT) = dens_zone*eint
             adv(i,j,k,UEDEN) = dens_zone*eint + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/dens_zone
             adv(i,j,k,UTEMP) = temp_zone
             adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

             ! for the next iteration
             pres_above = pres_zone

          end do
       end do
    enddo

  end subroutine hse_bc_zlo


  subroutine hse_bc_zhi(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc, density_only)

    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP
    use interpolate_module
    use network, only: nspec

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer bc(3,2,*)
    integer domlo(3), domhi(3)
    real(rt)         delta(3), xlo(3), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)
    logical, intent(in), optional :: density_only

    integer i,j,k

    logical :: just_density

    if (present(density_only)) then
       just_density = density_only
    else
       just_density = .false.
    endif

    if (just_density) then
       ! we cannot rely on any of the other state variables being
       ! initialized, so the only thing we can do here is extrapolate
       ! the density from the domain interior
       do k = domhi(3)+1, adv_h3
          do j = adv_l2, adv_h2
             do i = adv_l1, adv_h1
                adv(i,j,k,URHO) = adv(i,j,domhi(3),URHO)
             enddo
          enddo
       enddo
       return
    endif

    ! outflow, zeroing the normal velocity

    do j = adv_l2, adv_h2
       do i = adv_l1, adv_h1
          do k = domhi(3)+1, adv_h3
             adv(i,j,k,:) = adv(i,j,domhi(3),:)

             adv(i,j,k,UMZ) = max(ZERO, adv(i,j,k,UMZ))

             adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
          enddo
       enddo
    enddo
  end subroutine hse_bc_zhi

end module hse_bc_module

