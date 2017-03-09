module hse_bc_module

  use bl_types
  use bl_error_module
  use prob_params_module

  use amrex_fort_module, only : rt => c_real
  implicit none

  integer, parameter :: MAX_ITER = 100
  real(rt)        , parameter :: TOL = 1.e-8_rt
  
contains

  ! this expects a 2-d problem
  subroutine hse_bc_ylo(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc,density_only)

    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   UFS, UTEMP, const_grav
    use eos_module
    use network, only: nspec

    use amrex_fort_module, only : rt => c_real
    implicit none
    integer :: adv_l1, adv_l2, adv_h1, adv_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)
    logical, intent(in), optional :: density_only

    integer :: i, j, iter
    real(rt)         :: y, y0, slope
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
       do i=adv_l1,adv_h1
          slope = (adv(i,domlo(2)+1,URHO) - adv(i,domlo(2),URHO))/delta(2)
          y0 = problo(2) + dble(domlo(2) + HALF)*delta(2)

          do j=domlo(2)-1,adv_l2,-1
             y = problo(2) + dble(j + HALF)*delta(2)
             adv(i,j,URHO) = adv(i,domlo(2),URHO) + slope*(y - y0)
          enddo
       enddo

       return
    endif

    ! HSE integration to get density, pressure

    do i = adv_l1, adv_h1

       ! temperature and species held constant in BCs
       temp_zone = adv(i,domlo(2),UTEMP)
       X_zone(:) = adv(i,domlo(2),UFS:UFS-1+nspec)/adv(i,domlo(2),URHO)

       ! get pressure in zone above
       eos_state%rho = adv(i,domlo(2),URHO)
       eos_state%T = temp_zone
       eos_state%xn(:) = X_zone(:)

       call eos(eos_input_rt, eos_state)

       eint = eos_state%e
       pres_above = eos_state%p

       ! this do loop counts backwards since we want to work downward
       do j = domlo(2)-1, adv_l2, -1

          ! initial guesses
          dens_zone = adv(i,j+1,URHO)

          converged_hse = .FALSE.

          do iter = 1, MAX_ITER

             ! pressure needed from HSE
             p_want = pres_above - &
                  delta(2)*HALF*(dens_zone + adv(i,j+1,URHO))*const_grav

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
             drho = A/(dpdr + HALF*delta(2)*const_grav)

             dens_zone = max(0.9_rt*dens_zone, &
                             min(dens_zone + drho, 1.1_rt*dens_zone))

             ! convergence?
             if (abs(drho) < TOL*dens_zone) then
                converged_hse = .TRUE.
                exit
             endif

          enddo

          if (.not. converged_hse) call bl_error("ERROR: failure to converge in -Y BC")

          ! velocity
          if (zero_vels) then

             ! zero normal momentum causes pi waves to pass through
             adv(i,j,UMY) = ZERO

             ! zero transverse momentum
             adv(i,j,UMX) = ZERO
             adv(i,j,UMZ) = ZERO
          else

             ! zero gradient velocity
             adv(i,j,UMX) = dens_zone*(adv(i,domlo(2),UMX)/adv(i,domlo(2),URHO))
             adv(i,j,UMY) = dens_zone*(adv(i,domlo(2),UMY)/adv(i,domlo(2),URHO))
             adv(i,j,UMZ) = dens_zone*(adv(i,domlo(2),UMZ)/adv(i,domlo(2),URHO))
          endif
          
          adv(i,j,URHO) = dens_zone

          eos_state%rho = dens_zone
          eos_state%T = temp_zone
          eos_state%xn(:) = X_zone

          call eos(eos_input_rt, eos_state)

          pres_zone = eos_state%p
          eint = eos_state%e

          adv(i,j,UEINT) = dens_zone*eint
          adv(i,j,UEDEN) = dens_zone*eint + HALF*sum(adv(i,j,UMX:UMZ)**2)/dens_zone
          adv(i,j,UTEMP) = temp_zone
          adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

          ! for the next iteration
          pres_above = pres_zone

       end do
    end do

  end subroutine hse_bc_ylo


  subroutine hse_bc_yhi(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc, density_only)

    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
    use interpolate_module
    use eos_module
    use network, only: nspec

    use amrex_fort_module, only : rt => c_real
    implicit none
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)
    logical, intent(in), optional :: density_only

    integer i,j

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
       do j = domhi(2)+1, adv_h2
          do i = adv_l1, adv_h1
             adv(i,j,URHO) = rho_ambient
          enddo
       enddo

       return
    endif

    ! outflow, zeroing the normal velocity

    do i = adv_l1, adv_h1
       do j = domhi(2)+1, adv_h2
          adv(i,j,:) = adv(i,domhi(2),:)

          adv(i,j,URHO) = rho_ambient
          adv(i,j,UTEMP) = T_ambient
          adv(i,j,UFS:UFS-1+nspec) = rho_ambient*xn_ambient(:)
          adv(i,j,UEINT) = rho_ambient*e_ambient

          !adv(i,j,UMY) = max(ZERO, adv(i,j,UMY))

          adv(i,j,UEDEN) = adv(i,j,UEINT) + HALF*sum(adv(i,j,UMX:UMZ)**2)/adv(i,j,URHO)
       enddo
    enddo

  end subroutine hse_bc_yhi

end module hse_bc_module

