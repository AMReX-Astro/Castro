module bc_ext_fill_module

  use bl_constants_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UTEMP, const_grav, &
                                 hse_zero_vels, hse_interp_temp, hse_reflect_vels, &
                                 xl_ext, xr_ext, EXT_HSE, EXT_INTERP
  use interpolate_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  include 'AMReX_bc_types.fi'

  private

  public :: ext_fill, ext_denfill

contains

  ! this module contains different routines for filling the
  ! hydrodynamics boundary conditions

  ! NOTE: the hydrostatic boundary conditions here rely on
  ! constant gravity

  subroutine ext_fill(adv, adv_l1, adv_h1, &
                      domlo, domhi, delta, xlo, time, bc) &
                      bind(C, name="ext_fill")

    use prob_params_module, only : problo
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec
    use model_parser_module

    use amrex_fort_module, only : rt => amrex_real
    integer adv_l1, adv_h1
    integer bc(2,*)
    integer domlo(1), domhi(1)
    real(rt)         delta(1), xlo(1), time
    real(rt)         adv(adv_l1:adv_h1,NVAR)

    integer i, j, q, n, iter, m
    real(rt)         x
    real(rt)         :: dens_above, dens_base, temp_above
    real(rt)         :: pres_above, p_want, pres_zone, A
    real(rt)         :: drho, dpdr, temp_zone, eint, X_zone(nspec), dens_zone

    integer, parameter :: MAX_ITER = 100
    real(rt)        , parameter :: TOL = 1.e-8_rt
    logical :: converged_hse

    type (eos_t) :: eos_state

    do n = 1, NVAR


       ! XLO
       if (bc(1,n) == EXT_DIR .and. adv_l1 < domlo(1)) then

          if (xl_ext == EXT_HSE) then

             ! we will fill all the variables when we consider URHO
             if (n == URHO) then


                   ! we are integrating along a column at constant i.
                   ! Make sure that our starting state is well-defined
                   dens_above = adv(domlo(1),URHO)

                   ! sometimes, we might be working in a corner
                   ! where the ghost cells above us have not yet
                   ! been initialized.  In that case, take the info
                   ! from the initial model
                   if (dens_above == ZERO) then
                      x = problo(1) + delta(1)*(dble(domlo(1)) + HALF)

                      dens_above = interpolate(x,npts_model,model_r, &
                                               model_state(:,idens_model))

                      temp_above = interpolate(x,npts_model,model_r, &
                                              model_state(:,itemp_model))

                      do m = 1, nspec
                         X_zone(m) = interpolate(x,npts_model,model_r, &
                                                 model_state(:,ispec_model-1+m))
                      enddo

                   else
                      temp_above = adv(domlo(1),UTEMP)
                      X_zone(:) = adv(domlo(1),UFS:UFS-1+nspec)/dens_above
                   endif

                   ! keep track of the density at the base of the domain
                   dens_base = dens_above

                   ! get pressure in this zone (the initial above zone)
                   eos_state%rho = dens_above
                   eos_state%T = temp_above
                   eos_state%xn(:) = X_zone(:)

                   call eos(eos_input_rt, eos_state)

                   eint = eos_state%e
                   pres_above = eos_state%p

                   ! integrate downward
                   do j = domlo(1)-1, adv_l1, -1
                      x = problo(1) + delta(1)*(dble(j) + HALF)

                      ! HSE integration to get density, pressure

                      ! initial guesses
                      dens_zone = dens_above

                      ! temperature and species held constant in BCs
                      if (hse_interp_temp == 1) then
                         temp_zone = interpolate(x,npts_model,model_r, &
                                                 model_state(:,itemp_model))
                      else
                         temp_zone = temp_above
                      endif

                      converged_hse = .FALSE.


                      do iter = 1, MAX_ITER

                         ! pressure needed from HSE
                         p_want = pres_above - &
                              delta(1)*HALF*(dens_zone + dens_above)*const_grav

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
                         drho = A/(dpdr + HALF*delta(1)*const_grav)

                         dens_zone = max(0.9_rt*dens_zone, &
                              min(dens_zone + drho, 1.1_rt*dens_zone))

                         ! convergence?
                         if (abs(drho) < TOL*dens_zone) then
                            converged_hse = .TRUE.
                            exit
                         endif

                      enddo

                      if (.not. converged_hse) then
                         print *, "j, domlo(2): ", j, domlo(1)
                         print *, "p_want:    ", p_want
                         print *, "dens_zone: ", dens_zone
                         print *, "temp_zone: ", temp_zone
                         print *, "drho:      ", drho
                         print *, " "
                         print *, "column info: "
                         print *, "   dens: ", adv(j:domlo(1),URHO)
                         print *, "   temp: ", adv(j:domlo(1),UTEMP)
                         call bl_error("ERROR in bc_ext_fill_1d: failure to converge in -X BC")
                      endif


                      ! velocity
                      if (hse_zero_vels == 1) then

                         ! zero normal momentum causes pi waves to pass through
                         adv(j,UMX) = ZERO
                      else

                         if (hse_reflect_vels == 1) then
                            adv(j,UMX) = -dens_zone*(adv(domlo(1),UMX)/dens_base)
                         else
                            adv(j,UMX) = dens_zone*(adv(domlo(1),UMX)/dens_base)
                         endif
                      endif

                      eos_state%rho = dens_zone
                      eos_state%T = temp_zone
                      eos_state%xn(:) = X_zone
                      
                      call eos(eos_input_rt, eos_state)

                      pres_zone = eos_state%p
                      eint = eos_state%e

                      ! store the final state
                      adv(j,URHO) = dens_zone
                      adv(j,UEINT) = dens_zone*eint
                      adv(j,UEDEN) = dens_zone*eint + &
                           HALF*sum(adv(j,UMX:UMZ)**2)/dens_zone
                      adv(j,UTEMP) = temp_zone
                      adv(j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                      ! for the next zone
                      dens_above = dens_zone
                      pres_above = pres_zone

                   enddo

             endif  ! n == URHO

          elseif (xl_ext == EXT_INTERP) then

             do j = domlo(1)-1, adv_l1, -1
                x = problo(1) + delta(1)*(dble(j) + HALF)

                   ! set all the variables even though we're testing on URHO
                   if (n == URHO) then

                      dens_zone = interpolate(x,npts_model,model_r, &
                                              model_state(:,idens_model))

                      temp_zone = interpolate(x,npts_model,model_r, &
                                              model_state(:,itemp_model))

                      do q = 1, nspec
                         X_zone(q) = interpolate(x,npts_model,model_r, &
                                                 model_state(:,ispec_model-1+q))
                      enddo

                      ! extrap normal momentum
                      adv(j,UMX) = ZERO

                      eos_state%rho = dens_zone
                      eos_state%T = temp_zone
                      eos_state%xn(:) = X_zone

                      call eos(eos_input_rt, eos_state)

                      pres_zone = eos_state%p
                      eint = eos_state%e

                      adv(j,URHO) = dens_zone
                      adv(j,UEINT) = dens_zone*eint
                      adv(j,UEDEN) = dens_zone*eint + &
                           HALF*sum(adv(j,UMX:UMZ)**2)/dens_zone
                      adv(j,UTEMP) = temp_zone
                      adv(j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
                   endif

                enddo
          endif  ! xl_ext check


       endif


       ! YHI
       if (bc(2,n) == EXT_DIR .and. adv_h1 > domhi(1)) then

          if (xr_ext == EXT_HSE) then
             call bl_error("ERROR: HSE boundaries not implemented for +X")

          elseif (xr_ext == EXT_INTERP) then
             ! interpolate thermodynamics from initial model

             do j = domhi(1)+1, adv_h1
                x = problo(1) + delta(1)*(dble(j) + HALF)

                   ! set all the variables even though we're testing on URHO
                   if (n == URHO) then

                      dens_zone = interpolate(x,npts_model,model_r, &
                                              model_state(:,idens_model))

                      temp_zone = interpolate(x,npts_model,model_r, &
                                              model_state(:,itemp_model))

                      do q = 1, nspec
                         X_zone(q) = interpolate(x,npts_model,model_r, &
                                                 model_state(:,ispec_model-1+q))
                      enddo


                      ! extrap normal momentum
                      adv(j,UMX) = ZERO

                      eos_state%rho = dens_zone
                      eos_state%T = temp_zone
                      eos_state%xn(:) = X_zone

                      call eos(eos_input_rt, eos_state)

                      pres_zone = eos_state%p
                      eint = eos_state%e

                      adv(j,URHO) = dens_zone
                      adv(j,UEINT) = dens_zone*eint
                      adv(j,UEDEN) = dens_zone*eint + &
                           HALF*sum(adv(j,UMX:UMZ)**2)/dens_zone
                      adv(j,UTEMP) = temp_zone
                      adv(j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                   endif

                enddo
          endif  ! xr_ext check

       endif

    enddo

  end subroutine ext_fill


  subroutine ext_denfill(adv,adv_l1,adv_h1, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ext_denfill")

    use prob_params_module, only : problo
    use interpolate_module
    use model_parser_module
    use bl_error_module

    use amrex_fort_module, only : rt => amrex_real
    integer adv_l1,adv_h1
    integer bc(2,*)
    integer domlo(1), domhi(1)
    real(rt)         delta(1), xlo(1), time
    real(rt)         adv(adv_l1:adv_h1)

    integer i,j
    real(rt)         x

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)


    ! XLO
    if ( bc(1,1) == EXT_DIR .and. adv_l1 < domlo(1)) then
       do j = adv_l1, domlo(1)-1
          x = problo(1) + delta(1)*(dble(j) + HALF)
             adv(j) = interpolate(x,npts_model,model_r,model_state(:,idens_model))
       end do
    end if

    ! XHI
    if ( bc(2,1) == EXT_DIR .and. adv_h1 > domhi(1)) then
       do j = domhi(1)+1, adv_h1
          x = problo(1) + delta(1)*(dble(j)+ HALF)
             adv(j) = interpolate(x,npts_model,model_r,model_state(:,idens_model))
       end do
    end if

  end subroutine ext_denfill

end module bc_ext_fill_module
