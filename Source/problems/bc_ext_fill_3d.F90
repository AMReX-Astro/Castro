module bc_ext_fill_module

  use amrex_constants_module, only: ZERO, HALF
  use amrex_error_module
  use amrex_fort_module, only: rt => amrex_real
  use amrex_filcc_module, only: amrex_filccn
  use interpolate_module, only: interpolate_sub
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UTEMP, const_grav, &
                                 hse_zero_vels, hse_interp_temp, hse_reflect_vels, &
                                 xl_ext, xr_ext, yl_ext, yr_ext, zl_ext,zr_ext,EXT_HSE, EXT_INTERP

  implicit none

  include 'AMReX_bc_types.fi'

contains

  ! this module contains different routines for filling the
  ! hydrodynamics boundary conditions

  ! NOTE: the hydrostatic boundary conditions here rely on
  ! constant gravity

  subroutine ext_fill(adv, adv_l1, adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                      domlo, domhi, delta, xlo, time, bc) &
                      bind(C, name="ext_fill")

    use prob_params_module, only : problo
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec
    use model_parser_module, only: model_r, model_state, npts_model, idens_model, itemp_model, ispec_model
 
    integer, intent(in) :: adv_l1, adv_h1, adv_l2, adv_h2, adv_l3, adv_h3
    integer, intent(in) :: bc(3,2,*)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

    integer :: i, j, k, q, n, iter, m, koff

    real(rt) :: y, z
    real(rt) :: dens_above, dens_base, temp_above
    real(rt) :: pres_above, p_want, pres_zone, A
    real(rt) :: drho, dpdr, temp_zone, eint, X_zone(nspec), dens_zone

    integer, parameter :: MAX_ITER = 100
    real(rt), parameter :: TOL = 1.e-8_rt
    logical :: converged_hse

    type (eos_t) :: eos_state

    do n = 1, NVAR

#ifndef AMREX_USE_CUDA
       !XLO
       if (bc(1,1,n) == EXT_DIR .and. xl_ext == EXT_HSE .and. adv_l1 < domlo(1)) then
          call amrex_error("ERROR: HSE boundaries not implemented for -X")
       end if

       !YLO
       if (bc(2,1,n) == EXT_DIR .and. yl_ext == EXT_HSE .and. adv_l2 < domlo(2)) then
          call amrex_error("ERROR: HSE boundaries not implemented for -Y")
       end if

       ! XHI
       if (bc(1,2,n) == EXT_DIR .and. xr_ext == EXT_HSE .and. adv_h1 > domhi(1)) then
          call amrex_error("ERROR: HSE boundaries not implemented for +X")
       end if
       ! YHI
      if (bc(2,2,n) == EXT_DIR .and. yr_ext == EXT_HSE .and. adv_h2 > domhi(2)) then
         call amrex_error("ERROR: HSE boundaries not implemented for +Y")
      end if
#endif

       ! ZLO
       if (bc(3,1,n) == EXT_DIR .and. adv_l3 < domlo(3)) then

          if (zl_ext == EXT_HSE) then

             ! we will fill all the variables when we consider URHO
             if (n == URHO) then
                do j= adv_l2, adv_h2
                  do i = adv_l1, adv_h1
                   ! we are integrating along a column at constant i.
                   ! Make sure that our starting state is well-defined
                    dens_above = adv(i,j,domlo(3),URHO)

                   ! sometimes, we might be working in a corner
                   ! where the ghost cells above us have not yet
                   ! been initialized.  In that case, take the info
                   ! from the initial model
                    if (dens_above == ZERO) then
                      z = problo(3) + delta(3)*(dble(domlo(3)) + HALF)

                      call interpolate_sub(dens_above, z,npts_model,model_r, &
                                           model_state(:,idens_model))

                      call interpolate_sub(temp_above, z,npts_model,model_r, &
                                           model_state(:,itemp_model))

                      do m = 1, nspec
                         call interpolate_sub(X_zone(m), z,npts_model,model_r, &
                                              model_state(:,ispec_model-1+m))
                      enddo

                    else
                      temp_above = adv(i,j,domlo(3),UTEMP)
                      X_zone(:) = adv(i,j,domlo(3),UFS:UFS-1+nspec)/dens_above
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
                   do k = domlo(3)-1, adv_l3, -1
                      z = problo(3) + delta(3)*(dble(k) + HALF)

                      ! HSE integration to get density, pressure

                      ! initial guesses
                      dens_zone = dens_above

                      ! temperature and species held constant in BCs
                      if (hse_interp_temp == 1) then
                         call interpolate_sub(temp_zone, z,npts_model,model_r, &
                                              model_state(:,itemp_model))
                      else
                         temp_zone = temp_above
                      endif

                      converged_hse = .FALSE.


                      do iter = 1, MAX_ITER

                         ! pressure needed from HSE
                         p_want = pres_above - &
                              delta(3)*HALF*(dens_zone + dens_above)*const_grav

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

#ifndef AMREX_USE_CUDA
                      if (.not. converged_hse) then
                         print *, "i, j, k,domlo(3): ", i, j, k, domlo(3)
                         print *, "p_want:    ", p_want
                         print *, "dens_zone: ", dens_zone
                         print *, "temp_zone: ", temp_zone
                         print *, "drho:      ", drho
                         print *, " "
                         print *, "column info: "
                         print *, "   dens: ", adv(i,j,k:domlo(3),URHO)
                         print *, "   temp: ", adv(i,j,k:domlo(3),UTEMP)
                         call amrex_error("ERROR in bc_ext_fill_1d: failure to converge in -Z BC")
                      endif
#endif

                      ! velocity
                      if (hse_zero_vels == 1) then

                         ! zero normal momentum causes pi waves to pass through
                         adv(i,j,k,UMX) = ZERO
                         adv(i,j,k,UMY) = ZERO
                         adv(i,j,k,UMZ) = ZERO

                      else

                         if (hse_reflect_vels == 1) then
                            ! reflect normal, zero gradient for transverse
                            ! note: we need to match the corresponding
                            ! zone on the other side of the interface
                            koff = domlo(3)-k-1

                            adv(i,j,k,UMZ) = -dens_zone*(adv(i,j,domlo(3)+koff,UMZ)/adv(i,j,domlo(3)+koff,URHO))

                            adv(i,j,k,UMX) = -dens_zone*(adv(i,j,domlo(3),UMX)/dens_base)
                            adv(i,j,k,UMY) = -dens_zone*(adv(i,j,domlo(3),UMY)/dens_base)
                         else
                            ! zero gradient
                            adv(i,j,k,UMX) = dens_zone*(adv(i,j,domlo(1),UMX)/dens_base)
                            adv(i,j,k,UMY) = dens_zone*(adv(i,j,domlo(2),UMY)/dens_base)
                            adv(i,j,k,UMZ) = dens_zone*(adv(i,j,domlo(3),UMZ)/dens_base)
                         endif
                      endif
                      eos_state%rho = dens_zone
                      eos_state%T = temp_zone
                      eos_state%xn(:) = X_zone

                      call eos(eos_input_rt, eos_state)

                      pres_zone = eos_state%p
                      eint = eos_state%e

                      ! store the final state
                      adv(i,j,k,URHO) = dens_zone
                      adv(i,j,k,UEINT) = dens_zone*eint
                      adv(i,j,k,UEDEN) = dens_zone*eint + &
                           HALF*sum(adv(i,j,k,UMX:UMZ)**2)/dens_zone
                      adv(i,j,k,UTEMP) = temp_zone
                      adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                      ! for the next zone
                      dens_above = dens_zone
                      pres_above = pres_zone

                    end do
                  end do
                end do
             endif  ! n == URHO

          elseif (zl_ext == EXT_INTERP) then

            do k = domlo(3)-1, adv_l3, -1
              z = problo(3) + delta(3)*(dble(k)+HALF)

              do j = adv_l2, adv_h2

                do i = adv_l1, adv_h1
                   ! set all the variables even though we're testing on URHO
                  if (n == URHO) then

                     call interpolate_sub(dens_zone, z,npts_model,model_r, &
                                          model_state(:,idens_model))

                     call interpolate_sub(temp_zone, z,npts_model,model_r, &
                                          model_state(:,itemp_model))

                    do q = 1, nspec
                       call interpolate_sub(X_zone(q), z,npts_model,model_r, &
                                            model_state(:,ispec_model-1+q))
                    enddo

                    ! extrap normal momentum
                    adv(i,j,k,UMZ) = min(ZERO, adv(i,j,domlo(3),UMZ))

                      ! zero transverse momentum
                    adv(i,j,k,UMX) = ZERO
                    adv(i,j,k,UMY) = ZERO

                    eos_state%rho = dens_zone
                    eos_state%T = temp_zone
                    eos_state%xn(:) = X_zone

                    call eos(eos_input_rt, eos_state)

                    pres_zone = eos_state%p
                    eint = eos_state%e

                    adv(i,j,k,URHO) = dens_zone
                    adv(i,j,k,UEINT) = dens_zone*eint
                    adv(i,j,k,UEDEN) = dens_zone*eint + &
                          HALF*sum(adv(i,j,k,UMX:UMZ)**2)/dens_zone
                    adv(i,j,k,UTEMP) = temp_zone
                    adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
                  endif
                end do
              enddo
            end do
          endif  ! zl_ext check


       endif


       ! ZHI
       if (bc(3,2,n) == EXT_DIR .and. adv_h3 > domhi(3)) then

          if (zr_ext == EXT_HSE) then
#ifndef AMREX_USE_CUDA
             call amrex_error("ERROR: HSE boundaries not implemented for +Z")
#endif

          elseif (zr_ext == EXT_INTERP) then
             ! interpolate thermodynamics from initial model

            do k = domhi(3)+1, adv_h3
              z = problo(3) + delta(3)*(dble(k) + HALF)

              do j = adv_l2,adv_h2
                do i = adv_l1, adv_h1

                   ! set all the variables even though we're testing on URHO
                  if (n == URHO) then
                     call interpolate_sub(dens_zone, z,npts_model,model_r, &
                                          model_state(:,idens_model))

                     call interpolate_sub(temp_zone, z,npts_model,model_r, &
                                          model_state(:,itemp_model))

                    do q = 1, nspec
                       call interpolate_sub(X_zone(q), z,npts_model,model_r, &
                                            model_state(:,ispec_model-1+q))
                    enddo


                    ! extrap normal momentum
                    adv(i,j,k,UMZ) = max(ZERO, adv(i,j,domhi(3),UMZ))

                    ! zero transverse momentum
                    adv(i,j,k,UMX) = ZERO
                    adv(i,j,k,UMY) = ZERO

                    eos_state%rho = dens_zone
                    eos_state%T = temp_zone
                    eos_state%xn(:) = X_zone

                    call eos(eos_input_rt, eos_state)

                    pres_zone = eos_state%p
                    eint = eos_state%e

                    adv(i,j,k,URHO) = dens_zone
                    adv(i,j,k,UEINT) = dens_zone*eint
                    adv(i,j,k,UEDEN) = dens_zone*eint + &
                          HALF*sum(adv(i,j,k,UMX:UMZ)**2)/dens_zone
                    adv(i,j,k,UTEMP) = temp_zone
                    adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                  endif
                end do
              enddo
            enddo
          endif  ! zr_ext check

       endif

    enddo


  end subroutine ext_fill


  subroutine ext_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, adv_h3,&
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ext_denfill")

    use prob_params_module, only : problo
    use interpolate_module
    use model_parser_module
    use amrex_error_module

    implicit none

    integer, intent(in) :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer, intent(in) :: bc(3,2)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

    integer :: i, j, k
    real(rt) :: y, z

#ifndef AMREX_USE_CUDA
    ! XLO
    if ( bc(1,1) == EXT_DIR .and. adv_l1 < domlo(1)) then
       call amrex_error("We shoundn't be here (xlo denfill)")
    end if

    ! XHI
    if ( bc(1,2) == EXT_DIR .and. adv_h1 > domhi(1)) then
       call amrex_error("We shoundn't be here (xhi denfill)")
    endif

    ! YLO
    if ( bc(2,1) == EXT_DIR .and. adv_l2 < domlo(2)) then
       call amrex_error("We shoundn't be here (ylo denfill)")
    end if

    ! YHI
    if ( bc(2,2) == EXT_DIR .and. adv_h2 > domhi(2)) then
       call amrex_error("We shoundn't be here (yhi denfill)")
    endif
#endif

    ! ZLO
    if ( bc(3,1) == EXT_DIR .and. adv_l3 < domlo(3)) then
       do k = adv_l3, domlo(3)-1
          z = problo(3) + delta(3)*(dble(k) + HALF)
          do j = adv_l2, adv_h2
             do i = adv_l1, adv_h1
                call interpolate_sub(adv(i,j,k), z,npts_model,model_r,model_state(:,idens_model))
             end do
          end do
       end do
    end if

    ! ZHI
    if ( bc(3,2) == EXT_DIR .and. adv_h3 > domhi(3)) then
       do k = domhi(3)+1, adv_h3
          z = problo(3) + delta(3)*(dble(k)+ HALF)
          do j = adv_l2, adv_h2
             do i = adv_l1, adv_h1
                call interpolate_sub(adv(i,j,k), z,npts_model,model_r,model_state(:,idens_model))
             end do
          end do
       end do
    end if

  end subroutine ext_denfill

end module bc_ext_fill_module
