module bc_ext_fill_module

  use bl_constants_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UTEMP, const_grav, &
                                 hse_zero_vels, &
                                 xl_ext, xr_ext, yl_ext, yr_ext, EXT_HSE, EXT_INTERP
  use interpolate_module

  implicit none

  include 'bc_types.fi'

  private

  public :: ext_fill, ext_denfill

contains

  ! this module contains different routines for filling the
  ! hydrodynamics boundary conditions

  ! NOTE: the hydrostatic boundary conditions here rely on
  ! constant gravity

  subroutine ext_fill(adv, adv_l1, adv_l2, adv_h1, adv_h2, &
                      domlo, domhi, delta, xlo, time, bc) &
                      bind(C, name="ext_fill")

    use probdata_module
    use eos_module
    use network, only: nspec
    use model_parser_module
    
    integer adv_l1, adv_l2, adv_h1, adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    double precision delta(2), xlo(2), time
    double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer i, j, q, n, iter, MAX_ITER
    double precision y
    double precision pres_above, p_want, pres_zone, A
    double precision drho, dpdr, temp_zone, eint, X_zone(nspec), dens_zone
    double precision TOL
    logical converged_hse

    type (eos_t) :: eos_state

    MAX_ITER = 100
    TOL = 1.d-8

    do n = 1, NVAR

       ! XLO
       if (bc(1,1,n) == EXT_DIR .and. xl_ext == EXT_HSE .and. adv_l1 < domlo(1)) then
          call bl_error("ERROR: HSE boundaries not implemented for -X")
       end if

       ! XHI
       if (bc(1,2,n) == EXT_DIR .and. xr_ext == EXT_HSE .and. adv_h1 > domhi(1)) then
          call bl_error("ERROR: HSE boundaries not implemented for +X")
       end if

       ! YLO
       if (bc(2,1,n) == EXT_DIR .and. adv_l2 < domlo(2)) then
          
          if (yl_ext == EXT_HSE) then

             ! this do loop counts backwards since we want to work downward
             do j = domlo(2)-1, adv_l2, -1
                y = xlo(2) + delta(2)*(dble(j-adv_l2) + HALF)
                
                do i = adv_l1, adv_h1
                   
                   ! set all the variables even though we're testing on URHO
                   if (n == URHO) then
                      
                      ! HSE integration to get density, pressure
                      
                      ! initial guesses
                      dens_zone = adv(i,j+1,URHO)
                      
                      ! temperature and species held constant in BCs
                      temp_zone = adv(i,j+1,UTEMP)
                      X_zone(:) = adv(i,j+1,UFS:UFS-1+nspec)/adv(i,j+1,URHO)
                      
                      ! get pressure in zone above
                      eos_state%rho = adv(i,j+1,URHO)
                      eos_state%T = adv(i,j+1,UTEMP)
                      eos_state%xn(:) = adv(i,j+1,UFS:UFS-1+nspec)/adv(i,j+1,URHO)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      eint = eos_state%e
                      pres_above = eos_state%p
                      
                      converged_hse = .FALSE.
                      
                      do iter = 1, MAX_ITER
                         
                         ! pressure needed from HSE
                         p_want = pres_above - &
                              delta(2)*HALF*(dens_zone + adv(i,j+1,URHO))*const_grav
                         
                         ! pressure from EOS
                         eos_state%rho = dens_zone
                         eos_state%T = temp_zone
                         eos_state%xn(:) = X_zone
                         
                         call eos(eos_input_rt, eos_state)
                         
                         pres_zone = eos_state%p
                         dpdr = eos_state%dpdr
                         eint = eos_state%e
                         
                         ! Newton-Raphson - we want to zero A = p_want - p(rho)
                         A = p_want - pres_zone
                         drho = A/(dpdr + HALF*delta(2)*const_grav)
                         
                         dens_zone = max(0.9_dp_t*dens_zone, &
                              min(dens_zone + drho, 1.1_dp_t*dens_zone))
                         
                         ! convergence?
                         if (abs(drho) < TOL*dens_zone) then
                            converged_hse = .TRUE.
                            exit
                         endif
                         
                      enddo
                      
                      if (.not. converged_hse) call bl_error("ERROR: failure to converge in -Y BC")
                      
                   endif
                   
                   ! velocity
                   if (hse_zero_vels == 1) then
                      
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
                   
                   eos_state%rho = dens_zone
                   eos_state%T = temp_zone
                   eos_state%xn(:) = X_zone
                   
                   call eos(eos_input_rt, eos_state)
                   
                   pres_zone = eos_state%p
                   eint = eos_state%e
                   
                   adv(i,j,URHO) = dens_zone
                   adv(i,j,UEINT) = dens_zone*eint
                   adv(i,j,UEDEN) = dens_zone*eint + & 
                        HALF*(adv(i,j,UMX)**2 + adv(i,j,UMY)**2 + adv(i,j,UMZ)**2)/dens_zone
                   adv(i,j,UTEMP) = temp_zone
                   adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
                   
                enddo
             enddo
             
          elseif (yl_ext == EXT_INTERP) then

             do j = domlo(2)-1, adv_l2, -1
                y = xlo(2) + delta(2)*(dble(j-adv_l2) + HALF)
                
                do i = adv_l1, adv_h1
             
                   ! set all the variables even though we're testing on URHO
                   if (n == URHO) then

                      dens_zone = interpolate(y,npts_model,model_r, &
                                              model_state(:,idens_model)) 
                   
                      temp_zone = interpolate(y,npts_model,model_r, &
                                              model_state(:,itemp_model))

                      do q = 1, nspec
                         X_zone(q) = interpolate(y,npts_model,model_r, &
                                                 model_state(:,ispec_model-1+q))
                      enddo
                   
                      ! extrap normal momentum
                      adv(i,j,UMY) = min(ZERO, adv(i,domhi(2),UMY))
                   
                      ! zero transverse momentum
                      adv(i,j,UMX) = ZERO
                      adv(i,j,UMZ) = ZERO
                      
                      eos_state%rho = dens_zone
                      eos_state%T = temp_zone
                      eos_state%xn(:) = X_zone
                      
                      call eos(eos_input_rt, eos_state)
                      
                      pres_zone = eos_state%p
                      eint = eos_state%e
                      
                      adv(i,j,URHO) = dens_zone
                      adv(i,j,UEINT) = dens_zone*eint
                      adv(i,j,UEDEN) = dens_zone*eint + &
                           HALF*sum(adv(i,j,UMX:UMZ)**2)/dens_zone
                      adv(i,j,UTEMP) = temp_zone
                      adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
                   endif

                enddo
             enddo
          endif  ! yl_ext check

          
       endif


       ! YHI
       if (bc(2,2,n) == EXT_DIR .and. adv_h2 > domhi(2)) then

          if (yr_ext == EXT_HSE) then
             call bl_error("ERROR: HSE boundaries not implemented for +Y")

          elseif (yr_ext == EXT_INTERP) then
             ! interpolate thermodynamics from initial model
             
             do j = domhi(2)+1, adv_h2
                y = xlo(2) + delta(2)*(dble(j-adv_l2) + HALF)

                do i = adv_l1, adv_h1

                   ! set all the variables even though we're testing on URHO
                   if (n == URHO) then

                      dens_zone = interpolate(y,npts_model,model_r, &
                                              model_state(:,idens_model)) 

                      temp_zone = interpolate(y,npts_model,model_r, &
                                              model_state(:,itemp_model))

                      do q = 1, nspec
                         X_zone(q) = interpolate(y,npts_model,model_r, &
                                                 model_state(:,ispec_model-1+q))
                      enddo


                      ! extrap normal momentum
                      adv(i,j,UMY) = max(ZERO, adv(i,domhi(2),UMY))

                      ! zero transverse momentum
                      adv(i,j,UMX) = ZERO
                      adv(i,j,UMZ) = ZERO

                      eos_state%rho = dens_zone
                      eos_state%T = temp_zone
                      eos_state%xn(:) = X_zone

                      call eos(eos_input_rt, eos_state)

                      pres_zone = eos_state%p
                      eint = eos_state%e

                      adv(i,j,URHO) = dens_zone
                      adv(i,j,UEINT) = dens_zone*eint
                      adv(i,j,UEDEN) = dens_zone*eint + &
                           HALF*sum(adv(i,j,UMX:UMZ)**2)/dens_zone
                      adv(i,j,UTEMP) = temp_zone
                      adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                   endif

                enddo
             enddo
          endif  ! yr_ext check

       endif

    enddo

  end subroutine ext_fill


  subroutine ext_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ext_denfill")

    use probdata_module
    use interpolate_module
    use model_parser_module
    use bl_error_module

    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    double precision delta(2), xlo(2), time
    double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

    integer i,j
    double precision y

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    ! XLO
    if ( bc(1,1,1) == EXT_DIR .and. adv_l1 < domlo(1)) then
       call bl_error("We shoundn't be here (xlo denfill)")
    end if

    ! XHI
    if ( bc(1,2,1) == EXT_DIR .and. adv_h1 > domhi(1)) then
       call bl_error("We shoundn't be here (xlo denfill)")
    endif

    ! YLO
    if ( bc(2,1,1) == EXT_DIR .and. adv_l2 < domlo(2)) then
       do j = adv_l2, domlo(2)-1
          y = xlo(2) + delta(2)*(dble(j-adv_l2) + HALF)
          do i = adv_l1, adv_h1
             adv(i,j) = interpolate(y,npts_model,model_r,model_state(:,idens_model))
          end do
       end do
    end if

    ! YHI
    if ( bc(2,2,1) == EXT_DIR .and. adv_h2 > domhi(2)) then
       do j = domhi(2)+1, adv_h2
          y = xlo(2) + delta(2)*(dble(j-adv_l2)+ HALF)
          do i = adv_l1, adv_h1
             adv(i,j) = interpolate(y,npts_model,model_r,model_state(:,idens_model))
          end do
       end do
    end if

  end subroutine ext_denfill

end module bc_ext_fill_module
