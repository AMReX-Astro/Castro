!AUG10
module bc_fill_module
  use bc_ext_fill_module
  use bl_constants_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none
  include 'AMReX_bc_types.fi'
  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")
 
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, &
                                   UFS, UTEMP, const_grav, &
                                   hse_zero_vels, hse_interp_temp, hse_reflect_vels, &
                                   xl_ext,xr_ext,yl_ext, yr_ext, EXT_HSE, EXT_INTERP
    use amrex_fort_module, only : rt => amrex_real
    use interpolate_module
    use eos_module
    use network, only: nspec
    use model_parser_module
    use bl_error_module
    use eos_type_module

    include 'AMReX_bc_types.fi'

    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt) delta(2), xlo(2), time
    real(rt) adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)
    
    integer i,j,q,n
    real(rt) y
    real(rt) pres_above,p_want,pres_zone
    real(rt) temp_zone,X_zone(nspec),dens_zone
    real(rt) :: y_base, dens_base, slope
  
    type (eos_t) :: eos_state
    !need to fix
    yl_ext=EXT_HSE
    do n = 1,NVAR
    call filcc(adv(:,:,n),adv_l1,adv_l2,adv_h1,adv_h2, &
    domlo,domhi,delta,xlo,bc(:,:,n))
    enddo
    do n = 1, NVAR
         
       ! XLO
       if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then

          ! we are periodic in x -- we should never get here
          call bl_error("ERROR: invalid BC in Prob_2d.f90")
        
       end if
     
       ! XHI
       if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
          
          ! we are periodic in x -- we should never get here
          call bl_error("ERROR: invalid BC in Prob_2d.f90")
        
       end if

        if (xr_ext == EXT_HSE .or. xl_ext == EXT_HSE .or. xr_ext == EXT_INTERP .or. &
        xl_ext == EXT_INTERP) then
        call bl_error("ERROR: HSE boundaries not implemented for +,- X")
        end if



    enddo


    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
      if (yl_ext == EXT_HSE) then
      call ext_fill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                    domlo,domhi,delta,xlo,time,bc)
      else
       y_base = xlo(2) + delta(2)*(dble(domlo(2)-adv_l2) + HALF)
       do i=adv_l1,adv_h1

          dens_base = adv(i,domlo(2),URHO)

          ! density slope
          slope = (adv(i,domlo(2)+1,URHO) - adv(i,domlo(2),URHO))/delta(2)
          ! this do loop counts backwards since we want to work downward
          do j=domlo(2)-1,adv_l2,-1
             y = xlo(2) + delta(2)*(dble(j-adv_l2) + HALF)
             ! zero-gradient catch-all -- this will get the radiation
             ! energy
             adv(i,j,:) = adv(i,j+1,:)

             ! HSE integration to get temperature, pressure
                    
             ! density is linear from the last two zones
             dens_zone = dens_base + slope*(y - y_base)

             ! temperature guess and species held constant in BCs
             temp_zone = adv(i,j+1,UTEMP)

             X_zone(:) = adv(i,j+1,UFS:UFS-1+nspec)/adv(i,j+1,URHO)
             
             ! get pressure in zone above
             eos_state%rho = adv(i,j+1,URHO)
             eos_state%T = adv(i,j+1,UTEMP)
             eos_state%xn(:) = adv(i,j+1,UFS:UFS-1+nspec)/adv(i,j+1,URHO)
             
             call eos(eos_input_rt, eos_state)
             
             pres_above = eos_state%p

             ! pressure needed from HSE
             p_want = pres_above - &
                  delta(2)*HALF*(dens_zone + adv(i,j+1,URHO))*const_grav

             ! EOS with HSE pressure + linear density profile yields T, e, ...
             eos_state%rho = dens_zone
             eos_state%T = temp_zone   ! guess
             eos_state%xn(:) = X_zone(:)
             eos_state%p = p_want
 

            
             call eos(eos_input_rp, eos_state)

             ! velocity
             if (zero_vels) then
                
                ! zero normal momentum causes pi waves to pass through
                adv(i,j,UMY) = ZERO
                
                ! zero transverse momentum
                adv(i,j,UMX) = ZERO
             else
                
                ! zero gradient velocity
                adv(i,j,UMY) = min(ZERO,dens_zone*(adv(i,domlo(2),UMY)/adv(i,domlo(2),URHO)))
             endif
           
             adv(i,j,URHO) = dens_zone
             adv(i,j,UEINT) = dens_zone*eos_state%e
             adv(i,j,UEDEN) = dens_zone*eos_state%e + & 
                  HALF*(adv(i,j,UMX)**2.0_rt+adv(i,j,UMY)**2.0_rt)/dens_zone
             adv(i,j,UTEMP) = eos_state%T
             adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
             
          end do
       end do
      end if
    end if

  
  
    ! YHI
       if ( bc(2,2,URHO).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then

        if (yr_ext == EXT_HSE) then
        call bl_error("ERROR: HSE boundaries not implemented for +Y")
        elseif (yr_ext == EXT_INTERP) then
          call ext_fill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc)
        end if


          do j=domhi(2)+1,adv_h2
             y = xlo(2) + delta(2)*(dble(j-adv_l2) + HALF)

             adv(adv_l1:adv_h1,j,:) = adv(adv_l1:adv_h1,domhi(2),:)

             do i=adv_l1,adv_h1
                dens_zone=adv(i,domhi(2),URHO)
                temp_zone=adv(i,domhi(2),UTEMP)

                   do q = 1, nspec
                      X_zone(q) = interpolate(y,npts_model,model_r, &
                                              model_state(:,ispec_model-1+q))
                   enddo

                   ! extrap normal momentum
                   adv(i,j,UMY) = min(ZERO,adv(i,domhi(2),UMY))
                   
                   ! zero transverse momentum
                   adv(i,j,UMX) = ZERO
                   
                   eos_state%rho = dens_zone
                   eos_state%T = temp_zone
                   eos_state%xn(:) = X_zone

                   
                   call eos(eos_input_rt, eos_state)
                   
                   adv(i,j,URHO) = dens_zone
                   adv(i,j,UEINT) = dens_zone*eos_state%e
                   adv(i,j,UEDEN) = dens_zone*eos_state%e + &
                        HALF*(adv(i,j,UMX)**2.0_rt+adv(i,j,UMY)**2.0_rt)/dens_zone
                   adv(i,j,UTEMP) = temp_zone
                   adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
               

             end do
          end do
       end if


  end subroutine ca_hypfill

  
  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")
    
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, &
         UFS, UTEMP, const_grav
    use bl_error_module
    use interpolate_module
    use model_parser_module

    implicit none
    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt) delta(2), xlo(2), time
    real(rt) adv(adv_l1:adv_h1,adv_l2:adv_h2)

    integer i,j,q,n
    real(rt) y
    real(rt) :: y_base, dens_base, slope
    real(rt) TOL

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       call bl_error("We shoundn't be here (xlo denfill)")
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       call bl_error("We shoundn't be here (xlo denfill)")
    endif


  end subroutine ca_denfill

  
  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt) delta(2), xlo(2), time
    real(rt) grav(grav_l1:grav_h1,grav_l2:grav_h2)
    integer :: i, j

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)


  end subroutine ca_gravxfill


  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use meth_params_module, only: const_grav

    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt) delta(2), xlo(2), time
    real(rt) grav(grav_l1:grav_h1,grav_l2:grav_h2)
    integer :: i, j

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)


  end subroutine ca_gravyfill


  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use meth_params_module, only: const_grav

    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt) delta(2), xlo(2), time
    real(rt) grav(grav_l1:grav_h1,grav_l2:grav_h2)
    integer :: i, j

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

 
  end subroutine ca_gravzfill


  subroutine ca_reactfill(react,react_l1,react_l2, &
                          react_h1,react_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: react_l1,react_l2,react_h1,react_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt) delta(2), xlo(2), time
    real(rt) react(react_l1:react_h1,react_l2:react_h2)

    call filcc(react,react_l1,react_l2,react_h1,react_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_reactfill


  subroutine ca_radfill(rad,rad_l1,rad_l2, &
                        rad_h1,rad_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: rad_l1,rad_l2,rad_h1,rad_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt) delta(2), xlo(2), time
    real(rt) rad(rad_l1:rad_h1,rad_l2:rad_h2)

    integer :: j

    call filcc(rad,rad_l1,rad_l2,rad_h1,rad_h2,domlo,domhi,delta,xlo,bc)

    ! we are inflow at the lower boundary, so we need to take the appropriate
    ! action for the radiation here (during the hydro step)

    ! this do loop counts backwards since we want to work downward
    !YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. rad_l2.lt.domlo(2)) then
       do j=rad_l2,domlo(2)-1

          ! zero-gradient catch-all -- this will get the radiation
          ! energy
          rad(rad_l1:rad_h1,j) = rad(rad_l1:rad_h1,domlo(2))
       enddo
    endif

    !YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. rad_h2.gt.domhi(2)) then
       do j=domhi(2)+1,rad_h2

          ! zero-gradient catch-all -- this will get the radiation
          ! energy
          rad(rad_l1:rad_h1,j) = rad(rad_l1:rad_h1,domhi(2))
       enddo
    endif


  end subroutine ca_radfill


  subroutine ca_phigravfill(phi,phi_l1,phi_l2, &
                            phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_h1,phi_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt) :: delta(2), xlo(2), time
    real(rt) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
         domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill

end module bc_fill_module
