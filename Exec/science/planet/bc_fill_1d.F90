module bc_fill_module

  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                        domlo,domhi,delta,xlo,time,bc) bind(C)
 
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX,  UEDEN, UEINT, &
                                   UFS, UTEMP, const_grav
    use interpolate_module
    use eos_module
    use network, only: nspec
    use model_parser_module
    use bl_error_module
    use eos_type_module

    include 'AMReX_bc_types.fi'

    integer adv_l1,adv_h1
    integer bc(2,*)
    integer domlo(1), domhi(1)
    double precision delta(1), xlo(1), time
    double precision adv(adv_l1:adv_h1,NVAR)
    
    integer i,j,q,n
    double precision x
    double precision pres_above,p_want,pres_zone
    double precision temp_zone,X_zone(nspec),dens_zone
    double precision :: x_base, dens_base, slope
  
    type (eos_t) :: eos_state


    do n = 1,NVAR
       call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
                  domlo,domhi,delta,xlo,bc(1,n))
    enddo

    ! XLO -- HSE with linear density profile, T found via iteration
    ! we do all variables at once here
    if ( bc(1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then

       x_base = xlo(1) + delta(1)*(float(domlo(1)-adv_l1) + 0.5d0)

          dens_base = adv(domlo(1),URHO)

          ! density slope
          slope = (adv(domlo(1)+1,URHO) - adv(domlo(1),URHO))/delta(1)
          do j=domlo(1)-1,adv_l1,-1
             x = xlo(1) + delta(1)*(float(j-adv_l1) + 0.5d0)
             adv(j,:) = adv(j+1,:)

             ! HSE integration to get temperature, pressure
                    
             ! density is linear from the last two zones
             dens_zone = dens_base + slope*(x - x_base)

             ! temperature guess and species held constant in BCs
             temp_zone = adv(j+1,UTEMP)
             X_zone(:) = adv(j+1,UFS:UFS-1+nspec)/adv(j+1,URHO)
             
             ! get pressure in zone above
             eos_state%rho = adv(j+1,URHO)
             eos_state%T = adv(j+1,UTEMP)
             eos_state%xn(:) = adv(j+1,UFS:UFS-1+nspec)/adv(j+1,URHO)
             
             call eos(eos_input_rt, eos_state)
             
             pres_above = eos_state%p

             ! pressure needed from HSE
             p_want = pres_above - &
                  delta(1)*0.5d0*(dens_zone + adv(j+1,URHO))*const_grav

             ! EOS with HSE pressure + linear density profile yields T, e, ...
             eos_state%rho = dens_zone
             eos_state%T = temp_zone   ! guess
             eos_state%xn(:) = X_zone(:)
             eos_state%p = p_want
 

            
             call eos(eos_input_rp, eos_state)

            ! velocity
             if (zero_vels) then
                
                ! zero gradient velocity
                adv(j,UMX) = dens_zone*(adv(domlo(1),UMX)/adv(domlo(1),URHO))
             endif
           
             adv(j,URHO) = dens_zone
             adv(j,UEINT) = dens_zone*eos_state%e
             adv(j,UEDEN) = dens_zone*eos_state%e + &
                  0.5d0*(adv(j,UMX)**2)/dens_zone
             adv(j,UTEMP) = eos_state%T
             adv(j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
             
          end do
    end if

  
  
    ! YHI
    do n = 1, nvar
       if ( bc(2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
          
          do j=domhi(1)+1,adv_h1
             x = xlo(1) + delta(1)*(float(j-adv_l1) + 0.5d0)
             ! zero-gradient catch-all -- this will get the radiation
             ! energy
             adv(j,:) = adv(j-1,:)
             adv(j,UTEMP)=min(adv(j-1,UTEMP),temp_zone)

                ! set all the variables even though we're testing on URHO
                if (n .eq. URHO) then
                     
                   dens_zone = interpolate(x,npts_model,model_r, &
                                           model_state(:,idens_model))

                   temp_zone = interpolate(x,npts_model,model_r, &
                                           model_state(:,itemp_model))
                   do q = 1, nspec
                      X_zone(q) = interpolate(x,npts_model,model_r, &
                                             model_state(:,ispec_model-1+q))
                   enddo

                   ! extrap normal momentum
                   adv(j,UMX) = min(0.d0,adv(domhi(1),UMX))
                   
                   ! zero transverse momentum

                   eos_state%rho = dens_zone
                   eos_state%T = temp_zone
                   eos_state%xn(:) = X_zone

                   
                   call eos(eos_input_rt, eos_state)
                   
                   adv(j,URHO) = dens_zone
                   adv(j,UEINT) = dens_zone*eos_state%e
                   adv(j,UEDEN) = dens_zone*eos_state%e + &
                        0.5d0*(adv(j,UMX)**2)/dens_zone
                   adv(j,UTEMP) = temp_zone
                   adv(j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
               
                end if

          end do
       end if
     
    end do

  end subroutine ca_hypfill

  
  subroutine ca_denfill(adv,adv_l1,adv_h1, &
                        domlo,domhi,delta,xlo,time,bc) bind(C)
    
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, &
         UFS, UTEMP, const_grav
    use bl_error_module
    use interpolate_module
    use model_parser_module

    implicit none
    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_h1
    integer bc(2,*)
    integer domlo(1), domhi(1)
    double precision delta(1), xlo(1), time
    double precision adv(adv_l1:adv_h1)

    integer i,j,q,n
    double precision x
    double precision :: x_base, dens_base, slope
    double precision TOL

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)



   end subroutine ca_denfill

  
  subroutine ca_gravxfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_h1
    integer :: bc(2,*)
    integer :: domlo(1), domhi(1)
    double precision delta(1), xlo(1), time
    double precision grav(grav_l1:grav_h1)
    integer :: i, j

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

    ! our lower boundary is inflow, so we need to make sure the
    ! gravitational acceleration is set correctly there
    !     YLO

  end subroutine ca_gravxfill


  subroutine ca_gravyfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use meth_params_module, only: const_grav

    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_h1
    integer :: bc(2,*)
    integer :: domlo(1), domhi(1)
    double precision delta(1), xlo(1), time
    double precision grav(grav_l1:grav_h1)
    integer :: i, j

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

    ! our lower boundary is inflow, so we need to make sure the
    ! gravitational acceleration is set correctly there
    !     YLO

  end subroutine ca_gravyfill


  subroutine ca_gravzfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use meth_params_module, only: const_grav

    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_h1
    integer :: bc(2,*)
    integer :: domlo(1), domhi(1)
    double precision delta(1), xlo(1), time
    double precision grav(grav_l1:grav_h1)
    integer :: i, j

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

    ! our lower boundary is inflow, so we need to make sure the
    ! gravitational acceleration is set correctly there
    !     YLO

  end subroutine ca_gravzfill


  subroutine ca_reactfill(react,react_l1, &
                          react_h1,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: react_l1,react_h1
    integer :: bc(2,*)
    integer :: domlo(1), domhi(1)
    double precision delta(1), xlo(1), time
    double precision react(react_l1:react_h1)

    call filcc(react,react_l1,react_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_reactfill


  subroutine ca_radfill(rad,rad_l1, &
                        rad_h1,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: rad_l1,rad_h1
    integer :: bc(2,*)
    integer :: domlo(1), domhi(1)
    double precision delta(1), xlo(1), time
    double precision rad(rad_l1:rad_h1)

    integer :: j

    call filcc(rad,rad_l1,rad_h1,domlo,domhi,delta,xlo,bc)

    ! we are inflow at the lower boundary, so we need to take the appropriate
    ! action for the radiation here (during the hydro step)

    ! this do loop counts backwards since we want to work downward

  end subroutine ca_radfill

  
  subroutine ca_phigravfill(phi,phi_l1, &
                            phi_h1,domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_h1
    integer          :: bc(2,*)
    integer          :: domlo(1), domhi(1)
    double precision :: delta(1), xlo(1), time
    double precision :: phi(phi_l1:phi_h1)

    call filcc(phi,phi_l1,phi_h1, &
         domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill

end module bc_fill_module
