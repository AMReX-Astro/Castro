!Aug25
module bc_fill_module
  use bc_ext_fill_module
  use bl_constants_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none
  include 'AMReX_bc_types.fi'
  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")
 
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX,  UEDEN, UEINT, &
                                   UFS, UTEMP, const_grav, &
                                   hse_zero_vels, hse_interp_temp, hse_reflect_vels, &
                                   xl_ext, xr_ext, EXT_HSE, EXT_INTERP
    use amrex_fort_module, only : rt => amrex_real
    use interpolate_module
    use eos_module
    use network, only: nspec
    use model_parser_module
    use bl_error_module
    use eos_type_module

    include 'AMReX_bc_types.fi'

    integer adv_l1,adv_h1
    integer bc(1,2,*)
    integer domlo(1), domhi(1)
    real(rt) delta(1), xlo(1), time
    real(rt) adv(adv_l1:adv_h1,NVAR)
    
    integer i,j,q,n
    real(rt) x
    real(rt) pres_above,p_want,pres_zone
    real(rt) temp_zone,X_zone(nspec),dens_zone
    real(rt) :: x_base, dens_base, slope
  
    type (eos_t) :: eos_state
    !Need to fix
    xl_ext=EXT_HSE

    do n = 1,NVAR
    call filcc(adv(:,n),adv_l1,adv_h1, &
    domlo,domhi,delta,xlo,bc(:,:,n))
    enddo
    ! XLO -- HSE with linear density profile, T found via iteration
    ! we do all variables at once here
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
      if (xl_ext == EXT_HSE) then
      call ext_fill(adv,adv_l1,adv_h1, &
                    domlo,domhi,delta,xlo,time,bc)
      else

       x_base = xlo(1) + delta(1)*(dble(domlo(1)-adv_l1) + HALF)

          dens_base = adv(domlo(1),URHO)

          ! density slope
          slope = (adv(domlo(1)+1,URHO) - adv(domlo(1),URHO))/delta(1)
          do j=domlo(1)-1,adv_l1,-1
             x = xlo(1) + delta(1)*(dble(j-adv_l1) + HALF)
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
                  delta(1)*HALF*(dens_zone + adv(j+1,URHO))*const_grav

             ! EOS with HSE pressure + linear density profile yields T, e, ...
             eos_state%rho = dens_zone
             eos_state%T = temp_zone   ! guess
             eos_state%xn(:) = X_zone(:)
             eos_state%p = p_want

             call eos(eos_input_rp, eos_state)

            ! velocity
             if (zero_vels) then
                
                ! zero gradient velocity
                adv(j,UMX) = ZERO
             else
                adv(j,UMX) = min(ZERO,dens_zone*(adv(domlo(1),UMX)/adv(domlo(1),URHO)))
             endif

             adv(j,URHO) = dens_zone
             adv(j,UEINT) = dens_zone*eos_state%e
             adv(j,UEDEN) = dens_zone*eos_state%e + &
                  HALF*(adv(j,UMX)**2.0_rt)/dens_zone
             adv(j,UTEMP) = eos_state%T
             adv(j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
             
          end do
      end if
    end if

  
  
    ! XHI
       if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
        if (xr_ext == EXT_HSE) then
        call bl_error("ERROR: HSE boundaries not implemented for +X")
        elseif (xr_ext == EXT_INTERP) then
          call ext_fill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)
        end if
          do j=domhi(1)+1,adv_h1
             x = xlo(1) + delta(1)*(dble(j-adv_l1) + HALF)

             !Need to check
             adv(j,:) = adv(domhi(1),:)
             dens_zone=adv(domhi(1),URHO)
             temp_zone=adv(domhi(1),UTEMP)
                   do q = 1, nspec
                      X_zone(q) = interpolate(x,npts_model,model_r, &
                                             model_state(:,ispec_model-1+q))
                   enddo

                   ! extrap normal momentum
                   adv(j,UMX) = min(ZERO,adv(domhi(1),UMX))
                   
                   ! zero transverse momentum

                   eos_state%rho = dens_zone
                   eos_state%T = temp_zone
                   eos_state%xn(:) = X_zone

                   
                   call eos(eos_input_rt, eos_state)
                   
                   adv(j,URHO) = dens_zone
                   adv(j,UEINT) = dens_zone*eos_state%e
                   adv(j,UEDEN) = dens_zone*eos_state%e + &
                        HALF*(adv(j,UMX)**2.0_rt)/dens_zone
                   adv(j,UTEMP) = temp_zone
                   adv(j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
               
          end do
    end if
  end subroutine ca_hypfill

  
  subroutine ca_denfill(adv,adv_l1,adv_h1, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")
    
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, &
         UFS, UTEMP, const_grav
    use bl_error_module
    use interpolate_module
    use model_parser_module

    implicit none
    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_h1
    integer bc(1,2,*)
    integer domlo(1), domhi(1)
    real(rt) delta(1), xlo(1), time
    real(rt) adv(adv_l1:adv_h1,NVAR)

    integer i,j,q,n
    real(rt) x
    real(rt) :: x_base, dens_base, slope
    real(rt) TOL

    call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

   end subroutine ca_denfill

  
  subroutine ca_gravxfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravxfill")

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_h1
    integer :: bc(1,2,*)
    integer :: domlo(1), domhi(1)
    real(rt) delta(1), xlo(1), time
    real(rt) grav(grav_l1:grav_h1)
    integer :: i, j

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)


  end subroutine ca_gravxfill


  subroutine ca_gravyfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravyfill")

    use probdata_module
    use meth_params_module, only: const_grav

    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_h1
    integer :: bc(1,2,*)
    integer :: domlo(1), domhi(1)
    real(rt) delta(1), xlo(1), time
    real(rt) grav(grav_l1:grav_h1)
    integer :: i, j

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)


  end subroutine ca_gravyfill


  subroutine ca_gravzfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravzfill")

    use probdata_module
    use meth_params_module, only: const_grav

    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_h1
    integer :: bc(1,2,*)
    integer :: domlo(1), domhi(1)
    real(rt) delta(1), xlo(1), time
    real(rt) grav(grav_l1:grav_h1)
    integer :: i, j

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravzfill


  subroutine ca_reactfill(react,react_l1, &
                          react_h1,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_reactfill")

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: react_l1,react_h1
    integer :: bc(1,2,*)
    integer :: domlo(1), domhi(1)
    real(rt) delta(1), xlo(1), time
    real(rt) react(react_l1:react_h1)

    call filcc(react,react_l1,react_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_reactfill


  subroutine ca_radfill(rad,rad_l1, &
                        rad_h1,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_radfill")


    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: rad_l1,rad_h1
    integer :: bc(1,2,*)
    integer :: domlo(1), domhi(1)
    real(rt) delta(1), xlo(1), time
    real(rt) rad(rad_l1:rad_h1)

    integer :: j

    call filcc(rad,rad_l1,rad_h1,domlo,domhi,delta,xlo,bc)

    if ( bc(1,1,1).eq.EXT_DIR .and. rad_l1.lt.domlo(1)) then
       do j=rad_l1, domlo(1)-1

          rad(j) = rad(domlo(1))
       enddo
    endif

    if ( bc(1,2,1).eq.EXT_DIR .and. rad_h1.gt.domhi(1)) then
       do j = domhi(1)+1, rad_h1
          rad(j) = rad(domhi(1))
       end do
    end if


  end subroutine ca_radfill

  
  subroutine ca_phigravfill(phi,phi_l1, &
                            phi_h1,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_phigravfill")


    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt) :: delta(1), xlo(1), time
    real(rt) :: phi(phi_l1:phi_h1)

    call filcc(phi,phi_l1,phi_h1, &
         domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill

end module bc_fill_module
