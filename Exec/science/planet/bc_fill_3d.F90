!AUG10
module bc_fill_module
  use bc_ext_fill_module
  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none
  include 'AMReX_bc_types.fi'
  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")
 
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY,UMZ, UEDEN, UEINT, &
                                   UFS, UTEMP, const_grav, &
                                   hse_zero_vels, hse_interp_temp, hse_reflect_vels, &
                                   xl_ext,xr_ext,yl_ext, yr_ext,zl_ext,zr_ext, EXT_HSE, EXT_INTERP
    use amrex_fort_module, only : rt => amrex_real
    use interpolate_module
    use eos_module
    use network, only: nspec
    use model_parser_module
    use amrex_error_module
    use eos_type_module

    include 'AMReX_bc_types.fi'

    integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer bc(3,2,*)
    integer domlo(3), domhi(3)
    real(rt) delta(3), xlo(3), time
    real(rt) adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)
    
    integer i,j,k,q,n
    real(rt) z!x,y,z
    real(rt) pres_above,p_want,pres_zone
    real(rt) temp_zone,X_zone(nspec),dens_zone
    real(rt) :: z_base,dens_base, slope
  
    type (eos_t) :: eos_state
    !need to fix
    zl_ext=EXT_HSE

    do n = 1,NVAR
    call filcc(adv(:,:,:,n),adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,&
    domlo,domhi,delta,xlo,bc(:,:,n))
    enddo
    do n = 1, NVAR
         
       ! XLO
       if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then

          ! we are periodic in -x -- we should never get here
          call amrex_error("ERROR: invalid BC in Prob_3d.f90")
        
       end if

       ! YLO
       if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then

          ! we are periodic in -y -- we should never get here
          call amrex_error("ERROR: invalid BC in Prob_3d.f90")
        
       end if


        ! XHI
       if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
          
          ! we are periodic in x -- we should never get here
          call amrex_error("ERROR: invalid BC in Prob_3d.f90")
        
       end if

         ! YHI
       if ( bc(2,2,n).eq.EXT_DIR .and. adv_h2.lt.domhi(2)) then

          ! we are periodic in +y -- we should never get here
          call amrex_error("ERROR: invalid BC in Prob_3d.f90")
        
       end if



        if (xr_ext == EXT_HSE .or. xl_ext == EXT_HSE .or. xr_ext == EXT_INTERP .or. &
        xl_ext == EXT_INTERP ) then
        call amrex_error("ERROR: HSE boundaries not implemented for +,- X")
        end if

        if (yr_ext == EXT_HSE .or. yl_ext == EXT_HSE .or. yr_ext == EXT_INTERP .or. &
        yl_ext == EXT_INTERP ) then
        call amrex_error("ERROR: HSE boundaries not implemented for +,- Y")
        end if


    enddo


    if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
      if (zl_ext == EXT_HSE) then
        call ext_fill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                    domlo,domhi,delta,xlo,time,bc)
      else
       z_base = xlo(3) + delta(3)*(dble(domlo(3)-adv_l3) + HALF)

       do i=adv_l1,adv_h1
          do j=adv_l2,adv_h2
            dens_base = adv(i,j,domlo(3),URHO)

          ! density slope
            slope = (adv(i,j,domlo(3)+1,URHO) - adv(i,j,domlo(3),URHO))/delta(3)
          ! this do loop counts backwards since we want to work downward
            do k=domlo(3)-1,adv_l3,-1
                z = xlo(3) + delta(3)*(dble(k-adv_l3) + HALF)
             ! zero-gradient catch-all -- this will get the radiation
             ! energy
             adv(i,j,k,:) = adv(i,j,k+1,:)

             ! HSE integration to get temperature, pressure
                    
             ! density is linear from the last two zones
             dens_zone = dens_base + slope*(z - z_base)

             ! temperature guess and species held constant in BCs
             temp_zone = adv(i,j,k+1,UTEMP)

             X_zone(:) = adv(i,j,k+1,UFS:UFS-1+nspec)/adv(i,j,k+1,URHO)
             
             ! get pressure in zone above
             eos_state%rho = adv(i,j,k+1,URHO)
             eos_state%T = adv(i,j,k+1,UTEMP)
             eos_state%xn(:) = adv(i,j,k+1,UFS:UFS-1+nspec)/adv(i,j,k+1,URHO)
             
             call eos(eos_input_rt, eos_state)
             
             pres_above = eos_state%p

             ! pressure needed from HSE
             p_want = pres_above - &
                  delta(3)*HALF*(dens_zone + adv(i,j,k+1,URHO))*const_grav

             ! EOS with HSE pressure + linear density profile yields T, e, ...
             eos_state%rho = dens_zone
             eos_state%T = temp_zone   ! guess
             eos_state%xn(:) = X_zone(:)
             eos_state%p = p_want
 

            
             call eos(eos_input_rp, eos_state)

             ! velocity
             if (zero_vels) then
                
                ! zero normal momentum causes pi waves to pass through
                adv(i,j,k,UMZ) = ZERO
                
                ! zero transverse momentum
                adv(i,j,k,UMX) = ZERO
                adv(i,j,k,UMY) = ZERO


             else
                
                ! zero gradient velocity
                adv(i,j,k,UMZ) = min(ZERO,dens_zone*(adv(i,j,domlo(3),UMZ)/adv(i,j,domlo(3),URHO)))
             endif
           
             adv(i,j,k,URHO) = dens_zone
             adv(i,j,k,UEINT) = dens_zone*eos_state%e
             adv(i,j,k,UEDEN) = dens_zone*eos_state%e + &
                  HALF*(adv(i,j,k,UMX)**2+adv(i,j,k,UMY)**2+adv(i,j,k,UMZ)**2)/dens_zone
             adv(i,j,k,UTEMP) = eos_state%T
             adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
             
            end do
          end do
        end do
      end if
    end if

  
  
    ! YHI
        if ( bc(3,2,URHO).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then

        if (yr_ext == EXT_HSE) then
          call amrex_error("ERROR: HSE boundaries not implemented for +Y")
        elseif (yr_ext == EXT_INTERP) then
          call ext_fill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc)
        end if

        if (xr_ext == EXT_HSE) then
          call amrex_error("ERROR: HSE boundaries not implemented for +X")
        elseif (xr_ext == EXT_INTERP) then
          call ext_fill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc)
        end if

        do k=domhi(3)+1,adv_h3
          z = xlo(3) + delta(3)*(dble(k-adv_l3) + HALF)
          adv(adv_l1:adv_h1,adv_l2:adv_h2,k,:) = adv(adv_l1:adv_h1,adv_l2:adv_h2,domhi(3),:)
          do j=adv_l2,adv_h2
             do i=adv_l1,adv_h1
                dens_zone=adv(i,j,domhi(3),URHO)
                temp_zone=adv(i,j,domhi(3),UTEMP)

                do q = 1, nspec
                  X_zone(q) = interpolate(z,npts_model,model_r, &
                                      model_state(:,ispec_model-1+q))
                enddo

                   ! extrap normal momentum
                adv(i,j,k,UMZ) = min(ZERO,adv(i,j,domhi(3),UMZ))
                   
                   ! zero transverse momentum
                adv(i,j,k,UMX) = ZERO
                adv(i,j,k,UMY) = ZERO
                   
                eos_state%rho = dens_zone
                eos_state%T = temp_zone
                eos_state%xn(:) = X_zone

                   
                call eos(eos_input_rt, eos_state)
                   
                adv(i,j,k,URHO) = dens_zone
                adv(i,j,k,UEINT) = dens_zone*eos_state%e
                adv(i,j,k,UEDEN) = dens_zone*eos_state%e + &
                        HALF*(adv(i,j,k,UMX)**2+adv(i,j,k,UMY)**2+adv(i,j,k,UMZ)**2)/dens_zone
                adv(i,j,k,UTEMP) = temp_zone
                adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
               

              end do
            end do
          end do
       end if
     
  end subroutine ca_hypfill

  
  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")
    
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY,UMZ, UEDEN, UEINT, &
         UFS, UTEMP, const_grav
    use amrex_error_module
    use interpolate_module
    use model_parser_module

    implicit none
    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer bc(3,2,*)
    integer domlo(3), domhi(3)
    real(rt) delta(3), xlo(3), time
    real(rt) adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

    integer i,j,k,q,n
    real(rt) z
    real(rt) ::  z_base,dens_base, slope
    real(rt) TOL

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       call amrex_error("We shoundn't be here (xlo denfill)")
    end if

    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       call amrex_error("We shoundn't be here (ylo denfill)")
    end if


    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       call amrex_error("We shoundn't be here (xhi denfill)")
    endif
    !     YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
       call amrex_error("We shoundn't be here (yhi denfill)")
    endif


  end subroutine ca_denfill

  
  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    real(rt) delta(3), xlo(3), time
    real(rt) grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)
    integer :: i, j,k

    call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravxfill


  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use meth_params_module, only: const_grav

    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    real(rt) delta(3), xlo(3), time
    real(rt) grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)
    integer :: i, j,k

    call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc)


  end subroutine ca_gravyfill


  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use meth_params_module, only: const_grav

    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    real(rt) delta(3), xlo(3), time
    real(rt) grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)
    integer :: i, j,k

    call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc)


  end subroutine ca_gravzfill


  subroutine ca_reactfill(react,react_l1,react_l2, react_l3,&
                          react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: react_l1,react_l2,react_l3,react_h1,react_h2,react_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    real(rt) delta(3), xlo(3), time
    real(rt) react(react_l1:react_h1,react_l2:react_h2,react_l3:react_h3)

    call filcc(react,react_l1,react_l2,react_l3,react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,bc)

  end subroutine ca_reactfill


  subroutine ca_radfill(rad,rad_l1,rad_l2,rad_l3, &
                        rad_h1,rad_h2,rad_h3,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: rad_l1,rad_l2,rad_l3,rad_h1,rad_h2,rad_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    real(rt) delta(3), xlo(3), time
    real(rt) rad(rad_l1:rad_h1,rad_l2:rad_h2,rad_l3:rad_h3)

    integer :: j

    call filcc(rad,rad_l1,rad_l2,rad_l3,rad_h1,rad_h2,rad_h3,domlo,domhi,delta,xlo,bc)

    ! we are inflow at the lower boundary, so we need to take the appropriate
    ! action for the radiation here (during the hydro step)

    ! this do loop counts backwards since we want to work downward
    !YLO
    if ( bc(3,1,1).eq.EXT_DIR .and. rad_l3.lt.domlo(3)) then
       do j=rad_l3,domlo(3)-1

          ! zero-gradient catch-all -- this will get the radiation
          ! energy
          rad(rad_l1:rad_h1,rad_l2:rad_h2,j) = rad(rad_l1:rad_h1,rad_l2:rad_h2,domlo(3))
       enddo
    endif

    !YHI
    if ( bc(3,2,1).eq.EXT_DIR .and. rad_h3.gt.domhi(3)) then
       do j=domhi(3)+1,rad_h3

          ! zero-gradient catch-all -- this will get the radiation
          ! energy
          rad(rad_l1:rad_h1,rad_l2:rad_h2,j) = rad(rad_l1:rad_h1,rad_l2:rad_h2,domhi(3))
       enddo
    endif


  end subroutine ca_radfill


  subroutine ca_phigravfill(phi,phi_l1,phi_l2, phi_l3,&
                            phi_h1,phi_h2,phi_h3,domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
    integer          :: bc(3,2,*)
    integer          :: domlo(3), domhi(3)
    real(rt) :: delta(3), xlo(3), time
    real(rt) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

    call filcc(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
         domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill

end module bc_fill_module
