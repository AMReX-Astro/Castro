module bc_fill_module

  use amrex_constants_module, only : ZERO, HALF
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")

    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
    use interpolate_module
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_rt
    use network, only: nspec
    use model_parser_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer bc(3,2,*)
    integer domlo(3), domhi(3)
    real(rt)         delta(3), xlo(3), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

    integer i,j,k,n,iter,MAX_ITER,l
    real(rt)         z
    real(rt)         pres_above,pres_below,pres_want,pres_zone
    real(rt)         drho,dpdr,temperature,eint,pressure,species(3),density
    real(rt)         TOL
    logical converged_hse

    type (eos_t) :: eos_state

    MAX_ITER = 100
    TOL = 1.e-8_rt

    do n = 1,NVAR
       call filcc(adv(adv_l1,adv_l2,adv_l3,n),adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                  domlo,domhi,delta,xlo,bc(1,1,n))
    enddo



    !  ZLO
    if ( bc(3,1,n) == EXT_DIR .and. adv_l3 < domlo(3)) then

       ! this do loop counts backwards since we want to work downward
       do k = domlo(3)-1, adv_l3, -1
          z = xlo(3) + delta(3)*(dble(k-adv_l3) + 0.5e0_rt)

          do j = adv_l2, adv_h2
             do i = adv_l1, adv_h1

                ! set all the variables even though we're testing on URHO
                if (n .eq. URHO) then

                   density = interpolate(z,npts_model,model_r, &
                                         model_state(:,idens_model))
                   temperature = interpolate(z,npts_model,model_r, &
                                             model_state(:,itemp_model))
                   do l = 1, nspec
                      species(l) = interpolate(z,npts_model,model_r, &
                                               model_state(:,ispec_model-1+l))
                   enddo

                   ! extrap normal momentum causes pi=0 at boundary
                   !                     adv(i,j,UMY) = adv(i,domlo(2),UMY)

                   ! zero normal momentum causes pi waves to pass through
                   adv(i,j,k,UMZ) = ZERO

                   ! zero transverse momentum
                   adv(i,j,k,UMX) = ZERO
                   adv(i,j,k,UMY) = ZERO

                   eos_state%rho = density
                   eos_state%T = temperature
                   eos_state%xn(:) = species(:)

                   call eos(eos_input_rt, eos_state)

                   adv(i,j,k,URHO) = density
                   adv(i,j,k,UEINT) = density*eos_state%e
                   adv(i,j,k,UEDEN) = density*eos_state%e + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/density
                   adv(i,j,k,UTEMP) = temperature
                   adv(i,j,k,UFS:UFS-1+nspec) = density*species


                endif

             enddo
          enddo
       enddo
    endif

    ! ZHI
    if ( bc(3,2,n) == EXT_DIR .and. adv_h3 > domhi(3)) then

       do k = domhi(3)+1, adv_h3
          z = xlo(3) + delta(3)*(dble(k-adv_l3) + 0.5e0_rt)

          do j = adv_l2, adv_h2
             do i = adv_l1, adv_h1

                ! set all the variables even though we're testing on URHO
                if (n .eq. URHO) then

                   density = interpolate(z,npts_model,model_r, &
                                         model_state(:,idens_model))
                   temperature = interpolate(z,npts_model,model_r, &
                                             model_state(:,itemp_model))
                   do l = 1, nspec
                      species(l) = interpolate(z,npts_model,model_r, &
                                               model_state(:,ispec_model-1+l))
                   enddo

                   ! extrap normal momentum
                   adv(i,j,k,UMZ) = adv(i,j,domhi(3),UMZ)

                   ! zero transverse momentum
                   adv(i,j,k,UMX) = ZERO
                   adv(i,j,k,UMY) = ZERO

                   eos_state%rho = density
                   eos_state%T = temperature
                   eos_state%xn(:) = species(:)

                   call eos(eos_input_rt, eos_state)

                   adv(i,j,k,URHO) = density
                   adv(i,j,k,UEINT) = density*eos_state%e
                   adv(i,j,k,UEDEN) = density*eos_state%e + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/density
                   adv(i,j,k,UTEMP) = temperature
                   adv(i,j,k,UFS:UFS-1+nspec) = density*species(:)

                end if

             enddo
          enddo
       enddo

    endif

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    use probdata_module
    use interpolate_module
    use model_parser_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer bc(3,2,*)
    integer domlo(3), domhi(3)
    real(rt)         delta(3), xlo(3), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

    integer i,j,k
    real(rt) :: z

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,delta,xlo,bc)

    ! ZLO
    if ( bc(3,1,1) == EXT_DIR .and. adv_l3 < domlo(3)) then
       do k = adv_l3, domlo(3)-1
          z = xlo(3) + delta(3)*(dble(k-adv_l3)+ 0.5e0_rt)
          do j = adv_l2, adv_h2
             do i = adv_l1, adv_h1
                adv(i,j,k) = interpolate(z,npts_model,model_r, &
                                         model_state(:,idens_model))
             enddo
          enddo
       enddo
    endif

    ! ZHI
    if ( bc(3,2,1) == EXT_DIR .and. adv_h3 > domhi(3)) then
       do k = domhi(3)+1, adv_h3
          z = xlo(3) + delta(3)*(dble(k-adv_l3)+ 0.5e0_rt)
          do j = adv_l2, adv_h2
             do i = adv_l1, adv_h1
                adv(i,j,k) = interpolate(z,npts_model,model_r, &
                                         model_state(:,idens_model))
             enddo
          enddo
       enddo
    endif

  end subroutine ca_denfill

  ! ::: -----------------------------------------------------------

  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravxfill")

    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    real(rt)         delta(3), xlo(3), time
    real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)
    integer :: bc_temp(3,2)

    bc_temp(:,:) = bc(:,:)

    if (bc(3,1) == EXT_DIR .and. grav_l3 < domlo(3)) then
       bc_temp(3,1) = FOEXTRAP
    endif
    if (bc(3,2) == EXT_DIR .and. grav_h3 > domhi(3)) then
       bc_temp(3,2) = FOEXTRAP
    endif

    call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_gravxfill


  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravyfill")

    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    real(rt)         delta(3), xlo(3), time
    real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)
    integer :: bc_temp(3,2)

    bc_temp(:,:) = bc(:,:)

    if (bc(3,1) == EXT_DIR .and. grav_l3 < domlo(3)) then
       bc_temp(3,1) = FOEXTRAP
    endif
    if (bc(3,2) == EXT_DIR .and. grav_h3 > domhi(3)) then
       bc_temp(3,2) = FOEXTRAP
    endif

    call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravzfill")

    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    real(rt)         delta(3), xlo(3), time
    real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)
    integer :: bc_temp(3,2)

    bc_temp(:,:) = bc(:,:)

    if (bc(3,1) == EXT_DIR .and. grav_l3 < domlo(3)) then
       bc_temp(3,1) = FOEXTRAP
    endif
    if (bc(3,2) == EXT_DIR .and. grav_h3 > domhi(3)) then
       bc_temp(3,2) = FOEXTRAP
    endif

    call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_gravzfill



  subroutine ca_reactfill(react,react_l1,react_l2,react_l3, &
                          react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer :: react_l1,react_l2,react_l3,react_h1,react_h2,react_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    real(rt)         delta(3), xlo(3), time
    real(rt)         react(react_l1:react_h1,react_l2:react_h2,react_l3:react_h3)
    integer :: bc_temp(3,2)

    bc_temp(:,:) = bc(:,:)

    if (bc(3,1) == EXT_DIR .and. react_l3 < domlo(3)) then
       bc_temp(3,1) = FOEXTRAP
    endif
    if (bc(3,2) == EXT_DIR .and. react_h3 > domhi(3)) then
       bc_temp(3,2) = FOEXTRAP
    endif

    call filcc(react,react_l1,react_l2,react_l3,react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_reactfill



  subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_l3, &
                            phi_h1,phi_h2,phi_h3,domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
    integer          :: bc(3,2,*)
    integer          :: domlo(3), domhi(3)
    real(rt)         :: delta(3), xlo(3), time
    real(rt)         :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
    integer :: bc_temp(3,2)

    bc_temp(:,:) = bc(:,:)

    if (bc(3,1) == EXT_DIR .and. phi_l3 < domlo(3)) then
       bc_temp(3,1) = FOEXTRAP
    endif
    if (bc(3,2) == EXT_DIR .and. phi_h3 > domhi(3)) then
       bc_temp(3,2) = FOEXTRAP
    endif

    call filcc(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
               domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_phigravfill

end module bc_fill_module
