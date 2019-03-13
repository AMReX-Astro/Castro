module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_hypfill")

    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
    use interpolate_module
    use eos_module
    use network, only: nspec
    use amrex_constants_module, only : ZERO, HALF

    use model_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    include 'AMReX_bc_types.fi'
    
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer i,j,n

    real(rt)        , allocatable :: r_model(:), rho_model(:), T_model(:), &
         e_model(:), p_model(:)

    integer :: lo_model, hi_model

    ! we'll generate the initial model at the needed resolution
    call get_model_size(ymin, ymax, delta(2), lo_model, hi_model)

    allocate(  r_model(lo_model:hi_model))
    allocate(rho_model(lo_model:hi_model))
    allocate(  T_model(lo_model:hi_model))
    allocate(  e_model(lo_model:hi_model))
    allocate(  p_model(lo_model:hi_model))

    call get_model(ymin, ymax, delta(2), &
         pres_base, dens_base, do_isentropic, &
         xn_model, &
         r_model, rho_model, T_model, e_model, p_model, &
         lo_model, hi_model)


    do n=1,NVAR
       call filcc(adv(adv_l1,adv_l2,n),adv_l1,adv_l2,adv_h1,adv_h2, &
            domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    !        XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then

       do j=adv_l2,adv_h2
          do i=domlo(1)-1,adv_l1,-1

             ! zero transverse momentum
             adv(i,j,UMY) = ZERO
             adv(i,j,UMZ) = ZERO

             if (boundary_type .eq. 1) then
                ! extrapolate normal momentum
                ! enforces pi=0 at boundary
                adv(i,j,UMX) = adv(domlo(1),j,UMX)
             else
                ! zero normal momentum
                ! permits pi to pass through boundary
                adv(i,j,UMX) = ZERO
             end if

             adv(i,j,URHO) = rho_model(j)
             adv(i,j,UFS:UFS-1+nspec) = adv(i,j,URHO)*xn_model(:)
             adv(i,j,UEINT) = e_model(j)*adv(i,j,URHO)
             adv(i,j,UEDEN) = adv(i,j,UEINT) + HALF*sum(adv(i,j,UMX:UMZ)**2)/adv(i,j,URHO)
             adv(i,j,UTEMP) = T_model(j)

          end do
       end do

    end if

    !        XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then

       do j=adv_l2,adv_h2
          do i=domhi(1)+1,adv_h1

             ! zero transverse momentum
             adv(i,j,UMY) = ZERO
             adv(i,j,UMZ) = ZERO

             if (boundary_type .eq. 1) then
                ! extrapolate normal momentum
                ! enforces pi=0 at boundary
                adv(i,j,UMX) = adv(domhi(1),j,UMX)
             else
                ! zero normal momentum
                ! permits pi to pass through boundary
                adv(i,j,UMX) = ZERO
             end if

             adv(i,j,URHO) = rho_model(j)
             adv(i,j,UFS:UFS-1+nspec) = adv(i,j,URHO)*xn_model(:)
             adv(i,j,UEINT) = e_model(j)*adv(i,j,URHO)
             adv(i,j,UEDEN) = adv(i,j,UEINT) + HALF*sum(adv(i,j,UMX:UMZ)**2)/adv(i,j,URHO)
             adv(i,j,UTEMP) = T_model(j)

          end do
       end do

    end if


    !        YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       ! this do loop counts backwards since we want to work downward
       do j=domlo(2)-1,adv_l2,-1
          do i=adv_l1,adv_h1

             ! zero transverse momentum
             adv(i,j,UMX) = ZERO
             adv(i,j,UMZ) = ZERO

             if (boundary_type .eq. 1) then
                ! extrapolate normal momentum
                ! enforces pi=0 at boundary
                adv(i,j,UMY) = adv(i,domlo(2),UMY)
             else
                ! zero normal momentum
                ! permits pi to pass through boundary
                adv(i,j,UMY) = ZERO
             end if

             adv(i,j,URHO) = rho_model(j)
             adv(i,j,UFS:UFS-1+nspec) = adv(i,j,URHO)*xn_model(:)
             adv(i,j,UEINT) = e_model(j)*adv(i,j,URHO)
             adv(i,j,UEDEN) = adv(i,j,UEINT) + HALF*sum(adv(i,j,UMX:UMZ)**2)/adv(i,j,URHO)
             adv(i,j,UTEMP) = T_model(j)

          end do
       end do
    end if

    !        YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
       do j=domhi(2)+1,adv_h2
          do i=adv_l1,adv_h1

             ! zero transverse momentum
             adv(i,j,UMX) = ZERO
             adv(i,j,UMZ) = ZERO

             if (boundary_type .eq. 1) then
                ! extrapolate normal momentum
                ! enforces pi=0 at boundary
                adv(i,j,UMY) = adv(i,domhi(2),UMY)
             else
                ! zero normal momentum
                ! permits pi to pass through boundary
                adv(i,j,UMY) = ZERO
             end if

             adv(i,j,URHO) = rho_model(j)
             adv(i,j,UFS:UFS-1+nspec) = adv(i,j,URHO)*xn_model(:)
             adv(i,j,UEINT) = e_model(j)*adv(i,j,URHO)
             adv(i,j,UEDEN) = adv(i,j,UEINT) &
                  + HALF*sum(adv(i,j,UMX:UMZ)**2)/adv(i,j,URHO)
             adv(i,j,UTEMP) = T_model(j)

          end do
       end do
    end if

    deallocate(r_model, rho_model, T_model, p_model, e_model)

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    use probdata_module

    use model_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    include 'AMReX_bc_types.fi'
    
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2)

    integer i,j

    real(rt)        , allocatable :: r_model(:), rho_model(:), T_model(:), &
         e_model(:), p_model(:)

    integer :: lo_model, hi_model

    ! we'll generate the initial model at the needed resolution
    call get_model_size(ymin, ymax, delta(2), lo_model, hi_model)

    allocate(  r_model(lo_model:hi_model))
    allocate(rho_model(lo_model:hi_model))
    allocate(  T_model(lo_model:hi_model))
    allocate(  e_model(lo_model:hi_model))
    allocate(  p_model(lo_model:hi_model))

    call get_model(ymin, ymax, delta(2), &
         pres_base, dens_base, do_isentropic, &
         xn_model, &
         r_model, rho_model, T_model, e_model, p_model, &
         lo_model, hi_model)



    !     Note: this function should not be needed, technically, but is provided
    !     to filpatch because there are many times in the algorithm when just
    !     the density is needed.  We try to rig up the filling so that the same
    !     function is called here and in hypfill where all the states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       do j=adv_l2,adv_h2
          do i=adv_l1,domlo(1)-1

             adv(i,j) = rho_model(j)

          end do
       end do
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.lt.domhi(1)) then
       do j=adv_l2,adv_h2
          do i=domhi(1)+1,adv_h1

             adv(i,j) = rho_model(j)

          end do
       end do
    end if

    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       do j=adv_l2,domlo(2)-1
          do i=adv_l1,adv_h1

             adv(i,j) = rho_model(j)

          end do
       end do
    end if

    !     YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
       do j=domhi(2)+1,adv_h2
          do i=adv_l1,adv_h1

             adv(i,j) = rho_model(j)

          end do
       end do
    end if

    deallocate(r_model, rho_model, T_model, p_model, e_model)

  end subroutine ca_denfill



#ifdef GRAVITY
  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravxfill")

    use probdata_module
    use amrex_fort_module, only : rt => amrex_real
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravyfill")

    use probdata_module
    
    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravzfill")

    use probdata_module
    
    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravzfill



  subroutine ca_phigravfill(phi,phi_l1,phi_l2, &
                            phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_h1,phi_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
         domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill
#endif

end module bc_fill_module
