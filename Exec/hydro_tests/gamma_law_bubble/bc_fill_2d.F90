module bc_fill_module

  use amrex_fort_module, only : rt => c_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use prob_params_module, only: center
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP, const_grav
    use interpolate_module
    use eos_module
    use eos_type_module
    use network, only: nspec

    use amrex_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'
    
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer i,j,n
    real(rt)         y
    real(rt)         X_in(nspec)
    real(rt)         H

    integer npts_1d
    real(rt)        , allocatable :: pressure(:), density(:), temp(:), eint(:)
    real(rt)         const

    type (eos_t) :: eos_state

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute background state

    ! first make a 1D initial model for the entire domain
    npts_1d = (2.e0_rt*center(2)+1.e-8_rt) / delta(2)

    allocate(pressure(-5:npts_1d+4))
    allocate(density (-5:npts_1d+4))
    allocate(temp    (-5:npts_1d+4))
    allocate(eint    (-5:npts_1d+4))

    const = pres_base/dens_base**gamma_const

    pressure(0) = pres_base
    density(0)  = dens_base

    ! only initialize the first species
    X_in(1) = 1.e0_rt

    ! compute the pressure scale height (for an isothermal, ideal-gas
    ! atmosphere)
    H = pres_base / dens_base / abs(const_grav)

    do j=0,npts_1d+4

       ! initial guess
       temp(j) = 1000.e0_rt

       if (do_isentropic) then
          y = dble(j) * delta(2)
          density(j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*y/ &
               (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
       else
          y = (dble(j)+0.5e0_rt) * delta(2)
          density(j) = dens_base * exp(-y/H)
       end if

       if (j .gt. 0) then
          pressure(j) = pressure(j-1) - &
               delta(2) * 0.5e0_rt * (density(j)+density(j-1)) * abs(const_grav)
       end if

       eos_state%rho = density(j)
       eos_state%T = temp(j)
       eos_state%p = pressure(j)
       eos_state%xn(:) = X_in(:)

       call eos(eos_input_rp, eos_state)

       eint(j) = eos_state%e
       temp(j) = eos_state%T

    end do

    do j=-1,-5,-1

       ! initial guess
       temp(j) = 1000.e0_rt

       if (do_isentropic) then
          y = dble(j) * delta(2)
          density(j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*y/ &
               (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
       else
          y = (dble(j)+0.5e0_rt) * delta(2)
          density(j) = dens_base * exp(-y/H)
       end if

       pressure(j) = pressure(j+1) + &
            delta(2) * 0.5e0_rt * (density(j)+density(j+1)) * abs(const_grav)

       eos_state%rho = density(j)
       eos_state%T = temp(j)
       eos_state%p = pressure(j)
       eos_state%xn(:) = X_in(:)

       call eos(eos_input_rp, eos_state)

       eint(j) = eos_state%e
       temp(j) = eos_state%T

    end do

    ! end compute background state
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,NVAR
       call filcc(adv(adv_l1,adv_l2,n),adv_l1,adv_l2,adv_h1,adv_h2, &
            domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    !        XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then

       do j=adv_l2,adv_h2
          y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5e0_rt)
          do i=domlo(1)-1,adv_l1,-1

             ! zero transverse momentum
             adv(i,j,UMY) = 0.e0_rt

             if (boundary_type .eq. 1) then
                ! extrapolate normal momentum
                ! enforces pi=0 at boundary
                adv(i,j,UMX) = adv(domlo(1),j,UMX)
             else
                ! zero normal momentum
                ! permits pi to pass through boundary
                adv(i,j,UMX) = 0.e0_rt
             end if

             adv(i,j,URHO) = density(j)
             adv(i,j,UFS) = adv(i,j,URHO)
             adv(i,j,UEINT) = eint(j)*adv(i,j,URHO)
             adv(i,j,UEDEN) = adv(i,j,UEINT) &
                  + 0.5e0_rt*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/adv(i,j,URHO)
             adv(i,j,UTEMP) = temp(j)

          end do
       end do

    end if

    !        XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then

       do j=adv_l2,adv_h2
          y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5e0_rt)
          do i=domhi(1)+1,adv_h1

             ! zero transverse momentum
             adv(i,j,UMY) = 0.e0_rt

             if (boundary_type .eq. 1) then
                ! extrapolate normal momentum
                ! enforces pi=0 at boundary
                adv(i,j,UMX) = adv(domhi(1),j,UMX)
             else
                ! zero normal momentum
                ! permits pi to pass through boundary
                adv(i,j,UMX) = 0.e0_rt
             end if

             adv(i,j,URHO) = density(j)
             adv(i,j,UFS) = adv(i,j,URHO)
             adv(i,j,UEINT) = eint(j)*adv(i,j,URHO)
             adv(i,j,UEDEN) = adv(i,j,UEINT) &
                  + 0.5e0_rt*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/adv(i,j,URHO)
             adv(i,j,UTEMP) = temp(j)

          end do
       end do

    end if


    !        YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       ! this do loop counts backwards since we want to work downward
       do j=domlo(2)-1,adv_l2,-1
          y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5e0_rt)
          do i=adv_l1,adv_h1

             ! zero transverse momentum
             adv(i,j,UMX) = 0.e0_rt

             if (boundary_type .eq. 1) then
                ! extrapolate normal momentum
                ! enforces pi=0 at boundary
                adv(i,j,UMY) = adv(i,domlo(2),UMY)
             else
                ! zero normal momentum
                ! permits pi to pass through boundary
                adv(i,j,UMY) = 0.e0_rt
             end if

             adv(i,j,URHO) = density(j)
             adv(i,j,UFS) = adv(i,j,URHO)
             adv(i,j,UEINT) = eint(j)*adv(i,j,URHO)
             adv(i,j,UEDEN) = adv(i,j,UEINT) &
                  + 0.5e0_rt*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/adv(i,j,URHO)
             adv(i,j,UTEMP) = temp(j)

          end do
       end do
    end if

    !        YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
       do j=domhi(2)+1,adv_h2
          y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5e0_rt)
          do i=adv_l1,adv_h1

             ! zero transverse momentum
             adv(i,j,UMX) = 0.e0_rt

             if (boundary_type .eq. 1) then
                ! extrapolate normal momentum
                ! enforces pi=0 at boundary
                adv(i,j,UMY) = adv(i,domhi(2),UMY)
             else
                ! zero normal momentum
                ! permits pi to pass through boundary
                adv(i,j,UMY) = 0.e0_rt
             end if

             adv(i,j,URHO) = density(j)
             adv(i,j,UFS) = adv(i,j,URHO)
             adv(i,j,UEINT) = eint(j)*adv(i,j,URHO)
             adv(i,j,UEDEN) = adv(i,j,UEINT) &
                  + 0.5e0_rt*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/adv(i,j,URHO)
             adv(i,j,UTEMP) = temp(j)

          end do
       end do
    end if

    deallocate(pressure,density,temp,eint)

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use interpolate_module
    use eos_module, only: gamma_const
    use meth_params_module, only : const_grav

    use amrex_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'
    
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2)

    integer i,j
    real(rt)         y,H

    ! compute the pressure scale height (for an isothermal, ideal-gas
    ! atmosphere)
    H = pres_base / dens_base / abs(const_grav)

    !     Note: this function should not be needed, technically, but is provided
    !     to filpatch because there are many times in the algorithm when just
    !     the density is needed.  We try to rig up the filling so that the same
    !     function is called here and in hypfill where all the states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       do j=adv_l2,adv_h2
          do i=adv_l1,domlo(1)-1

             if (do_isentropic) then
                y = xlo(2) + delta(2)*float(j-adv_l2)
                adv(i,j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*y/ &
                     (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
             else
                y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5e0_rt)
                adv(i,j) = dens_base * exp(-y/H)
             end if

          end do
       end do
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.lt.domhi(1)) then
       do j=adv_l2,adv_h2
          do i=domhi(1)+1,adv_h1

             if (do_isentropic) then
                y = xlo(2) + delta(2)*float(j-adv_l2)
                adv(i,j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*y/ &
                     (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
             else
                y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5e0_rt)
                adv(i,j) = dens_base * exp(-y/H)
             end if

          end do
       end do
    end if

    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       do j=adv_l2,domlo(2)-1
          do i=adv_l1,adv_h1

             if (do_isentropic) then
                y = xlo(2) + delta(2)*float(j-adv_l2)
                adv(i,j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*y/ &
                     (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
             else
                y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5e0_rt)
                adv(i,j) = dens_base * exp(-y/H)
             end if

          end do
       end do
    end if

    !     YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
       do j=domhi(2)+1,adv_h2
          do i=adv_l1,adv_h1

             if (do_isentropic) then
                y = xlo(2) + delta(2)*float(j-adv_l2)
                adv(i,j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*y/ &
                     (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
             else
                y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5e0_rt)
                adv(i,j) = dens_base * exp(-y/H)
             end if

          end do
       end do
    end if

  end subroutine ca_denfill



#ifdef GRAVITY
  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    
    use amrex_fort_module, only : rt => c_real
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
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    
    use amrex_fort_module, only : rt => c_real
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
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module

    use amrex_fort_module, only : rt => c_real
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
                            phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    use amrex_fort_module, only : rt => c_real
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
