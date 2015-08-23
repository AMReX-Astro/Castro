subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                      domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  use prob_params_module, only: center
  use meth_params_module, only: NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
  use interpolate_module
  use eos_module
  use eos_type_module
  use network, only: nspec

  implicit none
  include 'bc_types.fi'

  integer          :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

  integer i,j,k,n
  double precision z
  double precision X_in(nspec)
  double precision H

  integer npts_1d
  double precision, allocatable :: pressure(:), density(:), temp(:), eint(:)
  double precision const

  type (eos_t) :: eos_state

  ! compute background state

  ! first make a 1D initial model for the entire domain
  npts_1d = (2.d0*center(3)+1.d-8) / delta(3)

  allocate(pressure(-5:npts_1d+4))
  allocate(density (-5:npts_1d+4))
  allocate(temp    (-5:npts_1d+4))
  allocate(eint    (-5:npts_1d+4))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  X_in(1) = 1.d0

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(gravity)

  do k=0,npts_1d+4

     ! initial guess
     temp(k) = 1000.d0

     if (do_isentropic) then
        z = dble(k) * delta(3)
        density(k) = dens_base*(gravity*dens_base*(gamma_const - 1.0)*z/ &
             (gamma_const*pres_base) + 1.d0)**(1.d0/(gamma_const - 1.d0))
     else
        z = (dble(k)+0.5d0) * delta(3)
        density(k) = dens_base * exp(-z/H)
     end if

     if (k .gt. 0) then
        pressure(k) = pressure(k-1) - &
             delta(3) * 0.5d0 * (density(k)+density(k-1)) * abs(gravity)
     end if

     eos_state%rho = density(k)
     eos_state%T = temp(k)
     eos_state%p = pressure(k)
     eos_state%xn(:) = X_in(:)

     call eos(eos_input_rp, eos_state)

     eint(k) = eos_state%e
     temp(k) = eos_state%T

  end do

  do k=-1,-5,-1

     ! initial guess
     temp(k) = 1000.d0

     if (do_isentropic) then
        z = dble(k) * delta(3)
        density(k) = dens_base*(gravity*dens_base*(gamma_const - 1.0)*z/ &
             (gamma_const*pres_base) + 1.d0)**(1.d0/(gamma_const - 1.d0))
     else
        z = (dble(k)+0.5d0) * delta(3)
        density(k) = dens_base * exp(-z/H)
     end if

     pressure(k) = pressure(k+1) + &
          delta(3) * 0.5d0 * (density(k)+density(k+1)) * abs(gravity)

     eos_state%rho = density(k)
     eos_state%T = temp(k)
     eos_state%p = pressure(k)
     eos_state%xn(:) = X_in(:)

     call eos(eos_input_rp, eos_state)

     eint(k) = eos_state%e
     temp(k) = eos_state%T

  end do

  ! end compute background state

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
                adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

  !  XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     call bl_error("should not have HSE at -x")
  end if

  !  XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     call bl_error("should not have HSE at +x")
  end if

  !  YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     call bl_error("should not have HSE at -x")
  end if

  ! YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     call bl_error("should not have HSE at +x")
  end if


  ! ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     ! this do loop counts backwards since we want to work downward
     do k=domlo(3)-1,adv_l3,-1
        z = xlo(3) + delta(3)*(float(k-adv_l3) + 0.5d0)

        do j = adv_l2, adv_h2
           do i = adv_l1, adv_h1

              ! zero transverse momentum
              adv(i,j,k,UMX:UMY) = 0.d0

              if (boundary_type .eq. 1) then
                 ! extrapolate normal momentum
                 ! enforces pi=0 at boundary
                 adv(i,j,k,UMZ) = adv(i,j,domlo(3),UMZ)
              else
                 ! zero normal momentum
                 ! permits pi to pass through boundary
                 adv(i,j,k,UMZ) = 0.d0
              end if

              adv(i,j,k,URHO) = density(k)
              adv(i,j,k,UFS) = adv(i,j,k,URHO)
              adv(i,j,k,UEINT) = eint(k)*adv(i,j,k,URHO)
              adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) &
                   + 0.5d0*(adv(i,j,k,UMX)**2 + &
                            adv(i,j,k,UMY)**2 + &
                            adv(i,j,k,UMZ)**2)/adv(i,j,k,URHO)
              adv(i,j,k,UTEMP) = temp(k)
           enddo
        end do

     end do
  end if

  ! ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     do k=domhi(3)+1,adv_h3
        z = xlo(3) + delta(3)*(float(k-adv_l3) + 0.5d0)
        do j = adv_l2, adv_h2
           do i = adv_l1, adv_h1

              ! zero transverse momentum
              adv(i,j,k,UMX:UMY) = 0.d0

              if (boundary_type .eq. 1) then
                 ! extrapolate normal momentum
                 ! enforces pi=0 at boundary
                 adv(i,j,k,UMZ) = adv(i,j,domhi(3),UMZ)
              else
                 ! zero normal momentum
                 ! permits pi to pass through boundary
                 adv(i,j,k,UMZ) = 0.d0
              end if

              adv(i,j,k,URHO) = density(k)
              adv(i,j,k,UFS) = adv(i,j,k,URHO)
              adv(i,j,k,UEINT) = eint(k)*adv(i,j,k,URHO)
              adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) &
                   + 0.5d0*(adv(i,j,k,UMX)**2 + &
                            adv(i,j,k,UMY)**2 + &
                            adv(i,j,k,UMZ)**2)/adv(i,j,k,URHO)
              adv(i,j,k,UTEMP) = temp(k)
           enddo
        end do

     end do
  end if

  deallocate(pressure,density,temp,eint)

end subroutine ca_hypfill

! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                      domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  use interpolate_module
  use eos_module, only: gamma_const

  implicit none
  include 'bc_types.fi'

  integer          :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  integer i,j,k
  double precision z,H

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(gravity)

  !     Note: this function should not be needed, technically, but is provided
  !     to filpatch because there are many times in the algorithm when just
  !     the density is needed.  We try to rig up the filling so that the same
  !     function is called here and in hypfill where all the states are filled.

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
             domlo,domhi,delta,xlo,bc)

  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     call bl_error("should not have HSE at -x")
  end if

  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.lt.domhi(1)) then
     call bl_error("should not have HSE at +x")
  end if

  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     call bl_error("should not have HSE at -y")
  end if

  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h1.lt.domhi(2)) then
     call bl_error("should not have HSE at +y")
  end if

  !     ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     do k=adv_l3,domlo(3)-1
        do j=adv_l2,adv_h2
           do i=adv_l1,adv_h1

              if (do_isentropic) then
                 z = xlo(3) + delta(3)*float(k-adv_l3)
                 adv(i,j,k) = dens_base*(gravity*dens_base*(gamma_const - 1.0)*z/ &
                      (gamma_const*pres_base) + 1.d0)**(1.d0/(gamma_const - 1.d0))
              else
                 z = xlo(3) + delta(3)*(float(k-adv_l3) + 0.5d0)
                 adv(i,j,k) = dens_base * exp(-z/H)
              end if
           enddo
        end do

     end do
  end if

  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     do k=domhi(3)+1,adv_h3
        do j=adv_l2,adv_h2
           do i=adv_l1,adv_h1

              if (do_isentropic) then
                 z = xlo(3) + delta(3)*float(k-adv_l3)
                 adv(i,j,k) = dens_base*(gravity*dens_base*(gamma_const - 1.0)*z/ &
                      (gamma_const*pres_base) + 1.d0)**(1.d0/(gamma_const - 1.d0))
              else
                 z = xlo(3) + delta(3)*(float(k-adv_l3) + 0.5d0)
                 adv(i,j,k) = dens_base * exp(-z/H)
              end if

           enddo
        end do
     end do
  end if

end subroutine ca_denfill

! ::: -----------------------------------------------------------


subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                       domlo,domhi,delta,xlo,time,bc)

  implicit none
  include 'bc_types.fi'

  integer          :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_gravxfill


subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                       domlo,domhi,delta,xlo,time,bc)

  implicit none
  include 'bc_types.fi'

  integer          :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_gravyfill


subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                       domlo,domhi,delta,xlo,time,bc)

  implicit none
  include 'bc_types.fi'

  integer          :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_gravzfill



subroutine ca_reactfill(react,react_l1,react_l2,react_l3, &
                        react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: react_l1,react_l2,react_l3,react_h1,react_h2,react_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: react(react_l1:react_h1,react_l2:react_h2,react_l3:react_h3)

  call filcc(react,react_l1,react_l2,react_l3,react_h1,react_h2,react_h3, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_reactfill



subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_l3, &
                          phi_h1,phi_h2,phi_h3,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

  call filcc(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_phigravfill



subroutine ca_radfill(rad,rad_l1,rad_l2,rad_l3, &
                      rad_h1,rad_h2,rad_h3,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: rad_l1,rad_l2,rad_l3,rad_h1,rad_h2,rad_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: rad(rad_l1:rad_h1,rad_l2:rad_h2,rad_l3:rad_h3)

  call filcc(rad,rad_l1,rad_l2,rad_l3,rad_h1,rad_h2,rad_h3, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_radfill
