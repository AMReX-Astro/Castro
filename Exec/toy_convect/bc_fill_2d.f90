! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                      domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP, const_grav
  use interpolate_module
  use eos_module
  use network, only: nspec
  use model_parser_module

  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

  integer i,j,q,n,iter
  double precision y
  double precision pres_above,p_want,pres_zone, A
  double precision drho,dpdr,temp_zone,eint,X_zone(nspec),dens_zone

  integer, parameter :: MAX_ITER = 100 
  double precision, parameter :: TOL = 1.d-8
  logical converged_hse

  type (eos_t) :: eos_state

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,n),adv_l1,adv_l2,adv_h1,adv_h2, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

  do n = 1, NVAR

     !        XLO
     if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then

        ! we are periodic in x -- we should never get here
        call bl_error("ERROR: invalid BC in Prob_2d.f90")

     end if

     !        XHI
     if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then

        ! we are periodic in x -- we should never get here
        call bl_error("ERROR: invalid BC in Prob_2d.f90")

     end if

     !        YLO
     if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then

        ! this do loop counts backwards since we want to work downward
        do j=domlo(2)-1,adv_l2,-1
           y = xlo(2) + delta(2)*(dble(j-adv_l2) + 0.5d0)

           do i=adv_l1,adv_h1

              ! set all the variables even though we're testing on URHO
              if (n .eq. URHO) then

                 if (interp_BC) then

                    dens_zone = interpolate(y,npts_model,model_r, &
                                            model_state(:,idens_model))

                    temp_zone = interpolate(y,npts_model,model_r, &
                                            model_state(:,itemp_model))

                    do q = 1, nspec
                       X_zone(q) = interpolate(y,npts_model,model_r, &
                                               model_state(:,ispec_model-1+q))
                    enddo

                 else

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
                            delta(2)*0.5d0*(dens_zone + adv(i,j+1,URHO))*const_grav

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
                       drho = A/(dpdr + 0.5*delta(2)*const_grav)

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
                 if (zero_vels) then

                    ! zero normal momentum causes pi waves to pass through
                    adv(i,j,UMY) = 0.d0

                    ! zero transverse momentum
                    adv(i,j,UMX) = 0.d0
                 else

                    ! zero gradient velocity
                    adv(i,j,UMX) = dens_zone*(adv(i,domlo(2),UMX)/adv(i,domlo(2),URHO))
                    adv(i,j,UMY) = dens_zone*(adv(i,domlo(2),UMY)/adv(i,domlo(2),URHO))
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
                      0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/dens_zone
                 adv(i,j,UTEMP) = temp_zone
                 adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

              end if

           end do
        end do
     end if

     !        YHI
     if ( bc(2,2,n).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then

        do j=domhi(2)+1,adv_h2
           y = xlo(2) + delta(2)*(dble(j-adv_l2) + 0.5d0)

           do i=adv_l1,adv_h1

              ! set all the variables even though we're testing on URHO
              if (n .eq. URHO) then

                 dens_zone = interpolate(y,npts_model,model_r, &
                                         model_state(:,idens_model))

                 temp_zone = interpolate(y,npts_model,model_r, &
                                         model_state(:,itemp_model))

                 do q = 1, nspec
                    X_zone(q) = interpolate(y,npts_model,model_r, &
                                            model_state(:,ispec_model-1+q))
                 enddo


                 ! extrap normal momentum
                 adv(i,j,UMY) = max(0.d0,adv(i,domhi(2),UMY))

                 ! zero transverse momentum
                 adv(i,j,UMX) = 0.d0

                 eos_state%rho = dens_zone
                 eos_state%T = temp_zone
                 eos_state%xn(:) = X_zone

                 call eos(eos_input_rt, eos_state)

                 pres_zone = eos_state%p
                 eint = eos_state%e

                 adv(i,j,URHO) = dens_zone
                 adv(i,j,UEINT) = dens_zone*eint
                 adv(i,j,UEDEN) = dens_zone*eint + &
                      0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/dens_zone
                 adv(i,j,UTEMP) = temp_zone
                 adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

              end if

           end do
        end do
     end if

  end do

end subroutine ca_hypfill

! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                      domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  use interpolate_module
  use model_parser_module
  use bl_error_module

  implicit none
  include 'bc_types.fi'
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

  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     call bl_error("We shoundn't be here (xlo denfill)")
  end if

  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     call bl_error("We shoundn't be here (xlo denfill)")
  endif


  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     do j=adv_l2,domlo(2)-1
        y = xlo(2) + delta(2)*(dble(j-adv_l2) + 0.5d0)
        do i=adv_l1,adv_h1
           adv(i,j) = interpolate(y,npts_model,model_r,model_state(:,idens_model))
        end do
     end do
  end if

  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     do j=domhi(2)+1,adv_h2
        y = xlo(2) + delta(2)*(dble(j-adv_l2)+ 0.5d0)
        do i=adv_l1,adv_h1
           adv(i,j) = interpolate(y,npts_model,model_r,model_state(:,idens_model))
        end do
     end do
  end if

end subroutine ca_denfill

subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                        domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  implicit none
  include 'bc_types.fi'

  integer :: grav_l1,grav_l2,grav_h1,grav_h2
  integer :: bc(2,2,*)
  integer :: domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)

  call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

end subroutine ca_gravxfill

! ::: -----------------------------------------------------------

subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                        domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  implicit none
  include 'bc_types.fi'

  integer :: grav_l1,grav_l2,grav_h1,grav_h2
  integer :: bc(2,2,*)
  integer :: domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)

  call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

end subroutine ca_gravyfill

! ::: -----------------------------------------------------------

subroutine ca_reactfill(react,react_l1,react_l2, &
                        react_h1,react_h2,domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  implicit none
  include 'bc_types.fi'

  integer :: react_l1,react_l2,react_h1,react_h2
  integer :: bc(2,2,*)
  integer :: domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision react(react_l1:react_h1,react_l2:react_h2)

  call filcc(react,react_l1,react_l2,react_h1,react_h2,domlo,domhi,delta,xlo,bc)

end subroutine ca_reactfill


subroutine ca_phigravfill(phi,phi_l1,phi_l2, &
                          phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: phi_l1,phi_l2,phi_h1,phi_h2
  integer          :: bc(2,2,*)
  integer          :: domlo(2), domhi(2)
  double precision :: delta(2), xlo(2), time
  double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

  call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_phigravfill

