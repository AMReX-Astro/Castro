subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                      domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
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

  integer i,j,n,iter,MAX_ITER,l
  double precision y
  double precision pres_above,pres_below,pres_want,pres_zone
  double precision drho,dpdr,temperature,eint,pressure,species(3),density
  double precision TOL
  logical converged_hse

  type (eos_t) :: eos_state

  MAX_ITER = 100
  TOL = 1.d-8

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,n),adv_l1,adv_l2,adv_h1,adv_h2, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

  do n = 1, NVAR

     !        XLO
     if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then

        do j=adv_l2,adv_h2
           y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5d0)
           do i=domlo(1)-1,adv_l1,-1

              ! set all the variables even though we're testing on URHO
              if (n .eq. URHO) then

                 density = interpolate(y,npts_model,model_r, &
                                       model_state(:,idens_model))
                 temperature = interpolate(y,npts_model,model_r, &
                                           model_state(:,itemp_model))
                 do l = 1, nspec
                    species(l) = interpolate(y,npts_model,model_r, &
                                             model_state(:,ispec_model-1+l))
                 enddo
  

                 ! extrap normal momentum
                 adv(i,j,UMX) = adv(domlo(1),j,UMX)

                 ! zero transverse momentum
                 adv(i,j,UMY) = 0.d0
                 adv(i,j,UMZ) = 0.d0
                 
                 eos_state%rho = density
                 eos_state%T = temperature
                 eos_state%xn(:) = species(:)

                 call eos(eos_input_rt, eos_state)

                 adv(i,j,URHO) = density
                 adv(i,j,UEINT) = density*eos_state%e
                 adv(i,j,UEDEN) = density*eos_state%e + 0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2+adv(i,j,UMZ)**2)/density
                 adv(i,j,UTEMP) = temperature
                 adv(i,j,UFS:UFS+2) = density*species

              end if

           end do
        end do

     end if

     !        XHI
     if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then

        do j=adv_l2,adv_h2
           y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5d0)
           do i=domhi(1)+1,adv_h1

              ! set all the variables even though we're testing on URHO
              if (n .eq. URHO) then

                 density = interpolate(y,npts_model,model_r, &
                                       model_state(:,idens_model))
                 temperature = interpolate(y,npts_model,model_r, &
                                           model_state(:,itemp_model))
                 do l = 1, nspec
                    species(l) = interpolate(y,npts_model,model_r, &
                                             model_state(:,ispec_model-1+l))
                 enddo


                 ! extrap normal momentum
                 adv(i,j,UMX) = adv(domhi(1),j,UMX)

                 ! zero transverse momentum
                 adv(i,j,UMY) = 0.d0
                 adv(i,j,UMZ) = 0.d0
                 
                 eos_state%rho = density
                 eos_state%T = temperature
                 eos_state%xn(:) = species(:)

                 call eos(eos_input_rt, eos_state)

                 adv(i,j,URHO) = density
                 adv(i,j,UEINT) = density*eos_state%e
                 adv(i,j,UEDEN) = density*eos_state%e + 0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2+adv(i,j,UMZ)**2)/density
                 adv(i,j,UTEMP) = temperature
                 adv(i,j,UFS:UFS+2) = density*species

              end if

           end do
        end do

     end if

     !        YLO
     if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
        ! this do loop counts backwards since we want to work downward
        do j=domlo(2)-1,adv_l2,-1
           y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5d0)
           do i=adv_l1,adv_h1

              ! set all the variables even though we're testing on URHO
              if (n .eq. URHO) then

                 density = interpolate(y,npts_model,model_r, &
                                       model_state(:,idens_model))
                 temperature = interpolate(y,npts_model,model_r, &
                                           model_state(:,itemp_model))
                 do l = 1, nspec
                    species(l) = interpolate(y,npts_model,model_r, &
                                             model_state(:,ispec_model-1+l))
                 enddo

                 ! extrap normal momentum causes pi=0 at boundary
                 !                     adv(i,j,UMY) = adv(i,domlo(2),UMY)

                 ! zero normal momentum causes pi waves to pass through
                 adv(i,j,UMY) = 0.d0

                 ! zero transverse momentum
                 adv(i,j,UMX) = 0.d0
                 adv(i,j,UMZ) = 0.d0
                 
                 eos_state%rho = density
                 eos_state%T = temperature
                 eos_state%xn(:) = species(:)

                 call eos(eos_input_rt, eos_state)

                 adv(i,j,URHO) = density
                 adv(i,j,UEINT) = density*eos_state%e
                 adv(i,j,UEDEN) = density*eos_state%e + 0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2+adv(i,j,UMZ)**2)/density
                 adv(i,j,UTEMP) = temperature
                 adv(i,j,UFS:UFS+2) = density*species


              end if

           end do
        end do
     end if

     !        YHI
     if ( bc(2,2,n).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
        do j=domhi(2)+1,adv_h2
           y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5d0)
           do i=adv_l1,adv_h1

              ! set all the variables even though we're testing on URHO
              if (n .eq. URHO) then

                 density = interpolate(y,npts_model,model_r, &
                                       model_state(:,idens_model))
                 temperature = interpolate(y,npts_model,model_r, &
                                           model_state(:,itemp_model))
                 do l = 1, nspec
                    species(l) = interpolate(y,npts_model,model_r, &
                                             model_state(:,ispec_model-1+l))
                 enddo
                 
                 ! extrap normal momentum
                 adv(i,j,UMY) = adv(i,domhi(2),UMY)

                 ! zero transverse momentum
                 adv(i,j,UMX) = 0.d0
                 adv(i,j,UMZ) = 0.d0
                 
                 eos_state%rho = density
                 eos_state%T = temperature
                 eos_state%xn(:) = species(:)

                 call eos(eos_input_rt, eos_state)

                 adv(i,j,URHO) = density
                 adv(i,j,UEINT) = density*eos_state%e
                 adv(i,j,UEDEN) = density*eos_state%e + 0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2+adv(i,j,UMZ)**2)/density
                 adv(i,j,UTEMP) = temperature
                 adv(i,j,UFS:UFS+2) = density*species

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
     do j=adv_l2,adv_h2
        y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5d0)
        do i=adv_l1,domlo(1)-1
           adv(i,j) = interpolate(y,npts_model,model_r, &
                                  model_state(:,idens_model))
        end do
     end do
  end if

  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     do j=adv_l2,adv_h2
        y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5d0)
        do i=domhi(1)+1,adv_h1
           adv(i,j) = interpolate(y,npts_model,model_r, &
                                  model_state(:,idens_model))
        end do
     end do
  end if

  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     do j=adv_l2,domlo(2)-1
        y = xlo(2) + delta(2)*(float(j-adv_l2)+ 0.5d0)
        do i=adv_l1,adv_h1
           adv(i,j) = interpolate(y,npts_model,model_r, &
                                  model_state(:,idens_model))
        end do
     end do
  end if

  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     do j=domhi(2)+1,adv_h2
        y = xlo(2) + delta(2)*(float(j-adv_l2)+ 0.5d0)
        do i=adv_l1,adv_h1
           adv(i,j) = interpolate(y,npts_model,model_r, &
                                  model_state(:,idens_model))
        end do
     end do
  end if

end subroutine ca_denfill

! ::: -----------------------------------------------------------

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

subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
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

end subroutine ca_gravzfill

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
