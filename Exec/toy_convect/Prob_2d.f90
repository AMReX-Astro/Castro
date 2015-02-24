subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use model_parser_module
  use bl_error_module

  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  double precision offset
  integer untin,i

  namelist /fortin/ model_name, apply_vel_field, &
       velpert_scale, velpert_amplitude, velpert_height_loc, num_vortices, &
       H_min, cutoff_density, interp_BC, zero_vels

  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  ! Build "probin" filename from C++ land --
  ! the name of file containing fortin namelist.

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do


  ! Namelist defaults
  apply_vel_field = .false.
  velpert_scale = 1.0d2
  velpert_amplitude = 1.0d2
  velpert_height_loc = 6.5d3
  num_vortices = 1
  H_min = 1.d-4
  cutoff_density = 50.d0
  interp_BC = .false.
  zero_vels = .false.

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)


  ! Read initial model
  call read_model_file(model_name)


  do i = 1, npts_model
     print *, i, model_r(i), model_state(i,idens_model)
  enddo

  ! velocity perturbation stuff
  offset = (probhi(1) - problo(1)) / (num_vortices)

  allocate(xloc_vortices(num_vortices))

  do i = 1, num_vortices
     xloc_vortices(i) = (dble(i-1) + 0.5d0) * offset + problo(1)
  enddo

end subroutine PROBINIT


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use probdata_module
  use interpolate_module
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec
  use model_parser_module

  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)

  double precision xdist,ydist,x,y,r,upert(2)
  integer i,j,n,vortex

  double precision temppres(state_l1:state_h1,state_l2:state_h2)

  type (eos_t) :: eos_state

  do j = lo(2), hi(2)
     y = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)

     do i = lo(1), hi(1)

        state(i,j,URHO)  = interpolate(y,npts_model,model_r, &
                                      model_state(:,idens_model))
        state(i,j,UTEMP) = interpolate(y,npts_model,model_r, &
                                       model_state(:,itemp_model))
        do n = 1, nspec
           state(i,j,UFS-1+n) = interpolate(y,npts_model,model_r, &
                                            model_state(:,ispec_model-1+n))
        enddo

     enddo
  enddo

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        eos_state%rho = state(i,j,URHO)
        eos_state%T = state(i,j,UTEMP)
        eos_state%xn(:) = state(i,j,UFS:)

        call eos(eos_input_rt, eos_state)

        state(i,j,UEINT) = eos_state%e
        temppres(i,j) = eos_state%p

     end do
  end do

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        state(i,j,UEDEN) = state(i,j,URHO) * state(i,j,UEINT)
        state(i,j,UEINT) = state(i,j,URHO) * state(i,j,UEINT)

        do n = 1,nspec
           state(i,j,UFS+n-1) = state(i,j,URHO) * state(i,j,UFS+n-1)
        end do

     enddo
  enddo

  ! Initial velocities = 0
  state(:,:,UMX:UMY) = 0.d0

  ! Now add the velocity perturbation
  if (apply_vel_field) then

     do j = lo(2), hi(2)
        y = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)
        ydist = y - velpert_height_loc

        do i = lo(1), hi(1)
           x = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

           upert = 0.d0

           ! loop over each vortex
           do vortex = 1, num_vortices

              xdist = x - xloc_vortices(vortex)

              r = sqrt(xdist**2 + ydist**2)

              upert(1) = upert(1) - (ydist/velpert_scale) * &
                   velpert_amplitude * exp( -r**2/(2.d0*velpert_scale**2)) &
                   * (-1.d0)**vortex

              upert(2) = upert(2) + (xdist/velpert_scale) * &
                   velpert_amplitude * exp(-r**2/(2.d0*velpert_scale**2)) &
                   * (-1.d0)**vortex
           enddo

           state(i,j,UMX) = state(i,j,URHO) * upert(1)
           state(i,j,UMY) = state(i,j,URHO) * upert(2)

        end do
     end do

  endif

end subroutine ca_initdata


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
