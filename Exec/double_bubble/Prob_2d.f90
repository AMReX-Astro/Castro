subroutine PROBINIT (init,name,namlen,problo,probhi)

  use prob_params_module, only: center
  use probdata_module
  use bl_error_module

  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ dens_base, pres_base, &
       pert_factor, y_pert_center, pert_width, &
       do_isentropic, &
       boundary_type, &
       single

  ! build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  do_isentropic = .false.
  single = .false.

  if (namlen .gt. maxlen) then
     call bl_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do


  ! model composition
  xn_model(:) = 0.0d0
  xn_model(1) = 1.0d0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! set local variable defaults
  center(1) = 0.5d0*(problo(1)+probhi(1))
  center(2) = 0.5d0*(problo(2)+probhi(2))

  ymin = problo(2)
  ymax = probhi(2)

  if (single) then
     left_bubble_x_center = problo(1) + 0.5d0*(probhi(1)-problo(1))
  else
     left_bubble_x_center = problo(1) + (probhi(1)-problo(1))/3.d0
     right_bubble_x_center = problo(1) + 2.d0*(probhi(1)-problo(1))/3.d0
  endif

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
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP, &
       const_grav
  use eos_module
  use eos_type_module
  use network

  use model_module

  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)

  integer i,j,npts_1d
  double precision z,xn(nspec),x,y,x1,y1,x2,y2,r1,r2,const

  double precision, allocatable :: r_model(:), rho_model(:), T_model(:), &
                                   e_model(:), p_model(:)

  integer :: lo_model, hi_model

  type (eos_t) :: eos_state

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


  if (.not. single) then
     x1 = left_bubble_x_center
     y1 = y_pert_center

     x2 = right_bubble_x_center
     y2 = y_pert_center

     do j=lo(2),hi(2)
        y = (dble(j)+0.5d0)*delta(2) + ymin

        do i=lo(1),hi(1)
           x = (dble(i)+0.5d0)*delta(1)

           r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / pert_width
           r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / pert_width

           eos_state % xn(:) = xn_model(:)
           state(i,j,UTEMP) = T_model(j)
           state(i,j,URHO) = rho_model(j)

           ! which bubble are we in? -- we want their rho perturbations to be the
           ! same so they have the same buoyancy
           if (r1 < 2.0d0) then
              state(i,j,URHO) = rho_model(j) * (1.d0 - (pert_factor * (1.d0 + tanh(2.d0-r1))))
              eos_state % xn(:) = 0.0d0
              eos_state % xn(2) = 1.0d0
           endif

           if (r2 < 2.0d0) then
              state(i,j,URHO) = rho_model(j) * (1.d0 - (pert_factor * (1.d0 + tanh(2.d0-r2))))
              eos_state % xn(:) = 0.0d0
              eos_state % xn(3) = 1.0d0
           endif

           eos_state % p = p_model(j)
           eos_state % rho = state(i,j,URHO)

           call eos(eos_input_rp, eos_state)

           state(i,j,UEINT) = eos_state % e
           state(i,j,UTEMP) = eos_state % T


           ! make state conservative
           state(i,j,UFS:UFS-1+nspec) = state(i,j,URHO)*eos_state % xn(:)
           state(i,j,UEINT) = state(i,j,URHO)*state(i,j,UEINT)

           ! assumes ke=0
           state(i,j,UEDEN) = state(i,j,UEINT)

           state(i,j,UMX:UMY) = 0.d0

        end do
     end do

  else

     x1 = left_bubble_x_center
     y1 = y_pert_center

     do j=lo(2),hi(2)
        y = (dble(j)+0.5d0)*delta(2) + ymin

        do i=lo(1),hi(1)
           x = (dble(i)+0.5d0)*delta(1)

           r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / pert_width

           eos_state % xn(:) = xn_model(:)
           state(i,j,UTEMP) = T_model(j)
           state(i,j,URHO) = rho_model(j)

           ! which bubble are we in? -- we want their rho perturbations to be the
           ! same so they have the same buoyancy
           if (r1 < 2.0d0) then
              state(i,j,URHO) = rho_model(j) * (1.d0 - (pert_factor * (1.d0 + tanh(2.d0-r1))))
              eos_state % xn(:) = 0.0d0
              eos_state % xn(2) = 1.0d0
           endif

           eos_state % p = p_model(j)
           eos_state % rho = state(i,j,URHO)

           call eos(eos_input_rp, eos_state)

           state(i,j,UEINT) = eos_state % e
           state(i,j,UTEMP) = eos_state % T


           ! make state conservative
           state(i,j,UFS:UFS-1+nspec) = state(i,j,URHO)*eos_state % xn(:)
           state(i,j,UEINT) = state(i,j,URHO)*state(i,j,UEINT)

           ! assumes ke=0
           state(i,j,UEDEN) = state(i,j,UEINT)

           state(i,j,UMX:UMY) = 0.d0

        end do
     end do

  endif

  deallocate(r_model, rho_model, T_model, p_model, e_model)

end subroutine ca_initdata

! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
  use interpolate_module
  use eos_module
  use network, only: nspec

  use model_module

  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

  integer i,j,n

  double precision, allocatable :: r_model(:), rho_model(:), T_model(:), &
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
           adv(i,j,UMY) = 0.d0

           if (boundary_type .eq. 1) then
              ! extrapolate normal momentum
              ! enforces pi=0 at boundary
              adv(i,j,UMX) = adv(domlo(1),j,UMX)
           else
              ! zero normal momentum
              ! permits pi to pass through boundary
              adv(i,j,UMX) = 0.d0
           end if

           adv(i,j,URHO) = rho_model(j)
           adv(i,j,UFS:UFS-1+nspec) = adv(i,j,URHO)*xn_model(:)
           adv(i,j,UEINT) = e_model(j)*adv(i,j,URHO)
           adv(i,j,UEDEN) = adv(i,j,UEINT) &
                + 0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/adv(i,j,URHO)
           adv(i,j,UTEMP) = T_model(j)

        end do
     end do

  end if

  !        XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then

     do j=adv_l2,adv_h2
        do i=domhi(1)+1,adv_h1

           ! zero transverse momentum
           adv(i,j,UMY) = 0.d0

           if (boundary_type .eq. 1) then
              ! extrapolate normal momentum
              ! enforces pi=0 at boundary
              adv(i,j,UMX) = adv(domhi(1),j,UMX)
           else
              ! zero normal momentum
              ! permits pi to pass through boundary
              adv(i,j,UMX) = 0.d0
           end if

           adv(i,j,URHO) = rho_model(j)
           adv(i,j,UFS:UFS-1+nspec) = adv(i,j,URHO)*xn_model(:)
           adv(i,j,UEINT) = e_model(j)*adv(i,j,URHO)
           adv(i,j,UEDEN) = adv(i,j,UEINT) &
                + 0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/adv(i,j,URHO)
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
           adv(i,j,UMX) = 0.d0

           if (boundary_type .eq. 1) then
              ! extrapolate normal momentum
              ! enforces pi=0 at boundary
              adv(i,j,UMY) = adv(i,domlo(2),UMY)
           else
              ! zero normal momentum
              ! permits pi to pass through boundary
              adv(i,j,UMY) = 0.d0
           end if

           adv(i,j,URHO) = rho_model(j)
           adv(i,j,UFS:UFS-1+nspec) = adv(i,j,URHO)*xn_model(:)
           adv(i,j,UEINT) = e_model(j)*adv(i,j,URHO)
           adv(i,j,UEDEN) = adv(i,j,UEINT) &
                + 0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/adv(i,j,URHO)
           adv(i,j,UTEMP) = T_model(j)

        end do
     end do
  end if

  !        YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     do j=domhi(2)+1,adv_h2
        do i=adv_l1,adv_h1

           ! zero transverse momentum
           adv(i,j,UMX) = 0.d0

           if (boundary_type .eq. 1) then
              ! extrapolate normal momentum
              ! enforces pi=0 at boundary
              adv(i,j,UMY) = adv(i,domhi(2),UMY)
           else
              ! zero normal momentum
              ! permits pi to pass through boundary
              adv(i,j,UMY) = 0.d0
           end if

           adv(i,j,URHO) = rho_model(j)
           adv(i,j,UFS:UFS-1+nspec) = adv(i,j,URHO)*xn_model(:)
           adv(i,j,UEINT) = e_model(j)*adv(i,j,URHO)
           adv(i,j,UEDEN) = adv(i,j,UEINT) &
                + 0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/adv(i,j,URHO)
           adv(i,j,UTEMP) = T_model(j)

        end do
     end do
  end if

  deallocate(r_model, rho_model, T_model, p_model, e_model)

end subroutine ca_hypfill

! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                      domlo,domhi,delta,xlo,time,bc)

  use probdata_module

  use model_module

  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

  integer i,j

  double precision, allocatable :: r_model(:), rho_model(:), T_model(:), &
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
