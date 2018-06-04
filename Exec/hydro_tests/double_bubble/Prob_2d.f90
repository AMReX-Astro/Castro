subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use prob_params_module, only: center
  use probdata_module
  use amrex_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(2), probhi(2)

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
     call amrex_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do


  ! model composition
  xn_model(:) = 0.0e0_rt
  xn_model(1) = 1.0e0_rt

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! set local variable defaults
  center(1) = 0.5e0_rt*(problo(1)+probhi(1))
  center(2) = 0.5e0_rt*(problo(2)+probhi(2))

  ymin = problo(2)
  ymax = probhi(2)

  if (single) then
     left_bubble_x_center = problo(1) + 0.5e0_rt*(probhi(1)-problo(1))
  else
     left_bubble_x_center = problo(1) + (probhi(1)-problo(1))/3.e0_rt
     right_bubble_x_center = problo(1) + 2.e0_rt*(probhi(1)-problo(1))/3.e0_rt
  endif

end subroutine amrex_probinit


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

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), time, delta(2)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)

  integer i,j,npts_1d
  real(rt)         z,xn(nspec),x,y,x1,y1,x2,y2,r1,r2,const

  real(rt)        , allocatable :: r_model(:), rho_model(:), T_model(:), &
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
        y = (dble(j)+0.5e0_rt)*delta(2) + ymin

        do i=lo(1),hi(1)
           x = (dble(i)+0.5e0_rt)*delta(1)

           r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / pert_width
           r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / pert_width

           eos_state % xn(:) = xn_model(:)
           state(i,j,UTEMP) = T_model(j)
           state(i,j,URHO) = rho_model(j)

           ! which bubble are we in? -- we want their rho perturbations to be the
           ! same so they have the same buoyancy
           if (r1 < 2.0e0_rt) then
              state(i,j,URHO) = rho_model(j) * (1.e0_rt - (pert_factor * (1.e0_rt + tanh(2.e0_rt-r1))))
              eos_state % xn(:) = 0.0e0_rt
              eos_state % xn(2) = 1.0e0_rt
           endif

           if (r2 < 2.0e0_rt) then
              state(i,j,URHO) = rho_model(j) * (1.e0_rt - (pert_factor * (1.e0_rt + tanh(2.e0_rt-r2))))
              eos_state % xn(:) = 0.0e0_rt
              eos_state % xn(3) = 1.0e0_rt
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

           state(i,j,UMX:UMY) = 0.e0_rt

        end do
     end do

  else

     x1 = left_bubble_x_center
     y1 = y_pert_center

     do j=lo(2),hi(2)
        y = (dble(j)+0.5e0_rt)*delta(2) + ymin

        do i=lo(1),hi(1)
           x = (dble(i)+0.5e0_rt)*delta(1)

           r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / pert_width

           eos_state % xn(:) = xn_model(:)
           state(i,j,UTEMP) = T_model(j)
           state(i,j,URHO) = rho_model(j)

           ! which bubble are we in? -- we want their rho perturbations to be the
           ! same so they have the same buoyancy
           if (r1 < 2.0e0_rt) then
              state(i,j,URHO) = rho_model(j) * (1.e0_rt - (pert_factor * (1.e0_rt + tanh(2.e0_rt-r1))))
              eos_state % xn(:) = 0.0e0_rt
              eos_state % xn(2) = 1.0e0_rt
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

           state(i,j,UMX:UMY) = 0.e0_rt

        end do
     end do

  endif

  deallocate(r_model, rho_model, T_model, p_model, e_model)

end subroutine ca_initdata


