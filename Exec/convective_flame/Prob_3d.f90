subroutine PROBINIT (init,name,namlen,problo,probhi)

  use bl_error_module
  use probdata_module
  use prob_params_module, only: center

  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(3), probhi(3)

  integer :: untin,i

  namelist /fortin/ pert_factor,dens_base,pres_base,y_pert_center, &
       pert_width,gravity,do_isentropic,boundary_type, &
       frac

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults
  frac = 0.5

  do_isentropic = .false.

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center(1) = (problo(1)+probhi(1))/2.0d0
  center(2) = (problo(2)+probhi(2))/2.0d0
  center(3) = (problo(3)+probhi(3))/2.0d0

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
                       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                       delta,xlo,xhi)
  use probdata_module
  use prob_params_module, only: center
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
  use eos_module
  use eos_type_module
  use network, only: nspec
  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision xlo(3), xhi(3), time, delta(3)
  double precision state(state_l1:state_h1, &
                         state_l2:state_h2, &
                         state_l3:state_h3,NVAR)

  integer :: i,j,k,npts_1d
  double precision :: H,z,xn(nspec),x,y,x1,y1,z1,r1,const
  double precision, allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state
  
  ! first make a 1D initial model for the entire domain
  npts_1d = (2.d0*center(3)+1.d-8) / delta(3)

  allocate(pressure(0:npts_1d-1))
  allocate(density (0:npts_1d-1))
  allocate(temp    (0:npts_1d-1))
  allocate(eint    (0:npts_1d-1))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  xn(:) = 0.0d0
  xn(1) = 1.d0

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(gravity)

  do k=0,npts_1d-1

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

     eos_state%p = pressure(k)
     eos_state%T = temp(k)
     eos_state%rho = density(k)
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)

     eint(k) = eos_state%e
     temp(k) = eos_state%T

  end do

  
  ! add an isobaric perturbation
  x1 = center(1)
  y1 = center(2)
  z1 = y_pert_center

  do k=lo(3),hi(3)
     z = (dble(k)+0.5d0)*delta(3)

     do j=lo(2),hi(2)
        y = (dble(j)+0.5d0)*delta(2)

        do i=lo(1),hi(1)
           x = (dble(i)+0.5d0)*delta(1)

           r1 = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2) / pert_width

           state(i,j,k,UTEMP) = temp(k) * (1.d0 + (pert_factor * (1.d0 + tanh(2.d0-r1))))
           state(i,j,k,UFS:UFS-1+nspec) = xn(:)

           eos_state%T = state(i,j,k,UTEMP)
           eos_state%rho = state(i,j,k,URHO)
           eos_state%p = pressure(k)
           eos_state%xn(:) = xn(:)

           call eos(eos_input_tp, eos_state)

           state(i,j,k,URHO) = eos_state%rho
           state(i,j,k,UEINT) = eos_state%e

           ! make state conservative
           state(i,j,k,UFS) = state(i,j,k,UFS)*state(i,j,k,URHO)
           state(i,j,k,UEINT) = state(i,j,k,UEINT)*state(i,j,k,URHO)

           ! assumes ke=0
           state(i,j,k,UEDEN) = state(i,j,k,UEINT)

           state(i,j,k,UMX:UMZ) = 0.d0
        enddo
     end do
  end do


  deallocate(pressure,density,temp,eint)

end subroutine ca_initdata

