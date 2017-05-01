subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use bl_error_module
  use probdata_module
  use prob_params_module, only: center

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ pert_factor,dens_base,pres_base,y_pert_center, &
       pert_width,do_isentropic,boundary_type, &
       frac

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults here
  frac = 0.5

  do_isentropic = .false.

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! set center variable in prob_params_module
  center(1) = frac*(problo(1)+probhi(1))
  center(2) = frac*(problo(2)+probhi(2))

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
  use prob_params_module, only: center
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP, const_grav
  use eos_module
  use eos_type_module
  
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), time, delta(2)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)

  integer i,j,npts_1d
  real(rt)         H,z,xn(1),x,y,x1,y1,r1,const
  real(rt)        , allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state
  
  ! first make a 1D initial model for the entire domain
  npts_1d = (2.e0_rt*center(2)+1.e-8_rt) / delta(2)

  allocate(pressure(0:npts_1d-1))
  allocate(density (0:npts_1d-1))
  allocate(temp    (0:npts_1d-1))
  allocate(eint    (0:npts_1d-1))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  xn(1) = 1.e0_rt

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(const_grav)

  do j=0,npts_1d-1

     ! initial guess
     temp(j) = 1000.e0_rt

     if (do_isentropic) then
        z = dble(j) * delta(2)
        density(j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*z/ &
             (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
     else
        z = (dble(j)+0.5e0_rt) * delta(2)
        density(j) = dens_base * exp(-z/H)
     end if

     if (j .gt. 0) then
        pressure(j) = pressure(j-1) - &
             delta(2) * 0.5e0_rt * (density(j)+density(j-1)) * abs(const_grav)
     end if

     eos_state%p = pressure(j)
     eos_state%T = temp(j)
     eos_state%rho = density(j)
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)

     eint(j) = eos_state%e
     temp(j) = eos_state%T

  end do

  
  ! add an isobaric perturbation
  x1 = center(1)
  y1 = y_pert_center

  do j=lo(2),hi(2)
     y = (dble(j)+0.5e0_rt)*delta(2)
     do i=lo(1),hi(1)
        x = (dble(i)+0.5e0_rt)*delta(1)

        r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / pert_width

        state(i,j,UTEMP) = temp(j) * (1.e0_rt + (pert_factor * (1.e0_rt + tanh(2.e0_rt-r1))))
        state(i,j,UFS) = 1.e0_rt

        eos_state%T = state(i,j,UTEMP)
        eos_state%rho = state(i,j,URHO)
        eos_state%p = pressure(j)
        eos_state%xn(:) = xn(:)

        call eos(eos_input_tp, eos_state)

        state(i,j,URHO) = eos_state%rho
        state(i,j,UEINT) = eos_state%e

        ! make state conservative
        state(i,j,UFS) = state(i,j,UFS)*state(i,j,URHO)
        state(i,j,UEINT) = state(i,j,UEINT)*state(i,j,URHO)

        ! assumes ke=0
        state(i,j,UEDEN) = state(i,j,UEINT)

        state(i,j,UMX:UMY) = 0.e0_rt

     end do
  end do


  deallocate(pressure,density,temp,eint)

end subroutine ca_initdata

