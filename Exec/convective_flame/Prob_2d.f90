subroutine PROBINIT (init,name,namlen,problo,probhi)

  use bl_types
  use bl_constants_module
  use bl_error_module
  use probdata_module
  use prob_params_module, only: center

  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ pert_factor, dens_base, pres_base, &
       x_pert_loc, pert_width, &
       cutoff_density, &
       pert_width, do_isentropic, boundary_type, &
       zero_vels

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults here
  zero_vels = .false.
  do_isentropic = .false.
  x_pert_loc = ONE
  pert_width = 0.1_dp_t
  pert_factor = ONE

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! set center variable in prob_params_module
  center(1) = HALF*(problo(1)+probhi(1))
  center(2) = HALF*(problo(2)+probhi(2))

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
  use prob_params_module, only: center, problo
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP, const_grav
  use eos_module
  use eos_type_module
  use network, only: nspec
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)

  integer i,j,npts_1d, j_floor
  double precision H,z,xn(nspec),x,y,x1,y1,r1,const
  double precision, allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state
  
  ! first make a 1D initial model for the entire domain
  npts_1d = (2.d0*center(2)+1.d-8) / delta(2)

  allocate(pressure(0:npts_1d-1))
  allocate(density (0:npts_1d-1))
  allocate(temp    (0:npts_1d-1))
  allocate(eint    (0:npts_1d-1))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  xn(:) = ZERO
  xn(1) = ONE

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(const_grav)

  j_floor = -1

  do j = 0, npts_1d-1

     ! initial guess
     temp(j) = 1000.d0

     if (do_isentropic) then
        z = dble(j) * delta(2)
        density(j) = dens_base*(const_grav*dens_base*(gamma_const - ONE)*z/ &
             (gamma_const*pres_base) + ONE)**(ONE/(gamma_const - ONE))
     else
        z = (dble(j)+HALF) * delta(2)
        density(j) = dens_base * exp(-z/H)
     end if

     if (density(j) < cutoff_density) then
        density(j) = cutoff_density
        temp(j) = temp(j-1)
        j_floor = j
        exit
     endif
     
     if (j > 0) then
        pressure(j) = pressure(j-1) - &
             delta(2) * HALF * (density(j)+density(j-1)) * abs(const_grav)
     end if

     if (pressure(j) < ZERO) then
        density(j) = cutoff_density
        temp(j) = temp(j-1)
        j_floor = j
        exit
     endif

     eos_state%p = pressure(j)
     eos_state%T = temp(j)
     eos_state%rho = density(j)
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)

     eint(j) = eos_state%e
     temp(j) = eos_state%T

  end do

  if (j_floor > 0) then
     eos_state%rho = density(j_floor)
     eos_state%T = temp(j_floor)
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)

     density(j_floor:) = eos_state%rho
     temp(j_floor:) = eos_state%T
     pressure(j_floor:) = eos_state%p
     eint(j_floor:) = eos_state%e
  endif

  
  ! add an isobaric perturbation
  do j=lo(2),hi(2)
     y = problo(2) + (dble(j)+HALF)*delta(2)

     do i=lo(1),hi(1)
        x = problo(1) + (dble(i)+HALF)*delta(1)

        if (density(j) > cutoff_density) then
           state(i,j,UTEMP) = temp(j) * (ONE + (pert_factor * (ONE + tanh((x_pert_loc-x)/pert_width)) ) )
        else
           state(i,j,UTEMP) = temp(j)
        endif

        state(i,j,UFS:UFS-1+nspec) = xn(:)

        eos_state%T = state(i,j,UTEMP)
        eos_state%rho = state(i,j,URHO)
        eos_state%p = pressure(j)
        eos_state%xn(:) = xn(:)

        call eos(eos_input_tp, eos_state)

        state(i,j,URHO) = eos_state%rho
        state(i,j,UEINT) = eos_state%e

        ! make state conservative
        state(i,j,UFS:UFS-1+nspec) = state(i,j,UFS:UFS-1+nspec)*state(i,j,URHO)
        state(i,j,UEINT) = state(i,j,UEINT)*state(i,j,URHO)

        ! assumes ke=0
        state(i,j,UEDEN) = state(i,j,UEINT)

        state(i,j,UMX:UMZ) = ZERO

     end do
  end do

  deallocate(pressure,density,temp,eint)

end subroutine ca_initdata

