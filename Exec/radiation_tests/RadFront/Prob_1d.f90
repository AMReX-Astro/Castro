
subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use eos_module
  use eos_type_module, only: eos_t, eos_input_rt
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(1), probhi(1)
  type(eos_t) :: eos_state

  integer untin,i

  namelist /fortin/rho_0, T_0

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !
  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     write(6,*) 'probin file name too long'
     stop
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  !default density and temp in domain
  rho_0 = 1.0e-5
  T_0 = 3.0e2


  !Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  eos_state % rho = rho_0
  eos_state % T   = T_0
  eos_state % xn  = 0.e0_rt
  eos_state % xn(1) = 1.e0_rt

  call eos(eos_input_rt, eos_state)

  rhoe_0 = rho_0 * eos_state % e

  !     domain extrema
  xmin = problo(1)
  xmax = probhi(1)

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
     state,state_l1,state_h1, &
     delta,xlo,xhi)
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
  use network, only : nspec
  use eos_type_module
  use eos_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: level, nscal
  integer :: lo(1), hi(1)
  integer :: state_l1,state_h1
  real(rt)         :: state(state_l1:state_h1,NVAR)
  real(rt)         :: time, delta(1)
  real(rt)         :: xlo(1), xhi(1)

  integer :: i
  real(rt)         :: c_v, eint
  type(eos_t)::eos_state

  do i = lo(1), hi(1)

     state(i,URHO) = rho_0
     state(i,UMX) = 0.e0_rt
     state(i,UEDEN) = rhoe_0
     state(i,UEINT) = rhoe_0

     ! set the composition to be all in the first species
     state(i,UFS:UFS-1+nspec) = 0.e0_rt
     state(i,UFS  ) = state(i,URHO)


     state(i,UTEMP) = T_0


  enddo

end subroutine ca_initdata


! :::
! ::: -----------------------------------------------------------
! :::
subroutine ca_initrad(level,time,lo,hi,nrad, &
     rad_state,rad_state_l1, &
     rad_state_h1, &
     delta,xlo,xhi)

  use probdata_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer :: level, nrad
  integer :: lo(1), hi(1)
  integer :: rad_state_l1,rad_state_h1
  real(rt)         :: xlo(1), xhi(1), time, delta(1)
  real(rt)         ::  rad_state(rad_state_l1:rad_state_h1, nrad)

  ! local variables
  integer :: i

  do i = lo(1), hi(1)
         rad_state(i,:) = 0.0
  enddo

end subroutine ca_initrad
