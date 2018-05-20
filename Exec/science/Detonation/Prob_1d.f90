subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module, only: T_l, T_r, dens, cfrac, idir, w_T, center_T, &
                             xn, ihe4, ic12, io16, smallx
  use network, only: network_species_index, nspec
  use bl_error_module, only: bl_error
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(1), probhi(1)

  integer :: untin,i

  namelist /fortin/ T_l, T_r, dens, cfrac, idir, w_T, center_T, smallx

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Set namelist defaults

  T_l = 1.e9_rt
  T_r = 5.e7_rt
  dens = 1.e8_rt
  smallx = 1.e-12_rt

  idir = 1                ! direction across which to jump
  cfrac = 0.5

  w_T = 5.e-4_rt          ! ratio of the width of temperature transition zone to the full domain
  center_T = 3.e-1_rt     ! central position parameter of teperature profile transition zone

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! get the species indices
  ihe4 = network_species_index("helium-4")
  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")

  if (ihe4 < 0 .or. ic12 < 0 .or. io16 < 0) then
     call bl_error("ERROR: species indices not found")
  endif

  ! make sure that the carbon fraction falls between 0 and 1
  if (cfrac > 1.e0_rt .or. cfrac < 0.e0_rt) then
     call bl_error("ERROR: cfrac must fall between 0 and 1")
  endif

  ! set the default mass fractions
  allocate(xn(nspec))

  xn(:) = smallx
  xn(ic12) = cfrac
  xn(ihe4) = 1.e0_rt - cfrac - (nspec - 1)*smallx

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
                       state,state_l1,state_h1,delta,xlo,xhi)

  use network, only: nspec
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rt
  use probdata_module, only: T_l, T_r, center_T, w_T, dens, xn
  use meth_params_module, only: NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFS
  use amrex_fort_module, only: rt => amrex_real
  use prob_params_module, only: problo, probhi

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(1), hi(1)
  integer,  intent(in   ) :: state_l1,state_h1
  real(rt), intent(inout) :: state(state_l1:state_h1,NVAR)
  real(rt), intent(in   ) :: time, delta(1)
  real(rt), intent(in   ) :: xlo(1), xhi(1)

  real(rt) :: sigma, width, c_T
  real(rt) :: xcen
  integer  :: i

  type (eos_t) :: eos_state

  width = w_T * (probhi(1) - problo(1))
  c_T = problo(1) + center_T * (probhi(1) - problo(1))

  do i = lo(1), hi(1)
     xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5e0_rt)

     state(i,URHO ) = dens

     sigma = 1.0 / (1.0 + exp(-(c_T - xcen)/ width))

     state(i,UTEMP) = T_l + (T_r - T_l) * (1 - sigma)

     state(i,UFS:UFS-1+nspec) = state(i,URHO)*xn(1:nspec)

     eos_state%rho = state(i,URHO)
     eos_state%T = state(i,UTEMP)
     eos_state%xn(:) = xn

     call eos(eos_input_rt, eos_state)

     state(i,UMX  ) = 0.e0_rt
     state(i,UEDEN) = state(i,URHO)*eos_state%e  ! if vel /= 0, then KE needs to be added
     state(i,UEINT) = state(i,URHO)*eos_state%e

  enddo

end subroutine ca_initdata
