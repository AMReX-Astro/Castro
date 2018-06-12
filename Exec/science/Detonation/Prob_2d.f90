subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module, only: T_l, T_r, dens, cfrac, idir, w_T, center_T, &
                             xn, ihe4, ic12, io16, smallx, vel
  use network, only: network_species_index, nspec
  use amrex_error_module, only: amrex_error
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(2), probhi(2)

  integer :: untin,i

  namelist /fortin/ T_l, T_r, dens, cfrac, idir, w_T, center_T, smallx, vel

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) call amrex_error("probin file name too long")

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

  w_T = 5.e-4_rt           ! ratio of the width of temperature transition zone to the full domain
  center_T = 3.e-1_rt      ! central position parameter of teperature profile transition zone

  vel = 0.e0_rt           ! infall velocity towards the transition point

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
     call amrex_error("ERROR: species indices not found")
  endif

  ! make sure that the carbon fraction falls between 0 and 1
  if (cfrac > 1.e0_rt .or. cfrac < 0.e0_rt) then
     call amrex_error("ERROR: cfrac must fall between 0 and 1")
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
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use network, only: nspec
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rt
  use probdata_module, only: T_l, T_r, center_T, w_T, dens, vel, xn
  use meth_params_module, only: NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
  use amrex_fort_module, only: rt => amrex_real
  use prob_params_module, only: problo, probhi

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(2), hi(2)
  integer,  intent(in   ) :: state_l1,state_l2,state_h1,state_h2
  real(rt), intent(inout) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  real(rt), intent(in   ) :: time, delta(2)
  real(rt), intent(in   ) :: xlo(2), xhi(2)

  real(rt) :: sigma, width, c_T
  real(rt) :: xcen
  integer  :: i, j

  type (eos_t) :: eos_state
  
  width = w_T * (probhi(1) - problo(1))
  c_T = problo(1) + center_T * (probhi(1) - problo(1))

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5e0_rt)

        state(i,j,URHO ) = dens

        sigma = 1.0 / (1.0 + exp(-(c_T - xcen)/ width))

        state(i,j,UTEMP) = T_l + (T_r - T_l) * (1 - sigma)

        state(i,j,UFS:UFS-1+nspec) = state(i,j,URHO)*xn(1:nspec)

        eos_state%rho = state(i,j,URHO)
        eos_state%T = state(i,j,UTEMP)
        eos_state%xn(:) = xn

        call eos(eos_input_rt, eos_state)

        state(i,j,UMX  ) = state(i,j,URHO) * (vel - 2 * vel * (1.0e0_rt - sigma))
        state(i,j,UMY  ) = 0.e0_rt
        state(i,j,UEINT) = state(i,j,URHO) * eos_state%e
        state(i,j,UEDEN) = state(i,j,UEINT) + 0.5e0_rt * sum(state(i,j,UMX:UMY)**2) / state(i,j,URHO)
     enddo
  enddo

end subroutine ca_initdata
