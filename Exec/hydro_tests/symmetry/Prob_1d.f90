subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use bl_constants_module
  use probdata_module
  use prob_params_module, only : center
  use bl_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt)         :: problo(1), probhi(1)

  integer :: untin,i

  namelist /fortin/ rho_ambient, rho_peak, t_ambient, sigma

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  rho_ambient = 1.e0_rt
  rho_peak = 2.e0_rt
  t_ambient = 1.e0_rt
  sigma = 0.1e0_rt

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), &
       form='formatted', status='old')
  read(untin, fortin)
  close(unit=untin)

  center(1) = HALF*(problo(1) + probhi(1))

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

  use bl_constants_module
  use probdata_module
  use eos_type_module
  use eos_module
  use meth_params_module , only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: level, nscal
  integer :: lo(1), hi(1)
  integer :: state_l1,state_h1
  real(rt)         :: xlo(1), xhi(1), time, delta(1)
  real(rt)         :: state(state_l1:state_h1,NVAR)

  real(rt)         :: xx

  integer :: i
  type(eos_t) :: eos_state

  do i = lo(1), hi(1)
     xx = problo(1) + delta(1)*dble(i+HALF)

     state(i,URHO) = rho_ambient + (rho_peak-rho_ambient)*exp(-(xx-center(1))**2/sigma**2)
     state(i,UMX:UMZ) = ZERO

     state(i,UFS) = state(i,URHO)

     state(i,UTEMP) = t_ambient

     eos_state % rho = state(i,URHO)
     eos_state % T = state(i,UTEMP)
     eos_state % xn(:) = state(i,UFS:UFS-1+nspec)/ eos_state % rho

     call eos(eos_input_rt, eos_state)

     ! assuming no KE
     state(i,UEDEN) = eos_state % rho * eos_state % e
     state(i,UEINT) = eos_state % rho * eos_state % e

  end do

end subroutine ca_initdata
