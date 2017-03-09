subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use bl_constants_module
  use probdata_module
  use prob_params_module, only: center
  use bl_error_module

  use amrex_fort_module, only : rt => c_real
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt)         :: problo(2), probhi(2)

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

  ! Set namelist defaults

  rho_ambient = 1.e0_rt
  rho_peak = 2.e0_rt
  t_ambient = 1.e0_rt
  sigma = 0.1e0_rt

  !     Read namelists
  open(newunit=untin, file=probin(1:namlen), &
       form='formatted',status='old')
  read(untin, fortin)
  close(unit=untin)

  !     Set explosion center
  center(1) = HALF*(problo(1)+probhi(1))
  center(2) = HALF*(problo(2)+probhi(2))

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

  use bl_constants_module
  use probdata_module
  use eos_type_module
  use eos_module
  use meth_params_module , only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module

  use amrex_fort_module, only : rt => c_real
  implicit none

  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  real(rt)         :: xlo(2), xhi(2), time, delta(2)
  real(rt)         :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

  real(rt)         :: xx, yy

  integer :: i, j
  type(eos_t) :: eos_state

  do j = lo(2), hi(2)
     yy = problo(2) + delta(2)*dble(j+HALF)

     do i = lo(1), hi(1)
        xx = problo(1) + delta(1)*dble(i+HALF)

        state(i,j,URHO) = rho_ambient + &
             (rho_peak-rho_ambient)*exp(-((xx-center(1))**2 + &
                                          (yy-center(2))**2)/sigma**2)
        state(i,j,UMX:UMZ) = ZERO

        state(i,j,UFS) = state(i,j,URHO)

        state(i,j,UTEMP) = t_ambient

        eos_state % rho = state(i,j,URHO)
        eos_state % T = state(i,j,UTEMP)
        eos_state % xn(:) = state(i,j,UFS:UFS-1+nspec)/ eos_state % rho

        call eos(eos_input_rt, eos_state)

        ! assuming no KE
        state(i,j,UEDEN) = eos_state % rho * eos_state % e
        state(i,j,UEINT) = eos_state % rho * eos_state % e

     enddo
  enddo

end subroutine ca_initdata
