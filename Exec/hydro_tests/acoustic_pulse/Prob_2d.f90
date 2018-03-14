subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use bl_constants_module
  use probdata_module
  use prob_params_module, only : center
  use bl_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(2), probhi(2)

  integer :: untin, i

  namelist /fortin/ rho0, drho0

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

  rho0 = 1.4_rt
  drho0 = 0.14_rt

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(unit=untin)

  ! set explosion center
  center(1) = HALF*(problo(1) + probhi(1))
  center(2) = HALF*(problo(2) + probhi(2))

  xn_zone(:) = ZERO
  xn_zone(1) = ONE

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
  use bl_constants_module, only: M_PI, FOUR3RD, ZERO, HALF, ONE
  use meth_params_module , only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS
  use prob_params_module, only : center, coord_type
  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec
  use extern_probin_module, only : eos_gamma

  implicit none

  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  real(rt) :: xlo(2), xhi(2), time, delta(2)
  real(rt) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

  real(rt) :: xmin,ymin
  real(rt) :: xx, yy
  real(rt) :: dist, p, eint

  integer :: i,j


  do j = lo(2), hi(2)
     yy = xlo(2) + delta(2)*dble(j-lo(2) + HALF)

     do i = lo(1), hi(1)
        xx = xlo(1) + delta(1)*dble(i-lo(1) + HALF)

        dist = sqrt((center(1)-xx)**2 + (center(2)-yy)**2)

        if (dist <= HALF) then
           state(i,j,URHO) = rho0 + drho0*exp(-16.d0*dist**2) * cos(M_PI*dist)**6
        else
           state(i,j,URHO) = rho0
        endif

        state(i,j,UMX:UMZ) = 0.e0_rt

        ! we are isentropic, so p = (dens/rho0)**Gamma_1
        p = (state(i,j,URHO)/rho0)**eos_gamma
        eint = p/(eos_gamma - ONE)

        state(i,j,UEDEN) = eint
        state(i,j,UEINT) = eint

        state(i,j,UFS:UFS-1+nspec) = state(i,j,URHO)*xn_zone(:)

     enddo
  enddo

end subroutine ca_initdata
