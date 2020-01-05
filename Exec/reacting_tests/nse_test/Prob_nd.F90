subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use amrex_constants_module
  use probdata_module, only: T_min, T_max, rho_ambient, width, xn, cfrac, ofrac
  use network, only: network_species_index, nspec
  use castro_error_module, only: castro_error
  use amrex_fort_module, only: rt => amrex_real
  use prob_params_module, only : center

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin,i

  integer, save :: ihe4, ic12, io16

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  real(rt) :: smallx

  call probdata_init(name, namlen)

  ! get the species indices
  ihe4 = network_species_index("helium-4")
  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")

  if (ihe4 < 0 .or. ic12 < 0 .or. io16 < 0) then
     call castro_error("ERROR: species indices not found")
  endif

  ! make sure that the carbon fraction falls between 0 and 1
  if (cfrac > 1.e0_rt .or. cfrac < 0.e0_rt) then
     call castro_error("ERROR: cfrac must fall between 0 and 1")
  endif

  ! make sure that the oxygen fraction falls between 0 and 1
  if (ofrac > 1.e0_rt .or. cfrac < 0.e0_rt) then
     call castro_error("ERROR: ofrac must fall between 0 and 1")
  endif

  ! make sure that the C/O fraction sums to no more than 1
  if (cfrac + ofrac > 1.e0_rt) then
     call castro_error("ERROR: cfrac + ofrac cannot exceed 1.")
  end if

  ! set the default mass fractions
  allocate(xn(nspec))

  xn(:) = smallx
  xn(ic12) = max(cfrac, smallx)
  xn(io16) = max(ofrac, smallx)
  xn(ihe4) = 1.e0_rt - cfrac - ofrac - (nspec - 2) * smallx

  center(:) = HALF*(problo(:) + probhi(:))

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
                       state,state_lo,state_hi, &
                       delta,xlo,xhi)

  use amrex_constants_module
  use network, only: nspec
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rt
  use probdata_module, only: T_min, T_max, rho_ambient, width, xn
  use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
  use amrex_fort_module, only: rt => amrex_real
  use prob_params_module, only: problo, probhi, center

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: state_lo(3), state_hi(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
  real(rt), intent(in   ) :: time, delta(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3)

  integer  :: i, j, k

  type (eos_t) :: eos_state

  real(rt) :: xcen, ycen, zcen, r
  real(rt) :: T

  do k = lo(3), hi(3)
     zcen = problo(3) + delta(3)*(dble(k) + HALF) - center(3)

     do j = lo(2), hi(2)
        xcen = problo(2) + delta(2)*(dble(j) + HALF) - center(2)

        do i = lo(1), hi(1)
           xcen = problo(1) + delta(1)*(dble(i) + HALF) - center(1)

           r = sqrt(xcen**2 + ycen**2 + zcen**2)

           T = T_min + (T_max - T_min)**exp(-(r/width)**2)

           state(i,j,k,URHO ) = rho_ambient
           state(i,j,k,URHO ) = T

           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO)*xn(1:nspec)

           eos_state % rho = state(i,j,k,URHO)
           eos_state % T = state(i,j,k,UTEMP)
           eos_state % xn(:) = xn(:)

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UMX:UMZ) = ZERO
           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * sum(state(i,j,k,UMX:UMZ)**2) / state(i,j,k,URHO)
        enddo
     enddo
  enddo

end subroutine ca_initdata
